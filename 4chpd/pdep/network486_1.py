species(
    label = 'O=[C]CC1[CH]C=CCCC1(2233)',
    structure = SMILES('O=[C]CC1[CH]C=CCCC1'),
    E0 = (138.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.288055,0.0260298,0.000201692,-3.15433e-07,1.33664e-10,16865,34.7321], Tmin=(100,'K'), Tmax=(912.077,'K')), NASAPolynomial(coeffs=[40.0235,-0.003371,1.18031e-05,-2.49925e-09,1.57314e-13,3591.2,-186.343], Tmin=(912.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(138.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(CCCJ=O) + radical(Allyl_S)"""),
)

species(
    label = 'CHPD(1)',
    structure = SMILES('C1C=CCCCC=1'),
    E0 = (72.8504,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3862.82,'J/mol'), sigma=(6.55596,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=603.36 K, Pc=31.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.815077,0.0542245,7.00763e-06,-4.80543e-08,2.21793e-11,8890.61,7.61919], Tmin=(100,'K'), Tmax=(987.481,'K')), NASAPolynomial(coeffs=[16.3155,0.0286278,-1.06044e-05,1.97595e-09,-1.42758e-13,4016.04,-76.1476], Tmin=(987.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.8504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = 'CH2CO(43)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=C1CC2[CH]C1[CH]CCC2(2695)',
    structure = SMILES('O=C1CC2[CH]C1[CH]CCC2'),
    E0 = (96.0236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.696066,0.10057,-7.80906e-05,3.15583e-08,-5.27614e-12,11719.8,4.72012], Tmin=(100,'K'), Tmax=(1381.59,'K')), NASAPolynomial(coeffs=[17.5469,0.0477523,-2.07461e-05,3.8874e-09,-2.69041e-13,6678.98,-89.1881], Tmin=(1381.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.0236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + polycyclic(PolycyclicRing) + radical(CCJCC=O) + radical(CCJCC=O)"""),
)

species(
    label = 'H(19)',
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
    label = 'O=[C]CC1=CC=CCCC1(2696)',
    structure = SMILES('O=[C]CC1=CC=CCCC1'),
    E0 = (94.4935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.883446,0.0905269,-4.7504e-05,-3.29259e-09,7.08375e-12,11555.1,20.493], Tmin=(100,'K'), Tmax=(1094.9,'K')), NASAPolynomial(coeffs=[22.2424,0.0367515,-1.59052e-05,3.0852e-09,-2.21896e-13,4650.21,-101.579], Tmin=(1094.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.4935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C=CC1[CH]C=CCCC1(2697)',
    structure = SMILES('O=C=CC1C=C[CH]CCC1'),
    E0 = (102.832,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2120,512.5,787.5,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.199291,0.027215,0.000197963,-3.15174e-07,1.34707e-10,12557.5,31.4317], Tmin=(100,'K'), Tmax=(910.393,'K')), NASAPolynomial(coeffs=[41.8256,-0.00845435,1.416e-05,-2.94532e-09,1.87999e-13,-1122.9,-198.992], Tmin=(910.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCd(CCO)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + ring(Cycloheptene) + radical(Allyl_S)"""),
)

species(
    label = 'O=[C]CC1C=C=CCCC1(2698)',
    structure = SMILES('O=[C]CC1C=C=CCCC1'),
    E0 = (218.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.110382,0.0785486,-3.97271e-05,4.38724e-09,1.61429e-12,26382.5,30.9035], Tmin=(100,'K'), Tmax=(1257.75,'K')), NASAPolynomial(coeffs=[15.7464,0.0443997,-1.8417e-05,3.38337e-09,-2.31803e-13,21106,-54.3518], Tmin=(1257.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]=O(360)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97567e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02477e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]1C=C[CH]CCC1(93)',
    structure = SMILES('[CH]1C=C[CH]CCC1'),
    E0 = (259.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17636,-0.0217516,0.000290736,-3.97461e-07,1.6188e-10,31302.2,19.7886], Tmin=(100,'K'), Tmax=(909.243,'K')), NASAPolynomial(coeffs=[37.6598,-0.0173238,1.86036e-05,-3.75571e-09,2.40281e-13,18214,-184.512], Tmin=(909.243,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(Allyl_S)"""),
)

species(
    label = 'O=[C]C[C]1CC=CCCC1(2699)',
    structure = SMILES('O=[C]C[C]1CC=CCCC1'),
    E0 = (182.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.192832,0.038488,0.000146435,-2.45522e-07,1.05814e-10,22188,38.0079], Tmin=(100,'K'), Tmax=(912.374,'K')), NASAPolynomial(coeffs=[34.6512,0.00547643,6.60991e-06,-1.52634e-09,9.5007e-14,10986.4,-152.003], Tmin=(912.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(Tertalkyl) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C[CH]C1[CH]C=CCCC1(2700)',
    structure = SMILES('[O]C=CC1[CH]C=CCCC1'),
    E0 = (123.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.268141,0.0125795,0.000269363,-4.05532e-07,1.70665e-10,14997.8,33.732], Tmin=(100,'K'), Tmax=(909.978,'K')), NASAPolynomial(coeffs=[48.2927,-0.017344,2.00359e-05,-4.07181e-09,2.60546e-13,-1243.77,-234.646], Tmin=(909.978,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]CC1C[C]=CCCC1(2701)',
    structure = SMILES('O=[C]CC1C[C]=CCCC1'),
    E0 = (235.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0625073,0.0357883,0.000168172,-2.77938e-07,1.19819e-10,28502.5,37.122], Tmin=(100,'K'), Tmax=(910.55,'K')), NASAPolynomial(coeffs=[38.7419,-0.00130186,1.04611e-05,-2.26401e-09,1.43982e-13,15952.3,-176.095], Tmin=(910.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(CCCJ=O) + radical(Cds_S)"""),
)

species(
    label = 'O=[C]CC1[CH]CCC=CC1(2702)',
    structure = SMILES('O=[C]CC1[CH]CCC=CC1'),
    E0 = (192.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.28411,0.0300065,0.000181376,-2.88645e-07,1.22653e-10,23287.6,38.4085], Tmin=(100,'K'), Tmax=(913.718,'K')), NASAPolynomial(coeffs=[37.8696,7.10418e-05,9.5492e-06,-2.05251e-09,1.27331e-13,10800.2,-170.275], Tmin=(913.718,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(Cs_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C][CH]C1CC=CCCC1(2703)',
    structure = SMILES('O=[C][CH]C1CC=CCCC1'),
    E0 = (165.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,1855,455,950,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.249188,0.0249442,0.000210042,-3.29255e-07,1.40273e-10,20045.6,35.7398], Tmin=(100,'K'), Tmax=(908.006,'K')), NASAPolynomial(coeffs=[41.6921,-0.0071154,1.43706e-05,-3.04243e-09,1.96323e-13,6315.04,-194.364], Tmin=(908.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CC[C]1[CH]C=CCCC1(2704)',
    structure = SMILES('O=CC[C]1C=C[CH]CCC1'),
    E0 = (111.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621238,0.0175934,0.000222975,-3.33725e-07,1.39115e-10,13617.5,31.6984], Tmin=(100,'K'), Tmax=(913.658,'K')), NASAPolynomial(coeffs=[38.3336,-0.0002384,1.04651e-05,-2.24013e-09,1.38528e-13,579.248,-180.476], Tmin=(913.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.771,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(Allyl_T) + radical(Allyl_S)"""),
)

species(
    label = 'O=[C]CC1CC=[C]CCC1(2705)',
    structure = SMILES('O=[C]CC1CC=[C]CCC1'),
    E0 = (235.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0625073,0.0357883,0.000168172,-2.77938e-07,1.19819e-10,28502.5,37.122], Tmin=(100,'K'), Tmax=(910.55,'K')), NASAPolynomial(coeffs=[38.7419,-0.00130186,1.04611e-05,-2.26401e-09,1.43982e-13,15952.3,-176.095], Tmin=(910.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(CCCJ=O) + radical(Cds_S)"""),
)

species(
    label = 'O=CCC1[CH]CC[CH]C=C1(2648)',
    structure = SMILES('O=CCC1[CH]C=CCC[CH]1'),
    E0 = (173.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.542445,0.021981,0.000204712,-3.10373e-07,1.29321e-10,21013.4,35.5552], Tmin=(100,'K'), Tmax=(918.217,'K')), NASAPolynomial(coeffs=[37.1845,0.00230002,8.25407e-06,-1.75492e-09,1.03539e-13,8384.95,-170.22], Tmin=(918.217,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(Cs_S) + radical(Allyl_S)"""),
)

species(
    label = 'O=[C]CC1C[CH]CC=CC1(2706)',
    structure = SMILES('O=[C]CC1C[CH]CC=CC1'),
    E0 = (192.026,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.369322,0.0314205,0.000169914,-2.72314e-07,1.16059e-10,23271.4,38.1111], Tmin=(100,'K'), Tmax=(911.25,'K')), NASAPolynomial(coeffs=[35.5373,0.00337335,8.1388e-06,-1.83042e-09,1.15154e-13,11617.2,-157.064], Tmin=(911.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(CCCJ=O) + radical(RCCJCC)"""),
)

species(
    label = 'O=[C]CC1CC=C[CH]CC1(2707)',
    structure = SMILES('O=[C]CC1CC=C[CH]CC1'),
    E0 = (138.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.288055,0.0260298,0.000201692,-3.15433e-07,1.33664e-10,16865,34.7321], Tmin=(100,'K'), Tmax=(912.077,'K')), NASAPolynomial(coeffs=[40.0235,-0.003371,1.18031e-05,-2.49925e-09,1.57314e-13,3591.2,-186.343], Tmin=(912.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(138.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(CCCJ=O) + radical(Allyl_S)"""),
)

species(
    label = 'O=CCC1[CH]C=CC[CH]C1(2708)',
    structure = SMILES('O=CCC1[CH]C=CC[CH]C1'),
    E0 = (173.178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.627209,0.0234017,0.000193218,-2.93986e-07,1.22694e-10,20997.3,35.2593], Tmin=(100,'K'), Tmax=(916.131,'K')), NASAPolynomial(coeffs=[34.8489,0.00560785,6.84055e-06,-1.53212e-09,9.1303e-14,9203.37,-156.989], Tmin=(916.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.178,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(RCCJCC)"""),
)

species(
    label = 'O=CCC1[CH][C]=CCCC1(2709)',
    structure = SMILES('O=CCC1[CH][C]=CCCC1'),
    E0 = (216.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.319709,0.0277785,0.00019144,-2.99555e-07,1.26427e-10,26228.4,34.2726], Tmin=(100,'K'), Tmax=(915.25,'K')), NASAPolynomial(coeffs=[38.0538,0.000931983,9.16334e-06,-1.96582e-09,1.20141e-13,13538.5,-176.022], Tmin=(915.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'O=CCC1[CH]C=C[CH]CC1(2710)',
    structure = SMILES('O=CCC1[CH]C=C[CH]CC1'),
    E0 = (119.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.546775,0.0180004,0.000225038,-3.37165e-07,1.40329e-10,14590.9,31.8773], Tmin=(100,'K'), Tmax=(916.299,'K')), NASAPolynomial(coeffs=[39.3338,-0.00113407,1.05034e-05,-2.20059e-09,1.33432e-13,1177.95,-186.261], Tmin=(916.299,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(Allyl_S)"""),
)

species(
    label = 'O=CCC1[CH]C=[C]CCC1(2711)',
    structure = SMILES('O=CCC1[CH]C=[C]CCC1'),
    E0 = (216.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.319709,0.0277785,0.00019144,-2.99555e-07,1.26427e-10,26228.4,34.2726], Tmin=(100,'K'), Tmax=(915.25,'K')), NASAPolynomial(coeffs=[38.0538,0.000931983,9.16334e-06,-1.96582e-09,1.20141e-13,13538.5,-176.022], Tmin=(915.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = 'O=[C]CC1CCCC2[CH]C21(2712)',
    structure = SMILES('O=[C]CC1CCCC2[CH]C21'),
    E0 = (271.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.344883,0.0676241,-8.7003e-06,-2.27206e-08,9.56405e-12,32822.2,33.2408], Tmin=(100,'K'), Tmax=(1170.26,'K')), NASAPolynomial(coeffs=[12.7356,0.05165,-2.20358e-05,4.13741e-09,-2.8831e-13,28115.9,-36.2025], Tmin=(1170.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + polycyclic(s2_3_6_ane) + radical(bicyclo[4.1.0]heptane-C3-7) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C1CC2[CH][CH]C1CCC2(2713)',
    structure = SMILES('O=C1CC2[CH][CH]C1CCC2'),
    E0 = (65.5504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.250525,0.10058,-8.41183e-05,3.99125e-08,-8.30648e-12,8030.97,2.98843], Tmin=(100,'K'), Tmax=(1089.99,'K')), NASAPolynomial(coeffs=[10.9411,0.0595098,-2.75993e-05,5.34412e-09,-3.77915e-13,5591.21,-51.9691], Tmin=(1089.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.5504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + polycyclic(PolycyclicRing) + radical(cycloheptane) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]CCC=CC=CC[C]=O(2714)',
    structure = SMILES('[CH2]CCC=CC=CC[C]=O'),
    E0 = (241.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,1855,455,950,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.969234,0.0944933,-6.31602e-05,1.67957e-08,-6.67739e-13,29223.5,41.7262], Tmin=(100,'K'), Tmax=(1236.09,'K')), NASAPolynomial(coeffs=[21.6849,0.0386643,-1.66238e-05,3.13749e-09,-2.19191e-13,22287.6,-77.7706], Tmin=(1236.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(RCCJ)"""),
)

species(
    label = 'O=CCC1=CC=CCCC1(2715)',
    structure = SMILES('O=CCC1=CC=CCCC1'),
    E0 = (-65.4671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.813019,0.0875653,-3.2346e-05,-1.79028e-08,1.15497e-11,-7685.11,19.4231], Tmin=(100,'K'), Tmax=(1092,'K')), NASAPolynomial(coeffs=[21.6298,0.041281,-1.81145e-05,3.53434e-09,-2.54866e-13,-14728.5,-100.633], Tmin=(1092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.4671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = 'O=C=CC1CC=CCCC1(2716)',
    structure = SMILES('O=C=CC1CC=CCCC1'),
    E0 = (-37.3446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.224076,0.0256091,0.000207716,-3.24885e-07,1.37762e-10,-4301.91,34.2987], Tmin=(100,'K'), Tmax=(912.146,'K')), NASAPolynomial(coeffs=[41.4628,-0.00560288,1.29801e-05,-2.7152e-09,1.71185e-13,-18049.8,-194.983], Tmin=(912.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-37.3446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + ring(Cycloheptene)"""),
)

species(
    label = 'O=CCC1C=C=CCCC1(2717)',
    structure = SMILES('O=CCC1C=C=CCCC1'),
    E0 = (58.0852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00376762,0.0752156,-2.35202e-05,-1.12481e-08,6.38554e-12,7140.67,29.6999], Tmin=(100,'K'), Tmax=(1190.94,'K')), NASAPolynomial(coeffs=[14.5188,0.0498821,-2.11393e-05,3.9476e-09,-2.73948e-13,2019.05,-49.8805], Tmin=(1190.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.0852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = 'O=[C]C[CH]C1C=CCCC1(2658)',
    structure = SMILES('O=[C]C[CH]C1C=CCCC1'),
    E0 = (176.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102406,0.0671597,1.06855e-05,-5.62124e-08,2.45172e-11,21382.4,36.3434], Tmin=(100,'K'), Tmax=(1016.24,'K')), NASAPolynomial(coeffs=[16.9063,0.043526,-1.71736e-05,3.22399e-09,-2.30202e-13,15772,-55.796], Tmin=(1016.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclohexene) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1[CH]CCCC=C1(2232)',
    structure = SMILES('O=[C]CC1[CH]CCCC=C1'),
    E0 = (193.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4431.72,'J/mol'), sigma=(7.51402,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=692.22 K, Pc=23.7 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281231,0.0297574,0.00018278,-2.90754e-07,1.23603e-10,23416.6,38.4535], Tmin=(100,'K'), Tmax=(913.276,'K')), NASAPolynomial(coeffs=[38.1056,-0.00043822,9.87456e-06,-2.11907e-09,1.31997e-13,10858.2,-171.526], Tmin=(913.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(CCCJ=O) + radical(Cs_S)"""),
)

species(
    label = 'C=C([O])C1[CH]C=CCCC1(2229)',
    structure = SMILES('C=C([O])C1[CH]C=CCCC1'),
    E0 = (113.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.260949,0.018523,0.000240441,-3.68747e-07,1.56178e-10,13859.3,34.0469], Tmin=(100,'K'), Tmax=(908.535,'K')), NASAPolynomial(coeffs=[45.2131,-0.0129624,1.76543e-05,-3.64977e-09,2.35042e-13,-1177.45,-216.31], Tmin=(908.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=C(C)OJ) + radical(Allyl_S)"""),
)

species(
    label = 'O=C1CC2CCCC=CC12(2236)',
    structure = SMILES('O=C1CC2CCCC=CC12'),
    E0 = (-192.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.412999,0.0997174,-7.75791e-05,3.21613e-08,-5.65384e-12,-22987.4,-1.40441], Tmin=(100,'K'), Tmax=(1291.99,'K')), NASAPolynomial(coeffs=[14.4432,0.0537229,-2.41796e-05,4.60725e-09,-3.22149e-13,-26826.2,-76.8828], Tmin=(1291.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(Cycloheptene) + ring(Cyclobutane)"""),
)

species(
    label = 'CO(28)',
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
    label = '[CH2]C1[CH]C=CCCC1(499)',
    structure = SMILES('[CH2]C1[CH]C=CCCC1'),
    E0 = (290.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32877,-0.00220167,0.000264302,-3.79443e-07,1.57582e-10,35078.2,27.6852], Tmin=(100,'K'), Tmax=(905.357,'K')), NASAPolynomial(coeffs=[39.594,-0.0122333,1.74407e-05,-3.6485e-09,2.37632e-13,21631.9,-189.113], Tmin=(905.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(Isobutyl)"""),
)

species(
    label = '[C]=O(353)',
    structure = SMILES('[C]=O'),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08916,0.00200416,-1.61661e-05,2.55058e-08,-1.16424e-11,52802.7,4.52505], Tmin=(100,'K'), Tmax=(856.11,'K')), NASAPolynomial(coeffs=[0.961625,0.00569045,-3.48044e-06,7.19202e-10,-5.08041e-14,53738.7,21.4663], Tmin=(856.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'O=[C]CC1C=CC=CCC1(2718)',
    structure = SMILES('O=[C]CC1C=CC=CCC1'),
    E0 = (94.4792,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17938,0.0998092,-8.17366e-05,3.42211e-08,-5.71623e-12,11561.4,21.8886], Tmin=(100,'K'), Tmax=(1434.82,'K')), NASAPolynomial(coeffs=[22.3949,0.0340882,-1.30297e-05,2.29731e-09,-1.53861e-13,4796.5,-100.355], Tmin=(1434.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.4792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1C=CC[CH]CC1(2641)',
    structure = SMILES('O=[C]CC1C=CC[CH]CC1'),
    E0 = (193.096,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.366351,0.0311724,0.000171314,-2.7442e-07,1.17008e-10,23400.5,38.1564], Tmin=(100,'K'), Tmax=(910.803,'K')), NASAPolynomial(coeffs=[35.7741,0.00286285,8.46487e-06,-1.89716e-09,1.19835e-13,11674.9,-158.319], Tmin=(910.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(CCCJ=O) + radical(RCCJCC)"""),
)

species(
    label = 'O=[C]CC1C=[C]CCCC1(2644)',
    structure = SMILES('O=[C]CC1C=[C]CCCC1'),
    E0 = (236.479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0592302,0.0355441,0.000169558,-2.80024e-07,1.20759e-10,28631.6,37.1684], Tmin=(100,'K'), Tmax=(910.142,'K')), NASAPolynomial(coeffs=[38.9796,-0.00181407,1.07882e-05,-2.33099e-09,1.48683e-13,16009.7,-177.356], Tmin=(910.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1C=CCC[CH]C1(2637)',
    structure = SMILES('O=[C]CC1C=CCC[CH]C1'),
    E0 = (193.096,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.366351,0.0311724,0.000171314,-2.7442e-07,1.17008e-10,23400.5,38.1564], Tmin=(100,'K'), Tmax=(910.803,'K')), NASAPolynomial(coeffs=[35.7741,0.00286285,8.46487e-06,-1.89716e-09,1.19835e-13,11674.9,-158.319], Tmin=(910.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(RCCJCC) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1[C]=CCCCC1(2640)',
    structure = SMILES('O=[C]CC1[C]=CCCCC1'),
    E0 = (236.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0592302,0.0355441,0.000169558,-2.80024e-07,1.20759e-10,28631.6,37.1684], Tmin=(100,'K'), Tmax=(910.142,'K')), NASAPolynomial(coeffs=[38.9796,-0.00181407,1.07882e-05,-2.33099e-09,1.48683e-13,16009.7,-177.356], Tmin=(910.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CCC1[C]=C[CH]CCC1(2719)',
    structure = SMILES('O=CCC1[C]=C[CH]CCC1'),
    E0 = (217.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.317501,0.0275212,0.000192875,-3.01706e-07,1.27396e-10,26357.5,34.3152], Tmin=(100,'K'), Tmax=(914.781,'K')), NASAPolynomial(coeffs=[38.2876,0.000426555,9.48644e-06,-2.03185e-09,1.24762e-13,13597.4,-177.26], Tmin=(914.781,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = 'O=[C]C[C]1C=CCCCC1(2636)',
    structure = SMILES('O=[C]C[C]1C=CCCCC1'),
    E0 = (130.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.362514,0.0256214,0.000199643,-3.12025e-07,1.32472e-10,15891.6,34.5532], Tmin=(100,'K'), Tmax=(909.276,'K')), NASAPolynomial(coeffs=[39.0289,-0.00248471,1.17701e-05,-2.54005e-09,1.62513e-13,2990.13,-180.59], Tmin=(909.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(CCCJ=O) + radical(Allyl_T)"""),
)

species(
    label = 'O=[C][CH]C1C=CCCCC1(2639)',
    structure = SMILES('O=[C][CH]C1C=CCCCC1'),
    E0 = (166.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.246399,0.0246937,0.000211452,-3.31376e-07,1.4123e-10,20174.6,35.7845], Tmin=(100,'K'), Tmax=(907.644,'K')), NASAPolynomial(coeffs=[41.9289,-0.007626,1.46967e-05,-3.10917e-09,2.01004e-13,6372.75,-195.619], Tmin=(907.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C1CC2CCC[CH][CH]C12(2720)',
    structure = SMILES('O=C1CC2CCC[CH][CH]C12'),
    E0 = (65.5504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.250525,0.10058,-8.41183e-05,3.99125e-08,-8.30648e-12,8030.97,2.98843], Tmin=(100,'K'), Tmax=(1089.99,'K')), NASAPolynomial(coeffs=[10.9411,0.0595098,-2.75993e-05,5.34412e-09,-3.77915e-13,5591.21,-51.9691], Tmin=(1089.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.5504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + polycyclic(PolycyclicRing) + radical(CCJCC=O) + radical(cycloheptane)"""),
)

species(
    label = '[CH2]CC(C=CC=C)C[C]=O(2721)',
    structure = SMILES('[CH2]CC(C=CC=C)C[C]=O'),
    E0 = (238.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.897539,0.101962,-8.74955e-05,3.98679e-08,-7.40362e-12,28847.3,41.4187], Tmin=(100,'K'), Tmax=(1278.49,'K')), NASAPolynomial(coeffs=[18.3204,0.041835,-1.69519e-05,3.08334e-09,-2.10713e-13,23933.3,-56.0186], Tmin=(1278.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C=CC1C=CCCCC1(2656)',
    structure = SMILES('O=C=CC1C=CCCCC1'),
    E0 = (-38.2802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.000117643,0.0324513,0.000188969,-3.06633e-07,1.3171e-10,-4408.1,33.2532], Tmin=(100,'K'), Tmax=(910.329,'K')), NASAPolynomial(coeffs=[41.7531,-0.00591857,1.31148e-05,-2.7636e-09,1.76369e-13,-18021.8,-197.279], Tmin=(910.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-38.2802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCd(CCO)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + ring(Cycloheptene)"""),
)

species(
    label = 'O=CCC1C=CC=CCC1(2722)',
    structure = SMILES('O=CCC1C=CC=CCC1'),
    E0 = (-65.4814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.97872,0.0954409,-6.22609e-05,1.48527e-08,4.34108e-13,-7684.58,20.3429], Tmin=(100,'K'), Tmax=(1145.42,'K')), NASAPolynomial(coeffs=[20.0633,0.0413014,-1.66928e-05,3.07423e-09,-2.13069e-13,-13773.8,-89.5686], Tmin=(1145.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.4814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = '[CH2]C1C=CC(C[C]=O)CC1(2723)',
    structure = SMILES('[CH2]C1C=CC(C[C]=O)CC1'),
    E0 = (173.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0161283,0.0707429,2.83809e-06,-5.08887e-08,2.37399e-11,21050.2,35.1138], Tmin=(100,'K'), Tmax=(992.119,'K')), NASAPolynomial(coeffs=[17.0731,0.0424489,-1.57761e-05,2.87273e-09,-2.02429e-13,15660.9,-57.2679], Tmin=(992.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclohexene) + radical(CCCJ=O) + radical(Isobutyl)"""),
)

species(
    label = 'O=C1CC2C=CC1CCC2(2235)',
    structure = SMILES('O=C1CC2C=CC1CCC2'),
    E0 = (-187.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.848449,0.10813,-0.000103672,5.53707e-08,-1.22597e-11,-22369.3,1.08808], Tmin=(100,'K'), Tmax=(1074.01,'K')), NASAPolynomial(coeffs=[15.1096,0.0486969,-2.06665e-05,3.84743e-09,-2.66626e-13,-25797.2,-77.0401], Tmin=(1074.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-187.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclohexane) + ring(Cycloheptene) + ring(Cyclohexene)"""),
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
    E0 = (138.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (156.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (311.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (325.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (445.653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (235.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (244.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (308.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (297.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (435.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (319.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (289.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (254.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (396.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (232.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (233.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (557.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (202.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (368.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (192.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (393.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (419.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (384.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (231.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (268.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (216.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (227.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (163.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (336.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (353.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (384.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (146.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (196.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (729.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (307.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (318.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (436.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (320.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (396.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (318.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (248.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (209.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (188.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (265.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (147.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (156.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (333.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (146.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['CHPD(1)', 'CH2CO(43)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=C1CC2[CH]C1[CH]CCC2(2695)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.2, Ea=(17.9912,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R6;doublebond_intra_HNd_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', 'O=[C]CC1=CC=CCCC1(2696)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(4.93712,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2570 used for Cds-CsCs_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'O=C=CC1[CH]C=CCCC1(2697)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.82e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;HJ] for rate rule [Cds-CsH_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(19)', 'O=[C]CC1C=C=CCCC1(2698)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CHPD(1)', '[CH2][C]=O(360)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00517591,'m^3/(mol*s)'), n=2.41374, Ea=(2.11244,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2CO(43)', '[CH]1C=C[CH]CCC1(93)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.568606,'m^3/(mol*s)'), n=1.93802, Ea=(45.6341,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ck;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=[C]C[C]1CC=CCCC1(2699)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.94e+08,'s^-1'), n=1.27, Ea=(125.938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 162 used for R2H_S;C_rad_out_Cs2;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_Cs2;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=C[CH]C1[CH]C=CCCC1(2700)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.74568e+08,'s^-1'), n=1.384, Ea=(159.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=[C]CC1C[C]=CCCC1(2701)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=[C]CC1[CH]CCC=CC1(2702)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 219 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=[C][CH]C1CC=CCCC1(2703)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=CC[C]1[CH]C=CCCC1(2704)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;XH_out] for rate rule [R3H_SS_Cs;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=[C]CC1CC=[C]CCC1(2705)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.96e+10,'s^-1'), n=0.83, Ea=(257.734,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 201 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=CCC1[CH]CC[CH]C=C1(2648)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.164999,'s^-1'), n=3.35583, Ea=(59.1095,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;XH_out] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=[C]CC1C[CH]CC=CC1(2706)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.502,'s^-1'), n=3.86, Ea=(41.6308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 332 used for R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=[C]CC1CC=C[CH]CC1(2707)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 334 used for R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=CCC1[CH]C=CC[CH]C1(2708)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.0508,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_H/NonDeC;XH_out] for rate rule [R5H_CCC_O;C_rad_out_H/NonDeC;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=CCC1[CH][C]=CCCC1(2709)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=CCC1[CH]C=C[CH]CC1(2710)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(120.747,'s^-1'), n=2.0375, Ea=(54.2874,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R6H_SSSSS;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=CCC1[CH]C=[C]CCC1(2711)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.44313e+09,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;XH_out] for rate rule [R6HJ_2;Cd_rad_out_Cs;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]1C=C[CH]CCC1(93)', '[CH2][C]=O(360)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=[C]CC1CCCC2[CH]C21(2712)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_cs] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=C1CC2[CH][CH]C1CCC2(2713)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.42658e+11,'s^-1'), n=0.0609759, Ea=(92.8718,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_CO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]CCC=CC=CC[C]=O(2714)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 49 used for R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=CCC1=CC=CCCC1(2715)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=C=CC1CC=CCCC1(2716)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=CCC1C=C=CCCC1(2717)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O=[C]C[CH]C1C=CCCC1(2658)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O=[C]CC1[CH]CCCC=C1(2232)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['C=C([O])C1[CH]C=CCCC1(2229)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=C1CC2CCCC=CC12(2236)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CO(28)', '[CH2]C1[CH]C=CCCC1(499)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1230.59,'m^3/(mol*s)'), n=1.0523, Ea=(25.6182,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [COm;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C1[CH]C=CCCC1(499)', '[C]=O(353)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(19)', 'O=[C]CC1C=CC=CCC1(2718)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.35e+08,'cm^3/(mol*s)'), n=1.64, Ea=(1.58992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2557 used for Cds-CsH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O=[C]CC1C=CC[CH]CC1(2641)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['O=[C]CC1C=[C]CCCC1(2644)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['O=[C]CC1C=CCC[CH]C1(2637)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 219 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=[C]CC1[C]=CCCCC1(2640)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.96e+10,'s^-1'), n=0.83, Ea=(257.734,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 201 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['O=CCC1[C]=C[CH]CCC1(2719)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=[C]C[C]1C=CCCCC1(2636)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.86e+10,'s^-1'), n=0.58, Ea=(109.579,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;C_rad_out_single;Cs_H_out_Cs2] for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 2.82842712475
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['O=[C][CH]C1C=CCCCC1(2639)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(376186,'s^-1'), n=1.54152, Ea=(43.2892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_SSMS;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=C1CC2CCC[CH][CH]C12(2720)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]CC(C=CC=C)C[C]=O(2721)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.2, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 46 used for R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction45',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=C=CC1C=CCCCC1(2656)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=CCC1C=CC=CCC1(2722)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_De] for rate rule [R6radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C1C=CC(C[C]=O)CC1(2723)'],
    products = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction48',
    reactants = ['O=[C]CC1[CH]C=CCCC1(2233)'],
    products = ['O=C1CC2C=CC1CCC2(2235)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSDS;Y_rad_out;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

network(
    label = '486',
    isomers = [
        'O=[C]CC1[CH]C=CCCC1(2233)',
    ],
    reactants = [
        ('CHPD(1)', 'CH2CO(43)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '486',
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

