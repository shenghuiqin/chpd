species(
    label = 'C=C([O])C1[CH]CCCC=C1(2228)',
    structure = SMILES('C=C([O])C1[CH]CCCC=C1'),
    E0 = (166.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0172044,0.0340656,0.000179494,-2.9441e-07,1.26911e-10,20218.7,36.1979], Tmin=(100,'K'), Tmax=(909.497,'K')), NASAPolynomial(coeffs=[40.6837,-0.00472382,1.24674e-05,-2.6529e-09,1.70052e-13,7028.59,-187.983], Tmin=(909.497,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=C(C)OJ) + radical(Cs_S)"""),
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
    label = '[CH2]C1([O])C2C=CCCCC21(2724)',
    structure = SMILES('[CH2]C1([O])C2C=CCCCC21'),
    E0 = (304.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.39474,0.0510012,7.4929e-05,-1.28362e-07,5.04515e-11,36806.8,28.7723], Tmin=(100,'K'), Tmax=(995.142,'K')), NASAPolynomial(coeffs=[19.783,0.0425157,-1.69586e-05,3.32117e-09,-2.47074e-13,29509.3,-81.9479], Tmin=(995.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Cyclopropane) + ring(Cycloheptene) + ring(Cyclopropane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
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
    label = 'C=C([O])C1C=CCCCC=1(2725)',
    structure = SMILES('C=C([O])C1C=CCCCC=1'),
    E0 = (67.6347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.799552,0.0917122,-5.5157e-05,5.27389e-09,4.7513e-12,8319.32,17.8904], Tmin=(100,'K'), Tmax=(1044.37,'K')), NASAPolynomial(coeffs=[19.7779,0.0382273,-1.4716e-05,2.68021e-09,-1.86577e-13,2639.96,-88.8902], Tmin=(1044.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.6347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1C=CCCC=C1(2726)',
    structure = SMILES('C=C([O])C1C=CCCC=C1'),
    E0 = (89.8004,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.319874,0.0232127,0.000209231,-3.29676e-07,1.41155e-10,10987.1,33.9279], Tmin=(100,'K'), Tmax=(905.857,'K')), NASAPolynomial(coeffs=[42.5265,-0.0118113,1.66102e-05,-3.47532e-09,2.2667e-13,-2869.14,-199.797], Tmin=(905.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.8004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ)"""),
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
    label = 'C=C([O])[C]1C=CCCCC1(2727)',
    structure = SMILES('C=C([O])C1=C[CH]CCCC1'),
    E0 = (92.2403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.179005,0.03098,0.000208934,-3.3931e-07,1.46698e-10,11301.6,32.2177], Tmin=(100,'K'), Tmax=(905.735,'K')), NASAPolynomial(coeffs=[46.1783,-0.0137481,1.80313e-05,-3.75969e-09,2.45576e-13,-3658.74,-223.069], Tmin=(905.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.2403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=C(C)OJ) + radical(Allyl_S)"""),
)

species(
    label = 'C=C([O])C1C=CCC[CH]C1(2728)',
    structure = SMILES('C=C([O])C1C=CCC[CH]C1'),
    E0 = (166.415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102615,0.0354759,0.000168053,-2.7812e-07,1.20342e-10,20202.5,35.8998], Tmin=(100,'K'), Tmax=(906.879,'K')), NASAPolynomial(coeffs=[38.3556,-0.0014285,1.1061e-05,-2.43175e-09,1.57952e-13,7843.76,-174.796], Tmin=(906.879,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(RCCJCC) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1[C]=CCCCC1(2729)',
    structure = SMILES('C=C([O])C1[C]=CCCCC1'),
    E0 = (209.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.204422,0.0398462,0.000166303,-2.83736e-07,1.24099e-10,25433.7,34.9115], Tmin=(100,'K'), Tmax=(906.426,'K')), NASAPolynomial(coeffs=[41.5618,-0.00610653,1.3385e-05,-2.86573e-09,1.86812e-13,12178.2,-193.836], Tmin=(906.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1C=CC[CH]CC1(2730)',
    structure = SMILES('C=C([O])C1C=CC[CH]CC1'),
    E0 = (166.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102615,0.0354759,0.000168053,-2.7812e-07,1.20342e-10,20202.5,35.8998], Tmin=(100,'K'), Tmax=(906.879,'K')), NASAPolynomial(coeffs=[38.3556,-0.0014285,1.1061e-05,-2.43175e-09,1.57952e-13,7843.76,-174.796], Tmin=(906.879,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(RCCJCC) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)[C]1[CH]CCCC=C1(2731)',
    structure = SMILES('C=C(O)C1[CH]CCC[CH]C=1'),
    E0 = (95.5482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.343664,0.0296797,0.000225641,-3.65331e-07,1.57968e-10,11710,30.2015], Tmin=(100,'K'), Tmax=(906.122,'K')), NASAPolynomial(coeffs=[49.9859,-0.0196376,2.11302e-05,-4.33266e-09,2.82601e-13,-4507.32,-246.806], Tmin=(906.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.5482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Allyl_S) + radical(Allyl_S)"""),
)

species(
    label = '[CH]=C(O)C1[CH]CCCC=C1(2732)',
    structure = SMILES('[CH]=C(O)C1[CH]CCCC=C1'),
    E0 = (275.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,3120,650,792.5,1650,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.418557,0.041496,0.000170686,-2.9475e-07,1.29408e-10,33381.2,36.6649], Tmin=(100,'K'), Tmax=(907.679,'K')), NASAPolynomial(coeffs=[44.6807,-0.0109465,1.55757e-05,-3.24771e-09,2.10709e-13,19167.2,-209.744], Tmin=(907.679,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Cs_S) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C1[C]=CCCC[CH]1(2733)',
    structure = SMILES('C=C(O)C1[C]=CCCC[CH]1'),
    E0 = (266.536,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.372784,0.0425205,0.000162698,-2.82964e-07,1.24353e-10,32264.6,36.5708], Tmin=(100,'K'), Tmax=(907.969,'K')), NASAPolynomial(coeffs=[43.2114,-0.00854698,1.42259e-05,-2.99101e-09,1.93777e-13,18540.4,-201.482], Tmin=(907.969,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Cds_S) + radical(Cs_S)"""),
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
    label = 'C=C([O])C1C=[C]CCCC1(2734)',
    structure = SMILES('C=C([O])C1C=[C]CCCC1'),
    E0 = (209.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.204422,0.0398462,0.000166303,-2.83736e-07,1.24099e-10,25433.7,34.9115], Tmin=(100,'K'), Tmax=(906.426,'K')), NASAPolynomial(coeffs=[41.5618,-0.00610653,1.3385e-05,-2.86573e-09,1.86812e-13,12178.2,-193.836], Tmin=(906.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([O])C1C=CCCCC1(2735)',
    structure = SMILES('[CH]=C([O])C1C=CCCCC1'),
    E0 = (219.053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,350,440,435,1725,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.250583,0.0388266,0.000174273,-2.95497e-07,1.29143e-10,26550.3,35.007], Tmin=(100,'K'), Tmax=(906.204,'K')), NASAPolynomial(coeffs=[43.0324,-0.00850816,1.4736e-05,-3.12273e-09,2.03769e-13,12804.6,-202.104], Tmin=(906.204,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)C1[CH][CH]CCC=C1(2736)',
    structure = SMILES('C=C(O)C1[CH][CH]CCC=C1'),
    E0 = (223.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0651157,0.0381421,0.000164478,-2.77393e-07,1.20617e-10,27033.4,37.5569], Tmin=(100,'K'), Tmax=(908.464,'K')), NASAPolynomial(coeffs=[40.0038,-0.00386647,1.19005e-05,-2.55668e-09,1.64887e-13,14206.4,-182.434], Tmin=(908.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Cs_S) + radical(RCCJCC)"""),
)

species(
    label = 'C=C(O)C1[CH]CCC[C]=C1(2737)',
    structure = SMILES('C=C(O)C1[CH]CCC[C]=C1'),
    E0 = (266.536,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.372784,0.0425205,0.000162698,-2.82964e-07,1.24353e-10,32264.6,36.5708], Tmin=(100,'K'), Tmax=(907.969,'K')), NASAPolynomial(coeffs=[43.2114,-0.00854698,1.42259e-05,-2.99101e-09,1.93777e-13,18540.4,-201.482], Tmin=(907.969,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Cds_S) + radical(Cs_S)"""),
)

species(
    label = 'C=C(O)C1[CH]C[CH]CC=C1(2738)',
    structure = SMILES('C=C(O)C1[CH]C[CH]CC=C1'),
    E0 = (223.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0651157,0.0381421,0.000164478,-2.77393e-07,1.20617e-10,27033.4,37.5569], Tmin=(100,'K'), Tmax=(908.464,'K')), NASAPolynomial(coeffs=[40.0038,-0.00386647,1.19005e-05,-2.55668e-09,1.64887e-13,14206.4,-182.434], Tmin=(908.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(RCCJCC) + radical(Cs_S)"""),
)

species(
    label = 'C=C(O)C1[CH]CC[CH]C=C1(2739)',
    structure = SMILES('C=C(O)C1[CH]C=CCC[CH]1'),
    E0 = (170.355,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.093447,0.0211862,0.00023688,-3.68041e-07,1.56464e-10,20690.2,35.7031], Tmin=(100,'K'), Tmax=(909.766,'K')), NASAPolynomial(coeffs=[46.8617,-0.0154009,1.84941e-05,-3.77477e-09,2.41982e-13,5185.06,-223.95], Tmin=(909.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Allyl_S) + radical(Cs_S)"""),
)

species(
    label = '[CH]1CCCC=CC1[C]1CO1(2740)',
    structure = SMILES('[CH]1CCCC=CC1[C]1CO1'),
    E0 = (313.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.429612,0.0230462,0.000213267,-3.34426e-07,1.44219e-10,37852.2,35.589], Tmin=(100,'K'), Tmax=(893.937,'K')), NASAPolynomial(coeffs=[39.9863,-0.0068595,1.66278e-05,-3.71042e-09,2.52955e-13,24902.6,-183.687], Tmin=(893.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + ring(Cycloheptene) + radical(C2CsJO) + radical(Cs_S)"""),
)

species(
    label = 'C=C([O])C1[CH]C2CCCC21(2741)',
    structure = SMILES('C=C([O])C1[CH]C2CCCC21'),
    E0 = (213.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.542728,0.0537709,4.90929e-05,-9.47359e-08,3.76959e-11,25775.7,30.7627], Tmin=(100,'K'), Tmax=(1000.37,'K')), NASAPolynomial(coeffs=[16.3828,0.044245,-1.73096e-05,3.28642e-09,-2.38069e-13,19914,-59.1203], Tmin=(1000.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + polycyclic(s2_4_5_ane) + radical(C=C(C)OJ) + radical(bicyclo[3.2.0]heptane-C5-6)"""),
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
    label = 'C=C1OC2[CH]C1[CH]CCC2(2660)',
    structure = SMILES('C=C1OC2[CH]C1[CH]CCC2'),
    E0 = (117.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39034,0.100314,-5.45459e-05,-5.82532e-09,1.03006e-11,14304,5.24383], Tmin=(100,'K'), Tmax=(1028.27,'K')), NASAPolynomial(coeffs=[24.5518,0.0360353,-1.42232e-05,2.68244e-09,-1.92424e-13,7032.07,-130.053], Tmin=(1028.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + polycyclic(PolycyclicRing) + radical(cycloheptane) + radical(CCJCO)"""),
)

species(
    label = 'C=C([O])C1C2[CH]CCCC21(2742)',
    structure = SMILES('C=C([O])C1C2[CH]CCCC21'),
    E0 = (190.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.383553,0.0596148,3.03961e-05,-7.42053e-08,3.01968e-11,23070.3,32.3027], Tmin=(100,'K'), Tmax=(1013.32,'K')), NASAPolynomial(coeffs=[15.9264,0.0452725,-1.7964e-05,3.39467e-09,-2.43633e-13,17506.7,-54.7977], Tmin=(1013.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + polycyclic(s2_3_6_ane) + radical(bicyclo[4.1.0]heptane-C6-2) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1OC2[CH]CCC[CH]C12(2743)',
    structure = SMILES('C=C1OC2[CH]CCC[CH]C12'),
    E0 = (117.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39034,0.100314,-5.45459e-05,-5.82532e-09,1.03006e-11,14304,5.24383], Tmin=(100,'K'), Tmax=(1028.27,'K')), NASAPolynomial(coeffs=[24.5518,0.0360353,-1.42232e-05,2.68244e-09,-1.92424e-13,7032.07,-130.053], Tmin=(1028.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + polycyclic(PolycyclicRing) + radical(cycloheptane) + radical(CCJCO)"""),
)

species(
    label = '[CH]=CCCCC=CC(=C)[O](2744)',
    structure = SMILES('[CH]=CCCCC=CC(=C)[O]'),
    E0 = (262.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,350,440,435,1725,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.41573,0.108155,-9.70168e-05,4.53325e-08,-8.4533e-12,31799.6,40.5109], Tmin=(100,'K'), Tmax=(1294.36,'K')), NASAPolynomial(coeffs=[21.8584,0.0362303,-1.36643e-05,2.4011e-09,-1.61247e-13,25774.6,-77.7783], Tmin=(1294.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC=CC(C=C)C(=C)[O](2745)',
    structure = SMILES('[CH2]CC=CC(C=C)C(=C)[O]'),
    E0 = (234.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.717512,0.0999942,-8.53296e-05,3.93628e-08,-7.4751e-12,28410,40.6406], Tmin=(100,'K'), Tmax=(1244.08,'K')), NASAPolynomial(coeffs=[16.6517,0.0441487,-1.79967e-05,3.2813e-09,-2.24517e-13,24088.3,-46.9492], Tmin=(1244.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)C1C=CCCCC=1(2746)',
    structure = SMILES('C=C(O)C1C=CCCCC=1'),
    E0 = (-70.1701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11826,0.0951546,-4.59345e-05,-1.37769e-08,1.34932e-11,-8239.82,17.5316], Tmin=(100,'K'), Tmax=(994.606,'K')), NASAPolynomial(coeffs=[23.0929,0.0355691,-1.30557e-05,2.38057e-09,-1.687e-13,-14924.8,-108.538], Tmin=(994.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-70.1701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = 'C=C(O)C1C=CCCC=C1(2747)',
    structure = SMILES('C=C(O)C1C=CCCC=C1'),
    E0 = (-48.0044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0434239,0.027142,0.00021697,-3.47193e-07,1.49446e-10,-5570.07,33.7312], Tmin=(100,'K'), Tmax=(906.125,'K')), NASAPolynomial(coeffs=[46.2604,-0.0151627,1.86626e-05,-3.86625e-09,2.52039e-13,-20616.2,-221.814], Tmin=(906.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.0044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene)"""),
)

species(
    label = 'C=C([O])[CH]C1C=CCCC1(2748)',
    structure = SMILES('[CH2]C([O])=CC1C=CCCC1'),
    E0 = (96.5702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.336718,0.077403,-9.22347e-06,-4.29466e-08,2.19182e-11,11786.8,32.4677], Tmin=(100,'K'), Tmax=(990.176,'K')), NASAPolynomial(coeffs=[19.1218,0.0399111,-1.4711e-05,2.68227e-09,-1.89774e-13,5917.78,-71.3939], Tmin=(990.176,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.5702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C1CCC=CC1C(=C)[O](2749)',
    structure = SMILES('[CH2]C1CCC=CC1C(=C)[O]'),
    E0 = (145.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.263772,0.0751695,-1.55739e-06,-5.25167e-08,2.59694e-11,17722.7,32.7631], Tmin=(100,'K'), Tmax=(971.645,'K')), NASAPolynomial(coeffs=[19.2171,0.0390023,-1.3696e-05,2.44924e-09,-1.72642e-13,11858.5,-71.3558], Tmin=(971.645,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclohexene) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1OC2CCCC=CC12(2683)',
    structure = SMILES('C=C1OC2CCCC=CC12'),
    E0 = (-140.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58416,0.105372,-6.66771e-05,5.55932e-09,6.55252e-12,-16646.7,-0.736304], Tmin=(100,'K'), Tmax=(1036.84,'K')), NASAPolynomial(coeffs=[25.0289,0.0358844,-1.41538e-05,2.654e-09,-1.89273e-13,-23949,-138.693], Tmin=(1036.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-140.206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(Cycloheptene) + ring(Oxetane)"""),
)

species(
    label = 'O(20)',
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
    label = 'C=[C]C1[CH]CCCC=C1(2750)',
    structure = SMILES('C=[C]C1[CH]CCCC=C1'),
    E0 = (481.095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.629505,0.023801,0.000185082,-2.87236e-07,1.21189e-10,58030.6,32.9145], Tmin=(100,'K'), Tmax=(913.184,'K')), NASAPolynomial(coeffs=[36.0083,-0.000236236,9.49649e-06,-2.03973e-09,1.26949e-13,46109.9,-164.446], Tmin=(913.184,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(481.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Cds_S) + radical(Cs_S)"""),
)

species(
    label = 'CC(=O)[C]1[CH]CCCC=C1(2751)',
    structure = SMILES('CC([O])=C1[CH]CCCC=C1'),
    E0 = (83.5133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.064899,0.0269012,0.000220927,-3.52424e-07,1.51554e-10,10249.3,32.335], Tmin=(100,'K'), Tmax=(905.965,'K')), NASAPolynomial(coeffs=[46.5479,-0.0149463,1.87515e-05,-3.89246e-09,2.53891e-13,-4925.17,-225.076], Tmin=(905.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.5133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(Cycloheptane) + radical(C=C(C)OJ) + radical(Allyl_S)"""),
)

species(
    label = 'CC(=O)C1[C]=CCCC[CH]1(2752)',
    structure = SMILES('CC(=O)C1[C]=CCCC[CH]1'),
    E0 = (246.407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202631,0.0328923,0.000172757,-2.79639e-07,1.19501e-10,29819.9,37.1174], Tmin=(100,'K'), Tmax=(912.897,'K')), NASAPolynomial(coeffs=[37.7257,0.000400959,9.38101e-06,-2.03228e-09,1.26904e-13,17471.9,-170.598], Tmin=(912.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCJCC=O) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)C1[CH][CH]CCC=C1(2753)',
    structure = SMILES('CC(=O)C1[CH][CH]CCC=C1'),
    E0 = (203.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.509697,0.0285211,0.000174512,-2.74036e-07,1.15751e-10,24588.8,38.1056], Tmin=(100,'K'), Tmax=(913.716,'K')), NASAPolynomial(coeffs=[34.5212,0.00507614,7.05865e-06,-1.59869e-09,9.80759e-14,13136.7,-151.567], Tmin=(913.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCJCC=O) + radical(RCCJCC)"""),
)

species(
    label = 'CC(=O)C1[CH]CCC[C]=C1(2754)',
    structure = SMILES('CC(=O)C1[CH]CCC[C]=C1'),
    E0 = (246.407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202631,0.0328923,0.000172757,-2.79639e-07,1.19501e-10,29819.9,37.1174], Tmin=(100,'K'), Tmax=(912.897,'K')), NASAPolynomial(coeffs=[37.7257,0.000400959,9.38101e-06,-2.03228e-09,1.26904e-13,17471.9,-170.598], Tmin=(912.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cds_S) + radical(CCJCC=O)"""),
)

species(
    label = 'CC(=O)C1[CH]C[CH]CC=C1(2755)',
    structure = SMILES('CC(=O)C1[CH]C[CH]CC=C1'),
    E0 = (203.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.509697,0.0285211,0.000174512,-2.74036e-07,1.15751e-10,24588.8,38.1056], Tmin=(100,'K'), Tmax=(913.716,'K')), NASAPolynomial(coeffs=[34.5212,0.00507614,7.05865e-06,-1.59869e-09,9.80759e-14,13136.7,-151.567], Tmin=(913.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCJCC=O) + radical(RCCJCC)"""),
)

species(
    label = 'CC(=O)C1[CH]CC[CH]C=C1(2756)',
    structure = SMILES('CC(=O)C1[CH]CC[CH]C=C1'),
    E0 = (149.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.428365,0.0231311,0.000206289,-3.17155e-07,1.33357e-10,18182.4,34.7268], Tmin=(100,'K'), Tmax=(914.235,'K')), NASAPolynomial(coeffs=[39.0083,-0.00166967,1.07238e-05,-2.26772e-09,1.40252e-13,5110.38,-180.851], Tmin=(914.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCJCC=O) + radical(Allyl_S)"""),
)

species(
    label = 'C=C1OC2CCC[CH][CH]C12(2692)',
    structure = SMILES('C=C1OC2CCC[CH][CH]C12'),
    E0 = (86.7115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.63492,0.1083,-8.73918e-05,3.52596e-08,-5.65566e-12,10644.8,5.99459], Tmin=(100,'K'), Tmax=(1489.8,'K')), NASAPolynomial(coeffs=[25.3537,0.0358382,-1.4434e-05,2.61206e-09,-1.77182e-13,2603.25,-134.969], Tmin=(1489.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.7115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + polycyclic(PolycyclicRing) + radical(cycloheptane) + radical(cycloheptane)"""),
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
    label = 'O=C1CC2[CH]CCC[CH]C12(2652)',
    structure = SMILES('O=C1CC2[CH]CCC[CH]C12'),
    E0 = (65.5504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.250525,0.10058,-8.41183e-05,3.99125e-08,-8.30648e-12,8030.97,2.98843], Tmin=(100,'K'), Tmax=(1089.99,'K')), NASAPolynomial(coeffs=[10.9411,0.0595098,-2.75993e-05,5.34412e-09,-3.77915e-13,5591.21,-51.9691], Tmin=(1089.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.5504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + polycyclic(PolycyclicRing) + radical(CCJCC=O) + radical(cycloheptane)"""),
)

species(
    label = 'CC(=O)C1C=CCCCC=1(2757)',
    structure = SMILES('CC(=O)C1C=CCCCC=1'),
    E0 = (-86.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.985283,0.0942121,-5.49751e-05,5.31141e-09,3.99898e-12,-10181.9,18.9939], Tmin=(100,'K'), Tmax=(1107.52,'K')), NASAPolynomial(coeffs=[21.0077,0.0403726,-1.67169e-05,3.14586e-09,-2.21739e-13,-16623,-96.4414], Tmin=(1107.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = 'CC(=O)C1C=CCCC=C1(2758)',
    structure = SMILES('CC(=O)C1C=CCCC=C1'),
    E0 = (-80.4018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.407484,0.0257988,0.000194378,-3.02143e-07,1.27391e-10,-9491.12,33.5994], Tmin=(100,'K'), Tmax=(914.394,'K')), NASAPolynomial(coeffs=[37.8922,0.000114834,9.65261e-06,-2.06594e-09,1.27311e-13,-22127.7,-175.501], Tmin=(914.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.4018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene)"""),
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
    label = '[CH2]C1=C[CH]CCCC=CO1(2759)',
    structure = SMILES('[CH2]C1=C[CH]CCCC=CO1'),
    E0 = (209.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.12147,0.0976277,-6.04345e-05,7.33902e-09,4.40308e-12,25351,22.1402], Tmin=(100,'K'), Tmax=(1056.22,'K')), NASAPolynomial(coeffs=[21.3796,0.0391382,-1.53228e-05,2.82055e-09,-1.97516e-13,19107.2,-94.7016], Tmin=(1056.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclononane) + radical(Allyl_S) + radical(C=C(O)CJ)"""),
)

species(
    label = 'CH2(T)(33)',
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
    label = 'O=[C]C1[CH]CCCC=C1(2760)',
    structure = SMILES('O=[C]C1[CH]CCCC=C1'),
    E0 = (222.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00257,0.011019,0.000214568,-3.22673e-07,1.35754e-10,26874.7,33.3384], Tmin=(100,'K'), Tmax=(909.01,'K')), NASAPolynomial(coeffs=[38.6467,-0.0114022,1.52188e-05,-3.13452e-09,2.01308e-13,14113.5,-177.231], Tmin=(909.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(CC(C)CJ=O) + radical(CCJCC=O)"""),
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
    E0 = (166.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (306.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (293.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (307.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (249.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (241.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (294.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (325.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (344.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (318.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (323.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (467.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (310.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (208.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (254.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (263.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (250.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (299.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (283.946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (231.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (419.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (397.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (297.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (291.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (445.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (192.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (215.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (331.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (261.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (229.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (191.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (326.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (412.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (174.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (724.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (283.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (290.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (229.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (279.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (263.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (209.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (291.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (445.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (233.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (229.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (191.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (438.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (174.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (522.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (603.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['CHPD(1)', 'CH2CO(43)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['[CH2]C1([O])C2C=CCCCC21(2724)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(139.863,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', 'C=C([O])C1C=CCCCC=1(2725)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.45e+08,'cm^3/(mol*s)'), n=1.64, Ea=(13.7654,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2607 used for Cds-CdCd_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdCd_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'C=C([O])C1C=CCCC=C1(2726)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.92e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CHPD(1)', '[CH2][C]=O(360)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0147295,'m^3/(mol*s)'), n=2.37207, Ea=(16.7338,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2CO(43)', '[CH]1C=C[CH]CCC1(93)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7051.36,'m^3/(mol*s)'), n=1.07323, Ea=(43.3386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C([O])[C]1C=CCCCC1(2727)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([O])C1C=CCC[CH]C1(2728)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C([O])C1[C]=CCCCC1(2729)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.37e+07,'s^-1'), n=1.41, Ea=(177.82,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 207 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cd_H_out_doubleC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C([O])C1C=CC[CH]CC1(2730)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.1e+09,'s^-1'), n=1.01, Ea=(152.298,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 369 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C(O)[C]1[CH]CCCC=C1(2731)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.19923e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(O)C1[CH]CCCC=C1(2732)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C(O)C1[C]=CCCC[CH]1(2733)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C([O])C1[CH]C=CCCC1(2229)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.502,'s^-1'), n=3.86, Ea=(41.6308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 332 used for R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C([O])C1C=[C]CCCC1(2734)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C([O])C1C=CCCCC1(2735)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C(O)C1[CH][CH]CCC=C1(2736)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R5HJ_1;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C(O)C1[CH]CCC[C]=C1(2737)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C(O)C1[CH]C[CH]CC=C1(2738)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(46.1,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R6HJ_2;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C(O)C1[CH]CC[CH]C=C1(2739)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(46.1,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;C_rad_out_H/NonDeC;XH_out] for rate rule [R6H_SMSSR;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]1C=C[CH]CCC1(93)', '[CH2][C]=O(360)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['[CH]1CCCC=CC1[C]1CO1(2740)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C([O])C1[CH]C2CCCC21(2741)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra_pri;radadd_intra_csHCs] for rate rule [Rn0c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['O=C1CC2CCC[CH][CH]C12(2720)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.82767e+08,'s^-1'), n=1.02667, Ea=(125.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C1OC2[CH]C1[CH]CCC2(2660)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.11e+10,'s^-1'), n=0.18, Ea=(278.654,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C([O])C1C2[CH]CCCC21(2742)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(26.3459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C1OC2[CH]CCC[CH]C12(2743)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=CCCCC=CC(=C)[O](2744)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.58366e+14,'s^-1'), n=-0.170131, Ea=(68.698,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;doublebond_intra_pri;radadd_intra_cdsingleH] + [R7;doublebond_intra_pri_HCd;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_HCd;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]CC=CC(C=C)C(=C)[O](2745)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.2, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 46 used for R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C(O)C1C=CCCCC=1(2746)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C(O)C1C=CCCC=C1(2747)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C([O])[CH]C1C=CCCC1(2748)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['[CH2]C1CCC=CC1C(=C)[O](2749)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C1OC2CCCC=CC12(2683)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(20)', 'C=[C]C1[CH]CCCC=C1(2750)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['CC(=O)[C]1[CH]CCCC=C1(2751)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.30996,'s^-1'), n=3.5644, Ea=(117.059,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CC(=O)C1[C]=CCCC[CH]1(2752)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CC(=O)C1[CH][CH]CCC=C1(2753)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(0.113548,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R5HJ_1;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CC(=O)C1[CH]CCC[C]=C1(2754)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CC(=O)C1[CH]C[CH]CC=C1(2755)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(138.3,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R6HJ_2;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['CC(=O)C1[CH]CC[CH]C=C1(2756)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R6H_RSSMS;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['C=C1OC2CCC[CH][CH]C12(2692)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.7342e+08,'s^-1'), n=0.920995, Ea=(125.348,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCs] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHCs]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['O=C1CC2[CH]C1[CH]CCC2(2695)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(8.11e+10,'s^-1'), n=0.18, Ea=(278.654,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H] for rate rule [Rn2c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['O=C1CC2[CH]CCC[CH]C12(2652)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2.54811e+07,'s^-1'), n=0.785425, Ea=(66.8554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['CC(=O)C1C=CCCCC=1(2757)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['CC(=O)C1C=CCCC=C1(2758)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O=[C]CC1[CH]CCCC=C1(2232)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    products = ['O=C1CC2CCCC=CC12(2236)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C1=C[CH]CCCC=CO1(2759)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction50',
    reactants = ['CH2(T)(33)', 'O=[C]C1[CH]CCCC=C1(2760)'],
    products = ['C=C([O])C1[CH]CCCC=C1(2228)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

network(
    label = '481',
    isomers = [
        'C=C([O])C1[CH]CCCC=C1(2228)',
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
    label = '481',
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

