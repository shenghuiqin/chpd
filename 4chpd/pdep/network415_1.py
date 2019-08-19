species(
    label = '[O]CC1[CH]CCCC=C1(2162)',
    structure = SMILES('[O]CC1[CH]CCCC=C1'),
    E0 = (205.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924269,0.0170756,0.000199121,-2.96867e-07,1.23168e-10,24895,34.6139], Tmin=(100,'K'), Tmax=(916.502,'K')), NASAPolynomial(coeffs=[34.1671,0.0030038,7.72676e-06,-1.67283e-09,1.0008e-13,13299.2,-152.883], Tmin=(916.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cs_S) + radical(CCOJ)"""),
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
    label = 'CH2O(35)',
    structure = SMILES('C=O'),
    E0 = (-119.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: FFCM1(-)"""),
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
    label = '[O]CC1C=CCCCC=1(2400)',
    structure = SMILES('[O]CC1C=CCCCC=1'),
    E0 = (107.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.172,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0748471,0.0842688,-5.96182e-05,2.16124e-08,-3.22534e-12,13119.5,16.7265], Tmin=(100,'K'), Tmax=(1534.22,'K')), NASAPolynomial(coeffs=[16.6218,0.0407378,-1.70584e-05,3.11894e-09,-2.11865e-13,7996.17,-70.9717], Tmin=(1534.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CCOJ)"""),
)

species(
    label = '[O]CC1C=CCCC=C1(2401)',
    structure = SMILES('[O]CC1C=CCCC=C1'),
    E0 = (122.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.172,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.908551,0.0186074,0.000189184,-2.87416e-07,1.20487e-10,14909.5,30.5843], Tmin=(100,'K'), Tmax=(912.847,'K')), NASAPolynomial(coeffs=[34.5398,-0.000364139,9.3744e-06,-2.01372e-09,1.25489e-13,3419.9,-157.902], Tmin=(912.847,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(CCOJ)"""),
)

species(
    label = 'O=CC1[CH]CCCC=C1(2402)',
    structure = SMILES('O=CC1[CH]CCCC=C1'),
    E0 = (63.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.172,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.944209,0.0123675,0.000216369,-3.23451e-07,1.35387e-10,7790.87,32.3908], Tmin=(100,'K'), Tmax=(911.024,'K')), NASAPolynomial(coeffs=[37.8727,-0.00700133,1.31863e-05,-2.74415e-09,1.73962e-13,-4862.47,-174.844], Tmin=(911.024,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2][O](354)',
    structure = SMILES('[CH2][O]'),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88409,-0.00363885,3.28543e-05,-4.13611e-08,1.59631e-11,23210.8,7.47983], Tmin=(100,'K'), Tmax=(933.06,'K')), NASAPolynomial(coeffs=[6.69335,0.000289989,8.61416e-07,-1.56351e-10,7.33778e-15,21991.3,-9.6043], Tmin=(933.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(H3COJ) + radical(CsJOH)"""),
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
    label = '[O]C[C]1C=CCCCC1(2378)',
    structure = SMILES('[O]CC1=C[CH]CCCC1'),
    E0 = (151.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.934286,0.0154598,0.000207472,-3.07541e-07,1.27506e-10,18396.4,31.4175], Tmin=(100,'K'), Tmax=(915.404,'K')), NASAPolynomial(coeffs=[34.7885,0.00223368,8.41447e-06,-1.81913e-09,1.10213e-13,6554.43,-159.745], Tmin=(915.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCOJ) + radical(Allyl_S)"""),
)

species(
    label = '[O]CC1C=CCC[CH]C1(2375)',
    structure = SMILES('[O]CC1C=CCC[CH]C1'),
    E0 = (205.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00896,0.0184968,0.000187627,-2.80485e-07,1.16546e-10,24878.9,34.3182], Tmin=(100,'K'), Tmax=(914.16,'K')), NASAPolynomial(coeffs=[31.8334,0.00630836,6.31511e-06,-1.45047e-09,8.7881e-14,14116.8,-139.664], Tmin=(914.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(RCCJCC) + radical(CCOJ)"""),
)

species(
    label = 'O[CH]C1[CH]CCCC=C1(2403)',
    structure = SMILES('O[CH]C1[CH]CCCC=C1'),
    E0 = (160.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,3615,1277.5,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.35188,0.0278464,0.000182803,-2.93171e-07,1.25654e-10,19456.3,35.3581], Tmin=(100,'K'), Tmax=(909.068,'K')), NASAPolynomial(coeffs=[39.2282,-0.00618141,1.28404e-05,-2.71107e-09,1.74139e-13,6725.87,-179.634], Tmin=(909.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cs_S) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]C1C=CCCCC1(2379)',
    structure = SMILES('[O][CH]C1C=CCCCC1'),
    E0 = (191.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.498585,0.0266058,0.000180935,-2.85128e-07,1.21029e-10,23197.1,33.1358], Tmin=(100,'K'), Tmax=(911.857,'K')), NASAPolynomial(coeffs=[36.6843,-0.000783885,9.9299e-06,-2.14028e-09,1.34649e-13,11137.3,-168.042], Tmin=(911.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O]CC1[C]=CCCCC1(2376)',
    structure = SMILES('[O]CC1[C]=CCCCC1'),
    E0 = (248.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.701595,0.0228718,0.000185857,-2.86067e-07,1.20285e-10,30110,33.3311], Tmin=(100,'K'), Tmax=(913.314,'K')), NASAPolynomial(coeffs=[35.0386,0.00163194,8.6382e-06,-1.88424e-09,1.16724e-13,18451.7,-158.699], Tmin=(913.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]CC1C=CC[CH]CC1(2373)',
    structure = SMILES('[O]CC1C=CC[CH]CC1'),
    E0 = (205.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00896,0.0184969,0.000187627,-2.80485e-07,1.16546e-10,24878.9,34.3182], Tmin=(100,'K'), Tmax=(914.16,'K')), NASAPolynomial(coeffs=[31.8334,0.00630833,6.31513e-06,-1.45047e-09,8.78813e-14,14116.8,-139.664], Tmin=(914.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCOJ) + radical(RCCJCC)"""),
)

species(
    label = 'OC[C]1[CH]CCCC=C1(2404)',
    structure = SMILES('OCC1[CH]CCC[CH]C=1'),
    E0 = (67.0434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.790387,0.0127388,0.000229593,-3.42275e-07,1.43092e-10,8233.09,29.9673], Tmin=(100,'K'), Tmax=(911.068,'K')), NASAPolynomial(coeffs=[39.4859,-0.00660546,1.35787e-05,-2.83665e-09,1.79686e-13,-5065.79,-187.402], Tmin=(911.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.0434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(Allyl_S)"""),
)

species(
    label = 'OCC1[C]=CCCC[CH]1(2405)',
    structure = SMILES('OCC1[C]=CCCC[CH]1'),
    E0 = (217.813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.55482,0.0241136,0.000187717,-2.94095e-07,1.24901e-10,26369.2,35.5536], Tmin=(100,'K'), Tmax=(910.413,'K')), NASAPolynomial(coeffs=[37.5809,-0.00376304,1.15472e-05,-2.45469e-09,1.56186e-13,14040.9,-170.282], Tmin=(910.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cds_S) + radical(Cs_S)"""),
)

species(
    label = '[O]CC1[CH]C=CCCC1(2163)',
    structure = SMILES('[O]CC1[CH]C=CCCC1'),
    E0 = (151.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4464.45,'J/mol'), sigma=(7.65423,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=697.34 K, Pc=22.59 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930807,0.013352,0.000218015,-3.21517e-07,1.33212e-10,18343.4,30.8934], Tmin=(100,'K'), Tmax=(915.053,'K')), NASAPolynomial(coeffs=[36.0839,7.2699e-05,9.65441e-06,-2.0528e-09,1.25379e-13,6032.6,-167.695], Tmin=(915.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCOJ) + radical(Allyl_S)"""),
)

species(
    label = '[O]CC1C=[C]CCCC1(2374)',
    structure = SMILES('[O]CC1C=[C]CCCC1'),
    E0 = (248.976,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.701859,0.0228685,0.00018587,-2.86084e-07,1.20293e-10,30110,33.3301], Tmin=(100,'K'), Tmax=(913.302,'K')), NASAPolynomial(coeffs=[35.0378,0.00163325,8.63743e-06,-1.88405e-09,1.16709e-13,18452,-158.695], Tmin=(913.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'OCC1[CH][CH]CCC=C1(2406)',
    structure = SMILES('OCC1[CH][CH]CCC=C1'),
    E0 = (174.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.861992,0.0197413,0.000189475,-2.88495e-07,1.21152e-10,21138.1,36.5414], Tmin=(100,'K'), Tmax=(911.079,'K')), NASAPolynomial(coeffs=[34.3754,0.000913992,9.2238e-06,-2.02085e-09,1.27337e-13,9706.13,-151.245], Tmin=(911.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(RCCJCC) + radical(Cs_S)"""),
)

species(
    label = 'OCC1[CH]CCC[C]=C1(2407)',
    structure = SMILES('OCC1[CH]CCC[C]=C1'),
    E0 = (217.813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.55482,0.0241136,0.000187717,-2.94095e-07,1.24901e-10,26369.2,35.5536], Tmin=(100,'K'), Tmax=(910.413,'K')), NASAPolynomial(coeffs=[37.5809,-0.00376304,1.15472e-05,-2.45469e-09,1.56186e-13,14040.9,-170.282], Tmin=(910.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cds_S) + radical(Cs_S)"""),
)

species(
    label = 'OCC1[CH]C[CH]CC=C1(2408)',
    structure = SMILES('OCC1[CH]C[CH]CC=C1'),
    E0 = (174.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.861994,0.0197413,0.000189476,-2.88495e-07,1.21152e-10,21138.1,36.5414], Tmin=(100,'K'), Tmax=(911.079,'K')), NASAPolynomial(coeffs=[34.3753,0.000914005,9.22379e-06,-2.02084e-09,1.27337e-13,9706.13,-151.245], Tmin=(911.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cs_S) + radical(RCCJCC)"""),
)

species(
    label = 'OCC1[CH]C=CCC[CH]1(2359)',
    structure = SMILES('OCC1[CH]C=CCC[CH]1'),
    E0 = (120.014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.783386,0.0146025,0.000219839,-3.29489e-07,1.37799e-10,14602.6,33.1181], Tmin=(100,'K'), Tmax=(912.323,'K')), NASAPolynomial(coeffs=[38.6258,-0.00532161,1.25631e-05,-2.62319e-09,1.64836e-13,1622.02,-179.275], Tmin=(912.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cs_S) + radical(Allyl_S)"""),
)

species(
    label = '[O]CC1[CH]C2CCCC12(2409)',
    structure = SMILES('[O]CC1[CH]C2CCCC12'),
    E0 = (250.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19083,0.0485721,2.82543e-05,-5.2021e-08,1.78701e-11,30261,27.7046], Tmin=(100,'K'), Tmax=(1134.41,'K')), NASAPolynomial(coeffs=[9.01015,0.054446,-2.37365e-05,4.52227e-09,-3.18421e-13,26334.9,-20.4904], Tmin=(1134.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + polycyclic(s2_4_5_ane) + radical(bicyclo[3.2.0]heptane-C5-6) + radical(CCOJ)"""),
)

species(
    label = '[CH]1CCCC2[CH]C1CO2(2317)',
    structure = SMILES('[CH]1CCCC2[CH]C1CO2'),
    E0 = (68.9059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.369015,0.0921238,-6.65303e-05,2.43736e-08,-3.67341e-12,8447.37,1.86455], Tmin=(100,'K'), Tmax=(1516.64,'K')), NASAPolynomial(coeffs=[17.7673,0.0442915,-1.92234e-05,3.57919e-09,-2.45739e-13,2946.05,-93.1866], Tmin=(1516.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.9059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + polycyclic(PolycyclicRing) + radical(cycloheptane) + radical(CCJCO)"""),
)

species(
    label = '[O]CC1C2[CH]CCCC21(2410)',
    structure = SMILES('[O]CC1C2[CH]CCCC21'),
    E0 = (228.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0385,0.0542372,1.06434e-05,-3.35544e-08,1.15264e-11,27555.4,29.2266], Tmin=(100,'K'), Tmax=(1204.89,'K')), NASAPolynomial(coeffs=[9.36063,0.0542245,-2.37192e-05,4.47997e-09,-3.11994e-13,23545.4,-20.7923], Tmin=(1204.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + polycyclic(s2_3_6_ane) + radical(bicyclo[4.1.0]heptane-C6-2) + radical(CCOJ)"""),
)

species(
    label = '[CH]1CCC[CH]C2OCC12(2394)',
    structure = SMILES('[CH]1CCC[CH]C2OCC12'),
    E0 = (68.9059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.369015,0.0921238,-6.65303e-05,2.43736e-08,-3.67341e-12,8447.37,1.86455], Tmin=(100,'K'), Tmax=(1516.64,'K')), NASAPolynomial(coeffs=[17.7673,0.0442915,-1.92234e-05,3.57919e-09,-2.45739e-13,2946.05,-93.1866], Tmin=(1516.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.9059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + polycyclic(PolycyclicRing) + radical(cycloheptane) + radical(CCJCO)"""),
)

species(
    label = '[CH]=CCCCC=CC[O](2411)',
    structure = SMILES('[CH]=CCCCC=CC[O]'),
    E0 = (323.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3120,650,792.5,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.627837,0.0815251,-5.89752e-05,2.48929e-08,-4.88319e-12,39030.4,36.2272], Tmin=(100,'K'), Tmax=(1094.19,'K')), NASAPolynomial(coeffs=[7.11427,0.0578128,-2.64685e-05,5.08728e-09,-3.57997e-13,37610.9,4.34999], Tmin=(1094.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC=CC(C=C)C[O](2412)',
    structure = SMILES('[CH2]CC=CC(C=C)C[O]'),
    E0 = (267.633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.155236,0.091699,-9.10643e-05,6.15009e-08,-1.90684e-11,32320.8,36.3024], Tmin=(100,'K'), Tmax=(750.749,'K')), NASAPolynomial(coeffs=[6.09556,0.0600527,-2.78426e-05,5.36679e-09,-3.77923e-13,31428.8,9.34602], Tmin=(750.749,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(RCCJ)"""),
)

species(
    label = 'OCC1C=CCCCC=1(2413)',
    structure = SMILES('OCC1C=CCCCC=1'),
    E0 = (-117.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.678438,0.089912,-5.79166e-05,1.13059e-08,1.99213e-12,-13997.7,18.0316], Tmin=(100,'K'), Tmax=(1079.95,'K')), NASAPolynomial(coeffs=[19.0597,0.0378257,-1.47689e-05,2.69458e-09,-1.86903e-13,-19486.8,-84.387], Tmin=(1079.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = 'O=CC1C=CCCCC1(2382)',
    structure = SMILES('O=CC1C=CCCCC1'),
    E0 = (-136.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.757104,0.017498,0.000208007,-3.14243e-07,1.31626e-10,-16246.1,30.5656], Tmin=(100,'K'), Tmax=(912.604,'K')), NASAPolynomial(coeffs=[37.2956,-0.00268358,1.1118e-05,-2.35266e-09,1.47244e-13,-28743.7,-174.304], Tmin=(912.604,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene)"""),
)

species(
    label = 'OCC1C=CCCC=C1(2414)',
    structure = SMILES('OCC1C=CCCC=C1'),
    E0 = (-103.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.565887,0.0211164,0.000202341,-3.13666e-07,1.33108e-10,-12219.3,30.9542], Tmin=(100,'K'), Tmax=(908.516,'K')), NASAPolynomial(coeffs=[39.1685,-0.0066738,1.3497e-05,-2.85035e-09,1.83399e-13,-25100.9,-183.868], Tmin=(908.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene)"""),
)

species(
    label = '[O]C[CH]C1C=CCCC1(2371)',
    structure = SMILES('[O]C[CH]C1C=CCCC1'),
    E0 = (188.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.747263,0.0544019,2.7584e-05,-6.35442e-08,2.48545e-11,22860.7,32.5005], Tmin=(100,'K'), Tmax=(1039.13,'K')), NASAPolynomial(coeffs=[13.2293,0.0465424,-1.90837e-05,3.61542e-09,-2.57656e-13,18096.9,-38.6374], Tmin=(1039.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1CCC=CC1C[O](2415)',
    structure = SMILES('[CH2]C1CCC=CC1C[O]'),
    E0 = (185.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.620793,0.0583683,1.78354e-05,-5.54307e-08,2.28244e-11,22400,31.2641], Tmin=(100,'K'), Tmax=(1020.23,'K')), NASAPolynomial(coeffs=[13.1709,0.0459517,-1.79967e-05,3.32692e-09,-2.34214e-13,17924.6,-38.9173], Tmin=(1020.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CCOJ) + radical(Isobutyl)"""),
)

species(
    label = 'C1=CC2COC2CCC1(2168)',
    structure = SMILES('C1=CC2COC2CCC1'),
    E0 = (-186.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.753896,0.0902954,-5.35548e-05,8.00983e-09,2.32605e-12,-22291.9,-0.998609], Tmin=(100,'K'), Tmax=(1143.88,'K')), NASAPolynomial(coeffs=[19.5814,0.0407752,-1.69293e-05,3.16442e-09,-2.21203e-13,-28356.6,-108.013], Tmin=(1143.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-186.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(Cycloheptene) + ring(Oxetane)"""),
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
    label = '[CH2]C1[CH]CCCC=C1(449)',
    structure = SMILES('[CH2]C1[CH]CCCC=C1'),
    E0 = (344.881,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32194,0.00152687,0.000245382,-3.54745e-07,1.47509e-10,41629.8,31.4066], Tmin=(100,'K'), Tmax=(906.009,'K')), NASAPolynomial(coeffs=[37.6731,-0.00929551,1.55092e-05,-3.26764e-09,2.1226e-13,28900.2,-174.278], Tmin=(906.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cs_S) + radical(Isobutyl)"""),
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
    E0 = (205.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (330.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (340.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (292.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (282.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (205.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (367.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (364.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (289.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (349.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (383.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (357.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (280.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (318.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (247.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (293.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (256.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (250.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (285.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (285.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (452.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (336.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (484.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (231.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (255.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (350.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (294.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (269.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (269.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (230.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (365.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (451.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (213.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (587.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['CHPD(1)', 'CH2O(35)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(19)', '[O]CC1C=CCCCC=1(2400)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.11e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.1294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2594 used for Cds-CdCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', '[O]CC1C=CCCC=C1(2401)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.92e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'O=CC1[CH]CCCC=C1(2402)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.6e+09,'cm^3/(mol*s)'), n=0.935, Ea=(17.4473,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2782 used for CO-CsH_O;HJ
Exact match found for rate rule [CO-CsH_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CHPD(1)', '[CH2][O](354)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0147295,'m^3/(mol*s)'), n=2.37207, Ea=(16.7338,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2O(35)', '[CH]1C=C[CH]CCC1(93)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(450.719,'m^3/(mol*s)'), n=0.996465, Ea=(65.4968,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-HH_O;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 61.9 to 65.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['[O]C[C]1C=CCCCC1(2378)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.41e+08,'s^-1'), n=1.52, Ea=(161.502,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['[O]CC1C=CCC[CH]C1(2375)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O[CH]C1[CH]CCCC=C1(2403)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4500,'s^-1'), n=2.62, Ea=(129.286,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 322 used for R2H_S;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['[O][CH]C1C=CCCCC1(2379)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.17661e+06,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['[O]CC1[C]=CCCCC1(2376)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.37e+07,'s^-1'), n=1.41, Ea=(177.82,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 207 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cd_H_out_doubleC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['[O]CC1C=CC[CH]CC1(2373)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.1e+09,'s^-1'), n=1.01, Ea=(152.298,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 369 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['OC[C]1[CH]CCCC=C1(2404)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['OCC1[C]=CCCC[CH]1(2405)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['[O]CC1[CH]C=CCCC1(2163)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.502,'s^-1'), n=3.86, Ea=(41.6308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 332 used for R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]CC1C=[C]CCCC1(2374)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['OCC1[CH][CH]CCC=C1(2406)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.19599e+09,'s^-1'), n=0.63, Ea=(50.8774,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['OCC1[CH]CCC[C]=C1(2407)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['OCC1[CH]C[CH]CC=C1(2408)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.12e+09,'s^-1'), n=0, Ea=(79.7052,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6HJ_3;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['OCC1[CH]C=CCC[CH]1(2359)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.12e+09,'s^-1'), n=0, Ea=(79.7052,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6H;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6H_RSSMS;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]1C=C[CH]CCC1(93)', '[CH2][O](354)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['[O]CC1[CH]C2CCCC12(2409)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra_pri;radadd_intra_csHCs] for rate rule [Rn0c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['[CH]1CCCC2[CH]C1CO2(2317)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.11e+10,'s^-1'), n=0.18, Ea=(278.654,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['[O]CC1C2[CH]CCCC21(2410)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(26.1785,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['[CH]1CCC[CH]C2OCC12(2394)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=CCCCC=CC[O](2411)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]CC=CC(C=C)C[O](2412)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.2, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 46 used for R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['OCC1C=CCCCC=1(2413)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['O=CC1C=CCCCC1(2382)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['OCC1C=CCCC=C1(2414)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['[O]C[CH]C1C=CCCC1(2371)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['[CH2]C1CCC=CC1C[O](2415)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]CC1[CH]CCCC=C1(2162)'],
    products = ['C1=CC2COC2CCC1(2168)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O(20)', '[CH2]C1[CH]CCCC=C1(449)'],
    products = ['[O]CC1[CH]CCCC=C1(2162)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/Cs;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '415',
    isomers = [
        '[O]CC1[CH]CCCC=C1(2162)',
    ],
    reactants = [
        ('CHPD(1)', 'CH2O(35)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '415',
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

