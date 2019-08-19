species(
    label = 'C=[C]OC1[CH]C=CCCC1(2231)',
    structure = SMILES('C=[C]OC1[CH]C=CCCC1'),
    E0 = (232.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.458419,0.0330573,0.000216988,-3.54342e-07,1.53009e-10,28156.2,35.011], Tmin=(100,'K'), Tmax=(909.767,'K')), NASAPolynomial(coeffs=[49.5939,-0.0169661,1.91029e-05,-3.887e-09,2.49715e-13,12012,-240.404], Tmin=(909.767,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO) + radical(C=CCJCO)"""),
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
    label = 'C=[C]OC1=CC=CCCC1(2661)',
    structure = SMILES('C=[C]OC1=CC=CCCC1'),
    E0 = (234.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07293,0.0979998,-7.50916e-05,2.68105e-08,-3.14509e-12,28397,24.1053], Tmin=(100,'K'), Tmax=(1174.7,'K')), NASAPolynomial(coeffs=[21.2038,0.0363026,-1.43875e-05,2.61954e-09,-1.8025e-13,22186.5,-91.1115], Tmin=(1174.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC1C=C=CCCC1(2662)',
    structure = SMILES('C=[C]OC1C=C=CCCC1'),
    E0 = (334.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.249958,0.0721395,6.03778e-06,-6.06535e-08,2.82589e-11,40393.2,31.5516], Tmin=(100,'K'), Tmax=(995.788,'K')), NASAPolynomial(coeffs=[21.2095,0.0355797,-1.36673e-05,2.60074e-09,-1.89972e-13,33658.1,-84.2453], Tmin=(995.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(C=CJO)"""),
)

species(
    label = 'C#COC1[CH]C=CCCC1(2663)',
    structure = SMILES('C#COC1[CH]C=CCCC1'),
    E0 = (201.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2175,525,750,770,3400,2100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.670284,0.0369619,0.00020769,-3.4934e-07,1.52456e-10,24506.2,30.0995], Tmin=(100,'K'), Tmax=(909.262,'K')), NASAPolynomial(coeffs=[52.0038,-0.0224314,2.13808e-05,-4.29828e-09,2.77282e-13,7803.54,-258.185], Tmin=(909.262,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtOs) + group(Ct-CtH) + ring(Cycloheptene) + radical(C=CCJCO)"""),
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
    label = 'C=[C]O[C]1CC=CCCC1(2664)',
    structure = SMILES('C=[C]O[C]1CC=CCCC1'),
    E0 = (315.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.331784,0.0375049,0.000184201,-3.09547e-07,1.34665e-10,38163.4,38.8548], Tmin=(100,'K'), Tmax=(909.561,'K')), NASAPolynomial(coeffs=[45.2815,-0.011767,1.59044e-05,-3.28155e-09,2.10973e-13,23606.3,-211.288], Tmin=(909.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO) + radical(C2CsJOC(O))"""),
)

species(
    label = 'C=[C]OC1C[C]=CCCC1(2665)',
    structure = SMILES('C=[C]OC1C[C]=CCCC1'),
    E0 = (353.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.267778,0.0336427,0.00020013,-3.29091e-07,1.42476e-10,42688.8,37.6013], Tmin=(100,'K'), Tmax=(908.327,'K')), NASAPolynomial(coeffs=[46.2578,-0.0134181,1.72178e-05,-3.55438e-09,2.297e-13,27726.1,-218.223], Tmin=(908.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = '[CH]=COC1[CH]C=CCCC1(2666)',
    structure = SMILES('[CH]=COC1[CH]C=CCCC1'),
    E0 = (239.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.74849,0.034887,0.000227056,-3.7568e-07,1.6337e-10,29055.2,33.3273], Tmin=(100,'K'), Tmax=(908.007,'K')), NASAPolynomial(coeffs=[54.0294,-0.0240182,2.30376e-05,-4.65018e-09,3.01358e-13,11588,-267.065], Tmin=(908.007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CCJCO) + radical(Cds_P)"""),
)

species(
    label = 'C=CO[C]1[CH]C=CCCC1(2667)',
    structure = SMILES('C=CO[C]1[CH]C=CCCC1'),
    E0 = (192.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.766747,0.039774,0.000203138,-3.44348e-07,1.50503e-10,23413.2,34.4869], Tmin=(100,'K'), Tmax=(909.34,'K')), NASAPolynomial(coeffs=[51.5838,-0.0199677,2.03745e-05,-4.12064e-09,2.65699e-13,6841.43,-251.868], Tmin=(909.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CCJCO) + radical(C2CsJOC(O))"""),
)

species(
    label = 'C=[C]OC1[CH]CCC=CC1(2668)',
    structure = SMILES('C=[C]OC1[CH]CCC=CC1'),
    E0 = (315.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0545929,0.0239936,0.000233005,-3.67235e-07,1.57076e-10,38122.7,38.8591], Tmin=(100,'K'), Tmax=(908.46,'K')), NASAPolynomial(coeffs=[48.0451,-0.0172711,1.95841e-05,-3.99965e-09,2.58142e-13,22346.8,-227.304], Tmin=(908.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(CCJCO) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC1CC=[C]CCC1(2669)',
    structure = SMILES('C=[C]OC1CC=[C]CCC1'),
    E0 = (353.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.267778,0.0336427,0.00020013,-3.29091e-07,1.42476e-10,42688.8,37.6013], Tmin=(100,'K'), Tmax=(908.327,'K')), NASAPolynomial(coeffs=[46.2578,-0.0134181,1.72178e-05,-3.55438e-09,2.297e-13,27726.1,-218.223], Tmin=(908.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC1C[CH]CC=CC1(2670)',
    structure = SMILES('C=[C]OC1C[CH]CC=CC1'),
    E0 = (309.805,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.039987,0.0292631,0.000201915,-3.23527e-07,1.38743e-10,37457.7,38.587], Tmin=(100,'K'), Tmax=(908.779,'K')), NASAPolynomial(coeffs=[43.05,-0.0087372,1.48922e-05,-3.12e-09,2.00806e-13,23392.2,-199.174], Tmin=(908.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO) + radical(RCCJCC)"""),
)

species(
    label = 'C=COC1[CH]C=CCC[CH]1(2671)',
    structure = SMILES('C=COC1[CH]C=CCC[CH]1'),
    E0 = (192.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.489152,0.0262577,0.000251961,-4.02063e-07,1.72927e-10,23372.5,34.4897], Tmin=(100,'K'), Tmax=(908.361,'K')), NASAPolynomial(coeffs=[54.3465,-0.02547,2.4053e-05,-4.8385e-09,3.12846e-13,5582.33,-267.879], Tmin=(908.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CCJCO) + radical(CCJCO)"""),
)

species(
    label = 'C=[C]OC1CC=C[CH]CC1(2672)',
    structure = SMILES('C=[C]OC1CC=C[CH]CC1'),
    E0 = (256.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0418923,0.0238804,0.000233662,-3.66599e-07,1.56325e-10,31051.3,35.2102], Tmin=(100,'K'), Tmax=(909.77,'K')), NASAPolynomial(coeffs=[47.5369,-0.0154828,1.85573e-05,-3.78902e-09,2.42982e-13,15366,-228.456], Tmin=(909.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO) + radical(Allyl_S)"""),
)

species(
    label = '[CH]=[C]OC1CC=CCCC1(2673)',
    structure = SMILES('[CH]=[C]OC1CC=CCCC1'),
    E0 = (362.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.313923,0.0326229,0.000208101,-3.40852e-07,1.47519e-10,43805.4,37.6966], Tmin=(100,'K'), Tmax=(908.068,'K')), NASAPolynomial(coeffs=[47.7281,-0.0158192,1.85686e-05,-3.81133e-09,2.46653e-13,28352.5,-226.49], Tmin=(908.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = 'C=COC1[CH]C=CC[CH]C1(2674)',
    structure = SMILES('C=COC1[CH]C=CC[CH]C1'),
    E0 = (186.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.394869,0.0315309,0.000220857,-3.58335e-07,1.54584e-10,22707.4,34.2187], Tmin=(100,'K'), Tmax=(908.643,'K')), NASAPolynomial(coeffs=[49.3521,-0.0169375,1.9362e-05,-3.95904e-09,2.55527e-13,6627.39,-239.753], Tmin=(908.643,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CCJCO) + radical(RCCJCC)"""),
)

species(
    label = 'C=COC1[CH][C]=CCCC1(2675)',
    structure = SMILES('C=COC1[CH][C]=CCCC1'),
    E0 = (230.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.702257,0.0359057,0.00021909,-3.63925e-07,1.58329e-10,27938.6,33.2317], Tmin=(100,'K'), Tmax=(908.233,'K')), NASAPolynomial(coeffs=[52.5589,-0.0216166,2.16866e-05,-4.39317e-09,2.844e-13,10961.6,-258.796], Tmin=(908.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2][CH]OC1C=CC=CCC1(2574)',
    structure = SMILES('[CH2][CH]OC1C=CC=CCC1'),
    E0 = (275.287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.86808,0.114385,-9.72616e-05,3.69515e-08,-4.10161e-12,33333.2,23.7665], Tmin=(100,'K'), Tmax=(1077.25,'K')), NASAPolynomial(coeffs=[25.4825,0.0332948,-1.28465e-05,2.34671e-09,-1.63674e-13,26253,-115.732], Tmin=(1077.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CCsJOCs) + radical(CJCO)"""),
)

species(
    label = 'C=COC1[CH]C=[C]CCC1(2676)',
    structure = SMILES('C=COC1[CH]C=[C]CCC1'),
    E0 = (230.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.702257,0.0359057,0.00021909,-3.63925e-07,1.58329e-10,27938.6,33.2317], Tmin=(100,'K'), Tmax=(908.233,'K')), NASAPolynomial(coeffs=[52.5589,-0.0216166,2.16866e-05,-4.39317e-09,2.844e-13,10961.6,-258.796], Tmin=(908.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Cds_S) + radical(C=CCJCO)"""),
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
    label = 'C=[C]OC1CCCC2[CH]C21(2677)',
    structure = SMILES('C=[C]OC1CCCC2[CH]C21'),
    E0 = (389.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0819532,0.0649527,2.37125e-05,-7.223e-08,3.04448e-11,47005.4,33.4618], Tmin=(100,'K'), Tmax=(1011.16,'K')), NASAPolynomial(coeffs=[18.263,0.0427307,-1.70488e-05,3.25256e-09,-2.35457e-13,40787.8,-67.0166], Tmin=(1011.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + polycyclic(s2_3_6_ane) + radical(bicyclo[4.1.0]heptane-C3-7) + radical(C=CJO)"""),
)

species(
    label = 'C=C1OC2[CH][CH]C1CCC2(2678)',
    structure = SMILES('C=C1OC2[CH][CH]C1CCC2'),
    E0 = (117.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39034,0.100314,-5.45459e-05,-5.82532e-09,1.03006e-11,14304,5.24383], Tmin=(100,'K'), Tmax=(1028.27,'K')), NASAPolynomial(coeffs=[24.5518,0.0360353,-1.42232e-05,2.68244e-09,-1.92424e-13,7032.07,-130.053], Tmin=(1028.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + polycyclic(PolycyclicRing) + radical(CCJCO) + radical(cycloheptane)"""),
)

species(
    label = '[CH2]CCC=CC=CO[C]=C(2679)',
    structure = SMILES('[CH2]CCC=CC=CO[C]=C'),
    E0 = (384.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20729,0.101562,-7.77096e-05,2.54865e-08,-1.65535e-12,46398.1,44.2454], Tmin=(100,'K'), Tmax=(1072.48,'K')), NASAPolynomial(coeffs=[21.2364,0.0369536,-1.40591e-05,2.52546e-09,-1.7372e-13,40485.7,-70.7235], Tmin=(1072.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(384.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = 'C=COC1=CC=CCCC1(2680)',
    structure = SMILES('C=COC1=CC=CCCC1'),
    E0 = (-5.25101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13342,0.0946423,-4.34477e-05,-1.68142e-08,1.45073e-11,-430.608,21.1842], Tmin=(100,'K'), Tmax=(1000.11,'K')), NASAPolynomial(coeffs=[23.6929,0.034856,-1.3034e-05,2.41182e-09,-1.72517e-13,-7372.26,-108.469], Tmin=(1000.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.25101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = 'C=COC1C=C=CCCC1(2681)',
    structure = SMILES('C=COC1C=C=CCCC1'),
    E0 = (94.6743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.447325,0.070259,3.32137e-05,-9.95829e-08,4.44205e-11,11571.6,29.13], Tmin=(100,'K'), Tmax=(969.936,'K')), NASAPolynomial(coeffs=[25.1354,0.0318094,-1.10217e-05,2.09596e-09,-1.58102e-13,3454.77,-109.77], Tmin=(969.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.6743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = 'C=[C]O[CH]C1C=CCCC1(2682)',
    structure = SMILES('C=[C]O[CH]C1C=CCCC1'),
    E0 = (291.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1685,370,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.696224,0.0784154,8.98199e-06,-7.58575e-08,3.66719e-11,35297.4,35.5424], Tmin=(100,'K'), Tmax=(965.807,'K')), NASAPolynomial(coeffs=[25.339,0.0307097,-1.03018e-05,1.90802e-09,-1.41829e-13,27464.4,-103.673], Tmin=(965.807,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexene) + radical(C=CJO) + radical(CCsJOC(O))"""),
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
    label = 'H2CC(T)(358)',
    structure = SMILES('[C]=C'),
    E0 = (600.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94832,-0.00181667,2.0713e-05,-2.36697e-08,8.4549e-12,72206.8,5.31546], Tmin=(100,'K'), Tmax=(953.658,'K')), NASAPolynomial(coeffs=[4.2027,0.00473426,-1.57307e-06,2.8591e-10,-2.08647e-14,71811.9,2.28404], Tmin=(953.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]C1C=C[CH]CCC1(1866)',
    structure = SMILES('[O]C1[CH]C=CCCC1'),
    E0 = (151.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.990608,0.00629563,0.000236283,-3.52609e-07,1.48254e-10,18423.7,25.6234], Tmin=(100,'K'), Tmax=(908.766,'K')), NASAPolynomial(coeffs=[42.035,-0.0181151,1.86721e-05,-3.77455e-09,2.42977e-13,4511.77,-203.964], Tmin=(908.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C=[C]OC1C=CC=CCC1(2684)',
    structure = SMILES('C=[C]OC1C=CC=CCC1'),
    E0 = (210.852,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17089,0.0916157,-2.95462e-05,-3.94338e-08,2.47348e-11,25565.7,22.0082], Tmin=(100,'K'), Tmax=(969.661,'K')), NASAPolynomial(coeffs=[27.3156,0.026144,-8.7673e-06,1.6269e-09,-1.21165e-13,17594.8,-127.16], Tmin=(969.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC1C=CC[CH]CC1(2685)',
    structure = SMILES('C=[C]OC1C=CC[CH]CC1'),
    E0 = (309.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0779669,0.0261862,0.000213773,-3.37788e-07,1.44161e-10,37417.9,39.3594], Tmin=(100,'K'), Tmax=(910.015,'K')), NASAPolynomial(coeffs=[44.1412,-0.0105914,1.57661e-05,-3.26216e-09,2.08724e-13,22901.4,-204.762], Tmin=(910.015,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO) + radical(RCCJCC)"""),
)

species(
    label = 'C=[C]OC1C=[C]CCCC1(2686)',
    structure = SMILES('C=[C]OC1C=[C]CCCC1'),
    E0 = (352.852,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.22952,0.0305623,0.000212001,-3.43369e-07,1.47901e-10,42649.1,38.3727], Tmin=(100,'K'), Tmax=(909.52,'K')), NASAPolynomial(coeffs=[47.348,-0.0152704,1.80907e-05,-3.69628e-09,2.37596e-13,27235.7,-223.805], Tmin=(909.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]OC1C=CCC[CH]C1(2687)',
    structure = SMILES('C=[C]OC1C=CCC[CH]C1'),
    E0 = (309.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0779669,0.0261862,0.000213773,-3.37788e-07,1.44161e-10,37417.9,39.3594], Tmin=(100,'K'), Tmax=(910.015,'K')), NASAPolynomial(coeffs=[44.1412,-0.0105914,1.57661e-05,-3.26216e-09,2.08724e-13,22901.4,-204.762], Tmin=(910.015,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(RCCJCC) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC1[C]=CCCCC1(2688)',
    structure = SMILES('C=[C]OC1[C]=CCCCC1'),
    E0 = (352.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.22952,0.0305623,0.000212001,-3.43369e-07,1.47901e-10,42649.1,38.3727], Tmin=(100,'K'), Tmax=(909.52,'K')), NASAPolynomial(coeffs=[47.348,-0.0152704,1.80907e-05,-3.69628e-09,2.37596e-13,27235.7,-223.805], Tmin=(909.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = 'C=COC1[C]=C[CH]CCC1(2689)',
    structure = SMILES('C=COC1[C]=C[CH]CCC1'),
    E0 = (254.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.247891,0.0236537,0.000247614,-3.90432e-07,1.67056e-10,30793.9,34.2036], Tmin=(100,'K'), Tmax=(909.318,'K')), NASAPolynomial(coeffs=[51.5931,-0.0219875,2.2015e-05,-4.43735e-09,2.85585e-13,13824.9,-252.437], Tmin=(909.318,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]OC1[CH]CCCC=C1(2230)',
    structure = SMILES('C=[C]OC1[CH]CCCC=C1'),
    E0 = (314.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0160396,0.0209096,0.00024489,-3.81533e-07,1.62511e-10,38082.9,39.6294], Tmin=(100,'K'), Tmax=(909.553,'K')), NASAPolynomial(coeffs=[49.1346,-0.0191221,2.04562e-05,-4.14136e-09,2.66022e-13,21856.8,-232.882], Tmin=(909.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO) + radical(CCJCO)"""),
)

species(
    label = 'C=[C]O[C]1C=CCCCC1(2690)',
    structure = SMILES('C=[C]OC1=C[CH]CCCC1'),
    E0 = (279.764,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.284414,0.0251365,0.000206117,-3.20862e-07,1.35675e-10,33834.5,37.3852], Tmin=(100,'K'), Tmax=(913.075,'K')), NASAPolynomial(coeffs=[40.5321,-0.00394043,1.19986e-05,-2.51973e-09,1.57651e-13,20346.9,-186.737], Tmin=(913.075,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO) + radical(Allyl_S)"""),
)

species(
    label = '[CH]=[C]OC1C=CCCCC1(2691)',
    structure = SMILES('[CH]=[C]OC1C=CCCCC1'),
    E0 = (362.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,3120,650,792.5,1650,1685,370,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.275487,0.0295404,0.000219979,-3.55141e-07,1.5295e-10,43765.7,38.4674], Tmin=(100,'K'), Tmax=(909.224,'K')), NASAPolynomial(coeffs=[48.8176,-0.0176704,1.94408e-05,-3.95306e-09,2.54535e-13,27862.5,-232.068], Tmin=(909.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH2]CC(C=CC=C)O[C]=C(2693)',
    structure = SMILES('[CH2]CC(C=CC=C)O[C]=C'),
    E0 = (354.718,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34317,0.0988882,-5.21524e-05,-1.33627e-08,1.4907e-11,42871.7,43.1835], Tmin=(100,'K'), Tmax=(979.237,'K')), NASAPolynomial(coeffs=[25.8229,0.0297881,-1.04381e-05,1.90039e-09,-1.36687e-13,35543.9,-97.5565], Tmin=(979.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(RCCJ)"""),
)

species(
    label = 'C=COC1C=CC=CCC1(2577)',
    structure = SMILES('C=COC1C=CC=CCC1'),
    E0 = (-28.8922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.37143,0.0897545,-2.33919e-06,-7.8573e-08,4.10728e-11,-3255.7,19.599], Tmin=(100,'K'), Tmax=(952.157,'K')), NASAPolynomial(coeffs=[31.3453,0.0222033,-6.02603e-06,1.09994e-09,-8.74819e-14,-12654.2,-153.273], Tmin=(952.157,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.8922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = '[CH2]C1C=CC(CC1)O[C]=C(2694)',
    structure = SMILES('[CH2]C1C=CC(CC1)O[C]=C'),
    E0 = (290.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,1685,370,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.263665,0.065303,4.67105e-05,-1.15751e-07,5.13389e-11,35065.9,36.1688], Tmin=(100,'K'), Tmax=(951.027,'K')), NASAPolynomial(coeffs=[25.0936,0.0295767,-8.80753e-06,1.58581e-09,-1.19984e-13,27035.4,-101.755], Tmin=(951.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexene) + radical(C=CJO) + radical(Isobutyl)"""),
)

species(
    label = 'C=C1OC2C=CC1CCC2(2234)',
    structure = SMILES('C=C1OC2C=CC1CCC2'),
    E0 = (-134.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00051,0.114111,-8.94032e-05,2.80191e-08,-1.62071e-12,-15966.1,0.69395], Tmin=(100,'K'), Tmax=(1137.44,'K')), NASAPolynomial(coeffs=[27.9985,0.0328986,-1.43287e-05,2.78714e-09,-2.00471e-13,-24361.4,-154.803], Tmin=(1137.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(oxepane) + ring(36dihydro2hpyran)"""),
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
    E0 = (232.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (250.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (457.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (562.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (428.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (232.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (437.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (553.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (345.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (374.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (442.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (489.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (351.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (332.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (674.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (406.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (308.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (444.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (336.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (473.525,'kJ/mol'),
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
    E0 = (477.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (325.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (411.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (325.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (271.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (451.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (240.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (752.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (424.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (434.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (553.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (437.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (489.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (355.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (356.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (698.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (507.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (290.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (381.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (250.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (449.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (239.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['CHPD(1)', 'CH2CO(43)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=C1OC2[CH]C1[CH]CCC2(2660)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.2, Ea=(17.9912,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R6;doublebond_intra_HNd_pri;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', 'C=[C]OC1=CC=CCCC1(2661)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'C=[C]OC1C=C=CCCC1(2662)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(19)', 'C#COC1[CH]C=CCCC1(2663)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2CO(43)', '[CH]1C=C[CH]CCC1(93)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(23.3993,'m^3/(mol*s)'), n=2.021, Ea=(33.8475,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd;CJ]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 31.4 to 33.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C]O[C]1CC=CCCC1(2664)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(20108.5,'s^-1'), n=2.606, Ea=(121.63,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_NDMustO;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]OC1C[C]=CCCC1(2665)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=COC1[CH]C=CCCC1(2666)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=CO[C]1[CH]C=CCCC1(2667)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]OC1[CH]CCC=CC1(2668)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 219 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=[C]OC1CC=[C]CCC1(2669)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.96e+10,'s^-1'), n=0.83, Ea=(257.734,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 201 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]OC1C[CH]CC=CC1(2670)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.502,'s^-1'), n=3.86, Ea=(41.6308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 332 used for R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=COC1[CH]C=CCC[CH]1(2671)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.03311e+07,'s^-1'), n=1.28625, Ea=(100.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] + [R4H_SSS_OCs;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R4H_SSS_OCs;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]OC1CC=C[CH]CC1(2672)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 334 used for R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C]OC1CC=CCCC1(2673)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=COC1[CH]C=CC[CH]C1(2674)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.0994e+09,'s^-1'), n=0.735, Ea=(76.7137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS_OCC;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_SSSS_OCC;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=COC1[CH][C]=CCCC1(2675)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_doubleC] for rate rule [R5HJ_1;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]OC1C=CC=CCC1(2574)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(46.1,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;C_rad_out_1H;XH_out] for rate rule [R6H_SSSSS;C_rad_out_H/Cd;Cd_H_out_doubleC]
Euclidian distance = 2.82842712475
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=COC1[CH]C=[C]CCC1(2676)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.24775e+08,'s^-1'), n=1.35208, Ea=(243.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_single;Cd_H_out_doubleC] + [RnH;Cd_rad_out_Cs;XH_out] for rate rule [R6HJ_2;Cd_rad_out_Cs;Cd_H_out_doubleC]
Euclidian distance = 2.82842712475
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]1C=C[CH]CCC1(93)', '[CH2][C]=O(360)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=[C]OC1CCCC2[CH]C21(2677)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_cs] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=C1OC2[CH][CH]C1CCC2(2678)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.42658e+11,'s^-1'), n=0.0609759, Ea=(92.8718,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_cddouble]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]CCC=CC=CO[C]=C(2679)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_HNd;radadd_intra_cs2H] for rate rule [R7_linear;doublebond_intra_pri_HNd_O;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=COC1=CC=CCCC1(2680)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=COC1C=C=CCCC1(2681)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C]O[CH]C1C=CCCC1(2682)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=C1OC2CCCC=CC12(2683)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H2CC(T)(358)', '[O]C1C=C[CH]CCC1(1866)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(19)', 'C=[C]OC1C=CC=CCC1(2684)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.35e+08,'cm^3/(mol*s)'), n=1.64, Ea=(1.58992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2557 used for Cds-CsH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C]OC1C=CC[CH]CC1(2685)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C]OC1C=[C]CCCC1(2686)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C]OC1C=CCC[CH]C1(2687)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 219 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=[C]OC1[C]=CCCCC1(2688)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.96e+10,'s^-1'), n=0.83, Ea=(257.734,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 201 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=COC1[C]=C[CH]CCC1(2689)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 2.82842712475
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=[C]OC1[CH]CCCC=C1(2230)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(0.502,'s^-1'), n=3.86, Ea=(41.6308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 332 used for R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=[C]O[C]1C=CCCCC1(2690)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R4H_SDS;C_rad_out_NDMustO;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=[C]OC1C=CCCCC1(2691)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R7HJ_1;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=C1OC2CCC[CH][CH]C12(2692)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(8.88848e+09,'s^-1'), n=0.76695, Ea=(57.9972,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra_cddouble] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_cddouble]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]CC(C=CC=C)O[C]=C(2693)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.2, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 46 used for R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=COC1C=CC=CCC1(2577)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad_De] for rate rule [R6radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C1C=CC(CC1)O[C]=C(2694)'],
    products = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=[C]OC1[CH]C=CCCC1(2231)'],
    products = ['C=C1OC2C=CC1CCC2(2234)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSDS;Y_rad_out;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

network(
    label = '484',
    isomers = [
        'C=[C]OC1[CH]C=CCCC1(2231)',
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
    label = '484',
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

