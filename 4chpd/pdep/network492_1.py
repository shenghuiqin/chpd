species(
    label = 'C=C([O])C1[CH]C=CC=CC1(2239)',
    structure = SMILES('C=C([O])C1C=C[CH]C=CC1'),
    E0 = (187.021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.554994,0.0115722,0.000248281,-3.78125e-07,1.60819e-10,22677.8,32.2663], Tmin=(100,'K'), Tmax=(903.379,'K')), NASAPolynomial(coeffs=[45.8044,-0.0203488,2.1608e-05,-4.45426e-09,2.92929e-13,7629.38,-219.478], Tmin=(903.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'C7H8(11)',
    structure = SMILES('C1=CC=CCC=C1'),
    E0 = (165.503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3872.77,'J/mol'), sigma=(6.33517,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=604.92 K, Pc=34.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.976321,0.0470646,2.43918e-05,-7.34711e-08,3.37246e-11,20032,17.1406], Tmin=(100,'K'), Tmax=(948.557,'K')), NASAPolynomial(coeffs=[19.6256,0.016921,-4.63468e-06,8.3158e-10,-6.48659e-14,14312.2,-83.3476], Tmin=(948.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
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
    label = '[CH2]C1([O])C2C=CCC1C=C2(3413)',
    structure = SMILES('[CH2]C1([O])C2C=CCC1C=C2'),
    E0 = (386.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0797305,0.06859,-2.58244e-06,-4.06537e-08,1.84468e-11,46594.2,25.6425], Tmin=(100,'K'), Tmax=(1053.74,'K')), NASAPolynomial(coeffs=[17.7673,0.0400945,-1.7033e-05,3.29443e-09,-2.37623e-13,40721,-70.7965], Tmin=(1053.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
)

species(
    label = 'C=C1OC2[CH]CC1[CH]C=C2(3405)',
    structure = SMILES('C=C1OC2[CH]C=CC1C[CH]2'),
    E0 = (183.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.73526,0.0742818,8.3824e-05,-2.07551e-07,9.67207e-11,22304.1,3.20233], Tmin=(100,'K'), Tmax=(930.135,'K')), NASAPolynomial(coeffs=[49.9576,-0.0147177,1.23772e-05,-2.26177e-09,1.29997e-13,6921.47,-273.439], Tmin=(930.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(Cycloheptane) + ring(Oxane) + radical(CCJCO) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C1OC2[CH]C=CCC1[CH]2(3372)',
    structure = SMILES('C=C1OC2[CH]C=CCC1[CH]2'),
    E0 = (177.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30001,0.0917984,-2.28613e-05,-5.00662e-08,2.90001e-11,21520.6,0.32265], Tmin=(100,'K'), Tmax=(972.97,'K')), NASAPolynomial(coeffs=[29.4205,0.0239001,-8.2142e-06,1.58503e-09,-1.21713e-13,12778.4,-161.249], Tmin=(972.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Tetrahydrofuran) + ring(Cycloheptene) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(C=CCJCO)"""),
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
    label = 'C=C([O])C1C=CC=C=CC1(3414)',
    structure = SMILES('C=C([O])C1C=CC=C=CC1'),
    E0 = (290.698,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.100981,0.0748226,-2.43679e-05,-2.3611e-08,1.4789e-11,35124.1,25.7921], Tmin=(100,'K'), Tmax=(999.535,'K')), NASAPolynomial(coeffs=[18.8196,0.032964,-1.23624e-05,2.2719e-09,-1.61228e-13,29650.4,-73.9406], Tmin=(999.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1C=C=CC=CC1(3415)',
    structure = SMILES('C=C([O])C1C=C=CC=CC1'),
    E0 = (290.698,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.100981,0.0748226,-2.43679e-05,-2.3611e-08,1.4789e-11,35124.1,25.7921], Tmin=(100,'K'), Tmax=(999.535,'K')), NASAPolynomial(coeffs=[18.8196,0.032964,-1.23624e-05,2.2719e-09,-1.61228e-13,29650.4,-73.9406], Tmin=(999.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]1C=C[CH]CC=C1(125)',
    structure = SMILES('[CH]1C=C[CH]CC=C1'),
    E0 = (333.186,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.70798,-0.040239,0.000339088,-4.54205e-07,1.84685e-10,40184,19.5412], Tmin=(100,'K'), Tmax=(904.636,'K')), NASAPolynomial(coeffs=[40.627,-0.029508,2.54911e-05,-5.11056e-09,3.33198e-13,26023.7,-199.942], Tmin=(904.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C([O])C1C=CC[C]=CC1(3416)',
    structure = SMILES('C=C([O])C1C=CC[C]=CC1'),
    E0 = (323.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.213655,0.0284513,0.000186167,-3.04367e-07,1.32028e-10,39052.3,34.2939], Tmin=(100,'K'), Tmax=(904.74,'K')), NASAPolynomial(coeffs=[42.3395,-0.0136262,1.69085e-05,-3.52213e-09,2.31004e-13,25529.3,-197.33], Tmin=(904.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1C=[C]CC=CC1(3417)',
    structure = SMILES('C=C([O])C1C=[C]CC=CC1'),
    E0 = (323.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.213655,0.0284513,0.000186167,-3.04367e-07,1.32028e-10,39052.3,34.2939], Tmin=(100,'K'), Tmax=(904.74,'K')), NASAPolynomial(coeffs=[42.3395,-0.0136262,1.69085e-05,-3.52213e-09,2.31004e-13,25529.3,-197.33], Tmin=(904.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C(O)[C]1C=C[CH]C=CC1(3418)',
    structure = SMILES('C=C(O)C1=C[CH]C=C[CH]C1'),
    E0 = (169.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.189035,0.011177,0.000274061,-4.22189e-07,1.80836e-10,20591.7,29.9504], Tmin=(100,'K'), Tmax=(902.101,'K')), NASAPolynomial(coeffs=[52.9572,-0.0318283,2.80214e-05,-5.68835e-09,3.75586e-13,3300.64,-262.258], Tmin=(902.101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH]=C(O)C1C=C[CH]C=CC1(3419)',
    structure = SMILES('[CH]=C(O)C1C=C[CH]C=CC1'),
    E0 = (296.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,350,440,435,1725,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.11791,0.0190181,0.000239421,-3.78403e-07,1.63294e-10,35840.4,32.738], Tmin=(100,'K'), Tmax=(902.132,'K')), NASAPolynomial(coeffs=[49.8096,-0.0265853,2.47243e-05,-5.05098e-09,3.33745e-13,19764.7,-241.284], Tmin=(902.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Cds_P) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C([O])C1C=CCC=[C]C1(3420)',
    structure = SMILES('C=C([O])C1C=CCC=[C]C1'),
    E0 = (323.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.213655,0.0284513,0.000186167,-3.04367e-07,1.32028e-10,39052.3,34.2939], Tmin=(100,'K'), Tmax=(904.74,'K')), NASAPolynomial(coeffs=[42.3395,-0.0136262,1.69085e-05,-3.52213e-09,2.31004e-13,25529.3,-197.33], Tmin=(904.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1[C]=CCC=CC1(3421)',
    structure = SMILES('C=C([O])C1[C]=CCC=CC1'),
    E0 = (323.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.213655,0.0284513,0.000186167,-3.04367e-07,1.32028e-10,39052.3,34.2939], Tmin=(100,'K'), Tmax=(904.74,'K')), NASAPolynomial(coeffs=[42.3395,-0.0136262,1.69085e-05,-3.52213e-09,2.31004e-13,25529.3,-197.33], Tmin=(904.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)C1[CH]C=C[CH]C=C1(3422)',
    structure = SMILES('C=C(O)C1[CH]C=C[CH]C=C1'),
    E0 = (190.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,350,440,435,1725,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.390068,0.0102755,0.000264974,-4.04124e-07,1.72079e-10,23086.2,30.251], Tmin=(100,'K'), Tmax=(903.921,'K')), NASAPolynomial(coeffs=[49.6123,-0.0262385,2.47071e-05,-5.02728e-09,3.29959e-13,6780.75,-243.216], Tmin=(903.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C(O)C1[C]=C[CH]C=CC1(3423)',
    structure = SMILES('C=C(O)C1[C]=C[CH]C=CC1'),
    E0 = (287.058,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,350,440,435,1725,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.163789,0.0200413,0.000231437,-3.66621e-07,1.5824e-10,34723.7,32.6436], Tmin=(100,'K'), Tmax=(902.193,'K')), NASAPolynomial(coeffs=[48.3395,-0.0241845,2.33738e-05,-4.7941e-09,3.16798e-13,19138.1,-233.018], Tmin=(902.193,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.058,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(Cds_S)"""),
)

species(
    label = 'C=C([O])[C]1C=CCC=CC1(3424)',
    structure = SMILES('C=C([O])C1=C[CH]CC=CC1'),
    E0 = (205.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.239643,0.0195778,0.000228826,-3.59982e-07,1.54646e-10,24920.2,31.5981], Tmin=(100,'K'), Tmax=(904.276,'K')), NASAPolynomial(coeffs=[46.9548,-0.0212659,2.15537e-05,-4.4158e-09,2.89743e-13,9692.76,-226.557], Tmin=(904.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ) + radical(Allyl_S)"""),
)

species(
    label = 'C=C([O])C1[CH]C=CCC=C1(2240)',
    structure = SMILES('C=C([O])C1[CH]C=CCC=C1'),
    E0 = (226.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.439383,0.0186917,0.000219684,-3.41847e-07,1.4586e-10,27414.8,31.9033], Tmin=(100,'K'), Tmax=(906.58,'K')), NASAPolynomial(coeffs=[43.6162,-0.015687,1.82457e-05,-3.75623e-09,2.44242e-13,13170.3,-207.55], Tmin=(906.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ) + radical(Allyl_S)"""),
)

species(
    label = 'C=C(O)C1C=C[CH]C=[C]C1(3425)',
    structure = SMILES('C=C(O)C1C=C[CH]C=[C]C1'),
    E0 = (287.058,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,350,440,435,1725,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.163789,0.0200413,0.000231437,-3.66621e-07,1.5824e-10,34723.7,32.6436], Tmin=(100,'K'), Tmax=(902.193,'K')), NASAPolynomial(coeffs=[48.3395,-0.0241845,2.33738e-05,-4.7941e-09,3.16798e-13,19138.1,-233.018], Tmin=(902.193,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.058,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C(O)C1C=[C][CH]C=CC1(3426)',
    structure = SMILES('C=C(O)C1C=[C]C=C[CH]C1'),
    E0 = (270.102,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,350,440,435,1725,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.42554,0.106596,-8.92882e-05,2.93136e-08,-5.76517e-13,32692.4,16.5395], Tmin=(100,'K'), Tmax=(989.533,'K')), NASAPolynomial(coeffs=[23.9361,0.0292671,-1.02545e-05,1.79401e-09,-1.23582e-13,26439.8,-111.781], Tmin=(989.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(Allyl_S)"""),
)

species(
    label = '[CH]=C([O])C1C=CCC=CC1(3427)',
    structure = SMILES('[CH]=C([O])C1C=CCC=CC1'),
    E0 = (332.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.168193,0.0274229,0.000194171,-3.16176e-07,1.37095e-10,40168.9,34.3869], Tmin=(100,'K'), Tmax=(904.561,'K')), NASAPolynomial(coeffs=[43.8083,-0.0160248,1.82577e-05,-3.77869e-09,2.47925e-13,26156.4,-205.589], Tmin=(904.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C1C=C[CH][C]=CC1(3428)',
    structure = SMILES('C=C(O)C1[CH]C=C[C]=CC1'),
    E0 = (270.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,350,440,435,1725,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13877,0.094479,-4.67675e-05,-2.05853e-08,1.86158e-11,32753.5,17.8964], Tmin=(100,'K'), Tmax=(951.785,'K')), NASAPolynomial(coeffs=[26.0491,0.0249152,-7.57723e-06,1.30415e-09,-9.3566e-14,25553.6,-122.561], Tmin=(951.785,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(Allyl_S)"""),
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
    label = '[CH]1C=CCC(C=C1)[C]1CO1(3429)',
    structure = SMILES('[CH]1C=CCC(C=C1)[C]1CO1'),
    E0 = (333.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.963014,0.000599357,0.000281925,-4.18046e-07,1.7813e-10,40311.5,31.6734], Tmin=(100,'K'), Tmax=(890.812,'K')), NASAPolynomial(coeffs=[45.1511,-0.0225596,2.58119e-05,-5.52207e-09,3.76687e-13,25485.1,-215.43], Tmin=(890.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + radical(C2CsJO) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C([O])C1C=CC2[CH]C2C1(3430)',
    structure = SMILES('C=C([O])C1C=CC2[CH]C2C1'),
    E0 = (326.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.261106,0.067177,-8.26803e-06,-3.18681e-08,1.53123e-11,39454.8,29.2944], Tmin=(100,'K'), Tmax=(1042.26,'K')), NASAPolynomial(coeffs=[15.7119,0.0402841,-1.62004e-05,3.03586e-09,-2.15029e-13,34474,-54.3296], Tmin=(1042.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_3_6_ene_1) + radical(cyclopropane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1CC=CC2[CH]C21(3431)',
    structure = SMILES('C=C([O])C1CC=CC2[CH]C21'),
    E0 = (327.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.544996,0.0551453,3.36878e-05,-8.05981e-08,3.37919e-11,39516,30.6585], Tmin=(100,'K'), Tmax=(988.885,'K')), NASAPolynomial(coeffs=[17.5799,0.0363266,-1.37416e-05,2.59609e-09,-1.89072e-13,33697.9,-63.7174], Tmin=(988.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_3_6_ene_1) + radical(cyclopropane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=C1CC2[CH][CH]C1CC=C2(3105)',
    structure = SMILES('O=C1CC2[CH][CH]C1CC=C2'),
    E0 = (204.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0865876,0.0907585,-6.96845e-05,2.77505e-08,-4.62098e-12,24708,3.17735], Tmin=(100,'K'), Tmax=(1367.94,'K')), NASAPolynomial(coeffs=[15.1181,0.0462981,-2.09317e-05,3.99063e-09,-2.78698e-13,20548.2,-74.9404], Tmin=(1367.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclohexane) + ring(Cycloheptene) + ring(Cyclohexane) + radical(Cs_S) + radical(CCJCC=O)"""),
)

species(
    label = 'C=C1OC2C=C[CH][CH]C1C2(3330)',
    structure = SMILES('C=C1OC2[CH]C=C[CH]C1C2'),
    E0 = (118.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29152,0.0917232,-2.22714e-05,-4.94657e-08,2.83301e-11,14449.4,-3.31051], Tmin=(100,'K'), Tmax=(978.507,'K')), NASAPolynomial(coeffs=[28.9857,0.0255659,-9.17131e-06,1.77938e-09,-1.35532e-13,5766.02,-162.816], Tmin=(978.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Tetrahydrofuran) + ring(Cycloheptene) + ring(Tetrahydrofuran) + radical(C=CCJCO) + radical(Allyl_S)"""),
)

species(
    label = 'C=C1OC2[CH][CH]C=CCC12(3432)',
    structure = SMILES('C=C1OC2[CH]C=C[CH]CC12'),
    E0 = (118.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29152,0.0917232,-2.22714e-05,-4.94657e-08,2.83301e-11,14449.4,-3.31051], Tmin=(100,'K'), Tmax=(978.507,'K')), NASAPolynomial(coeffs=[28.9857,0.0255659,-9.17131e-06,1.77938e-09,-1.35532e-13,5766.02,-162.816], Tmin=(978.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(Cycloheptene) + ring(Oxetane) + radical(C=CCJCO) + radical(Allyl_S)"""),
)

species(
    label = 'C=C(O)C1C=C=CC=CC1(3433)',
    structure = SMILES('C=C(O)C1C=C=CC=CC1'),
    E0 = (152.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.445196,0.0785446,-1.60076e-05,-4.17404e-08,2.32363e-11,18566.1,25.526], Tmin=(100,'K'), Tmax=(974.348,'K')), NASAPolynomial(coeffs=[22.3706,0.0299175,-1.04836e-05,1.92155e-09,-1.39199e-13,11982.1,-94.9251], Tmin=(974.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = 'C=C(O)C1C=CC=C=CC1(3434)',
    structure = SMILES('C=C(O)C1C=CC=C=CC1'),
    E0 = (152.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.445196,0.0785446,-1.60076e-05,-4.17404e-08,2.32363e-11,18566.1,25.526], Tmin=(100,'K'), Tmax=(974.348,'K')), NASAPolynomial(coeffs=[22.3706,0.0299175,-1.04836e-05,1.92155e-09,-1.39199e-13,11982.1,-94.9251], Tmin=(974.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = 'C=C1OC2C=CCC1C=C2(2247)',
    structure = SMILES('C=C1OC2C=CCC1C=C2'),
    E0 = (-19.2214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62822,0.074511,7.02964e-05,-1.84794e-07,8.57177e-11,-2063.94,3.43806], Tmin=(100,'K'), Tmax=(943.723,'K')), NASAPolynomial(coeffs=[48.0718,-0.0114926,8.86712e-06,-1.43773e-09,6.79552e-14,-16995.4,-262.866], Tmin=(943.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.2214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(Cycloheptane) + ring(36dihydro2hpyran)"""),
)

species(
    label = 'C=CC1[CH]C1C=CC(=C)[O](3435)',
    structure = SMILES('C=CC1[CH]C1C=CC(=C)[O]'),
    E0 = (388.169,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2950,3150,900,1000,1100,350,440,435,1725,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.311583,0.0759831,-1.30626e-05,-4.45423e-08,2.4552e-11,46858.2,36.1637], Tmin=(100,'K'), Tmax=(962.318,'K')), NASAPolynomial(coeffs=[21.9705,0.0285074,-9.426e-06,1.68556e-09,-1.21505e-13,40479.5,-81.3388], Tmin=(962.318,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(cyclopropane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC1[CH]C=CC=C([O])C1(3436)',
    structure = SMILES('C=CC1C=C[CH]C=C([O])C1'),
    E0 = (187.021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.554994,0.0115722,0.000248281,-3.78125e-07,1.60819e-10,22677.8,32.2663], Tmin=(100,'K'), Tmax=(903.379,'K')), NASAPolynomial(coeffs=[45.8044,-0.0203488,2.1608e-05,-4.45426e-09,2.92929e-13,7629.38,-219.478], Tmin=(903.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ) + radical(C=CCJC=C)"""),
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
    label = 'C=[C]C1C=C[CH]C=CC1(3437)',
    structure = SMILES('C=[C]C1C=C[CH]C=CC1'),
    E0 = (501.617,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (118.176,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16776,0.00130438,0.000253866,-3.70922e-07,1.55069e-10,60489.7,28.9811], Tmin=(100,'K'), Tmax=(906.005,'K')), NASAPolynomial(coeffs=[41.1178,-0.0158422,1.86261e-05,-3.83851e-09,2.49612e-13,46715.5,-195.878], Tmin=(906.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'O=C1CC2[CH]CC1[CH]C=C2(3062)',
    structure = SMILES('O=C1CC2[CH]C=CC1C[CH]2'),
    E0 = (143.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.179081,0.0908091,-7.11651e-05,3.02088e-08,-5.58391e-12,17358,-2.07354], Tmin=(100,'K'), Tmax=(1206.55,'K')), NASAPolynomial(coeffs=[11.2706,0.0540381,-2.54511e-05,4.94998e-09,-3.50247e-13,14681.5,-57.6663], Tmin=(1206.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.228,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclohexane) + ring(Cycloheptene) + ring(Cyclohexane) + radical(Cs_S) + radical(Allyl_S)"""),
)

species(
    label = 'O=C1CC2[CH]C1C[CH]C=C2(3088)',
    structure = SMILES('O=C1CC2[CH]C=CCC1[CH]2'),
    E0 = (149.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.580062,0.0927089,-6.96612e-05,2.63438e-08,-4.03895e-12,18178.6,1.26619], Tmin=(100,'K'), Tmax=(1519.71,'K')), NASAPolynomial(coeffs=[19.8458,0.0389469,-1.65971e-05,3.06589e-09,-2.09664e-13,11970.2,-105.826], Tmin=(1519.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclopentane) + ring(Cycloheptene) + ring(Cyclopentane) + radical(Allyl_S) + radical(CCJCC=O)"""),
)

species(
    label = 'CC(=O)[C]1C=C[CH]C=CC1(3438)',
    structure = SMILES('CC(=O)[C]1C=C[CH]C=CC1'),
    E0 = (143.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23357,-0.00665721,0.000293957,-4.17078e-07,1.7184e-10,17413.1,31.1523], Tmin=(100,'K'), Tmax=(910.388,'K')), NASAPolynomial(coeffs=[42.9996,-0.0148097,1.84639e-05,-3.76111e-09,2.39353e-13,2541.66,-206.334], Tmin=(910.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJ(C)C=O) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC(=O)C1[CH]C=C[CH]C=C1(3439)',
    structure = SMILES('CC(=O)C1C=C[CH][CH]C=C1'),
    E0 = (201.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.806355,0.0153196,0.000212392,-3.1926e-07,1.33402e-10,24440,29.9546], Tmin=(100,'K'), Tmax=(914.384,'K')), NASAPolynomial(coeffs=[38.0359,-0.00495435,1.17416e-05,-2.42905e-09,1.50543e-13,11670.7,-178.919], Tmin=(914.384,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(Allyl_S)"""),
)

species(
    label = 'CC(=O)C1[C]=C[CH]C=CC1(3440)',
    structure = SMILES('CC(=O)C1[C]=C[CH]C=CC1'),
    E0 = (261.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.750315,0.0142518,0.000221904,-3.35919e-07,1.41621e-10,31630.1,33.209], Tmin=(100,'K'), Tmax=(908.299,'K')), NASAPolynomial(coeffs=[40.1666,-0.00996386,1.52237e-05,-3.17221e-09,2.04303e-13,18308.3,-187.078], Tmin=(908.299,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC(=O)C1C=C[CH]C=[C]C1(3441)',
    structure = SMILES('CC(=O)C1C=C[CH]C=[C]C1'),
    E0 = (261.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.750344,0.0142515,0.000221905,-3.35921e-07,1.41622e-10,31630.1,33.2089], Tmin=(100,'K'), Tmax=(908.298,'K')), NASAPolynomial(coeffs=[40.1666,-0.00996372,1.52236e-05,-3.17219e-09,2.04301e-13,18308.3,-187.078], Tmin=(908.298,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)C1C=[C][CH]C=CC1(3442)',
    structure = SMILES('CC(=O)C1C=[C]C=C[CH]C1'),
    E0 = (244.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.477271,0.096461,-8.30847e-05,3.84846e-08,-7.36331e-12,29582.9,15.8108], Tmin=(100,'K'), Tmax=(1230.46,'K')), NASAPolynomial(coeffs=[15.8638,0.0433394,-1.83267e-05,3.39862e-09,-2.34705e-13,25561.5,-66.4143], Tmin=(1230.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(Allyl_S)"""),
)

species(
    label = 'CC(=O)C1C=C[CH][C]=CC1(3443)',
    structure = SMILES('CC(=O)C1C=C[CH][C]=CC1'),
    E0 = (261.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.750315,0.0142518,0.000221904,-3.35919e-07,1.41621e-10,31630.1,33.209], Tmin=(100,'K'), Tmax=(908.299,'K')), NASAPolynomial(coeffs=[40.1666,-0.00996386,1.52237e-05,-3.17221e-09,2.04303e-13,18308.3,-187.078], Tmin=(908.299,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C1OC2[CH][CH]C1CC=C2(3373)',
    structure = SMILES('C=C1OC2[CH][CH]C1CC=C2'),
    E0 = (261.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03992,0.0543729,0.000132938,-2.5524e-07,1.12908e-10,31645,9.35821], Tmin=(100,'K'), Tmax=(931.252,'K')), NASAPolynomial(coeffs=[49.214,-0.0164431,1.33841e-05,-2.40912e-09,1.36184e-13,15996,-263.275], Tmin=(931.252,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(Cycloheptane) + ring(Oxane) + radical(CCJCO) + radical(Cs_S)"""),
)

species(
    label = 'O=C1CC2C=C[CH][CH]C1C2(3124)',
    structure = SMILES('O=C1CC2[CH][CH]C=CC1C2'),
    E0 = (143.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.179081,0.0908091,-7.11651e-05,3.02088e-08,-5.58391e-12,17358,-2.07354], Tmin=(100,'K'), Tmax=(1206.55,'K')), NASAPolynomial(coeffs=[11.2706,0.0540381,-2.54511e-05,4.94998e-09,-3.50247e-13,14681.5,-57.6663], Tmin=(1206.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.228,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclopentane) + ring(Cycloheptene) + ring(Cyclopentane) + radical(Cs_S) + radical(Allyl_S)"""),
)

species(
    label = 'O=C1CC2[CH][CH]C=CCC12(3444)',
    structure = SMILES('O=C1CC2[CH]C=C[CH]CC12'),
    E0 = (90.9331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.531055,0.0921501,-6.7327e-05,2.45709e-08,-3.63973e-12,11105.6,-2.51232], Tmin=(100,'K'), Tmax=(1562.36,'K')), NASAPolynomial(coeffs=[19.8608,0.0399429,-1.72042e-05,3.18345e-09,-2.17468e-13,4733.62,-109.99], Tmin=(1562.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.9331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(Cycloheptene) + ring(Cyclobutane) + radical(Allyl_S) + radical(Allyl_S)"""),
)

species(
    label = 'CC(=O)C1C=C=CC=CC1(3445)',
    structure = SMILES('CC(=O)C1C=C=CC=CC1'),
    E0 = (127.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0564434,0.0735613,-2.73028e-05,-1.05688e-08,7.34695e-12,15476.2,26.4075], Tmin=(100,'K'), Tmax=(1109.64,'K')), NASAPolynomial(coeffs=[15.6411,0.0417745,-1.73071e-05,3.23633e-09,-2.26607e-13,10515.9,-57.1672], Tmin=(1109.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = 'CC(=O)C1C=CC=C=CC1(3446)',
    structure = SMILES('CC(=O)C1C=CC=C=CC1'),
    E0 = (127.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0564434,0.0735613,-2.73028e-05,-1.05688e-08,7.34695e-12,15476.2,26.4075], Tmin=(100,'K'), Tmax=(1109.64,'K')), NASAPolynomial(coeffs=[15.6411,0.0417745,-1.73071e-05,3.23633e-09,-2.26607e-13,10515.9,-57.1672], Tmin=(1109.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = 'O=[C]CC1[CH]C=CC=CC1(2245)',
    structure = SMILES('O=[C]CC1C=C[CH]C=CC1'),
    E0 = (213.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4440.56,'J/mol'), sigma=(7.29278,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=693.61 K, Pc=25.98 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.819605,0.00725946,0.000251568,-3.74444e-07,1.57485e-10,25875.7,34.5197], Tmin=(100,'K'), Tmax=(906.275,'K')), NASAPolynomial(coeffs=[43.2142,-0.0160425,1.90032e-05,-3.91762e-09,2.54642e-13,11464.2,-202.952], Tmin=(906.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(CCCJ=O) + radical(C=CCJC=C)"""),
)

species(
    label = 'O=C1CC2C=CCC1C=C2(2249)',
    structure = SMILES('O=C1CC2C=CCC1C=C2'),
    E0 = (-75.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.690531,0.108527,-0.000124819,8.23051e-08,-2.24267e-11,-8965.52,-0.999858], Tmin=(100,'K'), Tmax=(884.597,'K')), NASAPolynomial(coeffs=[12.8763,0.0471805,-2.07961e-05,3.91008e-09,-2.71329e-13,-11365.8,-64.7884], Tmin=(884.597,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclohexane) + ring(1,4-Cycloheptadiene) + ring(Cyclohexene)"""),
)

species(
    label = '[CH2]C1=CC=C[CH]C=CCO1(3447)',
    structure = SMILES('C=C1C=C[CH]C=C[CH]CO1'),
    E0 = (213.727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.2053,0.0900176,-1.90596e-05,-5.55043e-08,3.18895e-11,25915.1,15.64], Tmin=(100,'K'), Tmax=(954.55,'K')), NASAPolynomial(coeffs=[29.1396,0.0224636,-6.56836e-06,1.18792e-09,-9.10702e-14,17406.4,-143.569], Tmin=(954.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclononane) + radical(C=CCJC=C) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C1=CCC=C[CH]C=CO1(3448)',
    structure = SMILES('[CH2]C1=CCC=CC=C[CH]O1'),
    E0 = (257.588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.966483,0.0917852,-4.70997e-05,-1.28769e-08,1.36023e-11,31175,20.2994], Tmin=(100,'K'), Tmax=(987.788,'K')), NASAPolynomial(coeffs=[23.9973,0.0292129,-1.05712e-05,1.94531e-09,-1.39939e-13,24364.1,-109.342], Tmin=(987.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclononane) + radical(C=CCJ(O)C) + radical(C=C(O)CJ)"""),
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
    label = 'O=[C]C1C=C[CH]C=CC1(3449)',
    structure = SMILES('O=[C]C1C=C[CH]C=CC1'),
    E0 = (237.287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.549,-0.00760758,0.000263671,-3.78911e-07,1.57864e-10,28685,29.4346], Tmin=(100,'K'), Tmax=(905.224,'K')), NASAPolynomial(coeffs=[41.098,-0.0217846,2.10717e-05,-4.27686e-09,2.78907e-13,14945.6,-193.769], Tmin=(905.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(CC(C)CJ=O) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH2]C1([O])C2C=CC=CCC21(3450)',
    structure = SMILES('[CH2]C1([O])C2C=CC=CCC21'),
    E0 = (400.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.255213,0.11318,-0.000181617,1.97262e-07,-8.36439e-11,48313,11.8192], Tmin=(100,'K'), Tmax=(790.348,'K')), NASAPolynomial(coeffs=[-4.65938,0.0889336,-4.728e-05,9.44808e-09,-6.69903e-13,50462.6,41.2255], Tmin=(790.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Cyclopropane) + ring(1,3-Cycloheptadiene) + ring(Cyclopropane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
)

species(
    label = 'C=C([O])C1=CC=CC=CC1(3451)',
    structure = SMILES('C=C([O])C1=CC=CC=CC1'),
    E0 = (139.622,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18641,0.0946572,-4.88675e-05,-2.30803e-08,2.11012e-11,16997.2,27.7705], Tmin=(100,'K'), Tmax=(936.83,'K')), NASAPolynomial(coeffs=[28.1712,0.0180352,-4.2033e-06,6.55487e-10,-4.87247e-14,9358.35,-123.36], Tmin=(936.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])[C]1CC=CC=CC1(3452)',
    structure = SMILES('[CH2]C([O])=C1CC=CC=CC1'),
    E0 = (219.064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102495,0.0213271,0.000228736,-3.65965e-07,1.58808e-10,26548.3,33.7538], Tmin=(100,'K'), Tmax=(900.255,'K')), NASAPolynomial(coeffs=[49.052,-0.0261678,2.46233e-05,-5.05983e-09,3.36335e-13,10846.1,-235.516], Tmin=(900.255,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cycloheptane) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1C[C]=CC=CC1(3453)',
    structure = SMILES('C=C([O])C1C[C]=CC=CC1'),
    E0 = (306.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.960811,0.0942007,-6.03883e-05,3.28831e-09,7.59372e-12,37017.8,20.7422], Tmin=(100,'K'), Tmax=(991.921,'K')), NASAPolynomial(coeffs=[22.8107,0.0298496,-1.07246e-05,1.93435e-09,-1.36447e-13,30751.8,-101.563], Tmin=(991.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C([O])C1CC=[C]C=CC1(3454)',
    structure = SMILES('C=C([O])C1CC=[C]C=CC1'),
    E0 = (267.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00132,0.0960983,-6.45681e-05,6.79049e-09,6.80822e-12,32346.3,20.0102], Tmin=(100,'K'), Tmax=(975.067,'K')), NASAPolynomial(coeffs=[22.3952,0.0305439,-1.0526e-05,1.84173e-09,-1.27658e-13,26337.4,-99.6901], Tmin=(975.067,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C([O])C1CC=CC=CC1(3455)',
    structure = SMILES('[CH]=C([O])C1CC=CC=CC1'),
    E0 = (315.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.988977,0.0929775,-5.17603e-05,-9.21865e-09,1.28948e-11,38133.6,20.7727], Tmin=(100,'K'), Tmax=(974.773,'K')), NASAPolynomial(coeffs=[24.1476,0.0276725,-9.50189e-06,1.70744e-09,-1.21974e-13,31435.2,-109.077], Tmin=(974.773,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)C1[CH]C=[C]C=CC1(3456)',
    structure = SMILES('C=C(O)C1[CH]C=[C]C=CC1'),
    E0 = (270.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,350,440,435,1725,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13878,0.0944791,-4.67676e-05,-2.05852e-08,1.86158e-11,32753.5,17.8964], Tmin=(100,'K'), Tmax=(951.785,'K')), NASAPolynomial(coeffs=[26.0491,0.0249152,-7.57722e-06,1.30415e-09,-9.35659e-14,25553.6,-122.561], Tmin=(951.785,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Allyl_S) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C([O])C1C[CH]C2C=CC21(3457)',
    structure = SMILES('C=C([O])C1C[CH]C2C=CC21'),
    E0 = (343.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.696157,0.0502744,4.73807e-05,-9.44977e-08,3.85568e-11,41484.2,29.9248], Tmin=(100,'K'), Tmax=(987.042,'K')), NASAPolynomial(coeffs=[17.6568,0.0359857,-1.36438e-05,2.60284e-09,-1.91302e-13,35483.9,-65.1139], Tmin=(987.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(Cs_S)"""),
)

species(
    label = 'C=C([O])C1CC2[CH]C=CC21(3458)',
    structure = SMILES('C=C([O])C1CC2[CH]C=CC21'),
    E0 = (247.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.690679,0.0461126,6.98489e-05,-1.23888e-07,5.03957e-11,29924.3,27.1187], Tmin=(100,'K'), Tmax=(970.339,'K')), NASAPolynomial(coeffs=[19.9561,0.0327449,-1.15893e-05,2.21281e-09,-1.66412e-13,23076,-81.2681], Tmin=(970.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_4_5_ene_1) + radical(cyclopentene-allyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C=CC=CC=CC(=C)[O](3459)',
    structure = SMILES('[CH2]C=CC=CC=CC(=C)[O]'),
    E0 = (203.407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3001,3007,3013,3019,3025,975,980,985,990,995,1000,1300,1315,1330,1345,1360,1375,400,420,440,460,480,500,1630,1640,1650,1660,1670,1680,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07636,'amu*angstrom^2'), symmetry=1, barrier=(24.7475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07689,'amu*angstrom^2'), symmetry=1, barrier=(24.7599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07778,'amu*angstrom^2'), symmetry=1, barrier=(24.7804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0778,'amu*angstrom^2'), symmetry=1, barrier=(24.7807,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.991524,0.0917953,-4.29473e-05,-2.51293e-08,2.09787e-11,24660.5,35.2898], Tmin=(100,'K'), Tmax=(932.52,'K')), NASAPolynomial(coeffs=[25.593,0.0233372,-6.13875e-06,9.65456e-10,-6.75132e-14,17720.8,-101.733], Tmin=(932.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)C1=CC=CC=CC1(3460)',
    structure = SMILES('C=C(O)C1=CC=CC=CC1'),
    E0 = (1.81738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.54465,0.0985264,-4.09182e-05,-4.08639e-08,2.95008e-11,439.813,27.5557], Tmin=(100,'K'), Tmax=(931.416,'K')), NASAPolynomial(coeffs=[31.8794,0.0147277,-2.17642e-06,2.70593e-10,-2.38586e-14,-8377.95,-145.232], Tmin=(931.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.81738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'C=C([O])C1[CH]CC=CC=C1(2238)',
    structure = SMILES('C=C([O])C1[CH]CC=CC=C1'),
    E0 = (262.341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04339,0.10061,-8.92457e-05,4.06621e-08,-7.38355e-12,31742.4,20.7418], Tmin=(100,'K'), Tmax=(1326.26,'K')), NASAPolynomial(coeffs=[21.3094,0.0331933,-1.29969e-05,2.33409e-09,-1.58673e-13,25813.3,-93.409], Tmin=(1326.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])[CH]C1C=CC=CC1(3461)',
    structure = SMILES('[CH2]C([O])=CC1C=CC=CC1'),
    E0 = (206.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.227949,0.0766482,-2.09872e-05,-3.08845e-08,1.83348e-11,24968.3,30.8308], Tmin=(100,'K'), Tmax=(978.032,'K')), NASAPolynomial(coeffs=[19.6452,0.0329817,-1.17009e-05,2.1058e-09,-1.48907e-13,19282.2,-73.8005], Tmin=(978.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C1OC2C=CC=CCC12(3411)',
    structure = SMILES('C=C1OC2C=CC=CCC12'),
    E0 = (-44.1526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.91853,0.152334,-0.00026791,2.63979e-07,-1.00472e-10,-5118.78,-15.5205], Tmin=(100,'K'), Tmax=(815.564,'K')), NASAPolynomial(coeffs=[4.3282,0.0751582,-4.0373e-05,8.01697e-09,-5.63114e-13,-4589.96,-34.895], Tmin=(815.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.1526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(1,3-Cycloheptadiene) + ring(Oxetane)"""),
)

species(
    label = '[CH2]C1([O])C2C=CC=CC1C2(3462)',
    structure = SMILES('[CH2]C1([O])C2C=CC=CC1C2'),
    E0 = (269.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.37184,0.174649,-0.000347888,3.58083e-07,-1.36124e-10,32599.1,-14.8757], Tmin=(100,'K'), Tmax=(852.2,'K')), NASAPolynomial(coeffs=[-0.709148,0.0849308,-4.57886e-05,8.96276e-09,-6.18821e-13,35290.2,-5.17947], Tmin=(852.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(1,3-Cycloheptadiene) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
)

species(
    label = 'C=C([O])C1C=CC=CC=C1(3463)',
    structure = SMILES('C=C([O])C1C=CC=CC=C1'),
    E0 = (164.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.930506,0.0933263,-6.03638e-05,-2.68579e-10,1.02316e-11,20030.7,28.3633], Tmin=(100,'K'), Tmax=(960.374,'K')), NASAPolynomial(coeffs=[23.9698,0.024692,-7.95062e-06,1.37861e-09,-9.7214e-14,13630.4,-99.1814], Tmin=(960.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(164.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1C=CC=[C]CC1(3464)',
    structure = SMILES('C=C([O])C1C=CC=[C]CC1'),
    E0 = (305.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19627,0.105658,-0.000100276,4.91801e-08,-9.58357e-12,36954.4,19.2041], Tmin=(100,'K'), Tmax=(1243.31,'K')), NASAPolynomial(coeffs=[21.3316,0.0331805,-1.28352e-05,2.29407e-09,-1.55895e-13,31352.6,-94.3856], Tmin=(1243.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])[C]1C=CC=CCC1(3465)',
    structure = SMILES('C=C([O])C1=C[CH]C=CCC1'),
    E0 = (166.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.354442,0.0124675,0.000257393,-3.96227e-07,1.69594e-10,20183.2,31.964], Tmin=(100,'K'), Tmax=(901.444,'K')), NASAPolynomial(coeffs=[49.1487,-0.0259374,2.49216e-05,-5.11515e-09,3.38541e-13,4149.51,-238.517], Tmin=(901.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C([O])C1C=C[C]=CCC1(3466)',
    structure = SMILES('C=C([O])C1C=C[C]=CCC1'),
    E0 = (266.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24934,0.107725,-0.000105173,5.38173e-08,-1.09522e-11,32283.4,18.5159], Tmin=(100,'K'), Tmax=(1193.57,'K')), NASAPolynomial(coeffs=[20.6714,0.0342624,-1.2849e-05,2.24972e-09,-1.50995e-13,27050.7,-91.1175], Tmin=(1193.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C([O])C1[C]=CC=CCC1(3467)',
    structure = SMILES('C=C([O])C1[C]=CC=CCC1'),
    E0 = (305.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19627,0.105658,-0.000100276,4.91801e-08,-9.58357e-12,36954.4,19.2041], Tmin=(100,'K'), Tmax=(1243.31,'K')), NASAPolynomial(coeffs=[21.3316,0.0331805,-1.28352e-05,2.29407e-09,-1.55895e-13,31352.6,-94.3856], Tmin=(1243.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C([O])C1C=[C]C=CCC1(3468)',
    structure = SMILES('C=C([O])C1C=[C]C=CCC1'),
    E0 = (266.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24934,0.107725,-0.000105173,5.38173e-08,-1.09522e-11,32283.4,18.5159], Tmin=(100,'K'), Tmax=(1193.57,'K')), NASAPolynomial(coeffs=[20.6714,0.0342624,-1.2849e-05,2.24972e-09,-1.50995e-13,27050.7,-91.1175], Tmin=(1193.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C([O])C1C=CC=CCC1(3469)',
    structure = SMILES('[CH]=C([O])C1C=CC=CCC1'),
    E0 = (314.895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,350,440,435,1725,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46759,0.107233,-0.00010112,4.85174e-08,-9.15322e-12,38080.9,20.1114], Tmin=(100,'K'), Tmax=(1292.08,'K')), NASAPolynomial(coeffs=[23.6978,0.0293259,-1.06749e-05,1.85075e-09,-1.23787e-13,31577.8,-107.746], Tmin=(1292.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1[CH]C2C=CC2C1(3470)',
    structure = SMILES('C=C([O])C1[CH]C2C=CC2C1'),
    E0 = (343.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.696157,0.0502744,4.73807e-05,-9.44977e-08,3.85568e-11,41484.2,29.9248], Tmin=(100,'K'), Tmax=(987.042,'K')), NASAPolynomial(coeffs=[17.6568,0.0359857,-1.36438e-05,2.60284e-09,-1.91302e-13,35483.9,-65.1139], Tmin=(987.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(Cs_S)"""),
)

species(
    label = 'C=C(O)C1C=CC=CC=C1(3471)',
    structure = SMILES('C=C(O)C1C=CC=CC=C1'),
    E0 = (27.1507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28074,0.0971032,-5.21078e-05,-1.84154e-08,1.87653e-11,3472.97,28.1197], Tmin=(100,'K'), Tmax=(946.62,'K')), NASAPolynomial(coeffs=[27.6263,0.021472,-5.97396e-06,1.00555e-09,-7.33264e-14,-4083.99,-120.763], Tmin=(946.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.1507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'C=C([O])C1=CC=CC[CH]C1(3472)',
    structure = SMILES('C=C([O])C1=CC=CC[CH]C1'),
    E0 = (241.428,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28193,0.104429,-9.72158e-05,4.68042e-08,-8.9143e-12,29237.1,20.581], Tmin=(100,'K'), Tmax=(1276.54,'K')), NASAPolynomial(coeffs=[22.0477,0.0313277,-1.1319e-05,1.94563e-09,-1.29242e-13,23280.7,-97.6671], Tmin=(1276.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(RCCJCC) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1C=CC=CC1C(=C)[O](3473)',
    structure = SMILES('[CH2]C1C=CC=CC1C(=C)[O]'),
    E0 = (256.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.154912,0.0741253,-1.17509e-05,-4.28194e-08,2.34639e-11,31033.2,31.1609], Tmin=(100,'K'), Tmax=(956.119,'K')), NASAPolynomial(coeffs=[19.9857,0.03155,-1.03535e-05,1.8047e-09,-1.26992e-13,25276.5,-75.066], Tmin=(956.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cyclohexadiene) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = 'C=C1OC2C=CC=CC1C2(3359)',
    structure = SMILES('C=C1OC2C=CC=CC1C2'),
    E0 = (-44.7007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.18919,0.16425,-0.000309637,3.12605e-07,-1.1894e-10,-5180.64,-16.9347], Tmin=(100,'K'), Tmax=(835.378,'K')), NASAPolynomial(coeffs=[2.1848,0.0795493,-4.30676e-05,8.50997e-09,-5.93332e-13,-3686.77,-23.9347], Tmin=(835.378,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.7007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Tetrahydrofuran) + ring(1,3-Cycloheptadiene) + ring(Tetrahydrofuran)"""),
)

species(
    label = 'CC(=O)C1C=C[C]=C[CH]C1(3474)',
    structure = SMILES('CC(=O)C1C=C[C]=C[CH]C1'),
    E0 = (244.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.477271,0.096461,-8.30847e-05,3.84846e-08,-7.36331e-12,29582.9,15.8108], Tmin=(100,'K'), Tmax=(1230.46,'K')), NASAPolynomial(coeffs=[15.8638,0.0433394,-1.83267e-05,3.39862e-09,-2.34705e-13,25561.5,-66.4143], Tmin=(1230.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(Allyl_S)"""),
)

species(
    label = 'CC(=O)C1C=CC=CC=C1(3475)',
    structure = SMILES('CC(=O)C1C=CC=CC=C1'),
    E0 = (-5.24664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11893,0.0988779,-8.40889e-05,3.63453e-08,-6.23677e-12,-435.304,29.0429], Tmin=(100,'K'), Tmax=(1403.64,'K')), NASAPolynomial(coeffs=[22.5326,0.0314772,-1.20609e-05,2.13513e-09,-1.43629e-13,-7074.93,-93.0815], Tmin=(1403.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.24664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'O=C1CC2C=CC=CC1C2(3143)',
    structure = SMILES('O=C1CC2C=CC=CC1C2'),
    E0 = (-95.5158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.85059,0.17071,-0.000363423,3.93675e-07,-1.53891e-10,-11317.7,-15.4845], Tmin=(100,'K'), Tmax=(856.811,'K')), NASAPolynomial(coeffs=[-8.39068,0.0972592,-5.27943e-05,1.03393e-08,-7.12965e-13,-6380.17,37.3304], Tmin=(856.811,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.5158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclopentane) + ring(1,3-Cycloheptadiene) + ring(Cyclopentane)"""),
)

species(
    label = 'CC(=O)C1=CC=CC=CC1(3476)',
    structure = SMILES('CC(=O)C1=CC=CC=CC1'),
    E0 = (-10.4574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.617167,0.0862913,-4.65193e-05,-2.15602e-09,6.80161e-12,-1078.44,29.1566], Tmin=(100,'K'), Tmax=(1066.36,'K')), NASAPolynomial(coeffs=[20.5926,0.0351051,-1.44297e-05,2.73434e-09,-1.94743e-13,-7215.1,-82.0953], Tmin=(1066.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.4574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'O=C1CC2C=CC=CCC12(2253)',
    structure = SMILES('O=C1CC2C=CC=CCC12'),
    E0 = (-94.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.79144,0.167054,-0.00034929,3.76319e-07,-1.47037e-10,-11181.1,-15.0288], Tmin=(100,'K'), Tmax=(854.275,'K')), NASAPolynomial(coeffs=[-7.28853,0.0951242,-5.14965e-05,1.00923e-08,-6.96947e-13,-6678.07,31.4847], Tmin=(854.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(1,3-Cycloheptadiene) + ring(Cyclobutane)"""),
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
    E0 = (187.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (392.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (207.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (308.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (519.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (519.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (319.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (480.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (480.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (345.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (488.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (472.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (472.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (296.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (331.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (321.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (418.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (329.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (314.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (376.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (332.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (499.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (418.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (432.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (432.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (300.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (230.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (236.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (204.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (204.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (194.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (535.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (334.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (744.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (207.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (256.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (318.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (243.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (305.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (296.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (284.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (333.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (266.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (256.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (253.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (204.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (204.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (459.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (195.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (527.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (571.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (618.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (412.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (359.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (327.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (342.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (506.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (414.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (359.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (497.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (414.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (366.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (310.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (265.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (308.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (358.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (194.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (289.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (376.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (505.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (384.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (414.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (349.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (527.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (352.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (414.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (204.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (357.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS79',
    E0 = (416.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS80',
    E0 = (194.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS81',
    E0 = (387.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS82',
    E0 = (204.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS83',
    E0 = (194.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS84',
    E0 = (265.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS85',
    E0 = (194.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C7H8(11)', 'CH2CO(43)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['[CH2]C1([O])C2C=CCC1C=C2(3413)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.249e+08,'s^-1'), n=0.846, Ea=(205.301,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS_D;doublebond_intra_2H;radadd_intra_cs] for rate rule [R6_SMS_D;doublebond_intra_2H;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C1OC2[CH]CC1[CH]C=C2(3405)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(414000,'s^-1'), n=1.14, Ea=(20.92,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R8;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C1OC2[CH]C=CCC1[CH]2(3372)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.291e+11,'s^-1'), n=0.234, Ea=(121.888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra;radadd_intra_O] for rate rule [R9;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(19)', 'C=C([O])C1C=CC=C=CC1(3414)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(19)', 'C=C([O])C1C=C=CC=CC1(3415)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', 'CH2CO(43)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00376298,'m^3/(mol*s)'), n=2.445, Ea=(46.6462,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CsJ-CdCsH] for rate rule [Ck_O;CsJ-CdCsH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([O])C1C=CC[C]=CC1(3416)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C1C=[C]CC=CC1(3417)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C(O)[C]1C=C[CH]C=CC1(3418)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(923167,'s^-1'), n=2.18104, Ea=(158.471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_noH] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out_noH] for rate rule [R3H_SS_2Cd;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(O)C1C=C[CH]C=CC1(3419)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C([O])C1C=CCC=[C]C1(3420)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C([O])C1[C]=CCC=CC1(3421)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C(O)C1[CH]C=C[CH]C=C1(3422)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.56e+06,'s^-1'), n=1.75, Ea=(105.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SS(Cd)S;C_rad_out_H/Cd;XH_out] for rate rule [R4H_SS(Cd)S;C_rad_out_H/Cd;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C(O)C1[C]=C[CH]C=CC1(3423)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C([O])[C]1C=CCC=CC1(3424)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.30755e+06,'s^-1'), n=1.805, Ea=(115.478,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_SDS;C_rad_out_OneDe/Cs;Cs_H_out_H/Cd]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(285077,'s^-1'), n=1.88833, Ea=(191.836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] + [R4H_SDS;C_rad_out_single;Cs_H_out_H/Cd] + [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C(O)C1C=C[CH]C=[C]C1(3425)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);Y_rad_out;XH_out] for rate rule [R5H_CCC(Cd);Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C(O)C1C=[C][CH]C=CC1(3426)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C([O])C1C=CCC=CC1(3427)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R6H_RSSMS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C(O)C1C=C[CH][C]=CC1(3428)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(13085.7,'s^-1'), n=1.78967, Ea=(62.3074,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;XH_out] for rate rule [R6H_DSSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', '[CH2][C]=O(360)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.42928e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H/CdCs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['[CH]1C=CCC(C=C1)[C]1CO1(3429)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C([O])C1C=CC2[CH]C2C1(3430)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_csHDe] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C([O])C1CC=CC2[CH]C21(3431)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_csHDe] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['O=C1CC2[CH][CH]C1CC=C2(3105)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.68243e+09,'s^-1'), n=0.4695, Ea=(113.018,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SDS_D;doublebond_intra;radadd_intra_cs] for rate rule [R6_SDS_D;doublebond_intra;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C1OC2C=C[CH][CH]C1C2(3330)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.95208e+10,'s^-1'), n=0.228823, Ea=(43.0818,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_beta_long;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C1OC2[CH][CH]C=CCC12(3432)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C(O)C1C=C=CC=CC1(3433)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C(O)C1C=CC=C=CC1(3434)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_De] for rate rule [R6radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C1OC2C=CCC1C=C2(2247)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=CC1[CH]C1C=CC(=C)[O](3435)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=CC1[CH]C=CC=C([O])C1(3436)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O(20)', 'C=[C]C1C=C[CH]C=CC1(3437)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['O=C1CC2[CH]CC1[CH]C=C2(3062)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(414000,'s^-1'), n=1.14, Ea=(20.92,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R7plus;doublebond_intra_HNd_pri;radadd_intra_cs2H] for rate rule [R8;doublebond_intra_HNd_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['O=C1CC2[CH]C1C[CH]C=C2(3088)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.58148e+08,'s^-1'), n=0.696039, Ea=(69.5697,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra;radadd_intra_cs2H] for rate rule [R9;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['CC(=O)[C]1C=C[CH]C=CC1(3438)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(25964.4,'s^-1'), n=2.44333, Ea=(131.615,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['CC(=O)C1[CH]C=C[CH]C=C1(3439)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(42200,'s^-1'), n=1.93, Ea=(56.484,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 86 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CC(=O)C1[C]=C[CH]C=CC1(3440)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CC(=O)C1C=C[CH]C=[C]C1(3441)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(63.58,'s^-1','*|/',3.66), n=2.81162, Ea=(34.4385,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CC(O2d)CC;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_CC(O2d)CC;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CC(=O)C1C=[C][CH]C=CC1(3442)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(263079,'s^-1'), n=1.73643, Ea=(39.8993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CC(=O)C1C=C[CH][C]=CC1(3443)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.17493e+06,'s^-1'), n=1.01509, Ea=(71.5445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_DSSSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C1OC2[CH][CH]C1CC=C2(3373)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4.00901e+09,'s^-1'), n=0.463766, Ea=(79.7739,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHCd] for rate rule [R6_SMS_CO;carbonyl_intra;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['O=C1CC2C=C[CH][CH]C1C2(3124)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(6.80916e+06,'s^-1'), n=0.895, Ea=(69.1406,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H] for rate rule [Rn2c7_beta_long;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['O=C1CC2[CH][CH]C=CCC12(3444)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.54811e+07,'s^-1'), n=0.785425, Ea=(66.8554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['CC(=O)C1C=C=CC=CC1(3445)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['CC(=O)C1C=CC=C=CC1(3446)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_De] for rate rule [R6radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['O=C1CC2C=CCC1C=C2(2249)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;C_rad_out_2H;Cpri_rad_out_single] + [R6_SSSDS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R6_SSSDS;C_rad_out_2H;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C1=CC=C[CH]C=CCO1(3447)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C1=CCC=C[CH]C=CO1(3448)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction52',
    reactants = ['CH2(T)(33)', 'O=[C]C1C=C[CH]C=CC1(3449)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction53',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['[CH2]C1([O])C2C=CC=CCC21(3450)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(225.593,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction54',
    reactants = ['H(19)', 'C=C([O])C1=CC=CC=CC1(3451)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1.26e+08,'cm^3/(mol*s)'), n=1.63, Ea=(7.9496,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for Cds-CdCs_Cds-(Cd-Cd-Cd)H;HJ
Exact match found for rate rule [Cds-CdCs_Cds-(Cd-Cd-Cd)H;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C7H8(11)', '[CH2][C]=O(360)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.11106,'m^3/(mol*s)'), n=1.92886, Ea=(1.83261,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds-OneDeH;CJ] + [Cds-CsH_Cds-(Cd-Cd-Cd)H;YJ] for rate rule [Cds-CsH_Cds-(Cd-Cd-Cd)H;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction56',
    reactants = ['C=C([O])[C]1CC=CC=CC1(3452)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(22577.1,'s^-1'), n=2.7035, Ea=(123.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['C=C([O])C1C[C]=CC=CC1(3453)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['C=C([O])C1CC=[C]C=CC1(3454)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH]=C([O])C1CC=CC=CC1(3455)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(148400,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 126 used for R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/(Cd-Cd-Cd)
Exact match found for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['C=C(O)C1[CH]C=[C]C=CC1(3456)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1.61756e+09,'s^-1'), n=1.29078, Ea=(226.405,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C([O])C1C[CH]C2C=CC21(3457)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(5.65845e+18,'s^-1'), n=-1.57151, Ea=(227.049,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 131 used for Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd
Exact match found for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction62',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C([O])C1CC2[CH]C=CC21(3458)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(6.33004e+14,'s^-1'), n=-0.792922, Ea=(179.405,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 112 used for Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2]C=CC=CC=CC(=C)[O](3459)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(4.72071e+06,'s^-1'), n=0.958205, Ea=(106.629,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R7_SDSD_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R7_SDSD_D;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction64',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C(O)C1=CC=CC=CC1(3460)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction65',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(4.25e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;CdHC_2] for rate rule [1_3_pentadiene;CH(CJ)_1;CdHC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction66',
    reactants = ['C=C([O])[CH]C1C=CC=CC1(3461)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(7.65366e+10,'s^-1'), n=0.611, Ea=(152.377,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ-CdH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction67',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C1OC2C=CC=CCC12(3411)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction68',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['[CH2]C1([O])C2C=CC=CC1C2(3462)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1.62826e+06,'s^-1'), n=1.44767, Ea=(102.625,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_2H;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_2H;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction69',
    reactants = ['H(19)', 'C=C([O])C1C=CC=CC=C1(3463)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(187.105,'m^3/(mol*s)'), n=1.3603, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 20 used for Cds-CsH_Cds-(Cd-Cd-Cd)H;HJ
Exact match found for rate rule [Cds-CsH_Cds-(Cd-Cd-Cd)H;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -1.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction70',
    reactants = ['C=C([O])C1C=CC=[C]CC1(3464)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction71',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C([O])[C]1C=CC=CCC1(3465)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(4.28857e+08,'s^-1'), n=1.2, Ea=(197.485,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS_Cs;C_rad_out_H/OneDe;Cs_H_out_noH] + [R3H_SS_Cs;C_rad_out_H/(Cd-Cd-Cd);Cs_H_out] for rate rule [R3H_SS_Cs;C_rad_out_H/(Cd-Cd-Cd);Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction72',
    reactants = ['C=C([O])C1C=C[C]=CCC1(3466)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction73',
    reactants = ['C=C([O])C1[C]=CC=CCC1(3467)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction74',
    reactants = ['C=C([O])C1C=[C]C=CCC1(3468)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(2.10713e+09,'s^-1'), n=0.59575, Ea=(261.003,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SDS;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] + [R4H_SDS;Cd_rad_out_Cd;Cs_H_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction75',
    reactants = ['[CH]=C([O])C1C=CC=CCC1(3469)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(14044.1,'s^-1'), n=1.88125, Ea=(37.5252,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 3.74165738677
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction76',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C([O])C1[CH]C2C=CC2C1(3470)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(5.65845e+18,'s^-1'), n=-1.57151, Ea=(227.049,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 131 used for Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd
Exact match found for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction77',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C(O)C1C=CC=CC=C1(3471)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction78',
    reactants = ['C=C([O])C1=CC=CC[CH]C1(3472)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS78',
    kinetics = Arrhenius(A=(3.94565e+09,'s^-1'), n=0.909333, Ea=(116.555,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH(CJ)_1;unsaturated_end] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH(CJ)_1;Cd(C)C_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction79',
    reactants = ['[CH2]C1C=CC=CC1C(=C)[O](3473)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS79',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction80',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C1OC2C=CC=CC1C2(3359)'],
    transitionState = 'TS80',
    kinetics = Arrhenius(A=(1.24579e+11,'s^-1'), n=0.1555, Ea=(7.322,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R5_SSSS;C_rad_out_single;Ypri_rad_out] for rate rule [R5_SSSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction81',
    reactants = ['CC(=O)C1C=C[C]=C[CH]C1(3474)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS81',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6H_SMSSR;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction82',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['CC(=O)C1C=CC=CC=C1(3475)'],
    transitionState = 'TS82',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction83',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['O=C1CC2C=CC=CC1C2(3143)'],
    transitionState = 'TS83',
    kinetics = Arrhenius(A=(1.24579e+11,'s^-1'), n=0.1555, Ea=(7.322,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R5_SSSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R5_SSSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction84',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['CC(=O)C1=CC=CC=CC1(3476)'],
    transitionState = 'TS84',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction85',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['O=C1CC2C=CC=CCC12(2253)'],
    transitionState = 'TS85',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

network(
    label = '492',
    isomers = [
        'C=C([O])C1[CH]C=CC=CC1(2239)',
    ],
    reactants = [
        ('C7H8(11)', 'CH2CO(43)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '492',
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

