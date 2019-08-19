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
    label = '[CH2]C1([O])C2C=CCC=CC21(3477)',
    structure = SMILES('[CH2]C1([O])C2C=CCC=CC21'),
    E0 = (419.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.827793,0.039177,9.66512e-05,-1.51369e-07,5.92688e-11,50553.7,28.1339], Tmin=(100,'K'), Tmax=(980.562,'K')), NASAPolynomial(coeffs=[20.557,0.0348865,-1.33369e-05,2.65133e-09,-2.02587e-13,43021.7,-85.3382], Tmin=(980.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Cyclopropane) + ring(1,4-Cycloheptadiene) + ring(Cyclopropane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
)

species(
    label = 'C=C1OC2[CH]CC=CC1[CH]2(3401)',
    structure = SMILES('C=C1OC2[CH]CC=CC1[CH]2'),
    E0 = (259.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16084,0.0945834,-4.83993e-05,-1.4339e-08,1.44257e-11,31425,2.73497], Tmin=(100,'K'), Tmax=(995.59,'K')), NASAPolynomial(coeffs=[25.6793,0.0280848,-1.04908e-05,1.98146e-09,-1.44907e-13,24032,-136.924], Tmin=(995.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Tetrahydrofuran) + ring(Cycloheptene) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(CCJCO)"""),
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
    label = 'C=C([O])C1C=CCC=CC=1(3478)',
    structure = SMILES('C=C([O])C1C=CCC=CC=1'),
    E0 = (160.288,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.58682,0.0839792,-3.59616e-05,-2.21476e-08,1.6973e-11,19458.5,27.2251], Tmin=(100,'K'), Tmax=(966.278,'K')), NASAPolynomial(coeffs=[22.6556,0.0272367,-9.15148e-06,1.63023e-09,-1.16433e-13,13124.1,-93.6442], Tmin=(966.278,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1C=C=CCC=C1(3479)',
    structure = SMILES('C=C([O])C1C=C=CCC=C1'),
    E0 = (315.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.231406,0.0696703,-2.22791e-05,-1.6981e-08,1.03342e-11,38042.5,27.7018], Tmin=(100,'K'), Tmax=(1048.19,'K')), NASAPolynomial(coeffs=[15.6085,0.0377857,-1.4997e-05,2.77607e-09,-1.94841e-13,33346.8,-54.2294], Tmin=(1048.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(C=C(C)OJ)"""),
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
    label = 'C=C(O)[C]1[CH]C=CCC=C1(3480)',
    structure = SMILES('C=C(O)C1[CH]C=CC[CH]C=1'),
    E0 = (169.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.189035,0.011177,0.000274061,-4.22189e-07,1.80836e-10,20591.7,29.9504], Tmin=(100,'K'), Tmax=(902.101,'K')), NASAPolynomial(coeffs=[52.9572,-0.0318283,2.80214e-05,-5.68835e-09,3.75586e-13,3300.64,-262.258], Tmin=(902.101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH]=C(O)C1[CH]C=CCC=C1(3481)',
    structure = SMILES('[CH]=C(O)C1[CH]C=CCC=C1'),
    E0 = (335.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3120,650,792.5,1650,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00344269,0.0261238,0.000210873,-3.42189e-07,1.48361e-10,40577.3,32.3709], Tmin=(100,'K'), Tmax=(905.068,'K')), NASAPolynomial(coeffs=[47.616,-0.0219142,2.13566e-05,-4.35167e-09,2.84951e-13,25307.8,-229.326], Tmin=(905.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(Cds_P)"""),
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
    label = 'C=C(O)C1[C]=CCC=C[CH]1(3482)',
    structure = SMILES('C=C(O)C1[C]=CCC=C[CH]1'),
    E0 = (326.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0492445,0.027148,0.000202886,-3.30403e-07,1.43306e-10,39460.7,32.2768], Tmin=(100,'K'), Tmax=(905.237,'K')), NASAPolynomial(coeffs=[46.1463,-0.0195141,2.00065e-05,-4.09489e-09,2.68011e-13,24681.1,-221.062], Tmin=(905.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(Allyl_S)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4629.19,'J/mol'), sigma=(7.55705,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=723.07 K, Pc=24.34 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.554994,0.0115722,0.000248281,-3.78125e-07,1.60819e-10,22677.8,32.2663], Tmin=(100,'K'), Tmax=(903.379,'K')), NASAPolynomial(coeffs=[45.8044,-0.0203488,2.1608e-05,-4.45426e-09,2.92929e-13,7629.38,-219.478], Tmin=(903.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ) + radical(C=CCJC=C)"""),
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
    label = 'C=C(O)C1[CH]C=CC[C]=C1(3483)',
    structure = SMILES('C=C(O)C1[CH]C=CC[C]=C1'),
    E0 = (326.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0492445,0.027148,0.000202886,-3.30403e-07,1.43306e-10,39460.7,32.2768], Tmin=(100,'K'), Tmax=(905.237,'K')), NASAPolynomial(coeffs=[46.1463,-0.0195141,2.00065e-05,-4.09489e-09,2.68011e-13,24681.1,-221.062], Tmin=(905.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = 'C=C(O)C1[CH][C]=CCC=C1(3484)',
    structure = SMILES('C=C(O)C1[CH][C]=CCC=C1'),
    E0 = (326.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0492445,0.027148,0.000202886,-3.30403e-07,1.43306e-10,39460.7,32.2768], Tmin=(100,'K'), Tmax=(905.237,'K')), NASAPolynomial(coeffs=[46.1463,-0.0195141,2.00065e-05,-4.09489e-09,2.68011e-13,24681.1,-221.062], Tmin=(905.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(Cds_S)"""),
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
    label = 'C=C(O)C1[CH]C=[C]CC=C1(3485)',
    structure = SMILES('C=C(O)C1[CH]C=[C]CC=C1'),
    E0 = (326.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0492445,0.027148,0.000202886,-3.30403e-07,1.43306e-10,39460.7,32.2768], Tmin=(100,'K'), Tmax=(905.237,'K')), NASAPolynomial(coeffs=[46.1463,-0.0195141,2.00065e-05,-4.09489e-09,2.68011e-13,24681.1,-221.062], Tmin=(905.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]1C=CCC=CC1[C]1CO1(3486)',
    structure = SMILES('[CH]1C=CC(C=CC1)[C]1CO1'),
    E0 = (371.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612743,0.0193369,0.000212392,-3.33703e-07,1.44671e-10,44832.3,29.7397], Tmin=(100,'K'), Tmax=(890.801,'K')), NASAPolynomial(coeffs=[40.4663,-0.0128217,1.93542e-05,-4.24092e-09,2.90646e-13,31907.6,-190.614], Tmin=(890.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + radical(C2CsJO) + radical(Allyl_S)"""),
)

species(
    label = 'C=C([O])C1C=CCC2[CH]C21(3487)',
    structure = SMILES('C=C([O])C1C=CCC2[CH]C21'),
    E0 = (327.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.20857,0.0696554,-1.77112e-05,-2.02151e-08,1.07483e-11,39538.4,29.15], Tmin=(100,'K'), Tmax=(1076.33,'K')), NASAPolynomial(coeffs=[15.2153,0.0415092,-1.69827e-05,3.17794e-09,-2.2352e-13,34707.9,-51.7856], Tmin=(1076.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_3_6_ene_2) + radical(cyclopropane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1[CH]C2CC=CC21(3488)',
    structure = SMILES('C=C([O])C1[CH]C2CC=CC21'),
    E0 = (308.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.649363,0.0520033,4.24263e-05,-8.8877e-08,3.64317e-11,37201.4,30.7217], Tmin=(100,'K'), Tmax=(990.148,'K')), NASAPolynomial(coeffs=[17.3889,0.0368407,-1.40793e-05,2.67925e-09,-1.95974e-13,31314.9,-62.8573], Tmin=(990.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_4_5_ene_1) + radical(cyclobutane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=C1CC2C=CC[CH][CH]C12(3090)',
    structure = SMILES('O=C1CC2C=CC[CH][CH]C12'),
    E0 = (204.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.168272,0.0900197,-7.31462e-05,3.34315e-08,-6.6878e-12,24685.2,2.28471], Tmin=(100,'K'), Tmax=(1129.59,'K')), NASAPolynomial(coeffs=[10.5161,0.0533771,-2.44881e-05,4.71434e-09,-3.32162e-13,22347.4,-48.8985], Tmin=(1129.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(Cycloheptene) + ring(Cyclobutane) + radical(CCJCC=O) + radical(RCCJCC)"""),
)

species(
    label = 'C=C1OC2[CH]C1[CH]C=CC2(3311)',
    structure = SMILES('C=C1OC2[CH]C1[CH]C=CC2'),
    E0 = (201.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.879981,0.0825788,-6.03279e-06,-6.25309e-08,3.24071e-11,24415.6,0.5093], Tmin=(100,'K'), Tmax=(972.434,'K')), NASAPolynomial(coeffs=[27.3501,0.0254063,-8.77311e-06,1.68618e-09,-1.28709e-13,16138,-149.227], Tmin=(972.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Tetrahydrofuran) + ring(Cycloheptene) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(Allyl_S)"""),
)

species(
    label = 'C=C1OC2[CH][CH]C1C=CC2(3343)',
    structure = SMILES('C=C1OC2[CH][CH]C1C=CC2'),
    E0 = (260.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32342,0.0690856,8.01946e-05,-1.93086e-07,8.90825e-11,31621.8,7.08081], Tmin=(100,'K'), Tmax=(933.594,'K')), NASAPolynomial(coeffs=[45.7697,-0.00983117,9.59987e-06,-1.72223e-09,9.37102e-14,17474.6,-245.555], Tmin=(933.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(Cycloheptane) + ring(Oxane) + radical(Cs_S) + radical(CCJCO)"""),
)

species(
    label = 'C=C([O])C1C2[CH]CC=CC21(3489)',
    structure = SMILES('C=C([O])C1C2[CH]CC=CC21'),
    E0 = (295.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.571728,0.0533315,4.05881e-05,-8.88145e-08,3.68843e-11,35743.2,30.7715], Tmin=(100,'K'), Tmax=(986.316,'K')), NASAPolynomial(coeffs=[18.1665,0.0353844,-1.33417e-05,2.53815e-09,-1.86349e-13,29674.5,-67.0399], Tmin=(986.316,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_3_6_ene_1) + radical(Cs_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1OC2[CH]CC=C[CH]C12(3406)',
    structure = SMILES('C=C1OC2[CH]CC=C[CH]C12'),
    E0 = (201.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.879981,0.0825788,-6.03279e-06,-6.25309e-08,3.24071e-11,24415.6,0.5093], Tmin=(100,'K'), Tmax=(972.434,'K')), NASAPolynomial(coeffs=[27.3501,0.0254063,-8.77311e-06,1.68618e-09,-1.28709e-13,16138,-149.227], Tmin=(972.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(Cycloheptene) + ring(Oxetane) + radical(CCJCO) + radical(Allyl_S)"""),
)

species(
    label = '[CH]=CCC=CC=CC(=C)[O](3490)',
    structure = SMILES('[CH]=CCC=CC=CC(=C)[O]'),
    E0 = (363.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.787332,'amu*angstrom^2'), symmetry=1, barrier=(18.1023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786551,'amu*angstrom^2'), symmetry=1, barrier=(18.0844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786478,'amu*angstrom^2'), symmetry=1, barrier=(18.0827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.78689,'amu*angstrom^2'), symmetry=1, barrier=(18.0921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.15307,0.0998955,-7.63033e-05,1.6492e-08,4.1078e-12,43892.7,37.7363], Tmin=(100,'K'), Tmax=(970.421,'K')), NASAPolynomial(coeffs=[23.7051,0.0269013,-9.02565e-06,1.566e-09,-1.08645e-13,37680.5,-88.5927], Tmin=(970.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C1C=CCC=CC=1(3491)',
    structure = SMILES('C=C(O)C1C=CCC=CC=1'),
    E0 = (22.4829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.938062,0.087771,-2.7774e-05,-4.01783e-08,2.54432e-11,2900.8,26.985], Tmin=(100,'K'), Tmax=(952.459,'K')), NASAPolynomial(coeffs=[26.3041,0.0240293,-7.18183e-06,1.25876e-09,-9.26745e-14,-4586.78,-115.18], Tmin=(952.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.4829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'C=C(O)C1C=C=CCC=C1(3492)',
    structure = SMILES('C=C(O)C1C=C=CCC=C1'),
    E0 = (177.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.10098,0.0732799,-1.36681e-05,-3.52018e-08,1.87065e-11,21483.9,27.3917], Tmin=(100,'K'), Tmax=(1000.7,'K')), NASAPolynomial(coeffs=[18.9564,0.0350697,-1.33026e-05,2.46825e-09,-1.76278e-13,15768.8,-74.0609], Tmin=(1000.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = 'C=C([O])[CH]C1C=CCC=C1(3493)',
    structure = SMILES('[CH2]C([O])=CC1C=CCC=C1'),
    E0 = (216.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.134155,0.0760948,-2.4428e-05,-2.35944e-08,1.48701e-11,26160.1,33.6729], Tmin=(100,'K'), Tmax=(990.645,'K')), NASAPolynomial(coeffs=[18.0929,0.0355132,-1.29709e-05,2.33681e-09,-1.63721e-13,20928.8,-62.2672], Tmin=(990.645,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cyclohexadiene) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C1OC2C=CCC=CC12(3380)',
    structure = SMILES('C=C1OC2C=CCC=CC12'),
    E0 = (-27.2028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06654,0.0902156,-3.28022e-05,-3.16806e-08,2.06627e-11,-3070.46,-0.804642], Tmin=(100,'K'), Tmax=(990.842,'K')), NASAPolynomial(coeffs=[26.3696,0.0273892,-1.02552e-05,1.97193e-09,-1.46742e-13,-10860.4,-144.789], Tmin=(990.842,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.2028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(1,4-Cycloheptadiene) + ring(Oxetane)"""),
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
    label = 'C=[C]C1[CH]C=CCC=C1(3494)',
    structure = SMILES('C=[C]C1[CH]C=CCC=C1'),
    E0 = (541.005,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (118.176,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05295,0.00841259,0.00022532,-3.34727e-07,1.40156e-10,65226.7,28.6154], Tmin=(100,'K'), Tmax=(909.615,'K')), NASAPolynomial(coeffs=[38.9323,-0.0111848,1.52663e-05,-3.14105e-09,2.00971e-13,52255.2,-183.965], Tmin=(909.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(Cds_S)"""),
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
    label = 'C=C([O])C1C=[C]CCC=C1(3495)',
    structure = SMILES('C=C([O])C1C=[C]CCC=C1'),
    E0 = (327.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.29353,0.0277332,0.000184719,-3.00753e-07,1.30328e-10,39590.2,34.4962], Tmin=(100,'K'), Tmax=(904.047,'K')), NASAPolynomial(coeffs=[41.3222,-0.0122859,1.6318e-05,-3.42288e-09,2.25058e-13,26388.9,-191.288], Tmin=(904.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1[C]=CCCC=C1(3496)',
    structure = SMILES('C=C([O])C1[C]=CCCC=C1'),
    E0 = (327.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.29353,0.0277332,0.000184719,-3.00753e-07,1.30328e-10,39590.2,34.4962], Tmin=(100,'K'), Tmax=(904.047,'K')), NASAPolynomial(coeffs=[41.3222,-0.0122859,1.6318e-05,-3.42288e-09,2.25058e-13,26388.9,-191.288], Tmin=(904.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)C1[C]=C[CH]CC=C1(3497)',
    structure = SMILES('C=C(O)C1[C]=C[CH]CC=C1'),
    E0 = (330.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.128954,0.0264321,0.000201429,-3.26777e-07,1.416e-10,39998.6,32.4797], Tmin=(100,'K'), Tmax=(904.617,'K')), NASAPolynomial(coeffs=[45.1292,-0.0181743,1.94163e-05,-3.99571e-09,2.62072e-13,25540.5,-215.022], Tmin=(904.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = 'C=C([O])[C]1C=CCCC=C1(3498)',
    structure = SMILES('[CH2]C(=O)[C]1C=CCCC=C1'),
    E0 = (187.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.601005,0.0139643,0.000233561,-3.51334e-07,1.47068e-10,22755.9,31.145], Tmin=(100,'K'), Tmax=(914.325,'K')), NASAPolynomial(coeffs=[42.3051,-0.0106005,1.48458e-05,-3.00348e-09,1.87246e-13,8530.25,-202.407], Tmin=(914.325,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJ(C=O)C=C) + radical(C2JC=O)"""),
)

species(
    label = '[CH]=C([O])C1C=CCCC=C1(3499)',
    structure = SMILES('[CH]=C([O])C1C=CCCC=C1'),
    E0 = (336.897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,350,440,435,1725,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.24708,0.0267172,0.000192675,-3.12496e-07,1.35364e-10,40706.8,34.5927], Tmin=(100,'K'), Tmax=(903.928,'K')), NASAPolynomial(coeffs=[42.7937,-0.0146893,1.767e-05,-3.68012e-09,2.42036e-13,27014.8,-199.562], Tmin=(903.928,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1C=CC2C[CH]C21(3500)',
    structure = SMILES('C=C([O])C1C=CC2C[CH]C21'),
    E0 = (307.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.370148,0.0639798,6.62467e-07,-4.03932e-08,1.80553e-11,37140,29.3408], Tmin=(100,'K'), Tmax=(1039.36,'K')), NASAPolynomial(coeffs=[15.4979,0.0408366,-1.65599e-05,3.12415e-09,-2.22353e-13,32100.8,-53.3401], Tmin=(1039.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_4_5_ene_1) + radical(cyclobutane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1[CH]C2CC2C=C1(3501)',
    structure = SMILES('C=C([O])C1[CH]C2CC2C=C1'),
    E0 = (295.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.294212,0.0652948,-1.16464e-06,-4.0288e-08,1.84621e-11,35681.6,29.3841], Tmin=(100,'K'), Tmax=(1031.58,'K')), NASAPolynomial(coeffs=[16.2289,0.0394557,-1.58643e-05,2.99269e-09,-2.13512e-13,30481.3,-57.2578], Tmin=(1031.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_3_6_ene_1) + radical(C=C(C)OJ) + radical(Cs_S)"""),
)

species(
    label = 'O=C1CC2[CH][CH]C1C=CC2(3137)',
    structure = SMILES('O=C1CC2[CH][CH]C1C=CC2'),
    E0 = (202.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102498,0.0916567,-7.43453e-05,3.28651e-08,-6.27969e-12,24432.2,1.80633], Tmin=(100,'K'), Tmax=(1182.45,'K')), NASAPolynomial(coeffs=[11.7213,0.0523525,-2.44858e-05,4.75431e-09,-3.36357e-13,21684.5,-56.195], Tmin=(1182.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclohexane) + ring(Cycloheptene) + ring(Cyclohexane) + radical(Cs_S) + radical(CCJCC=O)"""),
)

species(
    label = 'C=C1OC2C=CC[CH][CH]C12(3388)',
    structure = SMILES('C=C1OC2C=CC[CH][CH]C12'),
    E0 = (249.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.817974,0.0894305,-4.73914e-05,-4.9687e-09,8.45057e-12,30136.2,4.93089], Tmin=(100,'K'), Tmax=(1053.21,'K')), NASAPolynomial(coeffs=[22.1461,0.0341165,-1.40473e-05,2.68483e-09,-1.92889e-13,23529.7,-115.448], Tmin=(1053.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(Cycloheptene) + ring(Oxetane) + radical(Cs_S) + radical(RCCJCC)"""),
)

species(
    label = '[CH]=CC(C=CC=C)C(=C)[O](3502)',
    structure = SMILES('[CH]=CC(C=CC=C)C(=C)[O]'),
    E0 = (387.07,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,3120,650,792.5,1650,350,440,435,1725,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.627202,'amu*angstrom^2'), symmetry=1, barrier=(14.4206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627229,'amu*angstrom^2'), symmetry=1, barrier=(14.4212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627022,'amu*angstrom^2'), symmetry=1, barrier=(14.4165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.626946,'amu*angstrom^2'), symmetry=1, barrier=(14.4147,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.15418,0.103631,-9.96427e-05,4.9657e-08,-9.7906e-12,46747.5,38.7428], Tmin=(100,'K'), Tmax=(1234.35,'K')), NASAPolynomial(coeffs=[21.2688,0.0309678,-1.13405e-05,1.96513e-09,-1.31247e-13,41212,-74.1561], Tmin=(1234.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C1C=CC(C=C1)C(=C)[O](3503)',
    structure = SMILES('[CH2]C1C=CC(C=C1)C(=C)[O]'),
    E0 = (264.782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,350,440,435,1725,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.392956,0.0863102,-5.98776e-05,1.70615e-08,-2.4134e-13,32012.9,32.782], Tmin=(100,'K'), Tmax=(1057.19,'K')), NASAPolynomial(coeffs=[16.5604,0.0381142,-1.4124e-05,2.47949e-09,-1.67675e-13,27537.1,-54.1664], Tmin=(1057.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cyclohexadiene) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = 'C=C1OC2C=CC1C=CC2(2248)',
    structure = SMILES('C=C1OC2C=CC1C=CC2'),
    E0 = (-9.69085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.54541,0.0713639,8.04295e-05,-1.97114e-07,9.0777e-11,-919.337,3.8196], Tmin=(100,'K'), Tmax=(939.958,'K')), NASAPolynomial(coeffs=[48.5647,-0.0132544,1.02017e-05,-1.72147e-09,8.81429e-14,-16021.8,-265.056], Tmin=(939.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.69085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(Cycloheptane) + ring(36dihydro2hpyran)"""),
)

species(
    label = 'O=C1CC2[CH]CC=CC1[CH]2(3058)',
    structure = SMILES('O=C1CC2[CH]CC=CC1[CH]2'),
    E0 = (202.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102498,0.0916567,-7.43453e-05,3.28651e-08,-6.27969e-12,24432.2,1.80633], Tmin=(100,'K'), Tmax=(1182.45,'K')), NASAPolynomial(coeffs=[11.7213,0.0523525,-2.44858e-05,4.75431e-09,-3.36357e-13,21684.5,-56.195], Tmin=(1182.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclopentane) + ring(Cycloheptene) + ring(Cyclopentane) + radical(CCJCC=O) + radical(Cs_S)"""),
)

species(
    label = 'CC(=O)[C]1C=C[CH]CC=C1(3504)',
    structure = SMILES('CC(=O)[C]1C=C[CH]CC=C1'),
    E0 = (152.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.97141,0.00554079,0.000250824,-3.6483e-07,1.50895e-10,18463.7,29.151], Tmin=(100,'K'), Tmax=(914.501,'K')), NASAPolynomial(coeffs=[40.377,-0.00845818,1.40377e-05,-2.85853e-09,1.77269e-13,4634.52,-193.641], Tmin=(914.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJ(C=O)C=C)"""),
)

species(
    label = 'CC(=O)C1[C]=CC[CH]C=C1(3505)',
    structure = SMILES('CC(=O)C1[C]=CC[CH]C=C1'),
    E0 = (298.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579649,0.025093,0.000178814,-2.8168e-07,1.19515e-10,36077.5,32.3486], Tmin=(100,'K'), Tmax=(913.086,'K')), NASAPolynomial(coeffs=[36.756,-0.00288835,1.04015e-05,-2.19427e-09,1.37251e-13,24031.1,-168.681], Tmin=(913.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)C1[C]=C[CH]CC=C1(3506)',
    structure = SMILES('CC(=O)C1[C]=C[CH]CC=C1'),
    E0 = (298.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579649,0.025093,0.000178814,-2.8168e-07,1.19515e-10,36077.5,32.3486], Tmin=(100,'K'), Tmax=(913.086,'K')), NASAPolynomial(coeffs=[36.756,-0.00288835,1.04015e-05,-2.19427e-09,1.37251e-13,24031.1,-168.681], Tmin=(913.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)C1C=[C]C[CH]C=C1(3507)',
    structure = SMILES('CC(=O)C1C=[C]C[CH]C=C1'),
    E0 = (298.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579649,0.025093,0.000178814,-2.8168e-07,1.19515e-10,36077.5,32.3486], Tmin=(100,'K'), Tmax=(913.086,'K')), NASAPolynomial(coeffs=[36.756,-0.00288835,1.04015e-05,-2.19427e-09,1.37251e-13,24031.1,-168.681], Tmin=(913.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = 'CC(=O)C1C=[C][CH]CC=C1(3508)',
    structure = SMILES('CC(=O)C1C=[C][CH]CC=C1'),
    E0 = (298.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579649,0.025093,0.000178814,-2.8168e-07,1.19515e-10,36077.5,32.3486], Tmin=(100,'K'), Tmax=(913.086,'K')), NASAPolynomial(coeffs=[36.756,-0.00288835,1.04015e-05,-2.19427e-09,1.37251e-13,24031.1,-168.681], Tmin=(913.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(Allyl_S)"""),
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
    label = 'O=C1CC2[CH]C1[CH]C=CC2(3112)',
    structure = SMILES('O=C1CC2[CH]C1C=C[CH]C2'),
    E0 = (148.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.287527,0.0924059,-7.07701e-05,2.78601e-08,-4.53602e-12,18026.3,-0.464518], Tmin=(100,'K'), Tmax=(1411.22,'K')), NASAPolynomial(coeffs=[16.6829,0.0443043,-1.96422e-05,3.70705e-09,-2.57253e-13,13236.5,-88.1827], Tmin=(1411.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclopentane) + ring(Cycloheptene) + ring(Cyclopentane) + radical(Allyl_S) + radical(CCJCC=O)"""),
)

species(
    label = 'O=C1CC2[CH]CC=C[CH]C12(3063)',
    structure = SMILES('O=C1CC2[CH]C[CH]C=CC12'),
    E0 = (143.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.179079,0.0908091,-7.11652e-05,3.02088e-08,-5.58393e-12,17358,-2.07353], Tmin=(100,'K'), Tmax=(1206.55,'K')), NASAPolynomial(coeffs=[11.2706,0.0540381,-2.54511e-05,4.94998e-09,-3.50248e-13,14681.5,-57.6663], Tmin=(1206.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.228,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(Cycloheptene) + ring(Cyclobutane) + radical(Allyl_S) + radical(Cs_S)"""),
)

species(
    label = 'CC(=O)C1C=C=CCC=C1(3509)',
    structure = SMILES('CC(=O)C1C=C=CCC=C1'),
    E0 = (144.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.24531,0.0725926,-3.57434e-05,4.9404e-09,6.69063e-13,17568.1,27.6729], Tmin=(100,'K'), Tmax=(1417.86,'K')), NASAPolynomial(coeffs=[15.756,0.0423171,-1.79776e-05,3.2937e-09,-2.23123e-13,11814.5,-57.3517], Tmin=(1417.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
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
    label = 'O=[C]CC1[CH]C=CCC=C1(2246)',
    structure = SMILES('O=[C]CC1C=C[CH]CC=C1'),
    E0 = (251.272,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.46551,0.0260414,0.000181888,-2.89927e-07,1.23965e-10,30396.7,32.5997], Tmin=(100,'K'), Tmax=(909.75,'K')), NASAPolynomial(coeffs=[38.5536,-0.00634623,1.25698e-05,-2.64223e-09,1.6908e-13,17876.7,-178.273], Tmin=(909.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C1CC2C=CC1C=CC2(2250)',
    structure = SMILES('O=C1CC2C=CC1C=CC2'),
    E0 = (-76.4957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.511618,0.10077,-9.94395e-05,5.4281e-08,-1.22191e-11,-9039.08,0.0667863], Tmin=(100,'K'), Tmax=(1060.33,'K')), NASAPolynomial(coeffs=[14.629,0.043654,-1.86403e-05,3.47997e-09,-2.4156e-13,-12249.9,-73.8649], Tmin=(1060.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.4957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclohexane) + ring(1,4-Cycloheptadiene) + ring(Cyclohexene)"""),
)

species(
    label = '[CH2]C1=CC=C[CH]CC=CO1(3510)',
    structure = SMILES('[CH2]C1=C[CH]C=CCC=CO1'),
    E0 = (288.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.488047,0.0779228,-8.25815e-06,-5.346e-08,2.84072e-11,34934.7,21.3983], Tmin=(100,'K'), Tmax=(960.081,'K')), NASAPolynomial(coeffs=[23.4503,0.0283678,-9.23428e-06,1.65724e-09,-1.20717e-13,28025.5,-105.16], Tmin=(960.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclononane) + radical(C=C(O)CJ) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH2]C1=CC=CC[CH]C=CO1(3511)',
    structure = SMILES('C=C1[CH]C=CCC=C[CH]O1'),
    E0 = (246.441,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.932401,0.0905196,-4.34574e-05,-1.27256e-08,1.19708e-11,29833.3,19.7116], Tmin=(100,'K'), Tmax=(1029.63,'K')), NASAPolynomial(coeffs=[23.375,0.0329641,-1.33312e-05,2.55276e-09,-1.84838e-13,22873,-107.76], Tmin=(1029.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.441,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclononane) + radical(C=CCJCO) + radical(C=CCJ(O)C)"""),
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
    label = 'O=[C]C1C=C[CH]CC=C1(3512)',
    structure = SMILES('O=[C]C1C=C[CH]CC=C1'),
    E0 = (274.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38017,0.00321245,0.000220651,-3.24751e-07,1.35785e-10,33132.4,28.5675], Tmin=(100,'K'), Tmax=(909.121,'K')), NASAPolynomial(coeffs=[37.6751,-0.0146882,1.62374e-05,-3.29605e-09,2.11617e-13,20673.5,-175.303], Tmin=(909.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'CC(=O)C1[CH]C=[C]CC=C1(3513)',
    structure = SMILES('CC(=O)C1[CH]C=[C]CC=C1'),
    E0 = (359.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621979,0.0214829,0.000192668,-3.00323e-07,1.27447e-10,43438.5,36.4952], Tmin=(100,'K'), Tmax=(910.786,'K')), NASAPolynomial(coeffs=[38.4946,-0.00710366,1.28958e-05,-2.68661e-09,1.70924e-13,30826.6,-174.043], Tmin=(910.786,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(CCJCC=O) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)C1C=CCC=CC=1(3514)',
    structure = SMILES('CC(=O)C1C=CCC=CC=1'),
    E0 = (6.39494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.761656,0.0864355,-3.60545e-05,-2.10874e-08,1.54745e-11,956.7,28.2842], Tmin=(100,'K'), Tmax=(1004.18,'K')), NASAPolynomial(coeffs=[23.3479,0.0302397,-1.16244e-05,2.20353e-09,-1.60286e-13,-5894.13,-98.1332], Tmin=(1004.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.39494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'O=C1CC2C=CCC=CC12(2252)',
    structure = SMILES('O=C1CC2C=CCC=CC12'),
    E0 = (-78.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.359431,0.0922855,-7.0779e-05,2.78227e-08,-4.49674e-12,-9223.97,-0.676679], Tmin=(100,'K'), Tmax=(1429.71,'K')), NASAPolynomial(coeffs=[17.4235,0.0425333,-1.85813e-05,3.48343e-09,-2.40802e-13,-14308.9,-92.8262], Tmin=(1429.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(1,4-Cycloheptadiene) + ring(Cyclobutane)"""),
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
    E0 = (226.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (421.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (268.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (381.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (542.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (330.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (315.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (390.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (480.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (464.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (383.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (527.829,'kJ/mol'),
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
    E0 = (370.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (367.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (418.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (376.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (359.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (478.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (307.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (503.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (493.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (457.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (471.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (357.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (351.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (505.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (319.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (299.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (275.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (304.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (251.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (378.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (234.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (784.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (391.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (387.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (485.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (477.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (375.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (644.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (482.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (357.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (300.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (339.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (275.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (428.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (244.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (424.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (233.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (295.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (358.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (342.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (342.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (331.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (338.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (307.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (505.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (293.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (260.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (244.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (496.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (234.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (602.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (560.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (655.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (521.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (304.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (234.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C7H8(11)', 'CH2CO(43)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['[CH2]C1([O])C2C=CCC=CC21(3477)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(194.673,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C1OC2[CH]CC=CC1[CH]2(3401)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.2, Ea=(42.5884,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R6;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'C=C([O])C1C=CCC=CC=1(3478)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34e+08,'cm^3/(mol*s)'), n=1.64, Ea=(9.6232,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2609 used for Cds-CdCd_Cds-CdH;HJ
Exact match found for rate rule [Cds-CdCd_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(19)', 'C=C([O])C1C=C=CCC=C1(3479)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C7H8(11)', '[CH2][C]=O(360)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00666416,'m^3/(mol*s)'), n=2.62825, Ea=(4.69677,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', 'CH2CO(43)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7051.36,'m^3/(mol*s)'), n=1.07323, Ea=(43.3386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C([O])[C]1C=CCC=CC1(3424)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.95861e+10,'s^-1'), n=0.595417, Ea=(163.653,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C1C=CCC=[C]C1(3420)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C([O])C1[C]=CCC=CC1(3421)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C(O)[C]1[CH]C=CCC=C1(3480)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.19923e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(O)C1[CH]C=CCC=C1(3481)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C([O])C1C=CC[C]=CC1(3416)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 200 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C(O)C1[C]=CCC=C[CH]1(3482)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C([O])C1C=[C]CC=CC1(3417)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(285077,'s^-1'), n=1.88833, Ea=(191.836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] + [R4H_SDS;C_rad_out_single;Cs_H_out_H/Cd] + [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C([O])C1C=CCC=CC1(3427)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C(O)C1[CH]C=CC[C]=C1(3483)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C(O)C1[CH][C]=CCC=C1(3484)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C(O)C1[CH]C=C[CH]C=C1(3422)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(34428.8,'s^-1'), n=1.95377, Ea=(81.3736,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R6H_RS(Cd)SMS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C(O)C1[CH]C=[C]CC=C1(3485)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.44313e+09,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;XH_out] for rate rule [R6HJ_2;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', '[CH2][C]=O(360)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['[CH]1C=CCC=CC1[C]1CO1(3486)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C([O])C1C=CCC2[CH]C21(3487)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_cs] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C([O])C1[CH]C2CC=CC21(3488)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn0c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['O=C1CC2C=CC[CH][CH]C12(3090)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.94431e+07,'s^-1'), n=1.16299, Ea=(125.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs] for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C1OC2[CH]C1[CH]C=CC2(3311)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(8.11e+10,'s^-1'), n=0.18, Ea=(278.654,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C1OC2[CH][CH]C1C=CC2(3343)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.42658e+11,'s^-1'), n=0.0609759, Ea=(92.8718,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C([O])C1C2[CH]CC=CC21(3489)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(72.9485,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C1OC2[CH]CC=C[CH]C12(3406)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=CCC=CC=CC(=C)[O](3490)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.58366e+14,'s^-1'), n=-0.170131, Ea=(68.698,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;doublebond_intra_pri;radadd_intra_cdsingleH] + [R7;doublebond_intra_pri_HCd;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_HCd;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C(O)C1C=CCC=CC=1(3491)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C(O)C1C=C=CCC=C1(3492)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C([O])[CH]C1C=CCC=C1(3493)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7.65366e+10,'s^-1'), n=0.611, Ea=(152.377,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ-CdH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C1OC2C=CCC=CC12(3380)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(20)', 'C=[C]C1[CH]C=CCC=C1(3494)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['[CH2]C1([O])C2C=CCC1C=C2(3413)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(5.249e+08,'s^-1'), n=0.846, Ea=(165.135,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS_D;doublebond_intra_2H;radadd_intra_cs] for rate rule [R6_SMS_D;doublebond_intra_2H;radadd_intra_csHNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(19)', 'C=C([O])C1C=CC=CC=C1(3463)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(9.22566,'m^3/(mol*s)'), n=2.04274, Ea=(10.5229,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 101 used for Cds-CdH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CdH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=C([O])C1C=[C]CCC=C1(3495)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C([O])C1[C]=CCCC=C1(3496)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 200 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C(O)C1[C]=C[CH]CC=C1(3497)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C([O])[C]1C=CCCC=C1(3498)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=C([O])C1C=CCCC=C1(3499)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(3.692e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R6H_RSSMS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C([O])C1C=CC2C[CH]C21(3500)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn0c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C([O])C1[CH]C2CC2C=C1(3501)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(73.8271,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['O=C1CC2[CH][CH]C1C=CC2(3137)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(3.68243e+09,'s^-1'), n=0.4695, Ea=(113.018,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SDS_D;doublebond_intra;radadd_intra_cs] for rate rule [R6_SDS_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C1OC2C=CC[CH][CH]C12(3388)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=CC(C=CC=C)C(=C)[O](3502)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C(O)C1C=CC=CC=C1(3471)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_De] for rate rule [R6radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C1C=CC(C=C1)C(=C)[O](3503)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.31121e+11,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['C=C1OC2C=CC1C=CC2(2248)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction52',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['O=C1CC2[CH]CC=CC1[CH]2(3058)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.58148e+08,'s^-1'), n=0.696039, Ea=(69.5697,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra;radadd_intra_cs2H] for rate rule [R9;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['CC(=O)[C]1C=C[CH]CC=C1(3504)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(25964.4,'s^-1'), n=2.44333, Ea=(131.615,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['CC(=O)C1[C]=CC[CH]C=C1(3505)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['CC(=O)C1[C]=C[CH]CC=C1(3506)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['CC(=O)C1C=[C]C[CH]C=C1(3507)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['CC(=O)C1C=[C][CH]CC=C1(3508)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(263079,'s^-1'), n=1.73643, Ea=(39.8993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['CC(=O)C1[CH]C=C[CH]C=C1(3439)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(25169.2,'s^-1'), n=1.925, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6H_RSSMS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['O=C1CC2[CH]C1[CH]C=CC2(3112)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(8.11e+10,'s^-1'), n=0.18, Ea=(278.654,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H] for rate rule [Rn2c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction60',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['O=C1CC2[CH]CC=C[CH]C12(3063)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(2.54811e+07,'s^-1'), n=0.785425, Ea=(66.8554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction61',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['CC(=O)C1C=C=CCC=C1(3509)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction62',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['CC(=O)C1C=CC=CC=C1(3475)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_De] for rate rule [R6radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O=[C]CC1[CH]C=CCC=C1(2246)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction64',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['O=C1CC2C=CC1C=CC2(2250)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;C_rad_out_2H;Cpri_rad_out_single] + [R6_SSSDS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R6_SSSDS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH2]C1=CC=C[CH]CC=CO1(3510)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2]C1=CC=CC[CH]C=CO1(3511)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction67',
    reactants = ['CH2(T)(33)', 'O=[C]C1C=C[CH]CC=C1(3512)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction68',
    reactants = ['CC(=O)C1[CH]C=[C]CC=C1(3513)'],
    products = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction69',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['CC(=O)C1C=CCC=CC=1(3514)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction70',
    reactants = ['C=C([O])C1[CH]C=CCC=C1(2240)'],
    products = ['O=C1CC2C=CCC=CC12(2252)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

network(
    label = '493',
    isomers = [
        'C=C([O])C1[CH]C=CCC=C1(2240)',
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
    label = '493',
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

