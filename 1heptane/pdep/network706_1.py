species(
    label = '[CH2]C(CC)C(=C)[O](1127)',
    structure = SMILES('[CH2]C(CC)C(=C)[O]'),
    E0 = (56.2744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,474.866,474.906,475.039],'cm^-1')),
        HinderedRotor(inertia=(0.144122,'amu*angstrom^2'), symmetry=1, barrier=(3.31364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0874971,'amu*angstrom^2'), symmetry=1, barrier=(14.0659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0873058,'amu*angstrom^2'), symmetry=1, barrier=(14.0655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0030326,'amu*angstrom^2'), symmetry=1, barrier=(14.0665,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.626666,0.0650109,-3.33787e-05,-5.40408e-09,7.94584e-12,6897.98,29.7691], Tmin=(100,'K'), Tmax=(953.666,'K')), NASAPolynomial(coeffs=[14.5328,0.0282512,-9.48282e-06,1.60529e-09,-1.08168e-13,3264.87,-41.8022], Tmin=(953.666,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.2744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = 'butene1(127)',
    structure = SMILES('C=CCC'),
    E0 = (-16.4325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,252.555],'cm^-1')),
        HinderedRotor(inertia=(0.178654,'amu*angstrom^2'), symmetry=1, barrier=(7.72883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185589,'amu*angstrom^2'), symmetry=1, barrier=(7.72103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58773,0.0232778,1.93412e-05,-3.55496e-08,1.36906e-11,-1918.73,14.5751], Tmin=(100,'K'), Tmax=(1007.28,'K')), NASAPolynomial(coeffs=[7.20517,0.0236362,-9.0315e-06,1.65393e-09,-1.16019e-13,-3797.34,-12.4426], Tmin=(1007.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.4325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene1""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2CO(28)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]C1([O])CC1CC(3006)',
    structure = SMILES('[CH2]C1([O])CC1CC'),
    E0 = (212.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.499535,0.0642237,-2.28759e-05,-1.81153e-08,1.20407e-11,25655.3,25.2632], Tmin=(100,'K'), Tmax=(995.023,'K')), NASAPolynomial(coeffs=[16.6167,0.0276219,-1.01939e-05,1.85965e-09,-1.31573e-13,21052.4,-59.4244], Tmin=(995.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
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
    label = 'C=C([O])C(=C)CC(3027)',
    structure = SMILES('C=C([O])C(=C)CC'),
    E0 = (-42.403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,285.652,285.652,285.652],'cm^-1')),
        HinderedRotor(inertia=(0.236387,'amu*angstrom^2'), symmetry=1, barrier=(13.6876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236388,'amu*angstrom^2'), symmetry=1, barrier=(13.6876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236388,'amu*angstrom^2'), symmetry=1, barrier=(13.6876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.463697,0.0681644,-4.34914e-05,3.06498e-09,5.28787e-12,-4963.98,24.9441], Tmin=(100,'K'), Tmax=(971.207,'K')), NASAPolynomial(coeffs=[16.3403,0.0240564,-8.23576e-06,1.426e-09,-9.78444e-14,-9051.53,-56.3544], Tmin=(971.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C2H5(29)',
    structure = SMILES('C[CH2]'),
    E0 = (107.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1190.6,1642.82,1642.96,3622.23,3622.39],'cm^-1')),
        HinderedRotor(inertia=(0.866817,'amu*angstrom^2'), symmetry=1, barrier=(19.9298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24186,-0.00356905,4.82667e-05,-5.85401e-08,2.25805e-11,12969,4.44704], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.32196,0.0123931,-4.39681e-06,7.0352e-10,-4.18435e-14,12175.9,0.171104], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(107.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H5""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=CC(=C)[O](2398)',
    structure = SMILES('C=CC(=C)[O]'),
    E0 = (17.8331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.914251,'amu*angstrom^2'), symmetry=1, barrier=(21.0204,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98668,0.0359126,-9.09545e-06,-2.04235e-08,1.25939e-11,2225.07,15.9658], Tmin=(100,'K'), Tmax=(920.927,'K')), NASAPolynomial(coeffs=[13.0167,0.00983409,-2.17548e-06,3.0681e-10,-2.11581e-14,-732.214,-41.3653], Tmin=(920.927,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.8331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C][O](173)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97566e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02476e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CsCJ=O) + radical(CJC=O)"""),
)

species(
    label = '[CH2][CH]CC(130)',
    structure = SMILES('[CH2][CH]CC'),
    E0 = (255.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2031.24],'cm^-1')),
        HinderedRotor(inertia=(0.244974,'amu*angstrom^2'), symmetry=1, barrier=(5.63244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00192352,'amu*angstrom^2'), symmetry=1, barrier=(5.63177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244928,'amu*angstrom^2'), symmetry=1, barrier=(5.63137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98997,0.0287412,-9.51469e-06,4.19232e-10,1.90526e-13,30780.1,16.8971], Tmin=(100,'K'), Tmax=(2154.56,'K')), NASAPolynomial(coeffs=[12.4231,0.0182241,-7.06316e-06,1.16769e-09,-7.11818e-14,25091.4,-39.6212], Tmin=(2154.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])[C](C)CC(3028)',
    structure = SMILES('[CH2]C([O])=C(C)CC'),
    E0 = (-9.43951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,325,375,415,465,420,450,1700,1750,312.245,316.31],'cm^-1')),
        HinderedRotor(inertia=(0.00174689,'amu*angstrom^2'), symmetry=1, barrier=(0.119926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114303,'amu*angstrom^2'), symmetry=1, barrier=(7.97227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112806,'amu*angstrom^2'), symmetry=1, barrier=(7.91707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109505,'amu*angstrom^2'), symmetry=1, barrier=(7.94687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.28084,0.0735632,-6.17138e-05,2.79013e-08,-5.09813e-12,-994.535,28.2407], Tmin=(100,'K'), Tmax=(1312.98,'K')), NASAPolynomial(coeffs=[15.0975,0.028424,-1.01446e-05,1.71688e-09,-1.12401e-13,-4885.29,-47.2754], Tmin=(1312.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.43951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C(C)[CH]C(3029)',
    structure = SMILES('C=C([O])C(C)[CH]C'),
    E0 = (45.7342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.722008,0.0629665,-3.49956e-05,3.37847e-09,2.54455e-12,5626.3,30.2045], Tmin=(100,'K'), Tmax=(1084.24,'K')), NASAPolynomial(coeffs=[13.4125,0.0312883,-1.2115e-05,2.1883e-09,-1.5047e-13,1984.52,-36.1496], Tmin=(1084.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.7342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S)"""),
)

species(
    label = '[CH2][C](CC)C(=C)O(3030)',
    structure = SMILES('[CH2]C(O)=C([CH2])CC'),
    E0 = (4.25487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.106981,0.0728252,-3.81657e-05,-1.18142e-08,1.26968e-11,663.585,27.6411], Tmin=(100,'K'), Tmax=(943.2,'K')), NASAPolynomial(coeffs=[19.6146,0.0214947,-6.46768e-06,1.0755e-09,-7.45944e-14,-4412.98,-72.7347], Tmin=(943.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.25487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(O)C([CH2])CC(1352)',
    structure = SMILES('[CH]=C(O)C([CH2])CC'),
    E0 = (165.566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,366.327],'cm^-1')),
        HinderedRotor(inertia=(0.00125627,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148565,'amu*angstrom^2'), symmetry=1, barrier=(14.1468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148547,'amu*angstrom^2'), symmetry=1, barrier=(14.1469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148563,'amu*angstrom^2'), symmetry=1, barrier=(14.1467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148562,'amu*angstrom^2'), symmetry=1, barrier=(14.1468,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.209793,0.072218,-4.14133e-05,-6.71457e-09,1.08318e-11,20059.7,30.1684], Tmin=(100,'K'), Tmax=(929.87,'K')), NASAPolynomial(coeffs=[18.4305,0.0221981,-6.47271e-06,1.03371e-09,-6.9441e-14,15445,-63.004], Tmin=(929.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]C)C(=C)O(3031)',
    structure = SMILES('[CH2]C([CH]C)C(=C)O'),
    E0 = (113.012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,426.041,426.047],'cm^-1')),
        HinderedRotor(inertia=(0.928781,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000928744,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118211,'amu*angstrom^2'), symmetry=1, barrier=(15.2261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00134104,'amu*angstrom^2'), symmetry=1, barrier=(15.2262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118221,'amu*angstrom^2'), symmetry=1, barrier=(15.2259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.456,0.067707,-3.70288e-05,-4.63305e-09,8.23312e-12,13729,31.4369], Tmin=(100,'K'), Tmax=(958.982,'K')), NASAPolynomial(coeffs=[16.2161,0.0257542,-8.60957e-06,1.47243e-09,-1.00578e-13,9612.63,-49.6384], Tmin=(958.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cs_S)"""),
)

species(
    label = '[CH2]CC(C)C(=C)[O](3032)',
    structure = SMILES('[CH2]CC(C)C(=C)[O]'),
    E0 = (56.4384,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,402.05,402.604,403.371],'cm^-1')),
        HinderedRotor(inertia=(0.00105728,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00103268,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0776163,'amu*angstrom^2'), symmetry=1, barrier=(9.02247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0774311,'amu*angstrom^2'), symmetry=1, barrier=(8.9854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61324,0.06444,-3.264e-05,-3.97159e-09,6.40779e-12,6918.72,29.8223], Tmin=(100,'K'), Tmax=(1006,'K')), NASAPolynomial(coeffs=[14.7265,0.0290099,-1.06567e-05,1.90102e-09,-1.31353e-13,3032.37,-43.5531], Tmin=(1006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.4384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C([O])C(C)CC(3033)',
    structure = SMILES('[CH]=C([O])C(C)CC'),
    E0 = (98.2882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,346.725,346.731],'cm^-1')),
        HinderedRotor(inertia=(0.00140223,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118928,'amu*angstrom^2'), symmetry=1, barrier=(10.1457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118925,'amu*angstrom^2'), symmetry=1, barrier=(10.1457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118926,'amu*angstrom^2'), symmetry=1, barrier=(10.1456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.498044,0.0672808,-3.90374e-05,1.42189e-09,4.80918e-12,11956,28.8523], Tmin=(100,'K'), Tmax=(1005.28,'K')), NASAPolynomial(coeffs=[15.1967,0.028431,-1.03678e-05,1.83934e-09,-1.26637e-13,8008.52,-47.0726], Tmin=(1005.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.2882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC([CH2])C(=C)O(3034)',
    structure = SMILES('[CH2]CC([CH2])C(=C)O'),
    E0 = (123.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,302.782,303.537],'cm^-1')),
        HinderedRotor(inertia=(0.00183093,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00182208,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00182547,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289806,'amu*angstrom^2'), symmetry=1, barrier=(18.9611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290091,'amu*angstrom^2'), symmetry=1, barrier=(18.9565,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.320745,0.0694321,-3.52334e-05,-1.17855e-08,1.22728e-11,15022.6,31.1534], Tmin=(100,'K'), Tmax=(933.388,'K')), NASAPolynomial(coeffs=[17.961,0.0227747,-6.75982e-06,1.09486e-09,-7.41088e-14,10468.9,-59.4878], Tmin=(933.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(CC)[C]1CO1(3035)',
    structure = SMILES('[CH2]C(CC)[C]1CO1'),
    E0 = (204.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.361062,0.0712031,-6.15402e-05,3.17709e-08,-6.5832e-12,24760.9,29.1601], Tmin=(100,'K'), Tmax=(1321.86,'K')), NASAPolynomial(coeffs=[11.4774,0.030902,-8.24748e-06,1.08024e-09,-5.76237e-14,22404.1,-25.3701], Tmin=(1321.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(Isobutyl)"""),
)

species(
    label = 'CCC1CC[C]1[O](3036)',
    structure = SMILES('CCC1CC[C]1[O]'),
    E0 = (191.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12429,0.0487826,1.30712e-05,-4.80321e-08,2.06891e-11,23183.2,25.2939], Tmin=(100,'K'), Tmax=(1004.34,'K')), NASAPolynomial(coeffs=[13.7691,0.03214,-1.24313e-05,2.32344e-09,-1.66149e-13,18942.7,-44.2307], Tmin=(1004.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C(O)C(=C)CC(3037)',
    structure = SMILES('C=C(O)C(=C)CC'),
    E0 = (-180.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.118881,0.0718772,-3.5012e-05,-1.53661e-08,1.39386e-11,-21522,24.6812], Tmin=(100,'K'), Tmax=(946.891,'K')), NASAPolynomial(coeffs=[19.9693,0.0208828,-6.28584e-06,1.05923e-09,-7.44793e-14,-26754.4,-77.7811], Tmin=(946.891,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-180.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CH2(S)(23)',
    structure = SMILES('[CH2]'),
    E0 = (419.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.36,2789.41,2993.36],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19195,-0.00230793,8.0509e-06,-6.60123e-09,1.95638e-12,50484.3,-0.754589], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.28556,0.00460255,-1.97412e-06,4.09548e-10,-3.34695e-14,50922.4,8.67684], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(419.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]C(C)C(=C)[O](1151)',
    structure = SMILES('[CH2]C(C)C(=C)[O]'),
    E0 = (80.0546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,348.303,348.547],'cm^-1')),
        HinderedRotor(inertia=(0.00138873,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00138742,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127038,'amu*angstrom^2'), symmetry=1, barrier=(10.9498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3675.75,'J/mol'), sigma=(6.3037,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=574.14 K, Pc=33.3 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36686,0.0490894,-1.53711e-05,-1.91732e-08,1.26765e-11,9731.3,24.8689], Tmin=(100,'K'), Tmax=(922.496,'K')), NASAPolynomial(coeffs=[13.2401,0.0204788,-6.04096e-06,9.61536e-10,-6.39148e-14,6567.47,-36.7297], Tmin=(922.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.0546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])[CH]CCC(2856)',
    structure = SMILES('[CH2]C([O])=CCCC'),
    E0 = (5.83526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,373.018,373.028,373.332],'cm^-1')),
        HinderedRotor(inertia=(0.00121155,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0813913,'amu*angstrom^2'), symmetry=1, barrier=(8.04061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0814427,'amu*angstrom^2'), symmetry=1, barrier=(8.04097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0813616,'amu*angstrom^2'), symmetry=1, barrier=(8.04126,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.525374,0.0678657,-4.37604e-05,7.97548e-09,2.2261e-12,834.377,28.8471], Tmin=(100,'K'), Tmax=(1013.82,'K')), NASAPolynomial(coeffs=[14.3435,0.0294284,-1.06846e-05,1.87219e-09,-1.27258e-13,-2793.92,-42.0832], Tmin=(1013.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.83526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C[CH]CC(1126)',
    structure = SMILES('C=C([O])C[CH]CC'),
    E0 = (53.6223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,180,180,1286.34,2707.76],'cm^-1')),
        HinderedRotor(inertia=(0.42994,'amu*angstrom^2'), symmetry=1, barrier=(13.8863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.601606,'amu*angstrom^2'), symmetry=1, barrier=(13.8321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0117751,'amu*angstrom^2'), symmetry=1, barrier=(13.9572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00266568,'amu*angstrom^2'), symmetry=1, barrier=(13.8472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3795.8,'J/mol'), sigma=(6.58551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.90 K, Pc=30.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1185,0.0634505,-4.53405e-05,1.82342e-08,-3.18071e-12,6552.72,29.0296], Tmin=(100,'K'), Tmax=(1288.34,'K')), NASAPolynomial(coeffs=[9.34835,0.0378989,-1.55913e-05,2.84021e-09,-1.9357e-13,4432.14,-12.7598], Tmin=(1288.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.6223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJCC)"""),
)

species(
    label = 'C=C1OCC1CC(3005)',
    structure = SMILES('C=C1OCC1CC'),
    E0 = (-110.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07693,0.0416437,5.49148e-05,-1.11919e-07,4.98343e-11,-13196.2,19.652], Tmin=(100,'K'), Tmax=(915.421,'K')), NASAPolynomial(coeffs=[20.0698,0.017643,-2.4188e-06,2.29639e-10,-1.79585e-14,-19145.2,-83.7995], Tmin=(915.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane)"""),
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
    label = 'C=C([O])[CH]CC(2377)',
    structure = SMILES('[CH2]C([O])=CCC'),
    E0 = (29.6155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,331.592,331.647],'cm^-1')),
        HinderedRotor(inertia=(0.00153279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110666,'amu*angstrom^2'), symmetry=1, barrier=(8.63422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11059,'amu*angstrom^2'), symmetry=1, barrier=(8.63483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2878,0.0517027,-2.50274e-05,-6.51759e-09,7.14957e-12,3666.72,23.8661], Tmin=(100,'K'), Tmax=(960.293,'K')), NASAPolynomial(coeffs=[12.8498,0.0219906,-7.43267e-06,1.27278e-09,-8.66509e-14,595.546,-35.8741], Tmin=(960.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.6155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C([C]=C)CC(3038)',
    structure = SMILES('[CH2]C([C]=C)CC'),
    E0 = (370.87,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,3169.9],'cm^-1')),
        HinderedRotor(inertia=(0.108848,'amu*angstrom^2'), symmetry=1, barrier=(2.50263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.61763,'amu*angstrom^2'), symmetry=1, barrier=(14.2005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0317755,'amu*angstrom^2'), symmetry=1, barrier=(14.2121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.5605,'amu*angstrom^2'), symmetry=1, barrier=(81.8629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15195,0.0557142,-3.08427e-05,5.12439e-09,1.11172e-12,44713.7,26.8011], Tmin=(100,'K'), Tmax=(1095.25,'K')), NASAPolynomial(coeffs=[10.5722,0.0315537,-1.17827e-05,2.06206e-09,-1.3843e-13,42035.8,-22.3083], Tmin=(1095.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](CC)C(C)=O(3039)',
    structure = SMILES('[CH2]C(CC)=C(C)[O]'),
    E0 = (-16.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.444542,0.0690553,-4.5e-05,9.21829e-09,1.5619e-12,-1891.55,27.2418], Tmin=(100,'K'), Tmax=(1046.03,'K')), NASAPolynomial(coeffs=[14.8174,0.0297777,-1.11668e-05,1.98939e-09,-1.36239e-13,-5756.48,-46.8472], Tmin=(1046.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C([CH]C)C(C)=O(3040)',
    structure = SMILES('[CH2]C([CH]C)C(C)=O'),
    E0 = (98.8978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0123,0.0687034,-5.64125e-05,2.6669e-08,-5.49419e-12,11999.7,28.912], Tmin=(100,'K'), Tmax=(1109.53,'K')), NASAPolynomial(coeffs=[8.98893,0.0399464,-1.7535e-05,3.3091e-09,-2.3069e-13,10229.6,-10.3997], Tmin=(1109.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.8978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]CC([CH2])C(C)=O(3041)',
    structure = SMILES('[CH2]CC([CH2])C(C)=O'),
    E0 = (104.242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,298.639,2005.93],'cm^-1')),
        HinderedRotor(inertia=(0.00238394,'amu*angstrom^2'), symmetry=1, barrier=(0.150538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0045427,'amu*angstrom^2'), symmetry=1, barrier=(0.28693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139076,'amu*angstrom^2'), symmetry=1, barrier=(8.89501,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143963,'amu*angstrom^2'), symmetry=1, barrier=(8.89497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140362,'amu*angstrom^2'), symmetry=1, barrier=(8.89398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82945,0.0747828,-7.52386e-05,4.72698e-08,-1.3021e-11,12647.2,28.8701], Tmin=(100,'K'), Tmax=(853.191,'K')), NASAPolynomial(coeffs=[7.59104,0.0430827,-1.95066e-05,3.72202e-09,-2.60752e-13,11493.4,-2.6771], Tmin=(853.191,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(RCCJ) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2][C]1OCC1CC(3042)',
    structure = SMILES('[CH2][C]1OCC1CC'),
    E0 = (201.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916238,0.066242,-4.29833e-05,4.34387e-09,6.86601e-12,24392.1,24.3641], Tmin=(100,'K'), Tmax=(786.994,'K')), NASAPolynomial(coeffs=[10.2445,0.0340793,-1.07475e-05,1.65842e-09,-1.02367e-13,22451.6,-21.4055], Tmin=(786.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=C(CC)C(C)=O(3043)',
    structure = SMILES('C=C(CC)C(C)=O'),
    E0 = (-192.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20402,0.0572791,-3.01321e-05,6.80353e-09,-5.29248e-13,-23046,25.7523], Tmin=(100,'K'), Tmax=(1895.04,'K')), NASAPolynomial(coeffs=[18.1059,0.027833,-1.17556e-05,2.07354e-09,-1.34116e-13,-30570.6,-69.5455], Tmin=(1895.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(CC)C[C]=O(1131)',
    structure = SMILES('[CH2]C(CC)C[C]=O'),
    E0 = (84.6847,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,319,1667.44],'cm^-1')),
        HinderedRotor(inertia=(0.00162915,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115981,'amu*angstrom^2'), symmetry=1, barrier=(8.56704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113744,'amu*angstrom^2'), symmetry=1, barrier=(8.54097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116247,'amu*angstrom^2'), symmetry=1, barrier=(8.56661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118006,'amu*angstrom^2'), symmetry=1, barrier=(8.52988,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.918469,0.0692796,-6.0227e-05,3.15738e-08,-7.15934e-12,10294.9,29.5123], Tmin=(100,'K'), Tmax=(1029.78,'K')), NASAPolynomial(coeffs=[8.74736,0.0388699,-1.59319e-05,2.89794e-09,-1.97739e-13,8682.49,-8.48737], Tmin=(1029.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.6847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(Isobutyl)"""),
)

species(
    label = 'CCC1CCC1=O(1132)',
    structure = SMILES('CCC1CCC1=O'),
    E0 = (-169.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59162,0.0385628,3.2717e-05,-6.18803e-08,2.40605e-11,-20256.8,23.7619], Tmin=(100,'K'), Tmax=(1013.01,'K')), NASAPolynomial(coeffs=[11.1232,0.0355431,-1.40706e-05,2.6448e-09,-1.88887e-13,-23964.1,-31.1124], Tmin=(1013.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone)"""),
)

species(
    label = '[CH2]C=C([CH2])OCC(3044)',
    structure = SMILES('[CH2]C=C([CH2])OCC'),
    E0 = (47.0917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,350,440,435,1725,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202883,0.070435,-3.32129e-05,-1.37502e-08,1.22472e-11,5812.31,27.0927], Tmin=(100,'K'), Tmax=(965.906,'K')), NASAPolynomial(coeffs=[18.7757,0.0238765,-8.05009e-06,1.41841e-09,-9.99588e-14,808.352,-69.1963], Tmin=(965.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.0917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=O)CC(2680)',
    structure = SMILES('[CH2]C([C]=O)CC'),
    E0 = (112.555,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,301.166],'cm^-1')),
        HinderedRotor(inertia=(0.00184205,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132548,'amu*angstrom^2'), symmetry=1, barrier=(8.73205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135339,'amu*angstrom^2'), symmetry=1, barrier=(8.73252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133137,'amu*angstrom^2'), symmetry=1, barrier=(8.73985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50737,0.057543,-5.11715e-05,2.66539e-08,-6.01188e-12,13624.7,24.2058], Tmin=(100,'K'), Tmax=(1029.09,'K')), NASAPolynomial(coeffs=[8.04044,0.0321494,-1.41578e-05,2.67556e-09,-1.86745e-13,12280.1,-7.49981], Tmin=(1029.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CC(C)CJ=O)"""),
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
    E0 = (56.2744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (212.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (169.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (152.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (163.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (238.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (193.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (203.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (213.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (357.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (197.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (139.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (142.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (166.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (415.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (287.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (191.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (119.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (499.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (216.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (216.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (64.5587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (411.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (613.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (173.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (142.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (156.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (202.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (119.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (330.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (64.5587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (360.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (494.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['butene1(127)', 'CH2CO(28)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['[CH2]C1([O])CC1CC(3006)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(155.929,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C([O])C(=C)CC(3027)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.88727,'m^3/(mol*s)'), n=1.9687, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 18 used for Cds-CdCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -8.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C2H5(29)', 'C=CC(=C)[O](2398)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1340,'cm^3/(mol*s)'), n=2.41, Ea=(27.1542,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 814 used for Cds-CdH_Cds-HH;CsJ-CsHH
Exact match found for rate rule [Cds-CdH_Cds-HH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['butene1(127)', 'C=[C][O](173)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]CC(130)', 'CH2CO(28)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3525.68,'m^3/(mol*s)'), n=1.07323, Ea=(43.3386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([O])[C](C)CC(3028)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['C=C([O])C(C)[CH]C(3029)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['[CH2][C](CC)C(=C)O(3030)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.19923e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(O)C([CH2])CC(1352)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH]C)C(=C)O(3031)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;O_H_out] for rate rule [R4H_SS(Cd)S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]CC(C)C(=C)[O](3032)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C([O])C(C)CC(3033)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CC([CH2])C(=C)O(3034)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);C_rad_out_2H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]CC(130)', 'C=[C][O](173)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['[CH2]C(CC)[C]1CO1(3035)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['CCC1CC[C]1[O](3036)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(135.515,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 133.4 to 135.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['C=C(O)C(=C)CC(3037)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH2(S)(23)', '[CH2]C(C)C(=C)[O](1151)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['C=C([O])[CH]CCC(2856)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['C=C([O])C[CH]CC(1126)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['C=C1OCC1CC(3005)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(19)', 'C=C([O])[CH]CC(2377)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(4)', '[CH2]C([C]=C)CC(3038)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['[CH2][C](CC)C(C)=O(3039)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.30996,'s^-1'), n=3.5644, Ea=(117.059,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['[CH2]C([CH]C)C(C)=O(3040)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]CC([CH2])C(C)=O(3041)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CC(O2d)CC;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['[CH2][C]1OCC1CC(3042)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(146.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['C=C(CC)C(C)=O(3043)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(CC)C[C]=O(1131)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(CC)C(=C)[O](1127)'],
    products = ['CCC1CCC1=O(1132)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C=C([CH2])OCC(3044)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2(19)', '[CH2]C([C]=O)CC(2680)'],
    products = ['[CH2]C(CC)C(=C)[O](1127)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '706',
    isomers = [
        '[CH2]C(CC)C(=C)[O](1127)',
    ],
    reactants = [
        ('butene1(127)', 'CH2CO(28)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '706',
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

