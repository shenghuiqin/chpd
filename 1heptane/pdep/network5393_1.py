species(
    label = '[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)',
    structure = SMILES('[CH2]C(CC)C([CH2])(O)C(O)=CC'),
    E0 = (-143.895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,1583.42,1600,2911.86,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15167,'amu*angstrom^2'), symmetry=1, barrier=(3.4872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15167,'amu*angstrom^2'), symmetry=1, barrier=(3.4872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15167,'amu*angstrom^2'), symmetry=1, barrier=(3.4872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15167,'amu*angstrom^2'), symmetry=1, barrier=(3.4872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15167,'amu*angstrom^2'), symmetry=1, barrier=(3.4872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15167,'amu*angstrom^2'), symmetry=1, barrier=(3.4872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15167,'amu*angstrom^2'), symmetry=1, barrier=(3.4872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15167,'amu*angstrom^2'), symmetry=1, barrier=(3.4872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15167,'amu*angstrom^2'), symmetry=1, barrier=(3.4872,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.00217,0.141975,-0.000141857,7.44976e-08,-1.54574e-11,-17044.2,49.9783], Tmin=(100,'K'), Tmax=(1177.76,'K')), NASAPolynomial(coeffs=[26.4357,0.041996,-1.45227e-05,2.42061e-09,-1.57758e-13,-23978.4,-96.859], Tmin=(1177.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(Isobutyl)"""),
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
    label = 'C=C(O)C(O)=CC(5562)',
    structure = SMILES('C=C(O)C(O)=CC'),
    E0 = (-369.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.87493,'amu*angstrom^2'), symmetry=1, barrier=(20.1164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.877997,'amu*angstrom^2'), symmetry=1, barrier=(20.1869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.874595,'amu*angstrom^2'), symmetry=1, barrier=(20.1087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.875946,'amu*angstrom^2'), symmetry=1, barrier=(20.1397,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4302.09,'J/mol'), sigma=(6.83849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.98 K, Pc=30.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.884099,0.0877395,-9.34163e-05,4.76671e-08,-9.0713e-12,-44294.8,25.8217], Tmin=(100,'K'), Tmax=(1473.72,'K')), NASAPolynomial(coeffs=[23.5689,0.00887539,-4.29798e-07,-1.49469e-10,1.60597e-14,-50145.5,-97.0299], Tmin=(1473.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-369.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(CC)C1(O)CC1(O)[CH]C(32642)',
    structure = SMILES('[CH2]C(CC)C1(O)CC1(O)[CH]C'),
    E0 = (-82.2269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35233,0.123624,-8.55693e-05,1.33525e-08,6.93465e-12,-9646.47,46.7347], Tmin=(100,'K'), Tmax=(966.027,'K')), NASAPolynomial(coeffs=[26.051,0.0418181,-1.41369e-05,2.42111e-09,-1.64909e-13,-16804.7,-97.9596], Tmin=(966.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-82.2269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CCJCO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1(O)C(CC)CC1(O)[CH]C(32643)',
    structure = SMILES('[CH2]C1(O)C(CC)CC1(O)[CH]C'),
    E0 = (-83.0596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42162,0.123812,-8.44928e-05,1.46141e-08,5.20014e-12,-9743.28,44.3623], Tmin=(100,'K'), Tmax=(1008.53,'K')), NASAPolynomial(coeffs=[26.3925,0.0436396,-1.59809e-05,2.85955e-09,-1.9855e-13,-17289.9,-103.494], Tmin=(1008.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.0596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(627.743,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CJC(C)2O) + radical(CCJCO)"""),
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
    label = '[CH2]C(O)(C(=C)CC)C(O)=CC(32644)',
    structure = SMILES('[CH2]C(O)(C(=C)CC)C(O)=CC'),
    E0 = (-239.194,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1600,1633.63,2862.12,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152558,'amu*angstrom^2'), symmetry=1, barrier=(3.50762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152558,'amu*angstrom^2'), symmetry=1, barrier=(3.50762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152558,'amu*angstrom^2'), symmetry=1, barrier=(3.50762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152558,'amu*angstrom^2'), symmetry=1, barrier=(3.50762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152558,'amu*angstrom^2'), symmetry=1, barrier=(3.50762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152558,'amu*angstrom^2'), symmetry=1, barrier=(3.50762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152558,'amu*angstrom^2'), symmetry=1, barrier=(3.50762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152558,'amu*angstrom^2'), symmetry=1, barrier=(3.50762,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (155.214,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52348,0.131601,-0.000120682,5.6949e-08,-1.0698e-11,-28523.6,45.653], Tmin=(100,'K'), Tmax=(1285.25,'K')), NASAPolynomial(coeffs=[26.1093,0.0424891,-1.66809e-05,3.00284e-09,-2.04664e-13,-35883.7,-99.6691], Tmin=(1285.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-239.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ)"""),
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
    label = '[CH2]C(CC)C(=C)O(3441)',
    structure = SMILES('[CH2]C(CC)C(=C)O'),
    E0 = (-81.5304,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,410.858,411.713],'cm^-1')),
        HinderedRotor(inertia=(0.000993573,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114963,'amu*angstrom^2'), symmetry=1, barrier=(13.8198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114624,'amu*angstrom^2'), symmetry=1, barrier=(13.813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115855,'amu*angstrom^2'), symmetry=1, barrier=(13.8347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114897,'amu*angstrom^2'), symmetry=1, barrier=(13.8184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277625,0.0687697,-2.50371e-05,-2.36995e-08,1.65623e-11,-9659.83,29.5215], Tmin=(100,'K'), Tmax=(936.351,'K')), NASAPolynomial(coeffs=[18.1998,0.025014,-7.49664e-06,1.23004e-09,-8.41033e-14,-14454.3,-63.4438], Tmin=(936.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.5304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = 'CC=[C]O(31159)',
    structure = SMILES('CC=[C]O'),
    E0 = (72.7982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.344811,'amu*angstrom^2'), symmetry=1, barrier=(7.92788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346062,'amu*angstrom^2'), symmetry=1, barrier=(7.95665,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.692,0.0334424,-4.69125e-05,4.57787e-08,-1.77014e-11,8798.16,15.2755], Tmin=(100,'K'), Tmax=(837.57,'K')), NASAPolynomial(coeffs=[2.44779,0.0239997,-1.10022e-05,2.07309e-09,-1.42159e-13,9211.19,18.6318], Tmin=(837.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.7982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(O)=C(O)[CH]C(4609)',
    structure = SMILES('[CH2]C(O)=C(O)[CH]C'),
    E0 = (-122.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,322.653],'cm^-1')),
        HinderedRotor(inertia=(0.00160507,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208419,'amu*angstrom^2'), symmetry=1, barrier=(15.514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20977,'amu*angstrom^2'), symmetry=1, barrier=(15.5164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209392,'amu*angstrom^2'), symmetry=1, barrier=(15.5175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207973,'amu*angstrom^2'), symmetry=1, barrier=(15.5093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.141743,0.0688338,-2.99459e-05,-3.072e-08,2.30897e-11,-14621,29.0351], Tmin=(100,'K'), Tmax=(905.863,'K')), NASAPolynomial(coeffs=[24.2058,0.00628117,1.26066e-06,-4.23747e-10,2.91719e-14,-20774,-94.5791], Tmin=(905.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(CCJCO)"""),
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
    label = '[CH2]C(O)(C=C)C(O)=CC(32267)',
    structure = SMILES('[CH2]C(O)(C=C)C(O)=CC'),
    E0 = (-177.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.161,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.21264,0.10183,-9.49263e-05,4.47376e-08,-8.28732e-12,-21149.1,37.2819], Tmin=(100,'K'), Tmax=(1314.41,'K')), NASAPolynomial(coeffs=[23.1597,0.0276596,-1.02824e-05,1.80588e-09,-1.21659e-13,-27556.1,-86.9633], Tmin=(1314.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'OH(5)',
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
    label = '[CH2]C(CC)C(=C)C(O)=CC(32645)',
    structure = SMILES('[CH2]C(CC)C(=C)C(O)=CC'),
    E0 = (-66.6835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.215,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96544,0.11801,-9.17593e-05,2.72284e-08,6.79697e-13,-7793.66,39.9927], Tmin=(100,'K'), Tmax=(979.652,'K')), NASAPolynomial(coeffs=[23.7165,0.0401544,-1.39003e-05,2.38363e-09,-1.61043e-13,-14121.4,-89.9939], Tmin=(979.652,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.6835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(O)([C](C)CC)C(O)=CC(32646)',
    structure = SMILES('[CH2]C(O)([C](C)CC)C(O)=CC'),
    E0 = (-196.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.90383,0.135343,-0.000123338,5.82807e-08,-1.09084e-11,-23359.4,48.0471], Tmin=(100,'K'), Tmax=(1297.65,'K')), NASAPolynomial(coeffs=[27.5533,0.0414584,-1.48123e-05,2.52563e-09,-1.66805e-13,-31263.9,-106.826], Tmin=(1297.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCJ(C)CO) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(CC)C(C)([O])C(O)=CC(32647)',
    structure = SMILES('[CH2]C(CC)C(C)([O])C(O)=CC'),
    E0 = (-128.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3615,1277.5,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,401.057,739.57,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155279,'amu*angstrom^2'), symmetry=1, barrier=(3.57017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155279,'amu*angstrom^2'), symmetry=1, barrier=(3.57017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155279,'amu*angstrom^2'), symmetry=1, barrier=(3.57017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155279,'amu*angstrom^2'), symmetry=1, barrier=(3.57017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155279,'amu*angstrom^2'), symmetry=1, barrier=(3.57017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155279,'amu*angstrom^2'), symmetry=1, barrier=(3.57017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155279,'amu*angstrom^2'), symmetry=1, barrier=(3.57017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155279,'amu*angstrom^2'), symmetry=1, barrier=(3.57017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.57166,0.130484,-0.000115276,5.3875e-08,-1.00813e-11,-15174.8,48.7778], Tmin=(100,'K'), Tmax=(1289.77,'K')), NASAPolynomial(coeffs=[24.8081,0.0455704,-1.65203e-05,2.82942e-09,-1.86956e-13,-22237.5,-90.2807], Tmin=(1289.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-128.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2][C](CC)C(C)(O)C(O)=CC(32648)',
    structure = SMILES('[CH2][C](CC)C(C)(O)C(O)=CC'),
    E0 = (-204.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45041,0.127529,-9.87215e-05,2.8871e-08,9.36939e-13,-24382.5,47.3173], Tmin=(100,'K'), Tmax=(979.575,'K')), NASAPolynomial(coeffs=[25.3267,0.0435056,-1.50794e-05,2.58711e-09,-1.74841e-13,-31235.1,-93.3181], Tmin=(979.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2]C(O)(C(O)=CC)C(C)[CH]C(32649)',
    structure = SMILES('[CH2]C(O)(C(O)=CC)C(C)[CH]C'),
    E0 = (-154.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.54751,0.135719,-0.000128871,6.45382e-08,-1.29512e-11,-18331.5,49.1233], Tmin=(100,'K'), Tmax=(1203.37,'K')), NASAPolynomial(coeffs=[24.1465,0.0469881,-1.82683e-05,3.26411e-09,-2.21507e-13,-24756.1,-84.6017], Tmin=(1203.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-154.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C([O])(C(O)=CC)C(C)CC(32650)',
    structure = SMILES('[CH2]C([O])(C(O)=CC)C(C)CC'),
    E0 = (-119.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42264,0.13138,-0.000116617,5.44837e-08,-1.02569e-11,-14178.2,47.3345], Tmin=(100,'K'), Tmax=(1273.38,'K')), NASAPolynomial(coeffs=[23.935,0.0485849,-1.90877e-05,3.42303e-09,-2.32362e-13,-20890.8,-86.1958], Tmin=(1273.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)2OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C([CH]C)C(C)(O)C(O)=CC(32651)',
    structure = SMILES('[CH2]C([CH]C)C(C)(O)C(O)=CC'),
    E0 = (-162.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.67836,0.134642,-0.000127026,6.34231e-08,-1.26088e-11,-19329.1,50.4989], Tmin=(100,'K'), Tmax=(1221.66,'K')), NASAPolynomial(coeffs=[25.0159,0.0439647,-1.56896e-05,2.66678e-09,-1.75726e-13,-26095.7,-88.6547], Tmin=(1221.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC(C)C([CH2])(O)C(O)=CC(32652)',
    structure = SMILES('[CH2]CC(C)C([CH2])(O)C(O)=CC'),
    E0 = (-143.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,1600,1641.66,2854.89,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153151,'amu*angstrom^2'), symmetry=1, barrier=(3.52123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153151,'amu*angstrom^2'), symmetry=1, barrier=(3.52123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153151,'amu*angstrom^2'), symmetry=1, barrier=(3.52123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153151,'amu*angstrom^2'), symmetry=1, barrier=(3.52123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153151,'amu*angstrom^2'), symmetry=1, barrier=(3.52123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153151,'amu*angstrom^2'), symmetry=1, barrier=(3.52123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153151,'amu*angstrom^2'), symmetry=1, barrier=(3.52123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153151,'amu*angstrom^2'), symmetry=1, barrier=(3.52123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153151,'amu*angstrom^2'), symmetry=1, barrier=(3.52123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.88879,0.139881,-0.000135646,6.86009e-08,-1.37572e-11,-17029,49.5784], Tmin=(100,'K'), Tmax=(1212.79,'K')), NASAPolynomial(coeffs=[26.4751,0.0430336,-1.58637e-05,2.75679e-09,-1.84359e-13,-24151.4,-97.7504], Tmin=(1212.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(CC)C(C)(O)C(=O)[CH]C(7005)',
    structure = SMILES('[CH2]C(CC)C(C)(O)C([O])=CC'),
    E0 = (-219.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.46728,0.131488,-0.000121867,6.08381e-08,-1.21867e-11,-26161.9,48.6849], Tmin=(100,'K'), Tmax=(1208.19,'K')), NASAPolynomial(coeffs=[23.1005,0.0468392,-1.67736e-05,2.84822e-09,-1.87272e-13,-32340,-79.5], Tmin=(1208.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-219.534,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(CC)C(C)(O)C(O)=[C]C(32653)',
    structure = SMILES('[CH2]C(CC)C(C)(O)C(O)=[C]C'),
    E0 = (-119.497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3580,3650,1210,1345,900,1100,1685,370,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1600,1644.72,2855.65,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153283,'amu*angstrom^2'), symmetry=1, barrier=(3.52429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153283,'amu*angstrom^2'), symmetry=1, barrier=(3.52429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153283,'amu*angstrom^2'), symmetry=1, barrier=(3.52429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153283,'amu*angstrom^2'), symmetry=1, barrier=(3.52429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153283,'amu*angstrom^2'), symmetry=1, barrier=(3.52429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153283,'amu*angstrom^2'), symmetry=1, barrier=(3.52429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153283,'amu*angstrom^2'), symmetry=1, barrier=(3.52429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153283,'amu*angstrom^2'), symmetry=1, barrier=(3.52429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153283,'amu*angstrom^2'), symmetry=1, barrier=(3.52429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.86392,0.140046,-0.000139169,7.31915e-08,-1.52612e-11,-14115.7,49.0802], Tmin=(100,'K'), Tmax=(1169.96,'K')), NASAPolynomial(coeffs=[25.4164,0.0433579,-1.52049e-05,2.55373e-09,-1.6709e-13,-20733,-91.7948], Tmin=(1169.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)(C(=O)[CH]C)C(C)CC(7006)',
    structure = SMILES('[CH2]C(O)(C([O])=CC)C(C)CC'),
    E0 = (-211.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.34136,0.132615,-0.000123854,6.20976e-08,-1.25777e-11,-25164.1,47.3275], Tmin=(100,'K'), Tmax=(1187.8,'K')), NASAPolynomial(coeffs=[22.2302,0.0498676,-1.93566e-05,3.44683e-09,-2.33174e-13,-31001.3,-75.4449], Tmin=(1187.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-211.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)(C(O)=[C]C)C(C)CC(32654)',
    structure = SMILES('[CH2]C(O)(C(O)=[C]C)C(C)CC'),
    E0 = (-111.135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3580,3650,1210,1345,900,1100,1685,370,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1600,1694.01,2808.05,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15455,'amu*angstrom^2'), symmetry=1, barrier=(3.55341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15455,'amu*angstrom^2'), symmetry=1, barrier=(3.55341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15455,'amu*angstrom^2'), symmetry=1, barrier=(3.55341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15455,'amu*angstrom^2'), symmetry=1, barrier=(3.55341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15455,'amu*angstrom^2'), symmetry=1, barrier=(3.55341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15455,'amu*angstrom^2'), symmetry=1, barrier=(3.55341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15455,'amu*angstrom^2'), symmetry=1, barrier=(3.55341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15455,'amu*angstrom^2'), symmetry=1, barrier=(3.55341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15455,'amu*angstrom^2'), symmetry=1, barrier=(3.55341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.74899,0.141288,-0.000141492,7.48107e-08,-1.57786e-11,-13117.4,47.7634], Tmin=(100,'K'), Tmax=(1150.93,'K')), NASAPolynomial(coeffs=[24.5634,0.0463648,-1.77788e-05,3.15074e-09,-2.12894e-13,-19404.3,-87.8423], Tmin=(1150.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Cds_S) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]CC([CH2])C(C)(O)C(O)=CC(32655)',
    structure = SMILES('[CH2]CC([CH2])C(C)(O)C(O)=CC'),
    E0 = (-152.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.70755,0.135181,-0.000121454,5.20371e-08,-7.102e-12,-18040.1,49.8313], Tmin=(100,'K'), Tmax=(994.957,'K')), NASAPolynomial(coeffs=[25.6785,0.0427613,-1.48383e-05,2.52068e-09,-1.68187e-13,-24762.8,-92.3694], Tmin=(994.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C(O)C(C)(O)C([CH2])CC(32656)',
    structure = SMILES('[CH2]C=C(O)C(C)(O)C([CH2])CC'),
    E0 = (-205.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.68294,0.131203,-9.97375e-05,2.27664e-08,4.97185e-12,-24501.8,48.2381], Tmin=(100,'K'), Tmax=(953.317,'K')), NASAPolynomial(coeffs=[27.846,0.0395656,-1.29157e-05,2.1671e-09,-1.46352e-13,-31979.3,-106.276], Tmin=(953.317,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-205.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C(O)C([CH2])(O)C(C)CC(32657)',
    structure = SMILES('[CH2]C=C(O)C([CH2])(O)C(C)CC'),
    E0 = (-197.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.88959,0.136261,-0.000115509,4.19515e-08,-3.09001e-12,-23489.6,48.0731], Tmin=(100,'K'), Tmax=(1003.15,'K')), NASAPolynomial(coeffs=[28.0468,0.0407984,-1.44742e-05,2.52565e-09,-1.72467e-13,-31099.9,-108.27], Tmin=(1003.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-197.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(CC)C1(O)CC(C)[C]1O(32658)',
    structure = SMILES('[CH2]C(CC)C1(O)CC(C)[C]1O'),
    E0 = (-100.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.2035,0.121904,-8.99517e-05,2.73978e-08,-7.36699e-13,-11909,45.0966], Tmin=(100,'K'), Tmax=(1038.84,'K')), NASAPolynomial(coeffs=[23.1941,0.0479705,-1.76479e-05,3.10559e-09,-2.11155e-13,-18473.2,-84.5961], Tmin=(1038.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(627.743,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1(O)[C](O)C(C)CC1CC(32659)',
    structure = SMILES('[CH2]C1(O)[C](O)C(C)CC1CC'),
    E0 = (-179.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.13493,0.118995,-7.96105e-05,1.56405e-08,3.53314e-12,-21313,42.4066], Tmin=(100,'K'), Tmax=(1024.74,'K')), NASAPolynomial(coeffs=[23.741,0.0472743,-1.74911e-05,3.11391e-09,-2.14118e-13,-28153.8,-90.564], Tmin=(1024.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(631.9,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CJC(C)2O) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(CC)C(C)(O)C(O)=CC(32660)',
    structure = SMILES('C=C(CC)C(C)(O)C(O)=CC'),
    E0 = (-452.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.53282,0.127196,-0.00010065,3.61186e-08,-3.75556e-12,-54190.5,44.8115], Tmin=(100,'K'), Tmax=(1104.69,'K')), NASAPolynomial(coeffs=[26.1506,0.0445262,-1.71702e-05,3.10437e-09,-2.13958e-13,-61820.7,-102.277], Tmin=(1104.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-452.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(C)C([CH2])(O)C(O)=CC(6746)',
    structure = SMILES('[CH2]C(C)C([CH2])(O)C(O)=CC'),
    E0 = (-123.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4755.42,'J/mol'), sigma=(8.06767,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=742.79 K, Pc=20.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10142,0.124193,-0.000117561,5.30125e-08,-7.67636e-12,-14620.4,44.5008], Tmin=(100,'K'), Tmax=(964.677,'K')), NASAPolynomial(coeffs=[24.0863,0.0359813,-1.20785e-05,2.00981e-09,-1.3266e-13,-20621,-85.8115], Tmin=(964.677,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-123.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(CC)C([CH2])(O)C(=C)O(32661)',
    structure = SMILES('[CH2]C(CC)C([CH2])(O)C(=C)O'),
    E0 = (-107.869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,1557.37,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150296,'amu*angstrom^2'), symmetry=1, barrier=(3.45561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150296,'amu*angstrom^2'), symmetry=1, barrier=(3.45561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150296,'amu*angstrom^2'), symmetry=1, barrier=(3.45561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150296,'amu*angstrom^2'), symmetry=1, barrier=(3.45561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150296,'amu*angstrom^2'), symmetry=1, barrier=(3.45561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150296,'amu*angstrom^2'), symmetry=1, barrier=(3.45561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150296,'amu*angstrom^2'), symmetry=1, barrier=(3.45561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150296,'amu*angstrom^2'), symmetry=1, barrier=(3.45561,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.13521,0.123804,-0.000112283,4.45743e-08,-3.92142e-12,-12742.7,45.0612], Tmin=(100,'K'), Tmax=(948.597,'K')), NASAPolynomial(coeffs=[25.1077,0.0344038,-1.12011e-05,1.84644e-09,-1.22178e-13,-19057.4,-90.9741], Tmin=(948.597,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-107.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(CC)CC(O)=C(O)[CH]C(7103)',
    structure = SMILES('[CH2]C(CC)CC(O)=C(O)[CH]C'),
    E0 = (-176.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.59054,0.126245,-8.25774e-05,1.89433e-09,1.31217e-11,-20957,49.2973], Tmin=(100,'K'), Tmax=(946.862,'K')), NASAPolynomial(coeffs=[29.3754,0.0362138,-1.12533e-05,1.87824e-09,-1.28821e-13,-29028,-113.828], Tmin=(946.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-176.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(Isobutyl) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C(CC)[C](O)CC(O)=CC(32662)',
    structure = SMILES('[CH2]C(CC)[C](O)CC(O)=CC'),
    E0 = (-155.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45188,0.139396,-0.000142452,7.90912e-08,-1.77305e-11,-18447.9,47.3651], Tmin=(100,'K'), Tmax=(1078.79,'K')), NASAPolynomial(coeffs=[21.2257,0.0516016,-2.03771e-05,3.65096e-09,-2.47652e-13,-23556.4,-68.6611], Tmin=(1078.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-155.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(O)(C[CH]CC)C(O)=CC(7224)',
    structure = SMILES('[CH2]C(O)(C[CH]CC)C(O)=CC'),
    E0 = (-152.172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,545.772,600.579,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158721,'amu*angstrom^2'), symmetry=1, barrier=(3.64932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158721,'amu*angstrom^2'), symmetry=1, barrier=(3.64932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158721,'amu*angstrom^2'), symmetry=1, barrier=(3.64932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158721,'amu*angstrom^2'), symmetry=1, barrier=(3.64932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158721,'amu*angstrom^2'), symmetry=1, barrier=(3.64932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158721,'amu*angstrom^2'), symmetry=1, barrier=(3.64932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158721,'amu*angstrom^2'), symmetry=1, barrier=(3.64932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158721,'amu*angstrom^2'), symmetry=1, barrier=(3.64932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158721,'amu*angstrom^2'), symmetry=1, barrier=(3.64932,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.18745,0.136202,-0.000138158,7.72943e-08,-1.76499e-11,-18079.2,48.1256], Tmin=(100,'K'), Tmax=(1053.43,'K')), NASAPolynomial(coeffs=[19.2013,0.0549864,-2.25138e-05,4.10833e-09,-2.8144e-13,-22585.6,-56.1762], Tmin=(1053.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C(O)([CH]CCC)C(O)=CC(32663)',
    structure = SMILES('[CH2]C(O)([CH]CCC)C(O)=CC'),
    E0 = (-146.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.85411,0.13774,-0.000131091,6.49346e-08,-1.27483e-11,-17389.6,50.4448], Tmin=(100,'K'), Tmax=(1238.36,'K')), NASAPolynomial(coeffs=[26.7079,0.0422526,-1.54286e-05,2.66823e-09,-1.77924e-13,-24711.2,-98.4947], Tmin=(1238.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(CCJCO)"""),
)

species(
    label = 'CC=C(O)C1(O)CCC1CC(32038)',
    structure = SMILES('CC=C(O)C1(O)CCC1CC'),
    E0 = (-410.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78224,0.101756,-1.50851e-05,-6.21985e-08,3.37494e-11,-49193.3,41.8175], Tmin=(100,'K'), Tmax=(961.725,'K')), NASAPolynomial(coeffs=[27.8692,0.0394303,-1.30184e-05,2.32155e-09,-1.66909e-13,-57717.6,-114.743], Tmin=(961.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-410.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(631.9,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C(CC)C([CH2])(O)C(=O)CC(6344)',
    structure = SMILES('[CH2]C(CC)C([CH2])(O)C(=O)CC'),
    E0 = (-144.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,466.279,682.674,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156583,'amu*angstrom^2'), symmetry=1, barrier=(3.60014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156583,'amu*angstrom^2'), symmetry=1, barrier=(3.60014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156583,'amu*angstrom^2'), symmetry=1, barrier=(3.60014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156583,'amu*angstrom^2'), symmetry=1, barrier=(3.60014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156583,'amu*angstrom^2'), symmetry=1, barrier=(3.60014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156583,'amu*angstrom^2'), symmetry=1, barrier=(3.60014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156583,'amu*angstrom^2'), symmetry=1, barrier=(3.60014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156583,'amu*angstrom^2'), symmetry=1, barrier=(3.60014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156583,'amu*angstrom^2'), symmetry=1, barrier=(3.60014,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4659.74,'J/mol'), sigma=(8.03739,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=727.84 K, Pc=20.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.69487,0.151224,-0.000177053,1.15042e-07,-3.02452e-11,-17129,46.2525], Tmin=(100,'K'), Tmax=(924.594,'K')), NASAPolynomial(coeffs=[19.1962,0.0565174,-2.34056e-05,4.25578e-09,-2.89521e-13,-21177,-57.6425], Tmin=(924.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-144.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(Isobutyl) + radical(CJC(C)(C=O)O)"""),
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
    label = '[CH2]C(CC)[C](O)C(O)=CC(32664)',
    structure = SMILES('[CH2]C(CC)C(O)=C(O)[CH]C'),
    E0 = (-154.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.85464,0.110094,-6.32477e-05,-1.38539e-08,1.87343e-11,-18397.1,44.4471], Tmin=(100,'K'), Tmax=(930.222,'K')), NASAPolynomial(coeffs=[28.3159,0.0279366,-7.48852e-06,1.16846e-09,-7.9945e-14,-26068.7,-109.991], Tmin=(930.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-154.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CCJCO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(O)([CH]CC)C(O)=CC(32283)',
    structure = SMILES('[CH2]C(O)([CH]CC)C(O)=CC'),
    E0 = (-122.948,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1600,1708.75,2790.8,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154726,'amu*angstrom^2'), symmetry=1, barrier=(3.55746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154726,'amu*angstrom^2'), symmetry=1, barrier=(3.55746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154726,'amu*angstrom^2'), symmetry=1, barrier=(3.55746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154726,'amu*angstrom^2'), symmetry=1, barrier=(3.55746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154726,'amu*angstrom^2'), symmetry=1, barrier=(3.55746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154726,'amu*angstrom^2'), symmetry=1, barrier=(3.55746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154726,'amu*angstrom^2'), symmetry=1, barrier=(3.55746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154726,'amu*angstrom^2'), symmetry=1, barrier=(3.55746,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07213,0.121402,-0.000112046,5.05475e-08,-8.11801e-12,-14558.1,45.3905], Tmin=(100,'K'), Tmax=(1043.08,'K')), NASAPolynomial(coeffs=[24.4754,0.0360055,-1.28373e-05,2.2205e-09,-1.4963e-13,-20988.9,-88.0842], Tmin=(1043.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(CCJCO)"""),
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
    E0 = (-143.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-82.2269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-83.0596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-27.4018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-98.6148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-0.914233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-119.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-42.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-31.7077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-41.2153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (10.4423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-25.9557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (2.96312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-60.2149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-57.2861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-60.8876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-60.8425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-75.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-101.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-78.0952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-91.1765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-88.2477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-88.2477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (132.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-18.3749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-58.1229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-80.4947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (296.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (311.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (13.4235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (13.4235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (16.0403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (16.0403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-135.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-0.781275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (226.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (258.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['butene1(127)', 'C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C(CC)C1(O)CC1(O)[CH]C(32642)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(61.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra_cs2H] for rate rule [R4_S_D;doublebond_intra_HNd_secNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 59.6 to 61.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C1(O)C(CC)CC1(O)[CH]C(32643)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.65e+06,'s^-1'), n=1.3, Ea=(60.8352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_HNd;radadd_intra_cs2H] for rate rule [R5_SS_D;doublebond_intra_HNd_secNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 58.7 to 60.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(O)(C(=C)CC)C(O)=CC(32644)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]CC(130)', 'C=C(O)C(O)=CC(5562)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-CdOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(CC)C(=C)O(3441)', 'CC=[C]O(31159)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00712612,'m^3/(mol*s)'), n=2.40979, Ea=(7.81798,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CdsJ] for rate rule [Cds-OsCs_Cds;CdsJ-O2s]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['butene1(127)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2H5(29)', '[CH2]C(O)(C=C)C(O)=CC(32267)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1020,'cm^3/(mol*s)'), n=2.41, Ea=(27.3634,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CsJ-CsHH] for rate rule [Cds-Cs\O2s/H_Cds-HH;CsJ-CsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)', '[CH2]C(CC)C(=C)C(O)=CC(32645)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.118912,'m^3/(mol*s)'), n=2.18935, Ea=(6.60377,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdCs_Cds-HH;YJ] for rate rule [Cds-CdCs_Cds-HH;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C(O)([C](C)CC)C(O)=CC(32646)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(CC)C(C)([O])C(O)=CC(32647)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2][C](CC)C(C)(O)C(O)=CC(32648)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C(O)(C(O)=CC)C(C)[CH]C(32649)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C([O])(C(O)=CC)C(C)CC(32650)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C([CH]C)C(C)(O)C(O)=CC(32651)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC(C)C([CH2])(O)C(O)=CC(32652)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C(CC)C(C)(O)C(=O)[CH]C(7005)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00963743,'s^-1'), n=3.795, Ea=(83.0524,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_2H;O_H_out] + [R4H_SS(Cd)S;C_rad_out_2H;XH_out] for rate rule [R4H_SS(Cd)S;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(CC)C(C)(O)C(O)=[C]C(32653)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C(O)(C(=O)[CH]C)C(C)CC(7006)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);C_rad_out_2H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)(C(O)=[C]C)C(C)CC(32654)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]CC([CH2])C(C)(O)C(O)=CC(32655)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 111 used for R5H_CCC;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C=C(O)C(C)(O)C([CH2])CC(32656)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(121000,'s^-1'), n=1.9, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 92 used for R5H_SSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C=C(O)C([CH2])(O)C(C)CC(32657)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8010,'s^-1'), n=1.94, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 93 used for R6H_RSSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R6H_RSSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH]CC(130)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C(CC)C1(O)CC(C)[C]1O(32658)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.82e+08,'s^-1'), n=0.91, Ea=(125.52,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra_secNd;radadd_intra_cs2H] for rate rule [R4_Cs_RR_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C1(O)[C](O)C(C)CC1CC(32659)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.95277e+09,'s^-1'), n=0.6, Ea=(85.772,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_secNd;radadd_intra_cs2H] for rate rule [R5_SS_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['C=C(CC)C(C)(O)C(O)=CC(32660)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(S)(23)', '[CH2]C(C)C([CH2])(O)C(O)=CC(6746)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(S)(23)', '[CH2]C(CC)C([CH2])(O)C(=C)O(32661)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.97e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C(CC)CC(O)=C(O)[CH]C(7103)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C(CC)[C](O)CC(O)=CC(32662)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C(O)(C[CH]CC)C(O)=CC(7224)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C(O)([CH]CCC)C(O)=CC(32663)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['CC=C(O)C1(O)CCC1CC(32038)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(=O)CC(6344)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(19)', '[CH2]C(CC)[C](O)C(O)=CC(32664)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction37',
    reactants = ['CH2(19)', '[CH2]C(O)([CH]CC)C(O)=CC(32283)'],
    products = ['[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '5393',
    isomers = [
        '[CH2]C(CC)C([CH2])(O)C(O)=CC(7019)',
    ],
    reactants = [
        ('butene1(127)', 'C=C(O)C(O)=CC(5562)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '5393',
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

