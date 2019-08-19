species(
    label = '[CH]=C(O)C([CH2])=CC(25042)',
    structure = SMILES('[CH]=C(O)C([CH2])=CC'),
    E0 = (205.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.937249,'amu*angstrom^2'), symmetry=1, barrier=(21.5492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937452,'amu*angstrom^2'), symmetry=1, barrier=(21.5539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.936533,'amu*angstrom^2'), symmetry=1, barrier=(21.5327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937209,'amu*angstrom^2'), symmetry=1, barrier=(21.5483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.538846,0.0833987,-8.4439e-05,4.19564e-08,-7.91826e-12,24834.1,26.6263], Tmin=(100,'K'), Tmax=(1435.68,'K')), NASAPolynomial(coeffs=[21.8884,0.0131594,-2.95159e-06,3.55353e-10,-1.90557e-14,19193.6,-86.8999], Tmin=(1435.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(Cds_P) + radical(Allyl_P)"""),
)

species(
    label = 'HCCOH(50)',
    structure = SMILES('C#CO'),
    E0 = (80.0402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3615,1277.5,1000,287.694,287.7,1588.42],'cm^-1')),
        HinderedRotor(inertia=(0.002037,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05541,0.0252003,-3.80822e-05,3.09891e-08,-9.898e-12,9768.72,12.2272], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.3751,0.00549429,-1.88137e-06,2.93804e-10,-1.71772e-14,8932.78,-8.24498], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(80.0402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CH3CHCCH2(18175)',
    structure = SMILES('C=C=CC'),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.759584,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74635,0.0218189,8.22353e-06,-2.14768e-08,8.55624e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.6,'K')), NASAPolynomial(coeffs=[6.82078,0.0192338,-7.45622e-06,1.36536e-09,-9.53195e-14,16028,-10.4333], Tmin=(1025.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C1([CH]C)C=C1O(26242)',
    structure = SMILES('[CH2]C1([CH]C)C=C1O'),
    E0 = (368.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.505468,0.0648724,-3.69532e-05,-4.4877e-09,7.78475e-12,44515.7,26.661], Tmin=(100,'K'), Tmax=(999.125,'K')), NASAPolynomial(coeffs=[18.2421,0.0194783,-7.25787e-06,1.35752e-09,-9.83105e-14,39693,-65.2905], Tmin=(999.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsCs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(Neopentyl) + radical(Cs_S)"""),
)

species(
    label = '[CH]=[C]O(172)',
    structure = SMILES('[CH]=[C]O'),
    E0 = (348.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,1058.91,1059.8],'cm^-1')),
        HinderedRotor(inertia=(0.315045,'amu*angstrom^2'), symmetry=1, barrier=(7.24351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22339,0.0198737,-3.18466e-05,3.10436e-08,-1.17427e-11,41917.3,13.7606], Tmin=(100,'K'), Tmax=(829.31,'K')), NASAPolynomial(coeffs=[3.78101,0.0111997,-5.33323e-06,1.02849e-09,-7.13893e-14,42030.6,12.4155], Tmin=(829.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = 'C=[C][CH]C(18176)',
    structure = SMILES('[CH2][C]=CC'),
    E0 = (361.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.352622,'amu*angstrom^2'), symmetry=1, barrier=(8.10748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828631,'amu*angstrom^2'), symmetry=1, barrier=(19.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12013e-13,43485.7,14.8304], Tmin=(100,'K'), Tmax=(2065.83,'K')), NASAPolynomial(coeffs=[10.7464,0.014324,-5.20136e-06,8.69079e-10,-5.48385e-14,40045.6,-31.3799], Tmin=(2065.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Allyl_P)"""),
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
    label = 'C#CC([CH2])=CC(24773)',
    structure = SMILES('C#CC([CH2])=CC'),
    E0 = (345.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,350,440,435,1725,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.26681,'amu*angstrom^2'), symmetry=1, barrier=(18.6855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812654,'amu*angstrom^2'), symmetry=1, barrier=(18.6845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00731,'amu*angstrom^2'), symmetry=1, barrier=(46.152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45204,0.0489282,-2.77372e-05,2.38041e-09,2.3133e-12,41678.9,19.6959], Tmin=(100,'K'), Tmax=(1060.46,'K')), NASAPolynomial(coeffs=[11.5674,0.0232035,-8.93177e-06,1.61098e-09,-1.10941e-13,38834.6,-32.9939], Tmin=(1060.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(O)C(C)=[C]C(26257)',
    structure = SMILES('[CH]=C(O)C(C)=[C]C'),
    E0 = (291.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,1685,370,3615,1277.5,1000,3120,650,792.5,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.746923,'amu*angstrom^2'), symmetry=1, barrier=(17.1732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744765,'amu*angstrom^2'), symmetry=1, barrier=(17.1236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744374,'amu*angstrom^2'), symmetry=1, barrier=(17.1146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.74412,'amu*angstrom^2'), symmetry=1, barrier=(17.1088,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0344997,0.083468,-9.37826e-05,5.42682e-08,-1.23022e-11,35187.4,24.7546], Tmin=(100,'K'), Tmax=(1083.6,'K')), NASAPolynomial(coeffs=[16.8757,0.0213014,-7.72836e-06,1.3256e-09,-8.79341e-14,31537.5,-57.8469], Tmin=(1083.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=C([O])C(=C)[CH]C(24154)',
    structure = SMILES('[CH2]C(=CC)C(=C)[O]'),
    E0 = (95.7161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,220.195,220.299],'cm^-1')),
        HinderedRotor(inertia=(0.498422,'amu*angstrom^2'), symmetry=1, barrier=(17.2584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.502508,'amu*angstrom^2'), symmetry=1, barrier=(17.2558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.50011,'amu*angstrom^2'), symmetry=1, barrier=(17.2478,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3895.05,'J/mol'), sigma=(6.46861,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=608.40 K, Pc=32.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.462933,0.0694052,-5.31685e-05,1.37995e-08,1.49892e-12,11646.8,24.1211], Tmin=(100,'K'), Tmax=(975.072,'K')), NASAPolynomial(coeffs=[16.3563,0.0219779,-7.54735e-06,1.29951e-09,-8.85431e-14,7702.57,-56.486], Tmin=(975.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.7161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=[C]C)C(=C)O(25043)',
    structure = SMILES('[CH2]C(=[C]C)C(=C)O'),
    E0 = (195.753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.874718,'amu*angstrom^2'), symmetry=1, barrier=(20.1115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.875747,'amu*angstrom^2'), symmetry=1, barrier=(20.1352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.875587,'amu*angstrom^2'), symmetry=1, barrier=(20.1315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.875535,'amu*angstrom^2'), symmetry=1, barrier=(20.1303,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.277594,0.0819422,-8.399e-05,4.30667e-08,-8.50469e-12,23708,25.7548], Tmin=(100,'K'), Tmax=(1311.91,'K')), NASAPolynomial(coeffs=[20.2459,0.0158985,-4.51255e-06,6.64232e-10,-4.04409e-14,18621.4,-77.6936], Tmin=(1311.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(O)[C](C)C=C(23546)',
    structure = SMILES('[CH]=C(O)C(C)=C[CH2]'),
    E0 = (171.564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(1.00873,'amu*angstrom^2'), symmetry=1, barrier=(23.1927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00912,'amu*angstrom^2'), symmetry=1, barrier=(23.2018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00929,'amu*angstrom^2'), symmetry=1, barrier=(23.2056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00881,'amu*angstrom^2'), symmetry=1, barrier=(23.1945,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.391669,0.0784946,-7.47318e-05,3.5532e-08,-6.42971e-12,20808,27.3802], Tmin=(100,'K'), Tmax=(1521.5,'K')), NASAPolynomial(coeffs=[20.2962,0.0151323,-3.41734e-06,4.07865e-10,-2.14373e-14,15551.4,-77.696], Tmin=(1521.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=C([O])C(C)=CC(25044)',
    structure = SMILES('[CH]=C([O])C(C)=CC'),
    E0 = (191.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3120,650,792.5,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.618456,'amu*angstrom^2'), symmetry=1, barrier=(14.2195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.616914,'amu*angstrom^2'), symmetry=1, barrier=(14.1841,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.617175,'amu*angstrom^2'), symmetry=1, barrier=(14.1901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.453956,0.0746521,-7.56419e-05,4.0921e-08,-8.84595e-12,23140.2,24.2767], Tmin=(100,'K'), Tmax=(1123.8,'K')), NASAPolynomial(coeffs=[14.386,0.025063,-9.45254e-06,1.65578e-09,-1.11012e-13,20008.8,-44.5633], Tmin=(1123.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2][C](C=C)C(=C)O(20212)',
    structure = SMILES('[CH2]C=C([CH2])C(=C)O'),
    E0 = (75.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.415839,0.0639851,-2.0497e-05,-3.20472e-08,2.09296e-11,9279.5,24.3486], Tmin=(100,'K'), Tmax=(924.187,'K')), NASAPolynomial(coeffs=[20.6449,0.0148656,-3.1546e-06,4.41571e-10,-3.13973e-14,3899.04,-80.5301], Tmin=(924.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=C(O)[C]1CC1C(26258)',
    structure = SMILES('[CH]=C(O)[C]1CC1C'),
    E0 = (266.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13948,0.0397375,4.8349e-05,-1.0329e-07,4.59097e-11,32190.4,24.5183], Tmin=(100,'K'), Tmax=(933.067,'K')), NASAPolynomial(coeffs=[21.5236,0.0115096,-1.37395e-06,1.86266e-10,-2.08928e-14,25811.2,-86.2108], Tmin=(933.067,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCJ(C)CO) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]1C(O)=CC1C(26259)',
    structure = SMILES('[CH2]C1=C(O)[CH]C1C'),
    E0 = (210.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.079,0.039958,5.13807e-05,-1.0906e-07,4.86416e-11,25445.3,24.6821], Tmin=(100,'K'), Tmax=(930.08,'K')), NASAPolynomial(coeffs=[22.5513,0.00984029,-4.0598e-07,-3.71801e-12,-7.97007e-15,18759.5,-91.8222], Tmin=(930.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(Allyl_P) + radical(CCJCO)"""),
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
    label = '[CH]=C(O)C([CH2])=C(26260)',
    structure = SMILES('[CH]=C(O)C([CH2])=C'),
    E0 = (241.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,3615,1277.5,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(1.13783,'amu*angstrom^2'), symmetry=1, barrier=(26.1611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13411,'amu*angstrom^2'), symmetry=1, barrier=(26.0754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13406,'amu*angstrom^2'), symmetry=1, barrier=(26.0744,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.15449,0.07084,-7.41487e-05,3.65742e-08,-6.65361e-12,29156.8,22.753], Tmin=(100,'K'), Tmax=(1578.78,'K')), NASAPolynomial(coeffs=[20.4858,0.00533979,6.30051e-07,-3.00534e-10,2.448e-14,24285.3,-81.0377], Tmin=(1578.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)C[C]=CC(23574)',
    structure = SMILES('[CH]=C(O)C[C]=CC'),
    E0 = (323.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.587373,'amu*angstrom^2'), symmetry=1, barrier=(13.5049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587943,'amu*angstrom^2'), symmetry=1, barrier=(13.518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587091,'amu*angstrom^2'), symmetry=1, barrier=(13.4984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587729,'amu*angstrom^2'), symmetry=1, barrier=(13.513,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3874.43,'J/mol'), sigma=(6.44487,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.18 K, Pc=32.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.267225,0.0729642,-6.9046e-05,3.34454e-08,-6.36613e-12,39002.5,28.9654], Tmin=(100,'K'), Tmax=(1282.76,'K')), NASAPolynomial(coeffs=[17.2925,0.0198743,-6.96474e-06,1.18075e-09,-7.7957e-14,34634.6,-57.4112], Tmin=(1282.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=C1CC=C1O(26221)',
    structure = SMILES('CC=C1CC=C1O'),
    E0 = (-30.1943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.921932,0.0481461,2.14223e-05,-7.37213e-08,3.50821e-11,-3502.71,20.6101], Tmin=(100,'K'), Tmax=(933.826,'K')), NASAPolynomial(coeffs=[20.6955,0.0132148,-2.41001e-06,3.64394e-10,-3.08305e-14,-9365.7,-85.051], Tmin=(933.826,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.1943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
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
    label = '[CH]=C(O)[C]=CC(26261)',
    structure = SMILES('[CH]C(O)=C=CC'),
    E0 = (266.803,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.06048,'amu*angstrom^2'), symmetry=1, barrier=(47.3744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05935,'amu*angstrom^2'), symmetry=1, barrier=(47.3486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05975,'amu*angstrom^2'), symmetry=1, barrier=(47.3577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0086,0.058279,-4.0558e-05,8.16587e-09,1.98278e-12,32203.3,21.7341], Tmin=(100,'K'), Tmax=(995.957,'K')), NASAPolynomial(coeffs=[13.8273,0.0217201,-7.97395e-06,1.40037e-09,-9.55978e-14,28909.7,-43.7725], Tmin=(995.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]=C(O)C(=C)C=C(23545)',
    structure = SMILES('[CH]=C(O)C(=C)C=C'),
    E0 = (197.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3120,650,792.5,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.0556,'amu*angstrom^2'), symmetry=1, barrier=(24.2704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05415,'amu*angstrom^2'), symmetry=1, barrier=(24.2371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05483,'amu*angstrom^2'), symmetry=1, barrier=(24.2527,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.69961,0.0610925,-3.12113e-05,-1.22996e-08,1.18587e-11,23942.9,22.5244], Tmin=(100,'K'), Tmax=(942.779,'K')), NASAPolynomial(coeffs=[17.9818,0.0160821,-4.64641e-06,7.70692e-10,-5.43801e-14,19425.9,-66.5071], Tmin=(942.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)C(=C)C[CH2](26262)',
    structure = SMILES('[CH]=C(O)C(=C)C[CH2]'),
    E0 = (272.135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.752542,'amu*angstrom^2'), symmetry=1, barrier=(17.3024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.752606,'amu*angstrom^2'), symmetry=1, barrier=(17.3039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.752677,'amu*angstrom^2'), symmetry=1, barrier=(17.3055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.752483,'amu*angstrom^2'), symmetry=1, barrier=(17.3011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.469202,0.0825342,-8.40349e-05,4.20404e-08,-7.9928e-12,32904.6,28.9876], Tmin=(100,'K'), Tmax=(1421.34,'K')), NASAPolynomial(coeffs=[21.5957,0.0130667,-2.94345e-06,3.56286e-10,-1.91652e-14,27376.8,-82.6015], Tmin=(1421.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)C(=[CH])CC(26263)',
    structure = SMILES('[CH]=C(O)C(=[CH])CC'),
    E0 = (313.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3115,3125,620,680,785,800,1600,1700,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.827228,'amu*angstrom^2'), symmetry=1, barrier=(19.0196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.826874,'amu*angstrom^2'), symmetry=1, barrier=(19.0115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.826986,'amu*angstrom^2'), symmetry=1, barrier=(19.014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.825566,'amu*angstrom^2'), symmetry=1, barrier=(18.9814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.538076,0.0848654,-8.88098e-05,4.55405e-08,-8.86816e-12,37939.8,27.8488], Tmin=(100,'K'), Tmax=(1382.35,'K')), NASAPolynomial(coeffs=[21.9446,0.0126918,-2.7706e-06,3.21684e-10,-1.66691e-14,32404,-85.4369], Tmin=(1382.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]C)C(=C)O(25054)',
    structure = SMILES('[CH]C(=CC)C(=C)O'),
    E0 = (177.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.03874,'amu*angstrom^2'), symmetry=1, barrier=(46.8745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03878,'amu*angstrom^2'), symmetry=1, barrier=(46.8756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03847,'amu*angstrom^2'), symmetry=1, barrier=(46.8685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0386,'amu*angstrom^2'), symmetry=1, barrier=(46.8714,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.113172,0.0756004,-5.30937e-05,8.4857e-09,4.24234e-12,21448.6,24.7039], Tmin=(100,'K'), Tmax=(964.048,'K')), NASAPolynomial(coeffs=[17.6468,0.0249289,-8.60479e-06,1.47631e-09,-1.00446e-13,17042,-64.5649], Tmin=(964.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([O])C(=C)CC(25055)',
    structure = SMILES('[CH]=C([O])C(=C)CC'),
    E0 = (204.693,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.618978,'amu*angstrom^2'), symmetry=1, barrier=(14.2315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.619923,'amu*angstrom^2'), symmetry=1, barrier=(14.2532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.61794,'amu*angstrom^2'), symmetry=1, barrier=(14.2077,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0782997,0.0752855,-7.23544e-05,3.57004e-08,-6.87572e-12,24769.4,26.7348], Tmin=(100,'K'), Tmax=(1311.32,'K')), NASAPolynomial(coeffs=[18.0314,0.0188387,-5.86005e-06,9.1606e-10,-5.75252e-14,20205.7,-64.1923], Tmin=(1311.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[CH][C]1CC=C1O(26264)',
    structure = SMILES('C[CH]C1C[CH]C=1O'),
    E0 = (209.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28869,0.0368128,5.26518e-05,-1.04492e-07,4.5414e-11,25282.3,24.1481], Tmin=(100,'K'), Tmax=(938.612,'K')), NASAPolynomial(coeffs=[20.3646,0.013427,-2.51881e-06,4.24528e-10,-3.80885e-14,19150.5,-80.262], Tmin=(938.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(Allyl_S) + radical(CCJCO)"""),
)

species(
    label = 'C=CC(=C)C(=C)O(20221)',
    structure = SMILES('C=CC(=C)C(=C)O'),
    E0 = (-49.0993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770347,0.0576068,-1.46885e-05,-2.95e-08,1.76938e-11,-5776.76,21.8673], Tmin=(100,'K'), Tmax=(946.952,'K')), NASAPolynomial(coeffs=[17.748,0.0189038,-5.67397e-06,9.67926e-10,-6.91199e-14,-10472.3,-66.9299], Tmin=(946.952,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.0993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH]=C(O)C(C)[C]=C(23547)',
    structure = SMILES('[CH]=C(O)C(C)[C]=C'),
    E0 = (325.559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.643511,'amu*angstrom^2'), symmetry=1, barrier=(14.7956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.643671,'amu*angstrom^2'), symmetry=1, barrier=(14.7993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.643554,'amu*angstrom^2'), symmetry=1, barrier=(14.7966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.643523,'amu*angstrom^2'), symmetry=1, barrier=(14.7959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3817.52,'J/mol'), sigma=(6.40826,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.29 K, Pc=32.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.26628,0.0784406,-8.30581e-05,4.57659e-08,-9.96864e-12,39293.4,26.746], Tmin=(100,'K'), Tmax=(1120.94,'K')), NASAPolynomial(coeffs=[15.8404,0.0228646,-8.68712e-06,1.53386e-09,-1.03542e-13,35801.9,-50.1679], Tmin=(1120.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=C1C(O)=CC1C(26241)',
    structure = SMILES('C=C1C(O)=CC1C'),
    E0 = (-27.7384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.665907,0.0565704,-2.59669e-06,-4.88844e-08,2.63459e-11,-3200.69,17.9232], Tmin=(100,'K'), Tmax=(932.479,'K')), NASAPolynomial(coeffs=[20.5008,0.0141465,-2.9775e-06,4.50216e-10,-3.45867e-14,-8754.51,-86.327], Tmin=(932.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.7384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = 'CHCH3(T)(95)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438699,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82363,-0.000909515,3.2138e-05,-3.7348e-08,1.3309e-11,41371.4,7.10948], Tmin=(100,'K'), Tmax=(960.812,'K')), NASAPolynomial(coeffs=[4.30487,0.00943069,-3.27559e-06,5.95121e-10,-4.27307e-14,40709.1,1.84202], Tmin=(960.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C(O)[C]=C(8482)',
    structure = SMILES('[CH]C(O)=C=C'),
    E0 = (302.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,3615,1277.5,1000,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.06322,'amu*angstrom^2'), symmetry=1, barrier=(47.4376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0627,'amu*angstrom^2'), symmetry=1, barrier=(47.4255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70935,0.0420531,-1.76635e-05,-1.33767e-08,1.01044e-11,36511.9,17.4136], Tmin=(100,'K'), Tmax=(934.463,'K')), NASAPolynomial(coeffs=[13.5123,0.01244,-3.69296e-06,6.01951e-10,-4.15644e-14,33393.1,-43.6136], Tmin=(934.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(302.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
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
    E0 = (225.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (368.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (525.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (463.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (385.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (453.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (397.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (608.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (354.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (299.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (238.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (709.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (298.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (343.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (660.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (493.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (213.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (648.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (409.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (385.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (459.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (608.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (289.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (343.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (229.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (420.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (212.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (646.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['HCCOH(50)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(20.6479,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 20.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['[CH2]C1([CH]C)C=C1O(26242)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(163.983,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_D;doublebond_intra_HNd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 162.2 to 164.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C]O(172)', 'CH3CHCCH2(18175)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['HCCOH(50)', 'C=[C][CH]C(18176)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;CJ] for rate rule [Ct-O_Ct;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OH(5)', 'C#CC([CH2])=CC(24773)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(38.799,'m^3/(mol*s)'), n=1.55533, Ea=(11.1484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Ct_Ct;OJ_pri] + [Ct-Cd_Ct-H;YJ] for rate rule [Ct-Cd_Ct-H;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C(O)C(C)=[C]C(26257)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['C=C([O])C(=C)[CH]C(24154)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['[CH2]C(=[C]C)C(=C)O(25043)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_single] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['[CH]=C(O)[C](C)C=C(23546)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['[CH]=C([O])C(C)=CC(25044)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.0148956,'s^-1'), n=3.795, Ea=(94.7676,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_2H;O_H_out] + [R4H_S(Cd)SS;C_rad_out_2H;XH_out] for rate rule [R4H_S(Cd)SS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['[CH2][C](C=C)C(=C)O(20212)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 286 used for R5H_DSMS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R5H_DSMS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]O(172)', 'C=[C][CH]C(18176)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['[CH]=C(O)[C]1CC1C(26258)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.68393e+12,'s^-1'), n=-0.105173, Ea=(93.5715,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe_HNd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['[CH2][C]1C(O)=CC1C(26259)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CH2(S)(23)', '[CH]=C(O)C([CH2])=C(26260)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(O)C[C]=CC(23574)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['CC=C1CC=C1O(26221)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH2(19)', '[CH]=C(O)[C]=CC(26261)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH]=C(O)C(=C)C=C(23545)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(O)C(=C)C[CH2](26262)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(O)C(=[CH])CC(26263)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['[CH]=C([CH]C)C(=C)O(25054)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.27529e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_singleH] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['[CH]=C([O])C(=C)CC(25055)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;O_H_out] for rate rule [R4H_S(Cd)SS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['C[CH][C]1CC=C1O(26264)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['C=CC(=C)C(=C)O(20221)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(O)C(C)[C]=C(23547)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C(O)C([CH2])=CC(25042)'],
    products = ['C=C1C(O)=CC1C(26241)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CHCH3(T)(95)', '[CH]=C(O)[C]=C(8482)'],
    products = ['[CH]=C(O)C([CH2])=CC(25042)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4314',
    isomers = [
        '[CH]=C(O)C([CH2])=CC(25042)',
    ],
    reactants = [
        ('HCCOH(50)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4314',
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

