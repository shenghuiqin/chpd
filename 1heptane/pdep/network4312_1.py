species(
    label = 'C=[C]C(C)C=[C]O(23559)',
    structure = SMILES('C=[C]C(C)C=[C]O'),
    E0 = (323.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.603266,'amu*angstrom^2'), symmetry=1, barrier=(13.8703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602863,'amu*angstrom^2'), symmetry=1, barrier=(13.861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60326,'amu*angstrom^2'), symmetry=1, barrier=(13.8701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602705,'amu*angstrom^2'), symmetry=1, barrier=(13.8574,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.528279,0.0735384,-7.46492e-05,4.05178e-08,-8.80462e-12,39092.4,28.3331], Tmin=(100,'K'), Tmax=(1116.7,'K')), NASAPolynomial(coeffs=[14.0405,0.0251369,-9.63296e-06,1.70258e-09,-1.14738e-13,36074.7,-38.3468], Tmin=(1116.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_S)"""),
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
    label = 'C=[C]C(C)C=C=O(24874)',
    structure = SMILES('C=[C]C(C)C=C=O'),
    E0 = (206.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,1685,370,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.570681,'amu*angstrom^2'), symmetry=1, barrier=(13.1211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570292,'amu*angstrom^2'), symmetry=1, barrier=(13.1121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570829,'amu*angstrom^2'), symmetry=1, barrier=(13.1245,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1467,0.0648707,-6.41104e-05,3.5778e-08,-8.34054e-12,24880.8,24.2417], Tmin=(100,'K'), Tmax=(1018.98,'K')), NASAPolynomial(coeffs=[9.95016,0.0303114,-1.3235e-05,2.49147e-09,-1.73589e-13,23086.8,-18.3951], Tmin=(1018.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCd(CCO)H) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=C=C(C)C=[C]O(26232)',
    structure = SMILES('C=C=C(C)C=[C]O'),
    E0 = (227.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.976981,'amu*angstrom^2'), symmetry=1, barrier=(22.4627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973515,'amu*angstrom^2'), symmetry=1, barrier=(22.383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973156,'amu*angstrom^2'), symmetry=1, barrier=(22.3748,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.176549,0.0764278,-7.74268e-05,3.85724e-08,-7.28371e-12,27568.9,27.4395], Tmin=(100,'K'), Tmax=(1446.11,'K')), NASAPolynomial(coeffs=[20.2704,0.0120227,-2.48111e-06,2.68985e-10,-1.31182e-14,22475.8,-75.9102], Tmin=(1446.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]C(C)C#CO(26233)',
    structure = SMILES('C=[C]C(C)C#CO'),
    E0 = (293.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,268.911,269.659],'cm^-1')),
        HinderedRotor(inertia=(0.480516,'amu*angstrom^2'), symmetry=1, barrier=(24.3846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287664,'amu*angstrom^2'), symmetry=1, barrier=(14.7792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288141,'amu*angstrom^2'), symmetry=1, barrier=(14.7791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287859,'amu*angstrom^2'), symmetry=1, barrier=(14.7698,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17037,0.0651082,-6.41419e-05,3.55245e-08,-8.2619e-12,35352.1,23.7203], Tmin=(100,'K'), Tmax=(1016.88,'K')), NASAPolynomial(coeffs=[9.80287,0.0311517,-1.40532e-05,2.68662e-09,-1.88799e-13,33596.5,-18.0711], Tmin=(1016.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CtH) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(C)C=[C]O(26234)',
    structure = SMILES('C#CC(C)C=[C]O'),
    E0 = (227.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2175,525,1380,1390,370,380,2900,435,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.923424,'amu*angstrom^2'), symmetry=1, barrier=(21.2313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.924962,'amu*angstrom^2'), symmetry=1, barrier=(21.2667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.926585,'amu*angstrom^2'), symmetry=1, barrier=(21.304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92646,'amu*angstrom^2'), symmetry=1, barrier=(21.3011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.111267,0.0764708,-8.13339e-05,4.36905e-08,-9.06803e-12,27531.6,26.6875], Tmin=(100,'K'), Tmax=(1226.76,'K')), NASAPolynomial(coeffs=[18.1998,0.0155697,-4.5189e-06,6.69692e-10,-4.07016e-14,23238.1,-63.6867], Tmin=(1226.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO)"""),
)

species(
    label = 'CH3(17)',
    structure = SMILES('[CH3]'),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=C=CC=[C]O(26235)',
    structure = SMILES('C=C=CC=[C]O'),
    E0 = (265.451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,1685,370,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.12687,'amu*angstrom^2'), symmetry=1, barrier=(25.9089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13022,'amu*angstrom^2'), symmetry=1, barrier=(25.9859,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23992,0.0497461,-2.38241e-05,-1.80422e-08,1.46169e-11,32036,21.7503], Tmin=(100,'K'), Tmax=(912.009,'K')), NASAPolynomial(coeffs=[17.9731,0.00578456,7.78995e-08,-1.32633e-10,8.74098e-15,27760,-64.1465], Tmin=(912.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO)"""),
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
    label = 'C=[C]C(C)[CH]C=O(24877)',
    structure = SMILES('C=[C]C(C)C=C[O]'),
    E0 = (225.692,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.757798,'amu*angstrom^2'), symmetry=1, barrier=(17.4233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756896,'amu*angstrom^2'), symmetry=1, barrier=(17.4025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756957,'amu*angstrom^2'), symmetry=1, barrier=(17.4039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.595106,0.0663933,-4.98608e-05,1.50176e-08,-3.13868e-13,27274.4,26.3748], Tmin=(100,'K'), Tmax=(1048.83,'K')), NASAPolynomial(coeffs=[15.3199,0.0240422,-9.03616e-06,1.61832e-09,-1.11427e-13,23426.2,-48.9858], Tmin=(1048.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C[C](C)C=[C]O(23557)',
    structure = SMILES('[CH2]C=C(C)C=[C]O'),
    E0 = (169.312,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.94182,'amu*angstrom^2'), symmetry=1, barrier=(21.6543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.942638,'amu*angstrom^2'), symmetry=1, barrier=(21.6731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.943764,'amu*angstrom^2'), symmetry=1, barrier=(21.699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.941309,'amu*angstrom^2'), symmetry=1, barrier=(21.6425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.639459,0.0611824,-2.20546e-05,-2.52954e-08,1.74848e-11,20496.3,27.01], Tmin=(100,'K'), Tmax=(923.008,'K')), NASAPolynomial(coeffs=[18.4153,0.0171899,-4.25839e-06,6.34914e-10,-4.33134e-14,15807.4,-64.9486], Tmin=(923.008,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C]C(C)[C]=CO(26236)',
    structure = SMILES('C=[C]C(C)[C]=CO'),
    E0 = (322.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.762066,'amu*angstrom^2'), symmetry=1, barrier=(17.5214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.762251,'amu*angstrom^2'), symmetry=1, barrier=(17.5256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.763211,'amu*angstrom^2'), symmetry=1, barrier=(17.5477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.76185,'amu*angstrom^2'), symmetry=1, barrier=(17.5164,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0478618,0.079199,-8.24936e-05,4.40184e-08,-9.16619e-12,38885.1,27.4009], Tmin=(100,'K'), Tmax=(1181.04,'K')), NASAPolynomial(coeffs=[17.7516,0.0192382,-6.33813e-06,1.02996e-09,-6.63433e-14,34703.4,-60.955], Tmin=(1181.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C)C=[C]O(23563)',
    structure = SMILES('[CH]=CC(C)C=[C]O'),
    E0 = (333.228,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.677179,'amu*angstrom^2'), symmetry=1, barrier=(15.5697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.677212,'amu*angstrom^2'), symmetry=1, barrier=(15.5704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.677669,'amu*angstrom^2'), symmetry=1, barrier=(15.5809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.676909,'amu*angstrom^2'), symmetry=1, barrier=(15.5635,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.30751,0.0745824,-7.39076e-05,3.81365e-08,-7.76593e-12,40216.6,29.0544], Tmin=(100,'K'), Tmax=(1199.18,'K')), NASAPolynomial(coeffs=[16.1898,0.0216047,-7.63978e-06,1.29557e-09,-8.54449e-14,36407.5,-50.4534], Tmin=(1199.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.228,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(C=C)C=[C]O(23341)',
    structure = SMILES('[CH2]C(C=C)C=[C]O'),
    E0 = (291.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,308.581,311.299],'cm^-1')),
        HinderedRotor(inertia=(0.19307,'amu*angstrom^2'), symmetry=1, barrier=(13.1808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192811,'amu*angstrom^2'), symmetry=1, barrier=(13.1621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233696,'amu*angstrom^2'), symmetry=1, barrier=(15.8598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190477,'amu*angstrom^2'), symmetry=1, barrier=(13.1799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3876.19,'J/mol'), sigma=(6.45323,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.45 K, Pc=32.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.516779,0.0713614,-6.49459e-05,2.71006e-08,-2.87808e-12,35155.2,29.6822], Tmin=(100,'K'), Tmax=(928.049,'K')), NASAPolynomial(coeffs=[15.0033,0.0222931,-7.24732e-06,1.17648e-09,-7.64266e-14,31890.5,-42.2273], Tmin=(928.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC(C)[C]=[C]O(23560)',
    structure = SMILES('C=CC(C)[C]=[C]O'),
    E0 = (323.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.603266,'amu*angstrom^2'), symmetry=1, barrier=(13.8703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602863,'amu*angstrom^2'), symmetry=1, barrier=(13.861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60326,'amu*angstrom^2'), symmetry=1, barrier=(13.8701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602705,'amu*angstrom^2'), symmetry=1, barrier=(13.8574,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.528279,0.0735384,-7.46492e-05,4.05178e-08,-8.80462e-12,39092.4,28.3331], Tmin=(100,'K'), Tmax=(1116.7,'K')), NASAPolynomial(coeffs=[14.0405,0.0251369,-9.63296e-06,1.70258e-09,-1.14738e-13,36074.7,-38.3468], Tmin=(1116.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=[C][C](C)C=CO(26237)',
    structure = SMILES('[CH2][C]=C(C)C=CO'),
    E0 = (167.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.40214,0.0639517,-1.96672e-05,-3.52572e-08,2.29674e-11,20278.5,25.2074], Tmin=(100,'K'), Tmax=(912.929,'K')), NASAPolynomial(coeffs=[21.3532,0.0125863,-1.70205e-06,1.35267e-10,-9.17201e-15,14768.2,-83.1883], Tmin=(912.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C([C]=C)C=CO(23562)',
    structure = SMILES('[CH2]C([C]=C)C=CO'),
    E0 = (289.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.767147,'amu*angstrom^2'), symmetry=1, barrier=(17.6382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.767557,'amu*angstrom^2'), symmetry=1, barrier=(17.6476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.767252,'amu*angstrom^2'), symmetry=1, barrier=(17.6406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.768349,'amu*angstrom^2'), symmetry=1, barrier=(17.6658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.196653,0.0797753,-8.24568e-05,4.32139e-08,-8.66164e-12,34957.9,29.5857], Tmin=(100,'K'), Tmax=(1347.83,'K')), NASAPolynomial(coeffs=[19.0451,0.0157889,-3.58751e-06,4.15161e-10,-2.05301e-14,30396.1,-66.6696], Tmin=(1347.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=[C]C(C)C=CO(26238)',
    structure = SMILES('[CH]=[C]C(C)C=CO'),
    E0 = (331.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.822409,'amu*angstrom^2'), symmetry=1, barrier=(18.9088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.82252,'amu*angstrom^2'), symmetry=1, barrier=(18.9114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.821516,'amu*angstrom^2'), symmetry=1, barrier=(18.8883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.822115,'amu*angstrom^2'), symmetry=1, barrier=(18.902,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.198843,0.0805077,-8.25127e-05,4.24391e-08,-8.40771e-12,40010.5,28.2183], Tmin=(100,'K'), Tmax=(1310.78,'K')), NASAPolynomial(coeffs=[19.8105,0.0158647,-4.43887e-06,6.456e-10,-3.89593e-14,35072.7,-72.556], Tmin=(1310.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(C)[CH][C]=O(20108)',
    structure = SMILES('C=CC(C)[CH][C]=O'),
    E0 = (172.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04789,0.0570797,-3.32893e-05,3.75577e-09,2.48006e-12,20876.4,27.5296], Tmin=(100,'K'), Tmax=(1044.1,'K')), NASAPolynomial(coeffs=[12.5453,0.0270081,-1.01651e-05,1.81075e-09,-1.23878e-13,17713.8,-32.0827], Tmin=(1044.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(CCJCHO)"""),
)

species(
    label = 'C=C=C(C)C=CO(26239)',
    structure = SMILES('C=C=C(C)C=CO'),
    E0 = (-11.8833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.193201,0.0679196,-2.73001e-05,-2.99055e-08,2.14416e-11,-1277.49,22.9789], Tmin=(100,'K'), Tmax=(920.068,'K')), NASAPolynomial(coeffs=[22.9648,0.0103741,-1.06624e-06,5.56426e-11,-5.44412e-15,-7222.38,-94.5191], Tmin=(920.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-11.8833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC(C)C=C=O(20117)',
    structure = SMILES('C=CC(C)C=C=O'),
    E0 = (-31.8108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.9074,0.0635453,-5.0945e-05,2.16043e-08,-3.75027e-12,-3711.03,24.6205], Tmin=(100,'K'), Tmax=(1355.76,'K')), NASAPolynomial(coeffs=[13.117,0.0275222,-1.10891e-05,2.00589e-09,-1.36338e-13,-7021.68,-37.9999], Tmin=(1355.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.8108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCd(CCO)H) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
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
    label = 'C=[C]CC=[C]O(26220)',
    structure = SMILES('C=[C]CC=[C]O'),
    E0 = (357.543,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,3615,1277.5,1000,3010,987.5,1337.5,450,1655,203.265,205.219],'cm^-1')),
        HinderedRotor(inertia=(0.488291,'amu*angstrom^2'), symmetry=1, barrier=(14.6756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.501044,'amu*angstrom^2'), symmetry=1, barrier=(14.6528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.46145,'amu*angstrom^2'), symmetry=1, barrier=(14.6894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37408,0.0502193,-3.25129e-05,5.47387e-10,5.20764e-12,43103.9,25.7099], Tmin=(100,'K'), Tmax=(951.123,'K')), NASAPolynomial(coeffs=[14.0925,0.0146303,-4.61399e-06,7.78019e-10,-5.36031e-14,39874.9,-39.2675], Tmin=(951.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = 'C=C(C)[CH]C=[C]O(26240)',
    structure = SMILES('[CH2]C(C)=CC=[C]O'),
    E0 = (167.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.6395,0.0611819,-2.20527e-05,-2.5298e-08,1.7486e-11,20320.2,27.1659], Tmin=(100,'K'), Tmax=(923.001,'K')), NASAPolynomial(coeffs=[18.4152,0.0171901,-4.25851e-06,6.34943e-10,-4.33158e-14,15631.3,-64.7919], Tmin=(923.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C([CH]C)C=[C]O(26141)',
    structure = SMILES('[CH2]C(C=[C]O)=CC'),
    E0 = (202.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.358423,0.0676042,-3.68312e-05,-1.2555e-08,1.33823e-11,24528.5,26.741], Tmin=(100,'K'), Tmax=(930.222,'K')), NASAPolynomial(coeffs=[19.9379,0.0152641,-3.79478e-06,5.79085e-10,-4.04404e-14,19507.7,-73.7099], Tmin=(930.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CJO)"""),
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
    label = 'H2CC(41)',
    structure = SMILES('[C]=C'),
    E0 = (401.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2480.69,'J/mol'), sigma=(4.48499,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=387.48 K, Pc=62.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28155,0.00697643,-2.38528e-06,-1.21078e-09,9.82042e-13,48319.2,5.92036], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.27807,0.00475623,-1.63007e-06,2.54623e-10,-1.4886e-14,48014,0.639979], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(401.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C[CH]C=[C]O(6301)',
    structure = SMILES('C[CH]C=[C]O'),
    E0 = (155.542,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,469.415],'cm^-1')),
        HinderedRotor(inertia=(0.107296,'amu*angstrom^2'), symmetry=1, barrier=(16.7767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107341,'amu*angstrom^2'), symmetry=1, barrier=(16.7766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107347,'amu*angstrom^2'), symmetry=1, barrier=(16.7768,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0197,0.0325941,8.31877e-06,-3.95085e-08,1.92255e-11,18788.8,19.0672], Tmin=(100,'K'), Tmax=(933.563,'K')), NASAPolynomial(coeffs=[13.3807,0.0117561,-2.93179e-06,4.69157e-10,-3.43545e-14,15454.4,-41.4596], Tmin=(933.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_S) + radical(C=CJO)"""),
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
    E0 = (323.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (448.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (456.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (520.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (452.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (434.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (514.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (463.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (486.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (432.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (556.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (438.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (475.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (499.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (474.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (368.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (436.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (475.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (709.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (387.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (363.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (777.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (468.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (418.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (332.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (556.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C(C)C=[C]O(23559)'],
    products = ['HCCOH(50)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'C=[C]C(C)C=C=O(24874)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.185e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""Estimated using template [Od_R;HJ] for rate rule [Od_Cdd;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C=C(C)C=[C]O(26232)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(77.9769,'m^3/(mol*s)'), n=1.64, Ea=(16.5268,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeCs_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=[C]C(C)C#CO(26233)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#CC(C)C=[C]O(26234)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH3(17)', 'C=C=CC=[C]O(26235)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00970044,'m^3/(mol*s)'), n=2.41, Ea=(33.193,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Ca;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]O(172)', 'CH3CHCCH2(18175)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00472174,'m^3/(mol*s)'), n=2.41, Ea=(20.1294,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HCCOH(50)', 'C=[C][CH]C(18176)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]C(C)C=[C]O(23559)'],
    products = ['C=[C]C(C)[CH]C=O(24877)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C[C](C)C=[C]O(23557)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.03e+09,'s^-1'), n=1.31, Ea=(263.174,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 156 used for R2H_S;C_rad_out_OneDe/Cs;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]C(C)[C]=CO(26236)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC(C)C=[C]O(23563)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]C(C)C=[C]O(23559)'],
    products = ['[CH2]C(C=C)C=[C]O(23341)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=CC(C)[C]=[C]O(23560)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.36833e+07,'s^-1'), n=1.41, Ea=(175.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cd_H_out_doubleC] for rate rule [R3H_SS_Cs;Cd_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(C)C=[C]O(23559)'],
    products = ['C=[C][C](C)C=CO(26237)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]C(C)C=[C]O(23559)'],
    products = ['[CH2]C([C]=C)C=CO(23562)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]C(C)C=CO(26238)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C]C(C)C=[C]O(23559)'],
    products = ['C=CC(C)[CH][C]=O(20108)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]O(172)', 'C=[C][CH]C(18176)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]C(C)C=[C]O(23559)'],
    products = ['C=C=C(C)C=CO(26239)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C]C(C)C=[C]O(23559)'],
    products = ['C=CC(C)C=C=O(20117)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(S)(23)', 'C=[C]CC=[C]O(26220)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=[C]C(C)C=[C]O(23559)'],
    products = ['C=C(C)[CH]C=[C]O(26240)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.95888e+10,'s^-1'), n=0.7315, Ea=(144.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HC)CJ;CJ;CH3] + [cCs(-HC)CJ;CdsJ;C] for rate rule [cCs(-HC)CJ;CdsJ;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]C(C)C=[C]O(23559)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C]C(C)C=[C]O(23559)'],
    products = ['C=C1C(O)=CC1C(26241)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriND_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H2CC(41)', 'C[CH]C=[C]O(6301)'],
    products = ['C=[C]C(C)C=[C]O(23559)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4312',
    isomers = [
        'C=[C]C(C)C=[C]O(23559)',
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
    label = '4312',
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

