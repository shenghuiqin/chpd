species(
    label = 'CC=[C]CC([O])O(21299)',
    structure = SMILES('CC=[C]CC([O])O'),
    E0 = (50.941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,280.451,282.67,2961.53],'cm^-1')),
        HinderedRotor(inertia=(0.111358,'amu*angstrom^2'), symmetry=1, barrier=(5.98255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110436,'amu*angstrom^2'), symmetry=1, barrier=(5.91327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239027,'amu*angstrom^2'), symmetry=1, barrier=(13.7021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.097224,'amu*angstrom^2'), symmetry=1, barrier=(5.97067,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10312,0.0704205,-9.28953e-05,8.23788e-08,-3.07053e-11,6224.63,29.478], Tmin=(100,'K'), Tmax=(781.226,'K')), NASAPolynomial(coeffs=[4.4893,0.0429091,-2.05377e-05,3.96225e-09,-2.76788e-13,6006.01,15.9647], Tmin=(781.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'HOCHO(59)',
    structure = SMILES('O=CO'),
    E0 = (-389.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1051.41,1051.94,1052.55,2787.25,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00861396,'amu*angstrom^2'), symmetry=1, barrier=(6.76824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4005.91,'J/mol'), sigma=(3.626,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89836,-0.00355878,3.55205e-05,-4.385e-08,1.71078e-11,-46778.6,7.34954], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.61383,0.00644964,-2.29083e-06,3.6716e-10,-2.18737e-14,-45330.3,0.847884], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-389.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HOCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'CC=C=CC([O])O(25362)',
    structure = SMILES('CC=C=CC([O])O'),
    E0 = (-23.6766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,540,610,2055,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.290758,'amu*angstrom^2'), symmetry=1, barrier=(6.68511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290617,'amu*angstrom^2'), symmetry=1, barrier=(6.68187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290238,'amu*angstrom^2'), symmetry=1, barrier=(6.67314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16189,0.069362,-9.59915e-05,8.64645e-08,-3.22245e-11,-2752.15,26.74], Tmin=(100,'K'), Tmax=(787.604,'K')), NASAPolynomial(coeffs=[4.65739,0.0404333,-1.96116e-05,3.79623e-09,-2.65308e-13,-2956.13,12.9114], Tmin=(787.604,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.6766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(CCOJ)"""),
)

species(
    label = 'CC=[C]CC(=O)O(25363)',
    structure = SMILES('CC=[C]CC(=O)O'),
    E0 = (-157.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47862,0.0501346,-2.51632e-05,3.17253e-09,4.56942e-13,-18788.7,24.3429], Tmin=(100,'K'), Tmax=(1556.11,'K')), NASAPolynomial(coeffs=[15.3258,0.025783,-1.2527e-05,2.40185e-09,-1.65164e-13,-24459.4,-52.9583], Tmin=(1556.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + radical(Cds_S)"""),
)

species(
    label = 'CC#CCC([O])O(25364)',
    structure = SMILES('CC#CCC([O])O'),
    E0 = (-26.9529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2250,500,550,3615,1277.5,1000,1380,1390,370,380,2900,435,248.989,249.022,2259.21],'cm^-1')),
        HinderedRotor(inertia=(0.217568,'amu*angstrom^2'), symmetry=1, barrier=(9.5757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0190073,'amu*angstrom^2'), symmetry=1, barrier=(68.8412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.56503,'amu*angstrom^2'), symmetry=1, barrier=(68.8399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217588,'amu*angstrom^2'), symmetry=1, barrier=(9.5757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32753,0.0619862,-7.03998e-05,5.25787e-08,-1.66681e-11,-3148.52,27.6597], Tmin=(100,'K'), Tmax=(841.181,'K')), NASAPolynomial(coeffs=[6.32886,0.0343299,-1.41748e-05,2.54351e-09,-1.70422e-13,-3852.87,5.2108], Tmin=(841.181,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.9529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCOJ)"""),
)

species(
    label = '[O][CH]O(171)',
    structure = SMILES('[O][CH]O'),
    E0 = (5.37818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,1989.68],'cm^-1')),
        HinderedRotor(inertia=(0.0937094,'amu*angstrom^2'), symmetry=1, barrier=(2.15456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40204,0.0209571,-4.96202e-05,5.95355e-08,-2.44634e-11,660.882,10.3084], Tmin=(100,'K'), Tmax=(868.739,'K')), NASAPolynomial(coeffs=[-0.343229,0.0179325,-9.40003e-06,1.81361e-09,-1.23835e-13,2076.48,32.2523], Tmin=(868.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.37818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCOJ)"""),
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
    label = 'C#CCC([O])O(25365)',
    structure = SMILES('C#CCC([O])O'),
    E0 = (15.2282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,3615,1277.5,1000,1380,1390,370,380,2900,435,268.756,1544.44],'cm^-1')),
        HinderedRotor(inertia=(0.129945,'amu*angstrom^2'), symmetry=1, barrier=(6.65677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12959,'amu*angstrom^2'), symmetry=1, barrier=(6.65721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21046,'amu*angstrom^2'), symmetry=1, barrier=(62.1657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71153,0.0542282,-7.26288e-05,5.93097e-08,-1.97479e-11,1910.19,22.8593], Tmin=(100,'K'), Tmax=(837.876,'K')), NASAPolynomial(coeffs=[6.37752,0.025989,-1.1397e-05,2.09477e-09,-1.41762e-13,1337.63,2.4232], Tmin=(837.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.2282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCOJ)"""),
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
    label = 'CC=[C]CC=O(24131)',
    structure = SMILES('CC=[C]CC=O'),
    E0 = (108.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,346.917],'cm^-1')),
        HinderedRotor(inertia=(0.0130479,'amu*angstrom^2'), symmetry=1, barrier=(13.4182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587155,'amu*angstrom^2'), symmetry=1, barrier=(13.4999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.583718,'amu*angstrom^2'), symmetry=1, barrier=(13.4208,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94517,0.0403263,-1.39653e-05,-3.21076e-09,1.9131e-12,13148.7,21.3521], Tmin=(100,'K'), Tmax=(1433.43,'K')), NASAPolynomial(coeffs=[11.7345,0.0259672,-1.24992e-05,2.41377e-09,-1.67731e-13,9010.97,-34.044], Tmin=(1433.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S)"""),
)

species(
    label = 'CC=C[CH]C([O])O(21301)',
    structure = SMILES('CC=C[CH]C([O])O'),
    E0 = (-69.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31194,0.0648238,-5.69032e-05,2.95007e-08,-6.85566e-12,-8325.11,25.4732], Tmin=(100,'K'), Tmax=(976.973,'K')), NASAPolynomial(coeffs=[7.3144,0.0402481,-1.91711e-05,3.75319e-09,-2.67086e-13,-9497.97,-3.34533], Tmin=(976.973,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = 'CC=[C]C[C](O)O(25366)',
    structure = SMILES('CC=[C]C[C](O)O'),
    E0 = (30.4823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.957178,0.0715501,-8.15516e-05,5.50027e-08,-1.56273e-11,3771.65,29.151], Tmin=(100,'K'), Tmax=(842.87,'K')), NASAPolynomial(coeffs=[8.7373,0.0346281,-1.58442e-05,3.03162e-09,-2.12425e-13,2460.12,-7.05364], Tmin=(842.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.4823,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Cs_P)"""),
)

species(
    label = 'C[C]=CCC([O])O(21297)',
    structure = SMILES('C[C]=CCC([O])O'),
    E0 = (50.941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,280.451,282.67,2961.53],'cm^-1')),
        HinderedRotor(inertia=(0.111358,'amu*angstrom^2'), symmetry=1, barrier=(5.98255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110436,'amu*angstrom^2'), symmetry=1, barrier=(5.91327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239027,'amu*angstrom^2'), symmetry=1, barrier=(13.7021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.097224,'amu*angstrom^2'), symmetry=1, barrier=(5.97067,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10312,0.0704205,-9.28953e-05,8.23788e-08,-3.07053e-11,6224.63,29.478], Tmin=(100,'K'), Tmax=(781.226,'K')), NASAPolynomial(coeffs=[4.4893,0.0429091,-2.05377e-05,3.96225e-09,-2.76788e-13,6006.01,15.9647], Tmin=(781.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'C[CH][CH]CC(=O)O(2686)',
    structure = SMILES('C[CH][CH]CC(=O)O'),
    E0 = (-122.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17012,0.0474662,-2.22047e-06,-6.36635e-08,5.79051e-11,-14641.3,27.12], Tmin=(100,'K'), Tmax=(471.07,'K')), NASAPolynomial(coeffs=[3.70061,0.0438014,-2.02627e-05,3.91947e-09,-2.77489e-13,-14889,19.7895], Tmin=(471.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCJCC=O) + radical(RCCJC)"""),
)

species(
    label = 'CC=[C][CH]C(O)O(25367)',
    structure = SMILES('CC=[C][CH]C(O)O'),
    E0 = (-57.8474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.778284,0.0739201,-7.60826e-05,4.32008e-08,-1.01704e-11,-6843.94,26.3008], Tmin=(100,'K'), Tmax=(1011.89,'K')), NASAPolynomial(coeffs=[11.231,0.0326006,-1.48319e-05,2.84688e-09,-2.0054e-13,-8959.36,-24.2513], Tmin=(1011.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.8474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C[CH]CC([O])O(21215)',
    structure = SMILES('[CH2]C=CCC([O])O'),
    E0 = (-35.4016,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,261.162,261.168,1562.56],'cm^-1')),
        HinderedRotor(inertia=(0.00247146,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17044,'amu*angstrom^2'), symmetry=1, barrier=(8.24966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00247143,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.521038,'amu*angstrom^2'), symmetry=1, barrier=(25.2195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4287.55,'J/mol'), sigma=(7.04051,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=669.70 K, Pc=27.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48626,0.0591337,-4.42185e-05,1.80161e-08,-3.2106e-12,-4170.77,27.9104], Tmin=(100,'K'), Tmax=(1241.16,'K')), NASAPolynomial(coeffs=[8.57453,0.0362899,-1.66108e-05,3.18718e-09,-2.23718e-13,-5930.31,-7.81778], Tmin=(1241.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.4016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCOJ)"""),
)

species(
    label = 'CC=CCC([O])[O](21303)',
    structure = SMILES('CC=CCC([O])[O]'),
    E0 = (38.8044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,1931.88,4000],'cm^-1')),
        HinderedRotor(inertia=(0.164297,'amu*angstrom^2'), symmetry=1, barrier=(3.77751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32723,'amu*angstrom^2'), symmetry=1, barrier=(30.5156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162596,'amu*angstrom^2'), symmetry=1, barrier=(3.73839,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22357,0.0442443,-1.94958e-05,2.88853e-09,-7.77708e-14,4667.55,21.4329], Tmin=(100,'K'), Tmax=(2597.69,'K')), NASAPolynomial(coeffs=[36.5856,0.0050177,-3.858e-06,6.75092e-10,-3.79648e-14,-16763,-179.255], Tmin=(2597.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.8044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = 'C[C]=[C]CC(O)O(25368)',
    structure = SMILES('C[C]=[C]CC(O)O'),
    E0 = (63.0777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1670,1700,300,440,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,300.873,4000],'cm^-1')),
        HinderedRotor(inertia=(0.247559,'amu*angstrom^2'), symmetry=1, barrier=(15.9027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159554,'amu*angstrom^2'), symmetry=1, barrier=(10.2496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159554,'amu*angstrom^2'), symmetry=1, barrier=(10.2496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159555,'amu*angstrom^2'), symmetry=1, barrier=(10.2495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159554,'amu*angstrom^2'), symmetry=1, barrier=(10.2495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.727082,0.0774908,-0.000104125,8.43308e-08,-2.82537e-11,7699.11,29.7506], Tmin=(100,'K'), Tmax=(782.844,'K')), NASAPolynomial(coeffs=[8.09132,0.0358188,-1.65292e-05,3.1357e-09,-2.16959e-13,6670.01,-3.18322], Tmin=(782.844,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.0777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]CC(O)O(21300)',
    structure = SMILES('[CH2]C=[C]CC(O)O'),
    E0 = (-23.2649,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,360.845,360.968],'cm^-1')),
        HinderedRotor(inertia=(0.000712306,'amu*angstrom^2'), symmetry=1, barrier=(8.08756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272923,'amu*angstrom^2'), symmetry=1, barrier=(25.2197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272724,'amu*angstrom^2'), symmetry=1, barrier=(25.2201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.087545,'amu*angstrom^2'), symmetry=1, barrier=(8.08654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171133,'amu*angstrom^2'), symmetry=1, barrier=(15.8138,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.884106,0.0688883,-6.51115e-05,3.33131e-08,-7.01493e-12,-2686.21,28.995], Tmin=(100,'K'), Tmax=(1128.36,'K')), NASAPolynomial(coeffs=[11.9831,0.0295427,-1.28071e-05,2.41031e-09,-1.68097e-13,-5190.96,-25.8918], Tmin=(1128.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.2649,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'CC=C=CC(O)O(25369)',
    structure = SMILES('CC=C=CC(O)O'),
    E0 = (-249.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05749,0.0688218,-7.08819e-05,4.24502e-08,-1.0786e-11,-29891.2,25.5767], Tmin=(100,'K'), Tmax=(932.282,'K')), NASAPolynomial(coeffs=[9.04462,0.0345523,-1.57435e-05,3.02094e-09,-2.12562e-13,-31380.4,-12.3965], Tmin=(932.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CC=CCC(=O)O(2690)',
    structure = SMILES('CC=CCC(=O)O'),
    E0 = (-394.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4171,0.0469088,-6.28232e-06,-1.70677e-08,7.07774e-12,-47388.6,24.071], Tmin=(100,'K'), Tmax=(1205.94,'K')), NASAPolynomial(coeffs=[12.4713,0.0321994,-1.52966e-05,3.01337e-09,-2.15088e-13,-51651.3,-37.9489], Tmin=(1205.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-394.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs)"""),
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
    label = 'C=[C]CC([O])O(25359)',
    structure = SMILES('C=[C]CC([O])O'),
    E0 = (86.9666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,3615,1277.5,1000,1380,1390,370,380,2900,435,238.441,1566.36,4000],'cm^-1')),
        HinderedRotor(inertia=(0.234323,'amu*angstrom^2'), symmetry=1, barrier=(9.45105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234426,'amu*angstrom^2'), symmetry=1, barrier=(9.45118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.764197,'amu*angstrom^2'), symmetry=1, barrier=(30.8136,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73203,0.0550369,-7.28715e-05,6.43307e-08,-2.39119e-11,10536.4,25.415], Tmin=(100,'K'), Tmax=(775.288,'K')), NASAPolynomial(coeffs=[4.59757,0.0329055,-1.58379e-05,3.06462e-09,-2.14506e-13,10312.8,13.7438], Tmin=(775.288,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.9666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C([CH]C)C([O])O(24164)',
    structure = SMILES('[CH2]C(=CC)C([O])O'),
    E0 = (-51.8111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47893,0.0634877,-4.95477e-05,5.76648e-10,2.02979e-11,-6148.44,25.87], Tmin=(100,'K'), Tmax=(533.4,'K')), NASAPolynomial(coeffs=[5.98219,0.0412814,-1.96197e-05,3.81552e-09,-2.69739e-13,-6793.36,5.43256], Tmin=(533.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.8111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCOJ)"""),
)

species(
    label = 'CC=C1CC(O)O1(25344)',
    structure = SMILES('CC=C1CC(O)O1'),
    E0 = (-291.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15747,0.043938,2.8486e-05,-7.9719e-08,3.75772e-11,-34939,20.7773], Tmin=(100,'K'), Tmax=(916.015,'K')), NASAPolynomial(coeffs=[19.1678,0.0138443,-1.74162e-06,1.44656e-10,-1.1714e-14,-40275.6,-75.6512], Tmin=(916.015,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-291.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(2methyleneoxetane)"""),
)

species(
    label = 'CC=[C]C[CH]OO(25370)',
    structure = SMILES('CC=[C]C[CH]OO'),
    E0 = (285.554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,211.41,2312.25],'cm^-1')),
        HinderedRotor(inertia=(0.343019,'amu*angstrom^2'), symmetry=1, barrier=(11.0408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29089,'amu*angstrom^2'), symmetry=1, barrier=(41.5071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28787,'amu*angstrom^2'), symmetry=1, barrier=(41.5035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.347265,'amu*angstrom^2'), symmetry=1, barrier=(11.0473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34726,'amu*angstrom^2'), symmetry=1, barrier=(11.0457,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35147,0.0682007,-5.07633e-05,-2.11486e-08,4.4851e-11,34429.7,25.6146], Tmin=(100,'K'), Tmax=(505.578,'K')), NASAPolynomial(coeffs=[6.97997,0.040632,-1.92958e-05,3.71904e-09,-2.60344e-13,33643.7,0.155328], Tmin=(505.578,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCsJOOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])O(3174)',
    structure = SMILES('[CH2]C([O])O'),
    E0 = (-19.5078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,1935.17],'cm^-1')),
        HinderedRotor(inertia=(0.200266,'amu*angstrom^2'), symmetry=1, barrier=(4.6045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200463,'amu*angstrom^2'), symmetry=1, barrier=(4.60904,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66061,0.0330653,-4.66659e-05,4.29179e-08,-1.60973e-11,-2301.52,17.2788], Tmin=(100,'K'), Tmax=(804.622,'K')), NASAPolynomial(coeffs=[3.98157,0.0198211,-9.52741e-06,1.83298e-09,-1.27458e-13,-2297.94,12.5363], Tmin=(804.622,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.5078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CJCO)"""),
)

species(
    label = '[C]=CC(24199)',
    structure = SMILES('[C]=CC'),
    E0 = (564.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.40488,'amu*angstrom^2'), symmetry=1, barrier=(9.30899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28451,0.0144722,-2.9151e-06,-2.04635e-09,7.91551e-13,67868.8,10.0416], Tmin=(100,'K'), Tmax=(1455.86,'K')), NASAPolynomial(coeffs=[5.67821,0.0121697,-4.94658e-06,9.00486e-10,-6.07653e-14,66718.9,-3.96135], Tmin=(1455.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CdCdJ2_triplet)"""),
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
    label = 'CC=[C]C[CH]O(24553)',
    structure = SMILES('CC=[C]C[CH]O'),
    E0 = (205.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,259.208,259.208],'cm^-1')),
        HinderedRotor(inertia=(0.152203,'amu*angstrom^2'), symmetry=1, barrier=(7.25684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152204,'amu*angstrom^2'), symmetry=1, barrier=(7.25684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152204,'amu*angstrom^2'), symmetry=1, barrier=(7.25684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152203,'amu*angstrom^2'), symmetry=1, barrier=(7.25684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15718,0.0667693,-8.34791e-05,6.49253e-08,-2.11989e-11,24846.1,24.8175], Tmin=(100,'K'), Tmax=(784.151,'K')), NASAPolynomial(coeffs=[7.30244,0.0331288,-1.47414e-05,2.75647e-09,-1.89473e-13,23952.9,-2.88592], Tmin=(784.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CCsJOH)"""),
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
    E0 = (50.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (196.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (101.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (203.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (167.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (50.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (181.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (231.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (208.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (208.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (284.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (192.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (126.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (212.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (152.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (240.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (144.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (366.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (129.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (139.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (506.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (220.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (59.2254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (378.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (544.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (448.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['HOCHO(59)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'CC=C=CC([O])O(25362)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.23e+08,'cm^3/(mol*s)'), n=1.64, Ea=(8.49352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2566 used for Cds-CsH_Ca;HJ
Exact match found for rate rule [Cds-CsH_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'CC=[C]CC(=O)O(25363)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'CC#CCC([O])O(25364)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.02e+09,'cm^3/(mol*s)'), n=1.64, Ea=(18.4933,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2702 used for Ct-Cs_Ct-Cs;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O][CH]O(171)', 'CH3CHCCH2(18175)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00429749,'m^3/(mol*s)'), n=2.45395, Ea=(16.6104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['HOCHO(59)', 'C=[C][CH]C(18176)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.87291,'m^3/(mol*s)'), n=1.39198, Ea=(79.0955,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-NdH_O;YJ] for rate rule [CO-NdH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 75.9 to 79.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3(17)', 'C#CCC([O])O(25365)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(178000,'cm^3/(mol*s)'), n=2.41, Ea=(30.1248,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2271 used for Ct-H_Ct-Cs;CsJ-HHH
Exact match found for rate rule [Ct-H_Ct-Cs;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OH(5)', 'CC=[C]CC=O(24131)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.04245,'cm^3/(mol*s)'), n=3.486, Ea=(94.7258,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [CO-CsH_O;OJ] for rate rule [CO-CsH_O;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['CC=C[CH]C([O])O(21301)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/O)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['CC=[C]C[C](O)O(25366)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[C]=CCC([O])O(21297)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['C[CH][CH]CC(=O)O(2686)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['CC=[C][CH]C(O)O(25367)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['C=C[CH]CC([O])O(21215)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['CC=CCC([O])[O](21303)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[C]=[C]CC(O)O(25368)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.68863e+10,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['[CH2]C=[C]CC(O)O(21300)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6Hall;O_rad_out;Cs_H_out_2H] for rate rule [R6HJ_3;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][CH]O(171)', 'C=[C][CH]C(18176)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['CC=C=CC(O)O(25369)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['CC=CCC(=O)O(2690)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(S)(23)', 'C=[C]CC([O])O(25359)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['CC=C1CC(O)O1(25344)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CC=[C]C[CH]OO(25370)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.263e+10,'s^-1'), n=0, Ea=(92.7451,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnOOH;C_rad_out_H/NonDeC] for rate rule [ROOH;C_rad_out_H/NonDeC]
Euclidian distance = 1.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([O])O(3174)', '[C]=CC(24199)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['O(4)', 'CC=[C]C[CH]O(24553)'],
    products = ['CC=[C]CC([O])O(21299)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '4225',
    isomers = [
        'CC=[C]CC([O])O(21299)',
    ],
    reactants = [
        ('HOCHO(59)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4225',
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

