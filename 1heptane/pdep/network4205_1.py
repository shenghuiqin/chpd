species(
    label = 'C=C([CH]C)C(C)[O](24158)',
    structure = SMILES('[CH2]C(=CC)C(C)[O]'),
    E0 = (117.451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,641.794],'cm^-1')),
        HinderedRotor(inertia=(0.0127877,'amu*angstrom^2'), symmetry=1, barrier=(3.73839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0561695,'amu*angstrom^2'), symmetry=1, barrier=(16.4169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714058,'amu*angstrom^2'), symmetry=1, barrier=(16.4176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714083,'amu*angstrom^2'), symmetry=1, barrier=(16.4182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.55884,0.0646756,-3.09708e-05,-4.41003e-09,5.77486e-12,14259.4,28.3506], Tmin=(100,'K'), Tmax=(1059.75,'K')), NASAPolynomial(coeffs=[15.0773,0.0306369,-1.21771e-05,2.2532e-09,-1.5797e-13,10016.4,-48.0356], Tmin=(1059.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = 'CH3CHO(37)',
    structure = SMILES('CC=O'),
    E0 = (-178.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1305.64,1305.66,1305.67,3976.84],'cm^-1')),
        HinderedRotor(inertia=(0.136163,'amu*angstrom^2'), symmetry=1, barrier=(3.13064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-178.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]C1([CH]C)OC1C(24643)',
    structure = SMILES('[CH2]C1([CH]C)OC1C'),
    E0 = (205.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.614771,0.0822994,-7.67629e-05,3.81299e-08,-7.19935e-12,24840.1,29.0311], Tmin=(100,'K'), Tmax=(1518.91,'K')), NASAPolynomial(coeffs=[17.1785,0.0219416,-3.82483e-06,2.65164e-10,-4.02535e-15,20992.1,-59.1223], Tmin=(1518.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CJC(C)OC) + radical(CCJCO)"""),
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
    label = '[CH2]C(=CC)C(C)=O(24669)',
    structure = SMILES('[CH2]C(=CC)C(C)=O'),
    E0 = (-49.5112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,1187.57],'cm^-1')),
        HinderedRotor(inertia=(0.186085,'amu*angstrom^2'), symmetry=1, barrier=(4.27846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190269,'amu*angstrom^2'), symmetry=1, barrier=(4.37465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188299,'amu*angstrom^2'), symmetry=1, barrier=(4.32937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.801933,'amu*angstrom^2'), symmetry=1, barrier=(18.438,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40859,0.0567725,-3.34456e-05,9.31711e-09,-1.04527e-12,-5862.37,23.3674], Tmin=(100,'K'), Tmax=(1959.11,'K')), NASAPolynomial(coeffs=[14.5087,0.0300264,-1.2968e-05,2.34898e-09,-1.56103e-13,-10995.4,-48.6429], Tmin=(1959.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.5112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + radical(C=C(C=O)CJ)"""),
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
    label = 'C=C([CH]C)C=O(18210)',
    structure = SMILES('[CH2]C(C=O)=CC'),
    E0 = (-4.35622,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.27086,'amu*angstrom^2'), symmetry=1, barrier=(6.2276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.365722,'amu*angstrom^2'), symmetry=1, barrier=(17.6859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.763853,'amu*angstrom^2'), symmetry=1, barrier=(17.5625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69953,0.0484287,-3.10729e-05,9.36972e-09,-1.12943e-12,-439.904,19.3609], Tmin=(100,'K'), Tmax=(1872.08,'K')), NASAPolynomial(coeffs=[14.043,0.0220553,-9.94153e-06,1.84473e-09,-1.24546e-13,-5061.55,-47.9293], Tmin=(1872.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.35622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C=O)CJ)"""),
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
    label = 'C[CH][O](176)',
    structure = SMILES('C[CH][O]'),
    E0 = (157.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1642.51],'cm^-1')),
        HinderedRotor(inertia=(0.123965,'amu*angstrom^2'), symmetry=1, barrier=(2.85019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65562,0.0114444,2.34936e-06,-4.83164e-09,1.17966e-12,18963.9,10.3625], Tmin=(100,'K'), Tmax=(1718.65,'K')), NASAPolynomial(coeffs=[6.06294,0.0136322,-6.35953e-06,1.18407e-09,-7.90642e-14,16985.9,-5.90233], Tmin=(1718.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(=CC)[C](C)O(24670)',
    structure = SMILES('[CH2]C([CH]C)=C(C)O'),
    E0 = (-13.5491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.339834,0.0670757,-2.60831e-05,-1.9353e-08,1.37499e-11,-1485.76,25.009], Tmin=(100,'K'), Tmax=(970.922,'K')), NASAPolynomial(coeffs=[18.1836,0.0246245,-8.48662e-06,1.51431e-09,-1.07271e-13,-6414.8,-68.0897], Tmin=(970.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.5491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(O)C(=C)[CH]C(24161)',
    structure = SMILES('[CH2]C(=CC)C([CH2])O'),
    E0 = (98.6792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.072416,0.0757481,-6.26083e-05,2.66347e-08,-4.52233e-12,12018.6,30.9304], Tmin=(100,'K'), Tmax=(1413.97,'K')), NASAPolynomial(coeffs=[17.696,0.0258918,-9.71777e-06,1.69712e-09,-1.13123e-13,7034.88,-60.1975], Tmin=(1413.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.6792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CJCO)"""),
)

species(
    label = 'CC=C(C)[C](C)[O](24671)',
    structure = SMILES('C[CH]C(C)=C(C)[O]'),
    E0 = (-27.2435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.497558,0.0678563,-4.90851e-05,1.86594e-08,-2.89264e-12,-3142.97,25.6763], Tmin=(100,'K'), Tmax=(1513.01,'K')), NASAPolynomial(coeffs=[14.9951,0.0295286,-1.10869e-05,1.91648e-09,-1.26151e-13,-7529.96,-50.2695], Tmin=(1513.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.2435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[C]=C(C)C(C)[O](24672)',
    structure = SMILES('C[C]=C(C)C(C)[O]'),
    E0 = (203.793,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,261.633,261.645],'cm^-1')),
        HinderedRotor(inertia=(0.00246254,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200512,'amu*angstrom^2'), symmetry=1, barrier=(9.74024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200524,'amu*angstrom^2'), symmetry=1, barrier=(9.74025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200515,'amu*angstrom^2'), symmetry=1, barrier=(9.74027,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.892379,0.0675312,-4.9747e-05,1.94965e-08,-3.21674e-12,24622.9,27.3405], Tmin=(100,'K'), Tmax=(1378.19,'K')), NASAPolynomial(coeffs=[11.7931,0.0358935,-1.53132e-05,2.84004e-09,-1.95312e-13,21618.2,-28.7458], Tmin=(1378.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([O])C(C)=CC(24673)',
    structure = SMILES('[CH2]C([O])C(C)=CC'),
    E0 = (177.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,329.505,330.207],'cm^-1')),
        HinderedRotor(inertia=(0.00155332,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141827,'amu*angstrom^2'), symmetry=1, barrier=(10.9472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142097,'amu*angstrom^2'), symmetry=1, barrier=(10.9497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1415,'amu*angstrom^2'), symmetry=1, barrier=(10.9444,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.506098,0.0711669,-5.49282e-05,2.21037e-08,-3.63544e-12,21483.5,28.8939], Tmin=(100,'K'), Tmax=(1423.37,'K')), NASAPolynomial(coeffs=[14.7888,0.0310294,-1.26299e-05,2.29236e-09,-1.55797e-13,17417.6,-45.0539], Tmin=(1423.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(=[C]C)C(C)O(24674)',
    structure = SMILES('[CH2]C(=[C]C)C(C)O'),
    E0 = (124.932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1685,370,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0253949,'amu*angstrom^2'), symmetry=1, barrier=(13.7224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.025473,'amu*angstrom^2'), symmetry=1, barrier=(13.7443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188368,'amu*angstrom^2'), symmetry=1, barrier=(4.33095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.597229,'amu*angstrom^2'), symmetry=1, barrier=(13.7315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.596684,'amu*angstrom^2'), symmetry=1, barrier=(13.7189,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.46359,0.072069,-5.73265e-05,2.39457e-08,-4.08223e-12,15157.7,29.3583], Tmin=(100,'K'), Tmax=(1381.1,'K')), NASAPolynomial(coeffs=[14.7652,0.0306482,-1.234e-05,2.23053e-09,-1.51463e-13,11207.3,-44.2566], Tmin=(1381.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[C](C)C(C)[O](19565)',
    structure = SMILES('[CH2]C=C(C)C(C)[O]'),
    E0 = (117.451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,641.794],'cm^-1')),
        HinderedRotor(inertia=(0.0127877,'amu*angstrom^2'), symmetry=1, barrier=(3.73839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0561695,'amu*angstrom^2'), symmetry=1, barrier=(16.4169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714058,'amu*angstrom^2'), symmetry=1, barrier=(16.4176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714083,'amu*angstrom^2'), symmetry=1, barrier=(16.4182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.55884,0.0646756,-3.09708e-05,-4.41003e-09,5.77486e-12,14259.4,28.3506], Tmin=(100,'K'), Tmax=(1059.75,'K')), NASAPolynomial(coeffs=[15.0773,0.0306369,-1.21771e-05,2.2532e-09,-1.5797e-13,10016.4,-48.0356], Tmin=(1059.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C](C=C)C(C)O(19569)',
    structure = SMILES('[CH2]C=C([CH2])C(C)O'),
    E0 = (38.5893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.435667,0.065741,-2.69863e-05,-1.41871e-08,1.06649e-11,4780.69,29.2628], Tmin=(100,'K'), Tmax=(1002.22,'K')), NASAPolynomial(coeffs=[16.9502,0.0271769,-1.0199e-05,1.87282e-09,-1.32756e-13,96.9977,-57.2989], Tmin=(1002.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.5893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'CC([O])[C]1CC1C(24675)',
    structure = SMILES('CC([O])[C]1CC1C'),
    E0 = (163.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32252,0.0393549,4.70874e-05,-8.87805e-08,3.65584e-11,19732,25.0581], Tmin=(100,'K'), Tmax=(967.935,'K')), NASAPolynomial(coeffs=[15.7936,0.0277774,-9.70425e-06,1.80725e-09,-1.33175e-13,14671.5,-55.9545], Tmin=(967.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2][C]1C(C)OC1C(24676)',
    structure = SMILES('[CH2][C]1C(C)OC1C'),
    E0 = (155.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51282,0.0374781,5.35367e-05,-1.02414e-07,4.56303e-11,18828.5,22.7583], Tmin=(100,'K'), Tmax=(887.312,'K')), NASAPolynomial(coeffs=[14.5291,0.0253663,-4.70757e-06,4.9132e-10,-2.71513e-14,14685.4,-48.8112], Tmin=(887.312,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(Isobutyl) + radical(CCJ(C)CO)"""),
)

species(
    label = 'CC=C(C)C(C)=O(24677)',
    structure = SMILES('CC=C(C)C(C)=O'),
    E0 = (-205.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83065,0.0541944,-2.67796e-05,5.50334e-09,-4.00969e-13,-24687.9,22.3534], Tmin=(100,'K'), Tmax=(2362.43,'K')), NASAPolynomial(coeffs=[25.8615,0.0195959,-8.67832e-06,1.48642e-09,-9.13574e-14,-37741.5,-117.838], Tmin=(2362.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-205.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH)"""),
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
    label = 'C=C([CH]C)C[O](24151)',
    structure = SMILES('[CH2]C(=CC)C[O]'),
    E0 = (155.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,350,440,435,1725,180,1015.72],'cm^-1')),
        HinderedRotor(inertia=(0.160412,'amu*angstrom^2'), symmetry=1, barrier=(3.68818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160666,'amu*angstrom^2'), symmetry=1, barrier=(3.69403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160864,'amu*angstrom^2'), symmetry=1, barrier=(3.69859,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3733.31,'J/mol'), sigma=(6.34131,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=583.13 K, Pc=33.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52369,0.0395386,1.72519e-05,-1.01149e-07,8.88762e-11,18701.4,20.0906], Tmin=(100,'K'), Tmax=(431.253,'K')), NASAPolynomial(coeffs=[3.36701,0.040821,-1.88765e-05,3.65575e-09,-2.59104e-13,18544,15.7498], Tmin=(431.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(=C)C(C)[O](24678)',
    structure = SMILES('[CH2]C(=C)C(C)[O]'),
    E0 = (153.476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,410.605,410.607],'cm^-1')),
        HinderedRotor(inertia=(0.000999952,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180847,'amu*angstrom^2'), symmetry=1, barrier=(21.6369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180856,'amu*angstrom^2'), symmetry=1, barrier=(21.6368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27027,0.0483764,-8.09619e-06,-2.54702e-08,1.34651e-11,18567.6,23.2955], Tmin=(100,'K'), Tmax=(993.121,'K')), NASAPolynomial(coeffs=[14.447,0.0218685,-8.18133e-06,1.52046e-09,-1.09278e-13,14640.3,-46.7791], Tmin=(993.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = 'CC=[C]CC(C)[O](19705)',
    structure = SMILES('CC=[C]CC(C)[O]'),
    E0 = (219.405,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,264.049,264.049,4000],'cm^-1')),
        HinderedRotor(inertia=(0.13516,'amu*angstrom^2'), symmetry=1, barrier=(6.6872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34652,'amu*angstrom^2'), symmetry=1, barrier=(17.1445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346521,'amu*angstrom^2'), symmetry=1, barrier=(17.1445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346523,'amu*angstrom^2'), symmetry=1, barrier=(17.1445,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998929,0.066046,-4.73771e-05,1.82447e-08,-2.98999e-12,26495.9,27.5358], Tmin=(100,'K'), Tmax=(1372.38,'K')), NASAPolynomial(coeffs=[10.8888,0.0372205,-1.58711e-05,2.93988e-09,-2.01976e-13,23781.4,-23.3077], Tmin=(1372.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'CC=C1COC1C(24668)',
    structure = SMILES('CC=C1COC1C'),
    E0 = (-68.7817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5534,0.036005,4.73255e-05,-8.14737e-08,3.1827e-11,-8168.72,23.989], Tmin=(100,'K'), Tmax=(997.491,'K')), NASAPolynomial(coeffs=[13.3299,0.0323515,-1.27011e-05,2.4349e-09,-1.77954e-13,-12685.7,-43.6615], Tmin=(997.491,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.7817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutane)"""),
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
    label = '[CH2]C([CH]C)=CC(24211)',
    structure = SMILES('[CH2]C([CH]C)=CC'),
    E0 = (201.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.147404,'amu*angstrom^2'), symmetry=1, barrier=(3.38911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0566179,'amu*angstrom^2'), symmetry=1, barrier=(22.7485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00845232,'amu*angstrom^2'), symmetry=1, barrier=(3.39199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52997,'amu*angstrom^2'), symmetry=1, barrier=(35.177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3144,0.0486072,-4.1552e-06,-2.37791e-08,1.10434e-11,24281.5,21.4566], Tmin=(100,'K'), Tmax=(1043.43,'K')), NASAPolynomial(coeffs=[11.482,0.0321073,-1.27487e-05,2.35693e-09,-1.65194e-13,20936.1,-33.8921], Tmin=(1043.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_S)"""),
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
    label = 'CC=[C]C(C)[O](24679)',
    structure = SMILES('CC=[C]C(C)[O]'),
    E0 = (242.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,242.5,242.716],'cm^-1')),
        HinderedRotor(inertia=(0.00288737,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262722,'amu*angstrom^2'), symmetry=1, barrier=(10.9494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261725,'amu*angstrom^2'), symmetry=1, barrier=(10.9048,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28711,0.0525351,-3.54933e-05,1.20898e-08,-1.668e-12,29311.4,25.1871], Tmin=(100,'K'), Tmax=(1678.08,'K')), NASAPolynomial(coeffs=[13.6691,0.0230202,-9.11045e-06,1.60835e-09,-1.06474e-13,25155.9,-40.9583], Tmin=(1678.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC(=C)C(C)[O](19563)',
    structure = SMILES('C=CC(=C)C(C)[O]'),
    E0 = (91.2392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,288.007,288.127,288.15],'cm^-1')),
        HinderedRotor(inertia=(0.296206,'amu*angstrom^2'), symmetry=1, barrier=(17.452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296908,'amu*angstrom^2'), symmetry=1, barrier=(17.4508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296833,'amu*angstrom^2'), symmetry=1, barrier=(17.4491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.638038,0.0595943,-1.33579e-05,-3.02915e-08,1.7183e-11,11107.5,26.65], Tmin=(100,'K'), Tmax=(972.872,'K')), NASAPolynomial(coeffs=[17.9001,0.0222804,-7.7235e-06,1.41035e-09,-1.021e-13,6155.79,-64.3413], Tmin=(972.872,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.2392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]CC(=C)C(C)[O](24680)',
    structure = SMILES('[CH2]CC(=C)C(C)[O]'),
    E0 = (184.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,1124.04,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0201206,'amu*angstrom^2'), symmetry=1, barrier=(18.0411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784658,'amu*angstrom^2'), symmetry=1, barrier=(18.0408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134503,'amu*angstrom^2'), symmetry=1, barrier=(18.0404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784606,'amu*angstrom^2'), symmetry=1, barrier=(18.0396,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.594582,0.0641919,-3.18119e-05,-2.82972e-09,5.1117e-12,22331.4,30.8349], Tmin=(100,'K'), Tmax=(1066.1,'K')), NASAPolynomial(coeffs=[14.8758,0.0303879,-1.20788e-05,2.23287e-09,-1.56321e-13,18162.3,-44.2494], Tmin=(1066.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJ) + radical(CC(C)OJ)"""),
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
    label = '[CH]=C(CC)C(C)[O](24681)',
    structure = SMILES('[CH]=C(CC)C(C)[O]'),
    E0 = (226.428,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,370.623,372.769],'cm^-1')),
        HinderedRotor(inertia=(0.00125742,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139188,'amu*angstrom^2'), symmetry=1, barrier=(12.782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129768,'amu*angstrom^2'), symmetry=1, barrier=(12.7465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132419,'amu*angstrom^2'), symmetry=1, barrier=(12.7182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.473251,0.0671033,-3.84455e-05,2.85252e-09,3.39859e-12,27368.9,29.8869], Tmin=(100,'K'), Tmax=(1071.99,'K')), NASAPolynomial(coeffs=[15.3854,0.029744,-1.17532e-05,2.1627e-09,-1.50909e-13,23121.2,-47.9919], Tmin=(1071.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([O])C(=C)CC(24682)',
    structure = SMILES('[CH2]C([O])C(=C)CC'),
    E0 = (190.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,345.029,345.035,345.063],'cm^-1')),
        HinderedRotor(inertia=(0.00141548,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00141544,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176276,'amu*angstrom^2'), symmetry=1, barrier=(14.8794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17607,'amu*angstrom^2'), symmetry=1, barrier=(14.8793,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.396501,0.0687613,-4.15727e-05,4.78082e-09,3.04216e-12,23101.2,30.3923], Tmin=(100,'K'), Tmax=(1063.33,'K')), NASAPolynomial(coeffs=[15.8668,0.0289465,-1.13369e-05,2.08092e-09,-1.45217e-13,18772.1,-50.079], Tmin=(1063.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH]=C([CH]C)C(C)O(24683)',
    structure = SMILES('[CH]C(=CC)C(C)O'),
    E0 = (106.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.241085,0.0726576,-4.93022e-05,1.71706e-08,-2.44319e-12,12925.5,30.5282], Tmin=(100,'K'), Tmax=(1625.21,'K')), NASAPolynomial(coeffs=[16.2561,0.0332406,-1.29214e-05,2.24683e-09,-1.47496e-13,7720.07,-54.5123], Tmin=(1625.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C[CH][C]1COC1C(24684)',
    structure = SMILES('C[CH][C]1COC1C'),
    E0 = (157.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67373,0.0369961,4.24224e-05,-7.95597e-08,3.37896e-11,19072.4,24.6825], Tmin=(100,'K'), Tmax=(917.853,'K')), NASAPolynomial(coeffs=[11.1378,0.0322692,-9.53152e-06,1.52295e-09,-1.01809e-13,15796.8,-28.5446], Tmin=(917.853,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.771,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCJ(C)CO) + radical(Cs_S)"""),
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
    label = 'C=CC(=C)C(C)O(19576)',
    structure = SMILES('C=CC(=C)C(C)O'),
    E0 = (-139.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.476309,0.0635763,-1.92267e-05,-2.52998e-08,1.5594e-11,-16593.1,27.4223], Tmin=(100,'K'), Tmax=(972.519,'K')), NASAPolynomial(coeffs=[17.9627,0.0240768,-8.31105e-06,1.49818e-09,-1.07109e-13,-21527.5,-64.3348], Tmin=(972.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-139.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=[C]C(C)C(C)[O](19568)',
    structure = SMILES('C=[C]C(C)C(C)[O]'),
    E0 = (223.678,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,528.76,546.798],'cm^-1')),
        HinderedRotor(inertia=(0.104883,'amu*angstrom^2'), symmetry=1, barrier=(2.41146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227495,'amu*angstrom^2'), symmetry=1, barrier=(12.5088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0637133,'amu*angstrom^2'), symmetry=1, barrier=(12.4928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0599492,'amu*angstrom^2'), symmetry=1, barrier=(12.5012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.463893,0.0686477,-4.96367e-05,1.85833e-08,-2.8272e-12,27037.1,29.6678], Tmin=(100,'K'), Tmax=(1539.81,'K')), NASAPolynomial(coeffs=[15.6202,0.0292757,-1.12826e-05,1.97779e-09,-1.31166e-13,22369.5,-49.995], Tmin=(1539.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C1C(C)OC1C(24642)',
    structure = SMILES('C=C1C(C)OC1C'),
    E0 = (-75.0725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44096,0.0318876,7.47466e-05,-1.20491e-07,4.8084e-11,-8915.04,22.9308], Tmin=(100,'K'), Tmax=(969.541,'K')), NASAPolynomial(coeffs=[17.8763,0.0249764,-8.77366e-06,1.72019e-09,-1.3285e-13,-14964.1,-70.6118], Tmin=(969.541,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.0725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
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
    label = 'C=[C]C(C)[O](3245)',
    structure = SMILES('C=[C]C(C)[O]'),
    E0 = (278.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,461.597,461.616],'cm^-1')),
        HinderedRotor(inertia=(0.0767265,'amu*angstrom^2'), symmetry=1, barrier=(11.6024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0767264,'amu*angstrom^2'), symmetry=1, barrier=(11.6031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05953,0.035711,-1.16803e-05,-8.94971e-09,5.51275e-12,33616.7,20.5932], Tmin=(100,'K'), Tmax=(1043.68,'K')), NASAPolynomial(coeffs=[10.303,0.0183923,-7.30569e-06,1.36087e-09,-9.60942e-14,31118.5,-23.2536], Tmin=(1043.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S)"""),
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
    E0 = (117.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (206.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (194.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (163.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (238.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (335.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (242.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (256.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (320.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (365.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (259.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (169.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (267.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (206.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (518.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (348.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (243.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (180.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (574.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (573.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (389.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (125.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (444.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (624.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (303.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (298.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (326.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (371.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (277.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (536.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (243.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (180.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (142.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (318.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (125.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (622.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['CH3CHO(37)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['[CH2]C1([CH]C)OC1C(24643)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(88.8598,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[CH2]C(=CC)C(C)=O(24669)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2826 used for CO-CdCs_O;HJ
Exact match found for rate rule [CO-CdCs_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH3(17)', 'C=C([CH]C)C=O(18210)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.61258,'m^3/(mol*s)'), n=1.485, Ea=(32.0285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-HHH] for rate rule [CO-CdH_O;CsJ-HHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH3CHO(37)', 'C=[C][CH]C(18176)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C[CH][O](176)', 'CH3CHCCH2(18175)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['[CH2]C(=CC)[C](C)O(24670)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['[CH2]C(O)C(=C)[CH]C(24161)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['CC=C(C)[C](C)[O](24671)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.4947e+07,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[C]=C(C)C(C)[O](24672)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([O])C(C)=CC(24673)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(32400,'s^-1'), n=2.04, Ea=(82.4248,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 89 used for R4H_SS(Cd)S;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SS(Cd)S;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(=[C]C)C(C)O(24674)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[C](C)C(C)[O](19565)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['[CH2][C](C=C)C(C)O(19569)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.45388e+06,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_SSMS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C[CH][O](176)', 'C=[C][CH]C(18176)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['CC([O])[C]1CC1C(24675)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['[CH2][C]1C(C)OC1C(24676)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['CC=C(C)C(C)=O(24677)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH2(S)(23)', 'C=C([CH]C)C[O](24151)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['CH2(S)(23)', '[CH2]C(=C)C(C)[O](24678)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC=[C]CC(C)[O](19705)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['CC=C1COC1C(24668)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O(4)', '[CH2]C([CH]C)=CC(24211)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4171.1,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(19)', 'CC=[C]C(C)[O](24679)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', 'C=CC(=C)C(C)[O](19563)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]CC(=C)C(C)[O](24680)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['[CH2][C](CC)C(C)=O(3039)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.18083e+09,'s^-1'), n=1.04667, Ea=(209.2,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C(CC)C(C)[O](24681)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C([O])C(=C)CC(24682)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R4H_SS(Cd)S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['[CH]=C([CH]C)C(C)O(24683)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.4207e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_singleH] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['C[CH][C]1COC1C(24684)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['C=C(CC)C(C)=O(3043)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['C=CC(=C)C(C)O(19576)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C]C(C)C(C)[O](19568)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C([CH]C)C(C)[O](24158)'],
    products = ['C=C1C(C)OC1C(24642)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CHCH3(T)(95)', 'C=[C]C(C)[O](3245)'],
    products = ['C=C([CH]C)C(C)[O](24158)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4205',
    isomers = [
        'C=C([CH]C)C(C)[O](24158)',
    ],
    reactants = [
        ('CH3CHO(37)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4205',
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

