species(
    label = '[CH2]OC([CH2])(O)C(O)=CC(5898)',
    structure = SMILES('[CH2]OC([CH2])(O)C(O)=CC'),
    E0 = (-233.751,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1600,1649.82,2843.12,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.70503,0.126885,-0.000146988,8.02368e-08,-1.63792e-11,-27854.7,39.4027], Tmin=(100,'K'), Tmax=(1336.11,'K')), NASAPolynomial(coeffs=[32.4133,0.00891226,-1.33163e-07,-2.29055e-10,2.22688e-14,-36093.3,-135.911], Tmin=(1336.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CsJOCH3) + radical(CJC(C)OC)"""),
)

species(
    label = 'CH2O(15)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]OC1(O)CC1(O)[CH]C(32228)',
    structure = SMILES('[CH2]OC1(O)CC1(O)[CH]C'),
    E0 = (-172.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10342,0.0868699,-2.07902e-05,-6.21261e-08,3.78884e-11,-20486.8,35.9628], Tmin=(100,'K'), Tmax=(921.718,'K')), NASAPolynomial(coeffs=[32.597,0.00683009,1.71533e-06,-4.69236e-10,2.65729e-14,-29511.8,-139.132], Tmin=(921.718,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + ring(Cyclopropane) + radical(CsJOCH3) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1(O)OCC1(O)[CH]C(32229)',
    structure = SMILES('[CH2]C1(O)OCC1(O)[CH]C'),
    E0 = (-172.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.42343,0.117811,-0.000123558,6.18534e-08,-1.11874e-11,-20441.8,42.0359], Tmin=(100,'K'), Tmax=(1665.8,'K')), NASAPolynomial(coeffs=[27.6047,0.00976697,3.93213e-06,-1.25512e-09,9.76292e-14,-26126,-109.524], Tmin=(1665.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CJC(C)OC) + radical(CCJCO)"""),
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
    label = '[CH2]OC(=C)O(3587)',
    structure = SMILES('[CH2]OC(=C)O'),
    E0 = (-113.088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3615,1277.5,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,385.624,385.986],'cm^-1')),
        HinderedRotor(inertia=(0.149807,'amu*angstrom^2'), symmetry=1, barrier=(15.8087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149921,'amu*angstrom^2'), symmetry=1, barrier=(15.8066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149522,'amu*angstrom^2'), symmetry=1, barrier=(15.8058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930845,0.0518764,-5.23343e-05,2.50764e-08,-4.45103e-12,-13477.5,20.1549], Tmin=(100,'K'), Tmax=(1612.24,'K')), NASAPolynomial(coeffs=[16.0622,0.00502412,-8.08492e-08,-1.12707e-10,1.03568e-14,-17146.4,-56.3195], Tmin=(1612.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-113.088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=COCJ)"""),
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
    label = '[CH2]OC(=C)C(O)=CC(32230)',
    structure = SMILES('[CH2]OC(=C)C(O)=CC'),
    E0 = (-159.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1086,0.0952661,-9.71941e-05,4.86991e-08,-9.30675e-12,-18948.7,30.7499], Tmin=(100,'K'), Tmax=(1393.23,'K')), NASAPolynomial(coeffs=[24.3204,0.0154245,-3.87597e-06,5.25264e-10,-3.07111e-14,-25371.1,-97.983], Tmin=(1393.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-159.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=COCJ)"""),
)

species(
    label = '[CH2]OC(C)([O])C(O)=CC(32231)',
    structure = SMILES('[CH2]OC(C)([O])C(O)=CC'),
    E0 = (-211.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,528.15,613.411,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160244,'amu*angstrom^2'), symmetry=1, barrier=(3.68432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160244,'amu*angstrom^2'), symmetry=1, barrier=(3.68432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160244,'amu*angstrom^2'), symmetry=1, barrier=(3.68432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160244,'amu*angstrom^2'), symmetry=1, barrier=(3.68432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160244,'amu*angstrom^2'), symmetry=1, barrier=(3.68432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160244,'amu*angstrom^2'), symmetry=1, barrier=(3.68432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8319,0.113065,-0.000123566,6.64674e-08,-1.36578e-11,-25180.3,37.3978], Tmin=(100,'K'), Tmax=(1271.12,'K')), NASAPolynomial(coeffs=[26.7598,0.0176048,-4.44152e-06,5.93945e-10,-3.40443e-14,-32005.8,-105.656], Tmin=(1271.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-211.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)OJ) + radical(CsJOCH3)"""),
)

species(
    label = '[CH2]C([O])(OC)C(O)=CC(32232)',
    structure = SMILES('[CH2]C([O])(OC)C(O)=CC'),
    E0 = (-188.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,455.473,689.38,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162194,'amu*angstrom^2'), symmetry=1, barrier=(3.72917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162194,'amu*angstrom^2'), symmetry=1, barrier=(3.72917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162194,'amu*angstrom^2'), symmetry=1, barrier=(3.72917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162194,'amu*angstrom^2'), symmetry=1, barrier=(3.72917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162194,'amu*angstrom^2'), symmetry=1, barrier=(3.72917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162194,'amu*angstrom^2'), symmetry=1, barrier=(3.72917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6797,0.119676,-0.000145815,8.91083e-08,-2.10421e-11,-22471.4,35.8128], Tmin=(100,'K'), Tmax=(1047.83,'K')), NASAPolynomial(coeffs=[23.5885,0.0232156,-7.72687e-06,1.25017e-09,-7.98968e-14,-27766.6,-87.272], Tmin=(1047.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)OJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]OC(C)(O)C(=O)[CH]C(5888)',
    structure = SMILES('[CH2]OC(C)(O)C([O])=CC'),
    E0 = (-306.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.1431,0.111811,-0.000116045,5.80873e-08,-1.09555e-11,-36617.1,39.1468], Tmin=(100,'K'), Tmax=(1452.31,'K')), NASAPolynomial(coeffs=[29.0454,0.0138526,-2.41526e-06,2.0993e-10,-8.38331e-15,-44404.5,-118.58], Tmin=(1452.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-306.457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CsJOCH3) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]OC(C)(O)C(O)=[C]C(32233)',
    structure = SMILES('[CH2]OC(C)(O)C(O)=[C]C'),
    E0 = (-206.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180,1600,1647.95,2844.48,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152437,'amu*angstrom^2'), symmetry=1, barrier=(3.50484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152437,'amu*angstrom^2'), symmetry=1, barrier=(3.50484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152437,'amu*angstrom^2'), symmetry=1, barrier=(3.50484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152437,'amu*angstrom^2'), symmetry=1, barrier=(3.50484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152437,'amu*angstrom^2'), symmetry=1, barrier=(3.50484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152437,'amu*angstrom^2'), symmetry=1, barrier=(3.50484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152437,'amu*angstrom^2'), symmetry=1, barrier=(3.50484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.48861,0.119835,-0.000131758,6.87271e-08,-1.3428e-11,-24573.3,39.3537], Tmin=(100,'K'), Tmax=(1403.25,'K')), NASAPolynomial(coeffs=[31.2858,0.0104959,-9.15755e-07,-6.88399e-11,1.05386e-14,-32765.9,-130.448], Tmin=(1403.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-206.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CsJOCH3)"""),
)

species(
    label = '[CH2]C(O)(OC)C(=O)[CH]C(5889)',
    structure = SMILES('[CH2]C(O)(OC)C([O])=CC'),
    E0 = (-283.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57128,0.113619,-0.000122129,6.05896e-08,-1.008e-11,-33926.7,36.0461], Tmin=(100,'K'), Tmax=(953.156,'K')), NASAPolynomial(coeffs=[25.2076,0.020665,-6.41565e-06,1.03819e-09,-6.86689e-14,-39914,-96.4908], Tmin=(953.156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-283.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C(O)(OC)C(O)=[C]C(32234)',
    structure = SMILES('[CH2]C(O)(OC)C(O)=[C]C'),
    E0 = (-183.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180,1600,1707.92,2791.7,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154521,'amu*angstrom^2'), symmetry=1, barrier=(3.55274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154521,'amu*angstrom^2'), symmetry=1, barrier=(3.55274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154521,'amu*angstrom^2'), symmetry=1, barrier=(3.55274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154521,'amu*angstrom^2'), symmetry=1, barrier=(3.55274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154521,'amu*angstrom^2'), symmetry=1, barrier=(3.55274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154521,'amu*angstrom^2'), symmetry=1, barrier=(3.55274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154521,'amu*angstrom^2'), symmetry=1, barrier=(3.55274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.943,0.121856,-0.000138157,7.10478e-08,-1.2228e-11,-21881.6,36.3537], Tmin=(100,'K'), Tmax=(935.554,'K')), NASAPolynomial(coeffs=[27.6492,0.0169893,-4.74258e-06,7.20327e-10,-4.66281e-14,-28366.3,-109.505], Tmin=(935.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CJC(C)OC) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C(O)C(C)(O)O[CH2](32235)',
    structure = SMILES('[CH2]C=C(O)C(C)(O)O[CH2]'),
    E0 = (-292.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.2405,0.121624,-0.000127922,6.24041e-08,-1.11837e-11,-34918.1,41.8847], Tmin=(100,'K'), Tmax=(1606.49,'K')), NASAPolynomial(coeffs=[33.6326,0.00620614,1.88748e-06,-6.11954e-10,4.62582e-14,-43719.1,-144.005], Tmin=(1606.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-292.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CsJOCH3) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C(O)C([CH2])(O)OC(32236)',
    structure = SMILES('[CH2]C=C(O)C([CH2])(O)OC'),
    E0 = (-270.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.78362,0.125006,-0.000140384,7.41687e-08,-1.46045e-11,-32223.3,39.1795], Tmin=(100,'K'), Tmax=(1404.57,'K')), NASAPolynomial(coeffs=[32.7341,0.00860128,2.23007e-07,-3.02177e-10,2.7042e-14,-40695.9,-138.882], Tmin=(1404.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-270.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CJC(C)OC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][O](167)',
    structure = SMILES('[CH2][O]'),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88409,-0.00363885,3.28543e-05,-4.13611e-08,1.59631e-11,23210.8,7.47983], Tmin=(100,'K'), Tmax=(933.06,'K')), NASAPolynomial(coeffs=[6.69335,0.000289989,8.61416e-07,-1.56351e-10,7.33778e-15,21991.3,-9.6043], Tmin=(933.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(H3COJ) + radical(CsJOH)"""),
)

species(
    label = '[CH2]OC1(O)CC(C)[C]1O(32237)',
    structure = SMILES('[CH2]OC1(O)CC(C)[C]1O'),
    E0 = (-190.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.870401,0.0842029,-2.21366e-05,-5.15043e-08,3.13989e-11,-22753,34.0202], Tmin=(100,'K'), Tmax=(936.744,'K')), NASAPolynomial(coeffs=[29.1036,0.014044,-2.39923e-06,3.56359e-10,-3.12861e-14,-30905.9,-122.172], Tmin=(936.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-190.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + ring(Cyclobutane) + radical(CsJOCH3) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1(O)OCC(C)[C]1O(32238)',
    structure = SMILES('[CH2]C1(O)OCC(C)[C]1O'),
    E0 = (-265.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.03506,0.111114,-0.000110947,5.39016e-08,-9.5385e-12,-31618,37.3318], Tmin=(100,'K'), Tmax=(1695.15,'K')), NASAPolynomial(coeffs=[25.5048,0.0138199,1.64814e-06,-8.02501e-10,6.6354e-14,-36990.7,-102.727], Tmin=(1695.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-265.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(C2CsJOH) + radical(CJC(C)OC)"""),
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
    label = '[CH2]OC([CH2])(O)C(=C)O(32239)',
    structure = SMILES('[CH2]OC([CH2])(O)C(=C)O'),
    E0 = (-197.726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,1600,1616.38,2871.27,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150531,'amu*angstrom^2'), symmetry=1, barrier=(3.46102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150531,'amu*angstrom^2'), symmetry=1, barrier=(3.46102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150531,'amu*angstrom^2'), symmetry=1, barrier=(3.46102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150531,'amu*angstrom^2'), symmetry=1, barrier=(3.46102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150531,'amu*angstrom^2'), symmetry=1, barrier=(3.46102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150531,'amu*angstrom^2'), symmetry=1, barrier=(3.46102,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.26824,0.113784,-0.000135115,7.31831e-08,-1.45426e-11,-23534.5,36.0287], Tmin=(100,'K'), Tmax=(1425.2,'K')), NASAPolynomial(coeffs=[31.4505,0.000505809,3.72952e-06,-9.42443e-10,7.00647e-14,-31252.4,-131.949], Tmin=(1425.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-197.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsJOCH3) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]O[C](O)CC(O)=CC(32240)',
    structure = SMILES('[CH2]O[C](O)CC(O)=CC'),
    E0 = (-213.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.751541,0.102779,-0.000109334,6.0316e-08,-1.32587e-11,-25481.6,37.6781], Tmin=(100,'K'), Tmax=(1103.94,'K')), NASAPolynomial(coeffs=[18.675,0.0323883,-1.36885e-05,2.55499e-09,-1.77853e-13,-29770.7,-57.9647], Tmin=(1103.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Cs_P) + radical(CsJOCH3)"""),
)

species(
    label = 'CC=C(O)C1(O)CCO1(31893)',
    structure = SMILES('CC=C(O)C1(O)CCO1'),
    E0 = (-492.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.688661,0.0758,1.48999e-05,-1.06056e-07,5.73465e-11,-59054.5,29.9502], Tmin=(100,'K'), Tmax=(884.566,'K')), NASAPolynomial(coeffs=[32.1826,0.0039247,6.60287e-06,-1.69152e-09,1.2172e-13,-67873.2,-141.579], Tmin=(884.566,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-492.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Oxetane)"""),
)

species(
    label = '[CH2]OC([CH2])(O)C(=O)CC(5498)',
    structure = SMILES('[CH2]OC([CH2])(O)C(=O)CC'),
    E0 = (-247.331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,470.296,666.997,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4372.41,'J/mol'), sigma=(7.2676,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=682.96 K, Pc=25.85 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59674,0.116279,-0.000132833,7.55921e-08,-1.6688e-11,-29539.6,37.3463], Tmin=(100,'K'), Tmax=(1116.22,'K')), NASAPolynomial(coeffs=[24.1455,0.0240295,-8.86446e-06,1.55052e-09,-1.0467e-13,-35286.3,-89.6749], Tmin=(1116.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-247.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-OdCsCs) + radical(CJC(C)OC) + radical(CsJOCH3)"""),
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
    label = '[CH2]C([O])(O)C(O)=CC(31647)',
    structure = SMILES('[CH2]C([O])(O)C(O)=CC'),
    E0 = (-206.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57438,0.103411,-0.000118869,6.43734e-08,-1.29463e-11,-24595.2,35.268], Tmin=(100,'K'), Tmax=(1387.58,'K')), NASAPolynomial(coeffs=[27.1597,0.00660338,8.89175e-07,-4.23158e-10,3.58034e-14,-31224,-107.922], Tmin=(1387.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-206.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)OJ) + radical(C=CC(O)2CJ)"""),
)

species(
    label = '[CH2]O[C](O)C(O)=CC(32241)',
    structure = SMILES('[CH2]O[C](O)C(O)=CC'),
    E0 = (-190.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.068858,0.0904638,-0.00010519,6.35191e-08,-1.52476e-11,-22784.8,32.3188], Tmin=(100,'K'), Tmax=(1015.19,'K')), NASAPolynomial(coeffs=[15.9616,0.0273021,-1.18666e-05,2.23509e-09,-1.5605e-13,-26039.6,-45.261], Tmin=(1015.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-190.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Cs_P) + radical(CsJOCH3)"""),
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
    E0 = (-233.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-172.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-172.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-213.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-32.4722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-118.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-72.5371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-91.5072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-150.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-162.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-200.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-150.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-178.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-178.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (70.0577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-108.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-147.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (222.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (32.2995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-225.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-105.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (175.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (190.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    products = ['CH2O(15)', 'C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    products = ['[CH2]OC1(O)CC1(O)[CH]C(32228)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(61.6907,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra_cs2H] for rate rule [R4_S_D;doublebond_intra_HNd_secNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 57.6 to 61.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    products = ['[CH2]C1(O)OCC1(O)[CH]C(32229)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.65e+06,'s^-1'), n=1.3, Ea=(61.2392,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_HNd;radadd_intra_cs2H] for rate rule [R5_SS_D;doublebond_intra_HNd_secNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 57.9 to 61.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2O(15)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2330,'cm^3/(mol*s)'), n=3.17, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-HH;YJ] for rate rule [Od_CO-HH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]OC(=C)O(3587)', 'CC=[C]O(31159)'],
    products = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00712612,'m^3/(mol*s)'), n=2.40979, Ea=(7.81798,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CdsJ] for rate rule [Cds-OsOs_Cds;CdsJ-O2s]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OH(5)', '[CH2]OC(=C)C(O)=CC(32230)'],
    products = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.812049,'m^3/(mol*s)'), n=2.01336, Ea=(12.3541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;YJ] for rate rule [Cds-CdOs_Cds;OJ_pri]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]OC(C)([O])C(O)=CC(32231)'],
    products = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C([O])(OC)C(O)=CC(32232)'],
    products = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.3e+09,'s^-1','*|/',2.51189), n=0.75, Ea=(97.0688,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R4H_SSS;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    products = ['[CH2]OC(C)(O)C(=O)[CH]C(5888)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00963743,'s^-1'), n=3.795, Ea=(83.0524,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_2H;O_H_out] + [R4H_SS(Cd)S;C_rad_out_2H;XH_out] for rate rule [R4H_SS(Cd)S;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]OC(C)(O)C(O)=[C]C(32233)'],
    products = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    products = ['[CH2]C(O)(OC)C(=O)[CH]C(5889)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.06724e+06,'s^-1'), n=1.18977, Ea=(32.9114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_2H;XH_out] for rate rule [R5H_SSSS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)(OC)C(O)=[C]C(32234)'],
    products = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    products = ['[CH2]C=C(O)C(C)(O)O[CH2](32235)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(121000,'s^-1'), n=1.9, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 92 used for R5H_SSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    products = ['[CH2]C=C(O)C([CH2])(O)OC(32236)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8010,'s^-1'), n=1.94, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 93 used for R6H_RSSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R6H_RSSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][O](167)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    products = ['[CH2]OC1(O)CC(C)[C]1O(32237)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.82e+08,'s^-1'), n=0.91, Ea=(125.52,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra_secNd;radadd_intra_cs2H] for rate rule [R4_Cs_RR_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    products = ['[CH2]C1(O)OCC(C)[C]1O(32238)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.95277e+09,'s^-1'), n=0.6, Ea=(85.772,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_secNd;radadd_intra_cs2H] for rate rule [R5_SS_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH2(S)(23)', '[CH2]OC([CH2])(O)C(=C)O(32239)'],
    products = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.97e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]O[C](O)CC(O)=CC(32240)'],
    products = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    products = ['CC=C(O)C1(O)CCO1(31893)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1713.16,'s^-1'), n=2.9, Ea=(127.847,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_Cs;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_Cs;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(19)', '[CH2]C([O])(O)C(O)=CC(31647)'],
    products = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(19)', '[CH2]O[C](O)C(O)=CC(32241)'],
    products = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '5235',
    isomers = [
        '[CH2]OC([CH2])(O)C(O)=CC(5898)',
    ],
    reactants = [
        ('CH2O(15)', 'C=C(O)C(O)=CC(5562)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '5235',
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

