species(
    label = '[CH2]OC(C)[C](O)C(=C)O(31890)',
    structure = SMILES('[CH2]OC(C)C(O)=C([CH2])O'),
    E0 = (-291.224,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42888,0.119365,-0.000131624,6.91502e-08,-1.36217e-11,-34775.8,39.7868], Tmin=(100,'K'), Tmax=(1390.47,'K')), NASAPolynomial(coeffs=[30.7325,0.0111721,-1.10403e-06,-4.70016e-11,9.6527e-15,-42760.7,-126.68], Tmin=(1390.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-291.224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(CsJOCC2)"""),
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
    label = '[CH2][C](O)C1(O)COC1C(32437)',
    structure = SMILES('[CH2][C](O)C1(O)COC1C'),
    E0 = (-155.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.717674,0.0998104,-9.42504e-05,3.51072e-08,6.84705e-13,-18469.5,32.6567], Tmin=(100,'K'), Tmax=(820.247,'K')), NASAPolynomial(coeffs=[19.134,0.0276109,-7.2202e-06,9.48398e-10,-5.20523e-14,-22554,-64.2287], Tmin=(820.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-155.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOH) + radical(CJCO)"""),
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
    label = '[CH2]OC(C)C(O)=C=C(32438)',
    structure = SMILES('[CH2]OC(C)C(O)=C=C'),
    E0 = (-58.9776,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.615914,0.0910937,-8.48472e-05,3.48133e-08,-4.07209e-12,-6918.03,31.5841], Tmin=(100,'K'), Tmax=(1015.12,'K')), NASAPolynomial(coeffs=[21.7262,0.0209039,-7.50267e-06,1.33769e-09,-9.31899e-14,-12373.6,-81.068], Tmin=(1015.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.9776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]C(O)=C(O)[C](C)OC(32439)',
    structure = SMILES('[CH2]C(O)=C(O)[C](C)OC'),
    E0 = (-294.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.43529,0.114272,-0.000134369,7.97522e-08,-1.83786e-11,-35172.3,36.747], Tmin=(100,'K'), Tmax=(1070.67,'K')), NASAPolynomial(coeffs=[22.5878,0.0245199,-8.62451e-06,1.45386e-09,-9.55752e-14,-40316.4,-80.7906], Tmin=(1070.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-294.105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH2]OC(C)C(O)=C(C)[O](32440)',
    structure = SMILES('[CH2]OC(C)C(O)=C(C)[O]'),
    E0 = (-312.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29948,0.106342,-0.00010651,4.93511e-08,-7.6049e-12,-37365.5,36.5414], Tmin=(100,'K'), Tmax=(995.496,'K')), NASAPolynomial(coeffs=[23.9281,0.0229046,-7.80228e-06,1.34021e-09,-9.12898e-14,-43276.7,-89.5156], Tmin=(995.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]C(O)=C(O)C([CH2])OC(32441)',
    structure = SMILES('[CH2]C(O)=C(O)C([CH2])OC'),
    E0 = (-264.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69873,0.116448,-0.000127589,6.35699e-08,-1.03948e-11,-31574.5,37.4989], Tmin=(100,'K'), Tmax=(940.485,'K')), NASAPolynomial(coeffs=[26.1958,0.0188687,-5.54739e-06,8.69802e-10,-5.69168e-14,-37752.7,-100.316], Tmin=(940.485,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-264.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]OC(C)C([O])=C(C)O(32442)',
    structure = SMILES('[CH2]OC(C)C([O])=C(C)O'),
    E0 = (-312.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29948,0.106342,-0.00010651,4.93511e-08,-7.6049e-12,-37365.5,36.5414], Tmin=(100,'K'), Tmax=(995.496,'K')), NASAPolynomial(coeffs=[23.9281,0.0229046,-7.80228e-06,1.34021e-09,-9.12898e-14,-43276.7,-89.5156], Tmin=(995.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]O[C](C)C(O)=C(C)O(32443)',
    structure = SMILES('[CH2]O[C](C)C(O)=C(C)O'),
    E0 = (-269.437,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,360,370,350,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.65484,0.119707,-0.000144959,8.74047e-08,-2.03924e-11,-32198,36.2637], Tmin=(100,'K'), Tmax=(1058.67,'K')), NASAPolynomial(coeffs=[23.8441,0.0233685,-8.4655e-06,1.4563e-09,-9.70402e-14,-37597.3,-88.2091], Tmin=(1058.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-269.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C2CsJOCs) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]C(O)=C([O])C(C)OC(32444)',
    structure = SMILES('[CH2]C(O)=C([O])C(C)OC'),
    E0 = (-337.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3566,0.10409,-0.00010665,5.49983e-08,-1.09845e-11,-40327.7,38.0222], Tmin=(100,'K'), Tmax=(1257.81,'K')), NASAPolynomial(coeffs=[24.0127,0.0218546,-6.72404e-06,1.05104e-09,-6.63859e-14,-46586.5,-89.6996], Tmin=(1257.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-337.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]OC([CH2])C(O)=C(C)O(32445)',
    structure = SMILES('[CH2]OC([CH2])C(O)=C(C)O'),
    E0 = (-239.63,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98239,0.122661,-0.000140997,7.50139e-08,-1.40853e-11,-28597.4,37.2441], Tmin=(100,'K'), Tmax=(955.611,'K')), NASAPolynomial(coeffs=[27.6581,0.0173647,-5.18406e-06,8.23815e-10,-5.43562e-14,-35119.6,-108.893], Tmin=(955.611,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-239.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CsJOCC2) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C([O])=C(O)C(C)OC(32446)',
    structure = SMILES('[CH2]C([O])=C(O)C(C)OC'),
    E0 = (-337.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3566,0.10409,-0.00010665,5.49983e-08,-1.09845e-11,-40327.7,38.0222], Tmin=(100,'K'), Tmax=(1257.81,'K')), NASAPolynomial(coeffs=[24.0127,0.0218546,-6.72404e-06,1.05104e-09,-6.63859e-14,-46586.5,-89.6996], Tmin=(1257.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-337.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]O[CH]C(567)',
    structure = SMILES('[CH2]O[CH]C'),
    E0 = (132.786,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1071.68,3135.02],'cm^-1')),
        HinderedRotor(inertia=(0.569812,'amu*angstrom^2'), symmetry=1, barrier=(13.5453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.589052,'amu*angstrom^2'), symmetry=1, barrier=(13.5435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.589118,'amu*angstrom^2'), symmetry=1, barrier=(13.545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07912,0.0410978,-4.15576e-05,2.34499e-08,-5.32469e-12,16040.6,15.834], Tmin=(100,'K'), Tmax=(1069.75,'K')), NASAPolynomial(coeffs=[8.97105,0.0153272,-5.42175e-06,9.29844e-10,-6.16823e-14,14566.1,-17.8801], Tmin=(1069.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CsJOCC) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C(O)=[C]O(31700)',
    structure = SMILES('[CH2]C(O)=[C]O'),
    E0 = (-18.5667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.841623,'amu*angstrom^2'), symmetry=1, barrier=(19.3506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842368,'amu*angstrom^2'), symmetry=1, barrier=(19.3677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.841193,'amu*angstrom^2'), symmetry=1, barrier=(19.3407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.351247,0.0579095,-6.81787e-05,3.58442e-08,-6.68351e-12,-2082.38,24.0014], Tmin=(100,'K'), Tmax=(1624.88,'K')), NASAPolynomial(coeffs=[17.4373,-0.00264896,4.80126e-06,-1.10446e-09,7.9234e-14,-5193.08,-59.2091], Tmin=(1624.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.5667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(O)CJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]OC(C)C1(O)C[C]1O(32447)',
    structure = SMILES('[CH2]OC(C)C1(O)C[C]1O'),
    E0 = (-157.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09815,0.108586,-0.000120829,6.91298e-08,-1.55687e-11,-18796.1,34.7232], Tmin=(100,'K'), Tmax=(1086.75,'K')), NASAPolynomial(coeffs=[20.4213,0.0293793,-1.15031e-05,2.0638e-09,-1.40617e-13,-23473.4,-70.886], Tmin=(1086.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + ring(Cyclopropane) + radical(CsJOCC2) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1(O)COC(C)[C]1O(32448)',
    structure = SMILES('[CH2]C1(O)COC(C)[C]1O'),
    E0 = (-233.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13987,0.100986,-0.000104393,5.74018e-08,-1.22087e-11,-27833.6,31.1535], Tmin=(100,'K'), Tmax=(1273.18,'K')), NASAPolynomial(coeffs=[19.7738,0.0265918,-6.50791e-06,7.86347e-10,-3.91977e-14,-32454.7,-72.0281], Tmin=(1273.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(C2CsJOH) + radical(CJC(C)2O)"""),
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
    label = '[CH2]OC[C](O)C(=C)O(32263)',
    structure = SMILES('[CH2]OCC(O)=C([CH2])O'),
    E0 = (-244.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,325,375,415,465,420,450,1700,1750,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10028,0.0957872,-0.000106641,5.73741e-08,-1.15814e-11,-29235.4,34.8607], Tmin=(100,'K'), Tmax=(1365.41,'K')), NASAPolynomial(coeffs=[24.2283,0.010625,-1.04184e-06,-6.45313e-11,1.18343e-14,-35130.4,-91.4818], Tmin=(1365.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CsJOCC) + radical(C=C(O)CJ)"""),
)

species(
    label = 'CC1OCCC(O)=C1O(31891)',
    structure = SMILES('CC1OCCC(O)=C1O'),
    E0 = (-591.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.359394,0.0678618,2.88859e-05,-1.09503e-07,5.43136e-11,-70979.5,26.9554], Tmin=(100,'K'), Tmax=(911.761,'K')), NASAPolynomial(coeffs=[29.5899,0.0104597,1.59778e-06,-5.46864e-10,3.42918e-14,-79516.2,-131.631], Tmin=(911.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-591.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(36dihydro2hpyran)"""),
)

species(
    label = '[CH2]OC(C)C(=O)C([CH2])O(32449)',
    structure = SMILES('[CH2]OC(C)C(=O)C([CH2])O'),
    E0 = (-207.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.89391,0.0976901,-9.18281e-05,4.3336e-08,-8.06362e-12,-24766.5,40.9001], Tmin=(100,'K'), Tmax=(1303.37,'K')), NASAPolynomial(coeffs=[21.974,0.0275091,-1.10592e-05,2.0231e-09,-1.39355e-13,-30727.5,-75.4831], Tmin=(1303.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-207.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-OdCsCs) + radical(CsJOCC2) + radical(CJCO)"""),
)

species(
    label = '[CH2]OC(C)C(O)C([CH2])=O(32450)',
    structure = SMILES('[CH2]OC(C)C(O)C(=C)[O]'),
    E0 = (-255.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.944082,0.103394,-0.000108078,5.80177e-08,-1.22954e-11,-30596.3,38.0074], Tmin=(100,'K'), Tmax=(1150.7,'K')), NASAPolynomial(coeffs=[20.3634,0.0293263,-1.15264e-05,2.07939e-09,-1.42282e-13,-35500,-67.7796], Tmin=(1150.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsJOCC2) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C(O)=C(O)C(C)[O](31703)',
    structure = SMILES('[CH2]C(O)=C(O)C(C)[O]'),
    E0 = (-265.194,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,323.209,323.211],'cm^-1')),
        HinderedRotor(inertia=(0.211459,'amu*angstrom^2'), symmetry=1, barrier=(15.6755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211453,'amu*angstrom^2'), symmetry=1, barrier=(15.6755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211456,'amu*angstrom^2'), symmetry=1, barrier=(15.6755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211456,'amu*angstrom^2'), symmetry=1, barrier=(15.6755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211446,'amu*angstrom^2'), symmetry=1, barrier=(15.6755,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.83479,0.0991158,-0.00010535,5.22413e-08,-9.5195e-12,-31660.3,36.7229], Tmin=(100,'K'), Tmax=(1580.23,'K')), NASAPolynomial(coeffs=[27.5949,0.00551477,1.63534e-06,-5.44642e-10,4.19531e-14,-38575.8,-111.176], Tmin=(1580.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-265.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]OC(C)C(O)=[C]O(32451)',
    structure = SMILES('[CH2]OC(C)C(O)=[C]O'),
    E0 = (-168.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1685,370,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47687,0.100916,-0.000112809,5.94443e-08,-1.16836e-11,-20064.5,36.8043], Tmin=(100,'K'), Tmax=(1405.77,'K')), NASAPolynomial(coeffs=[27.1424,0.00732768,2.17979e-08,-2.14511e-10,1.98141e-14,-26910,-106.743], Tmin=(1405.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(CsJOCC2) + radical(C=CJO)"""),
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
    label = '[CH2]OC(C)C(=O)C(=C)O(32452)',
    structure = SMILES('[CH2]OC(C)C(=O)C(=C)O'),
    E0 = (-324.019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.916109,0.0989098,-9.80773e-05,4.82214e-08,-9.28623e-12,-38785.7,34.625], Tmin=(100,'K'), Tmax=(1265.63,'K')), NASAPolynomial(coeffs=[22.5797,0.0246533,-1.00717e-05,1.86563e-09,-1.29759e-13,-44733.2,-84.2639], Tmin=(1265.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-324.019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]OC(C)=C(O)C(=C)O(32453)',
    structure = SMILES('[CH2]OC(C)=C(O)C(=C)O'),
    E0 = (-373.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,350,375,415,440,465,420,435,450,1700,1725,1750,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.38181,0.117288,-0.000131693,6.96113e-08,-1.36976e-11,-44703,35.3714], Tmin=(100,'K'), Tmax=(1411.04,'K')), NASAPolynomial(coeffs=[30.8745,0.00801803,4.07927e-07,-3.33147e-10,2.91528e-14,-52595.4,-131.232], Tmin=(1411.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-373.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=COCJ)"""),
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
    label = '[CH2]OC=C(O)C(=C)O(32454)',
    structure = SMILES('[CH2]OC=C(O)C(=C)O'),
    E0 = (-331.966,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52314,0.107422,-0.000122714,6.29643e-08,-1.16359e-11,-39660.3,33.8879], Tmin=(100,'K'), Tmax=(1597.43,'K')), NASAPolynomial(coeffs=[30.5916,-0.00260652,6.06001e-06,-1.40134e-09,1.00053e-13,-46781.2,-130.556], Tmin=(1597.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-331.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COCJ)"""),
)

species(
    label = '[CH2]OC(C)C([O])C(=C)O(32455)',
    structure = SMILES('[CH2]OC(C)C([O])C(=C)O'),
    E0 = (-163.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.54783,0.108019,-0.000110522,5.58954e-08,-1.09392e-11,-19435.4,38.4854], Tmin=(100,'K'), Tmax=(1258.44,'K')), NASAPolynomial(coeffs=[25.6793,0.021477,-7.36762e-06,1.24853e-09,-8.31168e-14,-26288.1,-99.1283], Tmin=(1258.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-163.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsJOCC2) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]O[C](C)C(O)C(=C)O(32456)',
    structure = SMILES('[CH2]O[C](C)C(O)C(=C)O'),
    E0 = (-213.011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,360,370,350,350,440,435,1725,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20955,0.115632,-0.000142338,9.04256e-08,-2.26354e-11,-25432.5,37.4125], Tmin=(100,'K'), Tmax=(979.467,'K')), NASAPolynomial(coeffs=[19.4722,0.0311703,-1.29893e-05,2.38512e-09,-1.63811e-13,-29484,-61.9361], Tmin=(979.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C2CsJOCs) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]OC([CH2])C(O)C(=C)O(32457)',
    structure = SMILES('[CH2]OC([CH2])C(O)C(=C)O'),
    E0 = (-183.205,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67218,0.120252,-0.000144509,8.63301e-08,-1.99818e-11,-21826.2,38.8721], Tmin=(100,'K'), Tmax=(1066.02,'K')), NASAPolynomial(coeffs=[23.9924,0.0239503,-8.99978e-06,1.58417e-09,-1.07077e-13,-27297.9,-86.585], Tmin=(1066.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsJOCC2) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH]=C(O)C(O)C(C)O[CH2](32458)',
    structure = SMILES('[CH]=C(O)C(O)C(C)O[CH2]'),
    E0 = (-146.619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5718,0.11311,-0.000124994,6.83956e-08,-1.44824e-11,-17425.5,39.1616], Tmin=(100,'K'), Tmax=(1164.47,'K')), NASAPolynomial(coeffs=[24.8879,0.0222197,-7.91403e-06,1.36648e-09,-9.18933e-14,-23587.7,-92.5198], Tmin=(1164.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsJOCC2) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)[C](O)C(C)OC(32459)',
    structure = SMILES('[CH]C(O)=C(O)C(C)OC'),
    E0 = (-263.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31297,0.103087,-8.05516e-05,1.9084e-08,3.4265e-12,-31432.9,36.781], Tmin=(100,'K'), Tmax=(974.152,'K')), NASAPolynomial(coeffs=[24.37,0.0272671,-9.44053e-06,1.65036e-09,-1.14536e-13,-37843,-93.67], Tmin=(974.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'CC1OCC[C](O)[C]1O(32460)',
    structure = SMILES('CC1OCC[C](O)[C]1O'),
    E0 = (-272.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.141659,0.0855602,-6.46506e-05,1.74274e-08,2.5643e-12,-32592.2,29.1475], Tmin=(100,'K'), Tmax=(875.962,'K')), NASAPolynomial(coeffs=[15.1146,0.035109,-1.1162e-05,1.76115e-09,-1.11646e-13,-36002.2,-46.6428], Tmin=(875.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxane) + radical(C2CsJOH) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(O)C(O)=C(C)OC(32461)',
    structure = SMILES('C=C(O)C(O)=C(C)OC'),
    E0 = (-565.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.53693,0.119814,-0.00013138,6.84844e-08,-1.33358e-11,-67789.9,33.4447], Tmin=(100,'K'), Tmax=(1422.41,'K')), NASAPolynomial(coeffs=[31.0633,0.0105141,-4.97962e-07,-1.79005e-10,1.9091e-14,-75850.2,-135.229], Tmin=(1422.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-565.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(O)C(=O)C(C)OC(32462)',
    structure = SMILES('C=C(O)C(=O)C(C)OC'),
    E0 = (-507.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.563301,0.0917219,-7.97154e-05,3.49176e-08,-6.101e-12,-60879.2,33.9442], Tmin=(100,'K'), Tmax=(1369.96,'K')), NASAPolynomial(coeffs=[20.3434,0.0306782,-1.28765e-05,2.39129e-09,-1.65319e-13,-66607.4,-73.4992], Tmin=(1369.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-507.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]O[CH]C(C)(O)C(=C)O(32463)',
    structure = SMILES('[CH2]O[CH]C(C)(O)C(=C)O'),
    E0 = (-223.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,497.804,639.398,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157098,'amu*angstrom^2'), symmetry=1, barrier=(3.612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157098,'amu*angstrom^2'), symmetry=1, barrier=(3.612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157098,'amu*angstrom^2'), symmetry=1, barrier=(3.612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157098,'amu*angstrom^2'), symmetry=1, barrier=(3.612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157098,'amu*angstrom^2'), symmetry=1, barrier=(3.612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157098,'amu*angstrom^2'), symmetry=1, barrier=(3.612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157098,'amu*angstrom^2'), symmetry=1, barrier=(3.612,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.91707,0.115839,-0.000131701,7.35133e-08,-1.55972e-11,-26662.9,40.3393], Tmin=(100,'K'), Tmax=(1250.3,'K')), NASAPolynomial(coeffs=[26.8149,0.0164475,-3.4952e-06,3.73834e-10,-1.71603e-14,-33263.6,-102.359], Tmin=(1250.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsJOCC) + radical(CCsJOCs)"""),
)

species(
    label = 'C=C(O)C1(O)COC1C(31896)',
    structure = SMILES('C=C(O)C1(O)COC1C'),
    E0 = (-464.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.339438,0.075596,-5.35978e-06,-7.03174e-08,4.11107e-11,-55742.5,30.8326], Tmin=(100,'K'), Tmax=(878.042,'K')), NASAPolynomial(coeffs=[25.7086,0.0143915,1.03769e-06,-6.44615e-10,5.26358e-14,-62531.7,-104.059], Tmin=(878.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-464.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Oxetane)"""),
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
    E0 = (-291.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-155.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-191.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-30.6057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-159.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-88.6485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-181.455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-132.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-119.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-258.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-183.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-210.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (75.8215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (114.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-60.0076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-205.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (175.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-283.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-148.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-151.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (116.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (212.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-86.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-150.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-169.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-0.782082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-81.6336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-114.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-10.6388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-130.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-265.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-227.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-273.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-66.2429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-282.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['CH2O(15)', 'C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['[CH2][C](O)C1(O)COC1C(32437)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.20551e+07,'s^-1'), n=1.225, Ea=(136.213,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 134.8 to 136.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2O(15)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.0064799,'m^3/(mol*s)'), n=2.4275, Ea=(50.8904,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_R;CsJ-CdCsH] for rate rule [Od_CO-HH;CsJ-CdCsH]
Euclidian distance = 3.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['OH(5)', '[CH2]OC(C)C(O)=C=C(32438)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(41610,'cm^3/(mol*s)'), n=2.487, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 214 used for Ca_Cds-HH;OJ_pri
Exact match found for rate rule [Ca_Cds-HH;OJ_pri]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -7.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['[CH2]C(O)=C(O)[C](C)OC(32439)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(25964.4,'s^-1'), n=2.44333, Ea=(131.615,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_2H;Cs_H_out_noH] for rate rule [R3H_SS_O;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['[CH2]OC(C)C(O)=C(C)[O](32440)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.4947e+07,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_2Cd;C_rad_out_2H;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(O)=C(O)C([CH2])OC(32441)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['[CH2]OC(C)C([O])=C(C)O(32442)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(493835,'s^-1'), n=1.76395, Ea=(159.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_2H;XH_out] for rate rule [R4H_SDS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]O[C](C)C(O)=C(C)O(32443)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_2H] for rate rule [R4H_SDS;C_rad_out_NDMustO;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['[CH2]C(O)=C([O])C(C)OC(32444)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.06724e+06,'s^-1'), n=1.18977, Ea=(32.9114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_2H;XH_out] for rate rule [R5H_SSSS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]OC([CH2])C(O)=C(C)O(32445)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(121000,'s^-1'), n=1.9, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 92 used for R5H_SSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['[CH2]C([O])=C(O)C(C)OC(32446)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(12584.6,'s^-1'), n=1.925, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;C_rad_out_2H;XH_out] for rate rule [R6H_RSSMS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][O](167)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.71464e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [C_rad/H/CdCs;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]O[CH]C(567)', '[CH2]C(O)=[C]O(31700)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/CsO] for rate rule [Cd_rad;C_rad/H/CsO]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['[CH2]OC(C)C1(O)C[C]1O(32447)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_NdNd;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['[CH2]C1(O)COC(C)[C]1O(32448)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.95277e+09,'s^-1'), n=0.6, Ea=(85.772,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_secNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH2(S)(23)', '[CH2]OC[C](O)C(=C)O(32263)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(873476,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H2/OneDeO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['CC1OCCC(O)=C1O(31891)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;C_rad_out_2H;Cpri_rad_out_2H] + [R6_SSSDS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R6_SSSDS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['[CH2]OC(C)C(=O)C([CH2])O(32449)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['[CH2]OC(C)C(O)C([CH2])=O(32450)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(19)', '[CH2]C(O)=C(O)C(C)[O](31703)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(19)', '[CH2]OC(C)C(O)=[C]O(32451)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2]OC(C)C(=O)C(=C)O(32452)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.39e+06,'cm^3/(mol*s)'), n=2.09, Ea=(25.5224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2835 used for Od_CO-CdCs;HJ
Exact match found for rate rule [Od_CO-CdCs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]OC(C)=C(O)C(=C)O(32453)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH3(17)', '[CH2]OC=C(O)C(=C)O(32454)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.0063345,'m^3/(mol*s)'), n=2.46822, Ea=(26.7748,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CsJ-HHH] for rate rule [Cds-OsH_Cds;CsJ-HHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]OC(C)C([O])C(=C)O(32455)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.52488e+09,'s^-1'), n=1.21745, Ea=(162.572,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cd] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]O[C](C)C(O)C(=C)O(32456)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(8.2e+09,'s^-1'), n=0.65, Ea=(131.378,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_OneDe] for rate rule [R2H_S;C_rad_out_NDMustO;Cs_H_out_Cd]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]OC([CH2])C(O)C(=C)O(32457)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.14e-16,'s^-1'), n=8.15, Ea=(68.6594,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_2H;Cs_H_out_OneDe] for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_Cd]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C(O)C(O)C(C)O[CH2](32458)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.04e+10,'s^-1'), n=0.59, Ea=(135.98,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_NonDe] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C(O)[C](O)C(C)OC(32459)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['CC1OCC[C](O)[C]1O(32460)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.46e+06,'s^-1'), n=1.02, Ea=(25.5224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['C=C(O)C(O)=C(C)OC(32461)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['C=C(O)C(=O)C(C)OC(32462)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radEndo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]O[CH]C(C)(O)C(=C)O(32463)'],
    products = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]OC(C)[C](O)C(=C)O(31890)'],
    products = ['C=C(O)C1(O)COC1C(31896)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_OneDe/O]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

network(
    label = '5232',
    isomers = [
        '[CH2]OC(C)[C](O)C(=C)O(31890)',
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
    label = '5232',
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

