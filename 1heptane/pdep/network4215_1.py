species(
    label = '[CH2]C(O)C(C)[C]=C(20928)',
    structure = SMILES('[CH2]C(O)C(C)[C]=C'),
    E0 = (204.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,320.28,3945.52],'cm^-1')),
        HinderedRotor(inertia=(0.141892,'amu*angstrom^2'), symmetry=1, barrier=(10.32,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.273729,'amu*angstrom^2'), symmetry=1, barrier=(19.9568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0730315,'amu*angstrom^2'), symmetry=1, barrier=(5.33058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27379,'amu*angstrom^2'), symmetry=1, barrier=(19.9575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.273906,'amu*angstrom^2'), symmetry=1, barrier=(19.9582,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.456708,0.074416,-6.43583e-05,3.00087e-08,-5.72844e-12,24775,30.5054], Tmin=(100,'K'), Tmax=(1243.63,'K')), NASAPolynomial(coeffs=[13.7861,0.0315434,-1.26476e-05,2.28835e-09,-1.55967e-13,21459.6,-36.7076], Tmin=(1243.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CJCO)"""),
)

species(
    label = 'CH2CHOH(42)',
    structure = SMILES('C=CO'),
    E0 = (-138.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.72808,'amu*angstrom^2'), symmetry=1, barrier=(39.7321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28758,0.0197013,1.96383e-06,-1.9439e-08,1.02617e-11,-16537.3,14.1333], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[7.49818,0.0103957,-3.66891e-06,5.85206e-10,-3.47374e-14,-18164.3,-13.8388], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-138.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]C(O)C(C)=C=C(25092)',
    structure = SMILES('[CH2]C(O)C(C)=C=C'),
    E0 = (123.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,314.025],'cm^-1')),
        HinderedRotor(inertia=(0.181909,'amu*angstrom^2'), symmetry=1, barrier=(12.7287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1819,'amu*angstrom^2'), symmetry=1, barrier=(12.7286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181827,'amu*angstrom^2'), symmetry=1, barrier=(12.7288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181832,'amu*angstrom^2'), symmetry=1, barrier=(12.7289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479977,0.0737079,-6.63395e-05,3.14653e-08,-6.0415e-12,15017.6,28.2305], Tmin=(100,'K'), Tmax=(1244.31,'K')), NASAPolynomial(coeffs=[14.6128,0.0282757,-1.1571e-05,2.12161e-09,-1.45879e-13,11500.5,-43.0413], Tmin=(1244.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJCO)"""),
)

species(
    label = 'C=[C]C(C)C(=C)O(25093)',
    structure = SMILES('C=[C]C(C)C(=C)O'),
    E0 = (78.4632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,191.069,191.255],'cm^-1')),
        HinderedRotor(inertia=(0.537604,'amu*angstrom^2'), symmetry=1, barrier=(13.9898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.529071,'amu*angstrom^2'), symmetry=1, barrier=(13.9882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.540397,'amu*angstrom^2'), symmetry=1, barrier=(13.9907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.538676,'amu*angstrom^2'), symmetry=1, barrier=(13.9898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.199685,0.076551,-7.19768e-05,3.53677e-08,-6.90749e-12,9579.69,26.5827], Tmin=(100,'K'), Tmax=(1242.4,'K')), NASAPolynomial(coeffs=[16.4084,0.0243651,-8.96989e-06,1.55807e-09,-1.04111e-13,5552.19,-55.1329], Tmin=(1242.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.4632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(C)C([CH2])O(25094)',
    structure = SMILES('C#CC(C)C([CH2])O'),
    E0 = (132.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2175,525,750,770,3400,2100,3000,3100,440,815,1455,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,319.957],'cm^-1')),
        HinderedRotor(inertia=(0.00164752,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150292,'amu*angstrom^2'), symmetry=1, barrier=(10.9141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150358,'amu*angstrom^2'), symmetry=1, barrier=(10.9141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.278249,'amu*angstrom^2'), symmetry=1, barrier=(20.1991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02477,'amu*angstrom^2'), symmetry=1, barrier=(74.4451,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277409,0.0758323,-7.29405e-05,3.78911e-08,-7.86261e-12,16016.8,28.4818], Tmin=(100,'K'), Tmax=(1172.07,'K')), NASAPolynomial(coeffs=[14.9339,0.0258122,-8.92407e-06,1.47821e-09,-9.56532e-14,12581.2,-44.554], Tmin=(1172.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]O(284)',
    structure = SMILES('[CH2][CH]O'),
    E0 = (143.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.217851,'amu*angstrom^2'), symmetry=1, barrier=(5.00882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0197382,'amu*angstrom^2'), symmetry=1, barrier=(14.867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84769,0.0229973,-1.85068e-05,7.77211e-09,-1.30725e-12,17300.5,13.0245], Tmin=(100,'K'), Tmax=(1418.34,'K')), NASAPolynomial(coeffs=[7.97636,0.00853337,-3.21015e-06,5.82141e-10,-3.99268e-14,15845.7,-13.5107], Tmin=(1418.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CJCO) + radical(CCsJOH)"""),
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
    label = '[CH2]C(O)C=C=C(25095)',
    structure = SMILES('[CH2]C(O)C=C=C'),
    E0 = (162.839,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,273.356],'cm^-1')),
        HinderedRotor(inertia=(0.274745,'amu*angstrom^2'), symmetry=1, barrier=(14.5808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275137,'amu*angstrom^2'), symmetry=1, barrier=(14.5806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274031,'amu*angstrom^2'), symmetry=1, barrier=(14.5823,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1135,0.0561689,-4.43731e-05,1.55638e-08,-1.44504e-12,19695.2,25.2009], Tmin=(100,'K'), Tmax=(1093.41,'K')), NASAPolynomial(coeffs=[13.899,0.0193221,-7.44185e-06,1.34896e-09,-9.32639e-14,16305.9,-40.3369], Tmin=(1093.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJCO)"""),
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
    label = 'C=[C]C(C)C=C(24201)',
    structure = SMILES('C=[C]C(C)C=C'),
    E0 = (293.022,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1685,370,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,272.335,275.352],'cm^-1')),
        HinderedRotor(inertia=(0.2164,'amu*angstrom^2'), symmetry=1, barrier=(11.4749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214712,'amu*angstrom^2'), symmetry=1, barrier=(11.4778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213609,'amu*angstrom^2'), symmetry=1, barrier=(11.4757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55688,0.0537055,-3.52508e-05,1.22891e-08,-1.83166e-12,35329.9,21.6469], Tmin=(100,'K'), Tmax=(1475.3,'K')), NASAPolynomial(coeffs=[9.43814,0.032337,-1.35245e-05,2.47131e-09,-1.67972e-13,33004.4,-19.4404], Tmin=(1475.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)[C](C)C=C(20926)',
    structure = SMILES('[CH2]C=C(C)C([CH2])O'),
    E0 = (98.6792,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.186455,'amu*angstrom^2'), symmetry=1, barrier=(4.28697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.742029,'amu*angstrom^2'), symmetry=1, barrier=(17.0607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0115517,'amu*angstrom^2'), symmetry=1, barrier=(4.26604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.742983,'amu*angstrom^2'), symmetry=1, barrier=(17.0826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.046655,'amu*angstrom^2'), symmetry=1, barrier=(17.0392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.072416,0.0757481,-6.26083e-05,2.66347e-08,-4.52233e-12,12018.6,30.9304], Tmin=(100,'K'), Tmax=(1413.97,'K')), NASAPolynomial(coeffs=[17.696,0.0258918,-9.71777e-06,1.69712e-09,-1.13123e-13,7034.88,-60.1975], Tmin=(1413.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.6792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CJCO) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C(C)[C](C)O(24689)',
    structure = SMILES('C=[C]C(C)[C](C)O'),
    E0 = (169.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.699461,0.0749737,-7.0828e-05,3.86555e-08,-8.92456e-12,20556.5,28.4332], Tmin=(100,'K'), Tmax=(1022.22,'K')), NASAPolynomial(coeffs=[10.2144,0.0377404,-1.61905e-05,3.02123e-09,-2.09409e-13,18611.3,-17.6796], Tmin=(1022.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C2CsJOH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C)C([CH2])O(20929)',
    structure = SMILES('[CH]=CC(C)C([CH2])O'),
    E0 = (214.161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,724.81],'cm^-1')),
        HinderedRotor(inertia=(0.0346757,'amu*angstrom^2'), symmetry=1, barrier=(12.9383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107107,'amu*angstrom^2'), symmetry=1, barrier=(2.46261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.562691,'amu*angstrom^2'), symmetry=1, barrier=(12.9374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.562639,'amu*angstrom^2'), symmetry=1, barrier=(12.9362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.562659,'amu*angstrom^2'), symmetry=1, barrier=(12.9366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.185085,0.0760002,-6.52471e-05,2.94128e-08,-5.32625e-12,25901.5,31.4133], Tmin=(100,'K'), Tmax=(1325.23,'K')), NASAPolynomial(coeffs=[16.2356,0.0275543,-1.04124e-05,1.82781e-09,-1.22458e-13,21647.4,-50.5409], Tmin=(1325.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CJCO)"""),
)

species(
    label = '[CH2][C](O)C(C)C=C(20927)',
    structure = SMILES('[CH2][C](O)C(C)C=C'),
    E0 = (143.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.452381,0.077155,-7.17033e-05,3.66094e-08,-7.6903e-12,17410.7,29.4736], Tmin=(100,'K'), Tmax=(1133.85,'K')), NASAPolynomial(coeffs=[12.8773,0.0333224,-1.37158e-05,2.51449e-09,-1.72775e-13,14593.1,-32.03], Tmin=(1133.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCO) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(O)C([CH2])C=C(19570)',
    structure = SMILES('[CH2]C(O)C([CH2])C=C'),
    E0 = (172.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,605.686,3168.84],'cm^-1')),
        HinderedRotor(inertia=(0.00452218,'amu*angstrom^2'), symmetry=1, barrier=(5.06236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.70963,'amu*angstrom^2'), symmetry=1, barrier=(16.3158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709266,'amu*angstrom^2'), symmetry=1, barrier=(16.3074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.708798,'amu*angstrom^2'), symmetry=1, barrier=(16.2967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.6003,'amu*angstrom^2'), symmetry=1, barrier=(82.778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3818.2,'J/mol'), sigma=(6.62498,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.39 K, Pc=29.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.167345,0.0754739,-6.57891e-05,3.08102e-08,-5.78924e-12,20849.9,32.8542], Tmin=(100,'K'), Tmax=(1285.7,'K')), NASAPolynomial(coeffs=[15.712,0.0271128,-9.36777e-06,1.55465e-09,-1.00661e-13,16852.7,-46.0463], Tmin=(1285.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CJCO)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3818.2,'J/mol'), sigma=(6.62498,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.39 K, Pc=29.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.463893,0.0686477,-4.96367e-05,1.85833e-08,-2.8272e-12,27037.1,29.6678], Tmin=(100,'K'), Tmax=(1539.81,'K')), NASAPolynomial(coeffs=[15.6202,0.0292757,-1.12826e-05,1.97779e-09,-1.31166e-13,22369.5,-49.995], Tmin=(1539.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=[C][C](C)C(C)O(24690)',
    structure = SMILES('[CH2][C]=C(C)C(C)O'),
    E0 = (124.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.46359,0.072069,-5.73265e-05,2.39457e-08,-4.08223e-12,15157.7,29.3583], Tmin=(100,'K'), Tmax=(1381.1,'K')), NASAPolynomial(coeffs=[14.7652,0.0306482,-1.234e-05,2.23053e-09,-1.51463e-13,11207.3,-44.2566], Tmin=(1381.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=C)C(C)O(19571)',
    structure = SMILES('[CH2]C([C]=C)C(C)O'),
    E0 = (198.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,303.126,303.259],'cm^-1')),
        HinderedRotor(inertia=(0.100895,'amu*angstrom^2'), symmetry=1, barrier=(6.59071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101021,'amu*angstrom^2'), symmetry=1, barrier=(6.58973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00182961,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00182565,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308549,'amu*angstrom^2'), symmetry=1, barrier=(20.1166,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.521918,0.0721561,-6.14939e-05,2.90921e-08,-5.66065e-12,23990.8,31.4188], Tmin=(100,'K'), Tmax=(1222.9,'K')), NASAPolynomial(coeffs=[12.871,0.0317625,-1.19466e-05,2.08085e-09,-1.38589e-13,20970.5,-30.6434], Tmin=(1222.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([O])C(C)C=C(19572)',
    structure = SMILES('[CH2]C([O])C(C)C=C'),
    E0 = (197.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,1197.87],'cm^-1')),
        HinderedRotor(inertia=(0.00173266,'amu*angstrom^2'), symmetry=1, barrier=(19.6727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121323,'amu*angstrom^2'), symmetry=1, barrier=(19.6532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.85577,'amu*angstrom^2'), symmetry=1, barrier=(19.6758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855519,'amu*angstrom^2'), symmetry=1, barrier=(19.6701,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.44992,0.0681862,-4.18449e-05,6.2521e-09,2.34773e-12,23881,29.8659], Tmin=(100,'K'), Tmax=(1066.24,'K')), NASAPolynomial(coeffs=[15.1529,0.0298393,-1.15483e-05,2.09651e-09,-1.4509e-13,19790.1,-46.4916], Tmin=(1066.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH]=[C]C(C)C(C)O(24691)',
    structure = SMILES('[CH]=[C]C(C)C(C)O'),
    E0 = (240.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3120,650,792.5,1650,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,318.916],'cm^-1')),
        HinderedRotor(inertia=(0.158729,'amu*angstrom^2'), symmetry=1, barrier=(11.4564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00165742,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158727,'amu*angstrom^2'), symmetry=1, barrier=(11.4564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158735,'amu*angstrom^2'), symmetry=1, barrier=(11.4564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158737,'amu*angstrom^2'), symmetry=1, barrier=(11.4564,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.549393,0.0725811,-6.06552e-05,2.73812e-08,-5.09026e-12,29041.9,29.942], Tmin=(100,'K'), Tmax=(1268.77,'K')), NASAPolynomial(coeffs=[13.3223,0.0323126,-1.30481e-05,2.36649e-09,-1.61361e-13,25800.7,-34.7206], Tmin=(1268.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C=C(C)C(C)O(24692)',
    structure = SMILES('C=C=C(C)C(C)O'),
    E0 = (-87.8048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.41144,0.0709935,-5.46124e-05,2.17759e-08,-3.51959e-12,-10424.8,27.8505], Tmin=(100,'K'), Tmax=(1458.78,'K')), NASAPolynomial(coeffs=[15.7288,0.028993,-1.14251e-05,2.03911e-09,-1.3718e-13,-14893.7,-51.8307], Tmin=(1458.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.8048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC(C)C(=C)O(20932)',
    structure = SMILES('C=CC(C)C(=C)O'),
    E0 = (-159.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281536,0.0714086,-4.54864e-05,4.24118e-09,4.67284e-12,-19025.8,25.8135], Tmin=(100,'K'), Tmax=(995.086,'K')), NASAPolynomial(coeffs=[16.9122,0.0259875,-9.32223e-06,1.65497e-09,-1.14707e-13,-23396.6,-59.6691], Tmin=(995.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-159.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(O)C[C]=C(25076)',
    structure = SMILES('[CH2]C(O)C[C]=C'),
    E0 = (236.658,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,372.181,372.223],'cm^-1')),
        HinderedRotor(inertia=(0.00121656,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105245,'amu*angstrom^2'), symmetry=1, barrier=(10.3435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105343,'amu*angstrom^2'), symmetry=1, barrier=(10.3426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10521,'amu*angstrom^2'), symmetry=1, barrier=(10.3419,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23604,0.0606913,-5.59104e-05,2.8476e-08,-6.00418e-12,28563,25.7111], Tmin=(100,'K'), Tmax=(1125.24,'K')), NASAPolynomial(coeffs=[10.6291,0.0273011,-1.13996e-05,2.10489e-09,-1.45187e-13,26449.2,-20.7131], Tmin=(1125.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)C(=C)[CH]C(24161)',
    structure = SMILES('[CH2]C(=CC)C([CH2])O'),
    E0 = (98.6792,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.186455,'amu*angstrom^2'), symmetry=1, barrier=(4.28697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.742029,'amu*angstrom^2'), symmetry=1, barrier=(17.0607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0115517,'amu*angstrom^2'), symmetry=1, barrier=(4.26604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.742983,'amu*angstrom^2'), symmetry=1, barrier=(17.0826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.046655,'amu*angstrom^2'), symmetry=1, barrier=(17.0392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3895.04,'J/mol'), sigma=(6.6861,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=608.40 K, Pc=29.57 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.072416,0.0757481,-6.26083e-05,2.66347e-08,-4.52233e-12,12018.6,30.9304], Tmin=(100,'K'), Tmax=(1413.97,'K')), NASAPolynomial(coeffs=[17.696,0.0258918,-9.71777e-06,1.69712e-09,-1.13123e-13,7034.88,-60.1975], Tmin=(1413.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.6792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)[CH]C(=C)C(25096)',
    structure = SMILES('[CH2]C(O)[CH]C(=C)C'),
    E0 = (76.6783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.117597,0.0811085,-7.30844e-05,3.46846e-08,-6.6616e-12,9365.58,27.2063], Tmin=(100,'K'), Tmax=(1244.46,'K')), NASAPolynomial(coeffs=[15.7138,0.0309784,-1.26605e-05,2.315e-09,-1.58867e-13,5483.82,-51.4473], Tmin=(1244.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.6783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(CJCO)"""),
)

species(
    label = 'C=[C]C(C)C[CH]O(20943)',
    structure = SMILES('C=[C]C(C)C[CH]O'),
    E0 = (186.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612223,0.0768467,-7.47792e-05,4.1832e-08,-9.83103e-12,22521.7,29.2167], Tmin=(100,'K'), Tmax=(1009.39,'K')), NASAPolynomial(coeffs=[10.6035,0.0372525,-1.59393e-05,2.96957e-09,-2.05626e-13,20504.7,-19.0788], Tmin=(1009.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCsJOH)"""),
)

species(
    label = 'C=C1CC(O)C1C(25097)',
    structure = SMILES('C=C1CC(O)C1C'),
    E0 = (-103.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.4513,-0.0108149,0.000133113,-1.2865e-07,2.97531e-11,-12839.6,-17.7409], Tmin=(100,'K'), Tmax=(1707.07,'K')), NASAPolynomial(coeffs=[73.7292,0.032536,-7.33583e-05,1.77402e-08,-1.31581e-12,-62364,-438.639], Tmin=(1707.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane)"""),
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
    label = '[CH2]C(O)[CH]C(2450)',
    structure = SMILES('[CH2]C(O)[CH]C'),
    E0 = (95.4345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,208.196],'cm^-1')),
        HinderedRotor(inertia=(0.00389407,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00388995,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123961,'amu*angstrom^2'), symmetry=1, barrier=(3.80763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123787,'amu*angstrom^2'), symmetry=1, barrier=(3.80541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60717,0.0498047,-3.98259e-05,1.7252e-08,-3.07178e-12,11566.5,22.2154], Tmin=(100,'K'), Tmax=(1321.87,'K')), NASAPolynomial(coeffs=[10.6575,0.0224181,-8.74868e-06,1.57867e-09,-1.07529e-13,9173.84,-23.9728], Tmin=(1321.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.4345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo library: CBS_QB3_1dHR + radical(CCJCO) + radical(CJCO)"""),
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
    label = 'C=[C]C(C)[CH]O(24544)',
    structure = SMILES('C=[C]C(C)[CH]O'),
    E0 = (210.038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,349.763,2545.98],'cm^-1')),
        HinderedRotor(inertia=(0.607579,'amu*angstrom^2'), symmetry=1, barrier=(13.9694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.607626,'amu*angstrom^2'), symmetry=1, barrier=(13.9705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608244,'amu*angstrom^2'), symmetry=1, barrier=(13.9847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0583731,'amu*angstrom^2'), symmetry=1, barrier=(13.9576,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23097,0.0624087,-6.22036e-05,3.54679e-08,-8.42595e-12,25360.2,24.7489], Tmin=(100,'K'), Tmax=(1004.1,'K')), NASAPolynomial(coeffs=[9.63878,0.028915,-1.21686e-05,2.24762e-09,-1.54856e-13,23671.8,-15.8484], Tmin=(1004.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOH) + radical(Cds_S)"""),
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
    E0 = (204.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (347.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (301.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (356.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (309.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (331.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (231.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (405.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (387.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (307.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (319.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (346.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (356.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (362.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (322.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (287.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (306.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (273.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (504.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (283.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (293.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (656.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (299.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (349.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (363.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (213.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (496.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (591.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['CH2CHOH(42)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', '[CH2]C(O)C(C)=C=C(25092)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.51e+07,'cm^3/(mol*s)'), n=1.64, Ea=(11.8407,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2579 used for Cds-CsCs_Ca;HJ
Exact match found for rate rule [Cds-CsCs_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=[C]C(C)C(=C)O(25093)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#CC(C)C([CH2])O(25094)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]O(284)', 'CH3CHCCH2(18175)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00472174,'m^3/(mol*s)'), n=2.41, Ea=(20.1294,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH3(17)', '[CH2]C(O)C=C=C(25095)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(10800,'cm^3/(mol*s)'), n=2.41, Ea=(32.1331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 597 used for Cds-CsH_Ca;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Ca;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2CHOH(42)', 'C=[C][CH]C(18176)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.30847,'m^3/(mol*s)'), n=2.06429, Ea=(8.71744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds_Cds;CJ] + [Cds-OsH_Cds;YJ] for rate rule [Cds-OsH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OH(5)', 'C=[C]C(C)C=C(24201)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00674586,'m^3/(mol*s)'), n=2.3625, Ea=(84.2031,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;OJ] for rate rule [Cds-CsH_Cds-HH;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['[CH2]C(O)[C](C)C=C(20926)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_Cs2]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['C=[C]C(C)[C](C)O(24689)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC(C)C([CH2])O(20929)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['[CH2][C](O)C(C)C=C(20927)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['[CH2]C(O)C([CH2])C=C(19570)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(C)C(C)[O](19568)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['C=[C][C](C)C(C)O(24690)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['[CH2]C([C]=C)C(C)O(19571)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['[CH2]C([O])C(C)C=C(19572)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]C(C)C(C)O(24691)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]O(284)', 'C=[C][CH]C(18176)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['C=C=C(C)C(C)O(24692)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['C=CC(C)C(=C)O(20932)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(S)(23)', '[CH2]C(O)C[C]=C(25076)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['[CH2]C(O)C(=C)[CH]C(24161)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['[CH2]C(O)[CH]C(=C)C(25096)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.95888e+10,'s^-1'), n=0.7315, Ea=(144.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HC)CJ;CJ;CH3] + [cCs(-HC)CJ;CdsJ;C] for rate rule [cCs(-HC)CJ;CdsJ;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['C=[C]C(C)C[CH]O(20943)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(O)C(C)[C]=C(20928)'],
    products = ['C=C1CC(O)C1C(25097)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H2CC(41)', '[CH2]C(O)[CH]C(2450)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(19)', 'C=[C]C(C)[CH]O(24544)'],
    products = ['[CH2]C(O)C(C)[C]=C(20928)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4215',
    isomers = [
        '[CH2]C(O)C(C)[C]=C(20928)',
    ],
    reactants = [
        ('CH2CHOH(42)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4215',
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

