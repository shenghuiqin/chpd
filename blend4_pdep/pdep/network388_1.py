species(
    label = 'S(252)(251)',
    structure = SMILES('CC(O)=C=CC=O'),
    E0 = (-191.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,540,610,2055,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.833405,'amu*angstrom^2'), symmetry=1, barrier=(19.1616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.82892,'amu*angstrom^2'), symmetry=1, barrier=(19.0585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828534,'amu*angstrom^2'), symmetry=1, barrier=(19.0496,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4176.29,'J/mol'), sigma=(6.51688,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=652.33 K, Pc=34.24 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.986472,0.0659047,-6.4425e-05,3.18942e-08,-6.35286e-12,-22904.6,21.9789], Tmin=(100,'K'), Tmax=(1200.1,'K')), NASAPolynomial(coeffs=[13.8493,0.0230323,-1.08392e-05,2.12683e-09,-1.51862e-13,-25992,-42.4231], Tmin=(1200.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-191.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
)

species(
    label = 'OH(5)(5)',
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
    label = 'C[C]=C=CC=O(6394)',
    structure = SMILES('CC#CC=C[O]'),
    E0 = (161.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.56284,'amu*angstrom^2'), symmetry=1, barrier=(35.9327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57121,'amu*angstrom^2'), symmetry=1, barrier=(36.1253,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80258,0.0332802,2.22527e-05,-6.26465e-08,2.95311e-11,19480.5,19.4671], Tmin=(100,'K'), Tmax=(920.512,'K')), NASAPolynomial(coeffs=[16.5497,0.00825504,-6.12625e-07,6.96751e-12,-3.27546e-15,15110.8,-59.4458], Tmin=(920.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CtHHH) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=COJ)"""),
)

species(
    label = 'H(3)(3)',
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
    label = 'CC(=O)[C]=CC=O(6096)',
    structure = SMILES('CC([O])=C=CC=O'),
    E0 = (-53.5394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.653958,'amu*angstrom^2'), symmetry=1, barrier=(15.0358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.652529,'amu*angstrom^2'), symmetry=1, barrier=(15.0029,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48531,0.0601767,-6.50487e-05,3.92374e-08,-9.97775e-12,-6352.84,21.7039], Tmin=(100,'K'), Tmax=(931.209,'K')), NASAPolynomial(coeffs=[8.83561,0.028604,-1.41916e-05,2.82842e-09,-2.03226e-13,-7721.79,-13.2332], Tmin=(931.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.5394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CH3(15)(16)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=CC=C=[C]O(6395)',
    structure = SMILES('[O]C=CC#CO'),
    E0 = (61.9116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2100,2250,500,550,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44972,'amu*angstrom^2'), symmetry=1, barrier=(33.3319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45335,'amu*angstrom^2'), symmetry=1, barrier=(33.4153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8124,0.0347446,3.22938e-06,-3.94368e-08,2.01747e-11,7537.22,19.3479], Tmin=(100,'K'), Tmax=(948.051,'K')), NASAPolynomial(coeffs=[17.0656,0.00367907,-2.9063e-07,7.7097e-11,-1.21406e-14,3149,-61.3163], Tmin=(948.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.9116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-CtH) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtOs) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(O)=C=CC=O(6396)',
    structure = SMILES('C=C(O)[C]=CC=O'),
    E0 = (-44.4643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1685,370,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.08335,'amu*angstrom^2'), symmetry=1, barrier=(24.9083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08324,'amu*angstrom^2'), symmetry=1, barrier=(24.9057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08275,'amu*angstrom^2'), symmetry=1, barrier=(24.8945,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.661659,0.0720618,-8.35572e-05,4.86205e-08,-1.10886e-11,-5226.36,21.8972], Tmin=(100,'K'), Tmax=(1074.06,'K')), NASAPolynomial(coeffs=[15.2881,0.0175901,-7.48361e-06,1.40174e-09,-9.79087e-14,-8368.29,-49.7119], Tmin=(1074.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.4643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=C)"""),
)

species(
    label = 'HCO(14)(15)',
    structure = SMILES('[CH]=O'),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=C=C(C)O(6397)',
    structure = SMILES('[CH]=C=C(C)O'),
    E0 = (86.6206,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(0.984925,'amu*angstrom^2'), symmetry=1, barrier=(22.6454,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10635,0.0503488,-5.11114e-05,2.5565e-08,-4.75667e-12,10533.6,19.0983], Tmin=(100,'K'), Tmax=(1545.92,'K')), NASAPolynomial(coeffs=[14.256,0.00720865,-4.07516e-07,-1.15074e-10,1.30724e-14,7557.26,-46.5464], Tmin=(1545.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.6206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ)"""),
)

species(
    label = 'CC(O)=C=[C]C=O(6398)',
    structure = SMILES('CC(O)=C=C=C[O]'),
    E0 = (5.39251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.19634,'amu*angstrom^2'), symmetry=1, barrier=(27.5061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19457,'amu*angstrom^2'), symmetry=1, barrier=(27.4654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.341136,0.0756378,-8.14065e-05,4.07995e-08,-7.55309e-12,821.889,25.803], Tmin=(100,'K'), Tmax=(1532.39,'K')), NASAPolynomial(coeffs=[22.2459,0.00415139,8.31734e-07,-3.13139e-10,2.45532e-14,-4629.67,-88.0078], Tmin=(1532.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.39251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'CC(O)=C=C[C]=O(6399)',
    structure = SMILES('CC(O)=C=C[C]=O'),
    E0 = (-36.6733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1855,455,950,540,610,2055,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.849569,'amu*angstrom^2'), symmetry=1, barrier=(19.5333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.849203,'amu*angstrom^2'), symmetry=1, barrier=(19.5249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.848549,'amu*angstrom^2'), symmetry=1, barrier=(19.5098,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11327,0.0653135,-7.04144e-05,3.86206e-08,-8.51992e-12,-4308.29,22.2071], Tmin=(100,'K'), Tmax=(1089.85,'K')), NASAPolynomial(coeffs=[12.9142,0.0220018,-1.08031e-05,2.15622e-09,-1.55423e-13,-6880.54,-35.7407], Tmin=(1089.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-36.6733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(CsCJ=O)"""),
)

species(
    label = 'CC([O])[C]=CC=O(6100)',
    structure = SMILES('CC([O])[C]=CC=O'),
    E0 = (155.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1685,370,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,339.848,340.532],'cm^-1')),
        HinderedRotor(inertia=(0.148104,'amu*angstrom^2'), symmetry=1, barrier=(12.0269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144957,'amu*angstrom^2'), symmetry=1, barrier=(12.0204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147659,'amu*angstrom^2'), symmetry=1, barrier=(12.0486,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69495,0.0553565,-4.39351e-05,1.76841e-08,-3.01337e-12,18767.3,24.9739], Tmin=(100,'K'), Tmax=(1312.37,'K')), NASAPolynomial(coeffs=[10.0004,0.0300422,-1.50016e-05,2.98625e-09,-2.13492e-13,16587.4,-17.3527], Tmin=(1312.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)[C]=CC=O(6400)',
    structure = SMILES('[CH2]C(O)[C]=CC=O'),
    E0 = (136.614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1685,370,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,351.201],'cm^-1')),
        HinderedRotor(inertia=(0.125612,'amu*angstrom^2'), symmetry=1, barrier=(10.9801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12546,'amu*angstrom^2'), symmetry=1, barrier=(10.9849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125299,'amu*angstrom^2'), symmetry=1, barrier=(10.9817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125445,'amu*angstrom^2'), symmetry=1, barrier=(10.9867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28861,0.0653089,-7.11008e-05,4.26404e-08,-1.07621e-11,16523.8,27.2814], Tmin=(100,'K'), Tmax=(937.534,'K')), NASAPolynomial(coeffs=[9.42849,0.0305799,-1.55362e-05,3.12905e-09,-2.26125e-13,14997.5,-11.4638], Tmin=(937.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = 'CC(O)=C=[C]C[O](6401)',
    structure = SMILES('C[C](O)C#CC[O]'),
    E0 = (183.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,360,370,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2250,500,550,2750,2850,1437.5,1250,1305,750,350,186.269,186.368,1386.05],'cm^-1')),
        HinderedRotor(inertia=(0.247197,'amu*angstrom^2'), symmetry=1, barrier=(6.09796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247449,'amu*angstrom^2'), symmetry=1, barrier=(6.09564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00447117,'amu*angstrom^2'), symmetry=1, barrier=(6.09566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.5961,'amu*angstrom^2'), symmetry=1, barrier=(63.7238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20723,0.0659955,-8.88176e-05,7.17049e-08,-2.30887e-11,22182.7,28.4473], Tmin=(100,'K'), Tmax=(893.986,'K')), NASAPolynomial(coeffs=[6.92393,0.0304159,-1.23387e-05,2.15905e-09,-1.41061e-13,21560.2,3.74343], Tmin=(893.986,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C2CsJOH) + radical(CCOJ)"""),
)

species(
    label = 'C[C](O)C=[C]C=O(6107)',
    structure = SMILES('CC(O)=C[C]=C[O]'),
    E0 = (-24.9975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17014,'amu*angstrom^2'), symmetry=1, barrier=(26.9037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16965,'amu*angstrom^2'), symmetry=1, barrier=(26.8925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16976,'amu*angstrom^2'), symmetry=1, barrier=(26.8951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.928288,0.0811938,-8.62323e-05,4.28041e-08,-7.74497e-12,-2805.62,28.4291], Tmin=(100,'K'), Tmax=(1623.21,'K')), NASAPolynomial(coeffs=[22.5678,0.00443465,2.12776e-06,-6.43719e-10,4.90632e-14,-7948.99,-88.6548], Tmin=(1623.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.9975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = 'CC(O)=[C]C[C]=O(6402)',
    structure = SMILES('CC(O)=[C]C[C]=O'),
    E0 = (54.0795,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1855,455,950,3615,1277.5,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,291.409],'cm^-1')),
        HinderedRotor(inertia=(0.221922,'amu*angstrom^2'), symmetry=1, barrier=(13.3735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221918,'amu*angstrom^2'), symmetry=1, barrier=(13.3735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221931,'amu*angstrom^2'), symmetry=1, barrier=(13.3735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221915,'amu*angstrom^2'), symmetry=1, barrier=(13.3735,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808721,0.0631018,-5.69119e-05,2.49064e-08,-4.28248e-12,6624.97,26.2834], Tmin=(100,'K'), Tmax=(1402.58,'K')), NASAPolynomial(coeffs=[16.9855,0.0169681,-7.57439e-06,1.4559e-09,-1.02641e-13,2087.07,-57.2331], Tmin=(1402.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.0795,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'CC([O])=CC=C[O](5355)',
    structure = SMILES('CC([O])=CC=C[O]'),
    E0 = (-86.1882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.912938,'amu*angstrom^2'), symmetry=1, barrier=(20.9902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.915391,'amu*angstrom^2'), symmetry=1, barrier=(21.0466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.756235,0.0562742,-1.51561e-05,-3.6247e-08,2.26938e-11,-10235.2,24.2724], Tmin=(100,'K'), Tmax=(917.521,'K')), NASAPolynomial(coeffs=[21.0957,0.00729401,3.01004e-08,-1.33448e-10,7.27657e-15,-15638.3,-81.2071], Tmin=(917.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.1882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)=C[CH]C=O(6403)',
    structure = SMILES('C=C(O)[CH]C=C[O]'),
    E0 = (-76.2221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.35339,'amu*angstrom^2'), symmetry=1, barrier=(31.117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35269,'amu*angstrom^2'), symmetry=1, barrier=(31.1009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35312,'amu*angstrom^2'), symmetry=1, barrier=(31.1109,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.400843,0.0594199,-7.08619e-06,-5.37493e-08,3.0632e-11,-9019.51,24.2686], Tmin=(100,'K'), Tmax=(924.396,'K')), NASAPolynomial(coeffs=[25.3551,0.00260988,2.06445e-06,-4.65225e-10,2.62503e-14,-15819.3,-105.985], Tmin=(924.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.2221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(C=COJ)"""),
)

species(
    label = 'C[C](O)C=C[C]=O(6109)',
    structure = SMILES('CC(O)=C[CH][C]=O'),
    E0 = (-21.6336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1855,455,950,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,375.437],'cm^-1')),
        HinderedRotor(inertia=(0.191575,'amu*angstrom^2'), symmetry=1, barrier=(19.1619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191575,'amu*angstrom^2'), symmetry=1, barrier=(19.1619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191575,'amu*angstrom^2'), symmetry=1, barrier=(19.1619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191575,'amu*angstrom^2'), symmetry=1, barrier=(19.1619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.943088,0.0514379,-8.03709e-06,-3.43369e-08,1.82981e-11,-2477.79,25.7054], Tmin=(100,'K'), Tmax=(990.28,'K')), NASAPolynomial(coeffs=[19.613,0.0127286,-4.99884e-06,1.04521e-09,-8.29423e-14,-7975.14,-73.2705], Tmin=(990.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.6336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(CCCJ=O)"""),
)

species(
    label = 'CC([O])=[C]CC=O(6404)',
    structure = SMILES('CC([O])=[C]CC=O'),
    E0 = (31.9237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,350,440,435,1725,353.709,353.715],'cm^-1')),
        HinderedRotor(inertia=(0.00134743,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122239,'amu*angstrom^2'), symmetry=1, barrier=(10.8533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122252,'amu*angstrom^2'), symmetry=1, barrier=(10.8531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54122,0.0527991,-3.7905e-05,1.30063e-08,-1.80161e-12,3928.46,24.3284], Tmin=(100,'K'), Tmax=(1644,'K')), NASAPolynomial(coeffs=[13.6166,0.0234186,-1.10978e-05,2.13555e-09,-1.48499e-13,-41.8947,-39.9311], Tmin=(1644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.9237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)=[C]CC=O(6405)',
    structure = SMILES('[CH2]C(O)=[C]CC=O'),
    E0 = (53.0355,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,350,440,435,1725,3615,1277.5,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,277.581],'cm^-1')),
        HinderedRotor(inertia=(0.330374,'amu*angstrom^2'), symmetry=1, barrier=(18.0947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00219251,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329484,'amu*angstrom^2'), symmetry=1, barrier=(18.0948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.330383,'amu*angstrom^2'), symmetry=1, barrier=(18.0955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.755888,0.0617819,-4.9183e-05,1.55815e-08,-8.70302e-13,6503.59,26.3414], Tmin=(100,'K'), Tmax=(1120.64,'K')), NASAPolynomial(coeffs=[17.0284,0.0171408,-7.42234e-06,1.44191e-09,-1.03839e-13,2012.42,-57.7836], Tmin=(1120.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.0355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(C=C(O)CJ)"""),
)

species(
    label = 'CC(O)[C]=C[C]=O(6406)',
    structure = SMILES('CC(O)[C]=C[C]=O'),
    E0 = (85.6467,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.589614,'amu*angstrom^2'), symmetry=1, barrier=(13.5564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.589732,'amu*angstrom^2'), symmetry=1, barrier=(13.5591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.589523,'amu*angstrom^2'), symmetry=1, barrier=(13.5543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.589926,'amu*angstrom^2'), symmetry=1, barrier=(13.5636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.938233,0.0734944,-9.53377e-05,6.89363e-08,-2.05817e-11,10405.7,27.4286], Tmin=(100,'K'), Tmax=(808.472,'K')), NASAPolynomial(coeffs=[9.63505,0.0304642,-1.54984e-05,3.09815e-09,-2.22077e-13,8999.52,-12.6791], Tmin=(808.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(C=CCJ=O)"""),
)

species(
    label = 'CC(=O)[C]=CC[O](6105)',
    structure = SMILES('CC([O])=C=CC[O]'),
    E0 = (106.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,180,180,3259.98],'cm^-1')),
        HinderedRotor(inertia=(0.418015,'amu*angstrom^2'), symmetry=1, barrier=(9.61098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.418255,'amu*angstrom^2'), symmetry=1, barrier=(9.6165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29486,0.0649342,-8.72606e-05,7.4083e-08,-2.59591e-11,12900.2,25.7155], Tmin=(100,'K'), Tmax=(808.554,'K')), NASAPolynomial(coeffs=[5.85469,0.0344025,-1.58267e-05,2.98776e-09,-2.05742e-13,12423.5,6.29782], Tmin=(808.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)=C=CC[O](6407)',
    structure = SMILES('C=C(O)[C]=CC[O]'),
    E0 = (115.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1685,370,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.564809,'amu*angstrom^2'), symmetry=1, barrier=(12.9861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.564643,'amu*angstrom^2'), symmetry=1, barrier=(12.9823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.9629,'amu*angstrom^2'), symmetry=1, barrier=(91.1148,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.75812,0.0731697,-9.16994e-05,6.31982e-08,-1.74356e-11,14014.6,24.8969], Tmin=(100,'K'), Tmax=(886.799,'K')), NASAPolynomial(coeffs=[11.4675,0.0248632,-9.98885e-06,1.7701e-09,-1.1798e-13,12115.2,-25.4831], Tmin=(886.799,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CJC=C)"""),
)

species(
    label = 'CC([O])=[C]C=CO(5357)',
    structure = SMILES('CC([O])=[C]C=CO'),
    E0 = (-28.6553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03859,'amu*angstrom^2'), symmetry=1, barrier=(23.8792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04236,'amu*angstrom^2'), symmetry=1, barrier=(23.9658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04283,'amu*angstrom^2'), symmetry=1, barrier=(23.9768,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.794495,0.082309,-9.15949e-05,4.78442e-08,-9.10871e-12,-3253.9,28.024], Tmin=(100,'K'), Tmax=(1537.86,'K')), NASAPolynomial(coeffs=[22.2025,0.0040977,2.63382e-06,-7.82469e-10,6.06784e-14,-8151.85,-85.7487], Tmin=(1537.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.6553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C(O)=C=C[CH]O(6408)',
    structure = SMILES('[CH2]C(O)=[C]C=CO'),
    E0 = (-7.54347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.15532,'amu*angstrom^2'), symmetry=1, barrier=(26.563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15194,'amu*angstrom^2'), symmetry=1, barrier=(26.4854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15365,'amu*angstrom^2'), symmetry=1, barrier=(26.5248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15138,'amu*angstrom^2'), symmetry=1, barrier=(26.4724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.71522,0.0928755,-0.000108217,5.68818e-08,-1.06769e-11,-673.167,30.5215], Tmin=(100,'K'), Tmax=(1599.7,'K')), NASAPolynomial(coeffs=[25.5084,-0.00234577,6.52757e-06,-1.54657e-09,1.12258e-13,-5909.3,-102.749], Tmin=(1599.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.54347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(O)CJ) + radical(C=CJC=C)"""),
)

species(
    label = 'CH2(S)(21)(22)',
    structure = SMILES('[CH2]'),
    E0 = (418.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1358.21,2621.43,3089.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19331,-0.00233105,8.15676e-06,-6.62986e-09,1.93233e-12,50366.2,-0.746734], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.13502,0.00289594,-8.16668e-07,1.13573e-10,-6.36263e-15,50504.1,4.06031], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(418.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=CC=C=CO(6409)',
    structure = SMILES('O=CC=C=CO'),
    E0 = (-149.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(1.15186,'amu*angstrom^2'), symmetry=1, barrier=(26.4835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1532,'amu*angstrom^2'), symmetry=1, barrier=(26.5143,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25519,0.05143,-4.02864e-05,6.80829e-09,3.07965e-12,-17880.3,19.0076], Tmin=(100,'K'), Tmax=(998.283,'K')), NASAPolynomial(coeffs=[16.745,0.00799049,-3.00262e-06,6.00153e-10,-4.62935e-14,-21901,-60.3438], Tmin=(998.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CC(O)=C[C]C=O(6412)',
    structure = SMILES('CC(O)=C[C]C=O'),
    E0 = (79.4141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.779071,'amu*angstrom^2'), symmetry=1, barrier=(17.9124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777529,'amu*angstrom^2'), symmetry=1, barrier=(17.8769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.55766,'amu*angstrom^2'), symmetry=1, barrier=(58.8055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57371,'amu*angstrom^2'), symmetry=1, barrier=(36.1826,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02913,0.0671072,-7.02132e-05,4.01889e-08,-9.44086e-12,9656.86,34.1533], Tmin=(100,'K'), Tmax=(1020.45,'K')), NASAPolynomial(coeffs=[11.0908,0.0276666,-1.22373e-05,2.31247e-09,-1.61439e-13,7603.39,-14.5922], Tmin=(1020.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.4141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'CO(10)(11)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(762.44,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35308e-07,1.51269e-10,-9.88872e-15,-14292.7,6.51157], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=C=C(C)O(6410)',
    structure = SMILES('C=C=C(C)O'),
    E0 = (-67.8562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(0.902443,'amu*angstrom^2'), symmetry=1, barrier=(20.7489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900607,'amu*angstrom^2'), symmetry=1, barrier=(20.7067,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69403,0.0423221,-1.88022e-05,-1.16794e-08,9.29029e-12,-8070.53,16.3593], Tmin=(100,'K'), Tmax=(945.589,'K')), NASAPolynomial(coeffs=[13.6968,0.0121323,-3.56472e-06,5.98969e-10,-4.24187e-14,-11260.7,-45.7418], Tmin=(945.589,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.8562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'S(247)(246)',
    structure = SMILES('CC(=O)C=CC=O'),
    E0 = (-249.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,338.299],'cm^-1')),
        HinderedRotor(inertia=(0.102894,'amu*angstrom^2'), symmetry=1, barrier=(8.35729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1029,'amu*angstrom^2'), symmetry=1, barrier=(8.35732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102901,'amu*angstrom^2'), symmetry=1, barrier=(8.35737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3947.01,'J/mol'), sigma=(6.07703,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=616.51 K, Pc=39.91 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72359,0.0430426,-2.26393e-05,4.45245e-09,-2.82809e-13,-30034.5,19.0288], Tmin=(100,'K'), Tmax=(2443.83,'K')), NASAPolynomial(coeffs=[29.6812,0.00674437,-5.16283e-06,9.95166e-10,-6.31701e-14,-45547.2,-139.895], Tmin=(2443.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'CC(O)=C=C=CO(6411)',
    structure = SMILES('CC(O)=C=C=CO'),
    E0 = (-136.07,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(1.10911,'amu*angstrom^2'), symmetry=1, barrier=(25.5006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1103,'amu*angstrom^2'), symmetry=1, barrier=(25.528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1107,'amu*angstrom^2'), symmetry=1, barrier=(25.5371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00863,0.0856659,-9.57221e-05,4.91098e-08,-9.18995e-12,-16164.1,26.7861], Tmin=(100,'K'), Tmax=(1548.43,'K')), NASAPolynomial(coeffs=[24.6736,0.00158147,2.91807e-06,-7.58233e-10,5.60698e-14,-21990.7,-101.478], Tmin=(1548.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    label = 'Ar(8)',
    structure = SMILES('[Ar]'),
    E0 = (-6.19426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,-745,4.3663], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,-745,4.3663], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.19426,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: BurkeH2O2"""),
)

transitionState(
    label = 'TS1',
    E0 = (189.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (158.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (198.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (167.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (119.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (217.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (175.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (183.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (164.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (206.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-2.13613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (81.9323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-7.94108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (2.02494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (41.7666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (56.8969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (78.0088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (110.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (131.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (140.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-3.68206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (17.4298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (269.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (136.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (118.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-109.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (7.79664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction29',
    reactants = ['OH(5)(5)', 'C[C]=C=CC=O(6394)'],
    products = ['S(252)(251)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.11458e+06,'m^3/(mol*s)'), n=0.373849, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;Cd_allenic] + [O_pri_rad;Y_rad] for rate rule [O_pri_rad;Cd_allenic]
Euclidian distance = 2.0
family: R_Recombination
Ea raised from -5.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)(3)', 'CC(=O)[C]=CC=O(6096)'],
    products = ['S(252)(251)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH3(15)(16)', 'O=CC=C=[C]O(6395)'],
    products = ['S(252)(251)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.705e+09,'cm^3/(mol*s)'), n=1.07, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 48 used for C_methyl;Cd_allenic
Exact match found for rate rule [Cd_allenic;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -9.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)(3)', '[CH2]C(O)=C=CC=O(6396)'],
    products = ['S(252)(251)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.156e+12,'cm^3/(mol*s)'), n=0.461, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 15 used for C_rad/H2/Cd;H_rad
Exact match found for rate rule [C_rad/H2/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['HCO(14)(15)', '[CH]=C=C(C)O(6397)'],
    products = ['S(252)(251)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.81e+13,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [Cd_rad;CO_pri_rad] for rate rule [Cd_allenic;CO_pri_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)(3)', 'CC(O)=C=[C]C=O(6398)'],
    products = ['S(252)(251)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)(3)', 'CC(O)=C=C[C]=O(6399)'],
    products = ['S(252)(251)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_sec_rad;Y_rad] for rate rule [CO_rad/OneDe;H_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CC([O])[C]=CC=O(6100)'],
    products = ['S(252)(251)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(O)[C]=CC=O(6400)'],
    products = ['S(252)(251)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CC(O)=C=[C]C[O](6401)'],
    products = ['S(252)(251)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C[C](O)C=[C]C=O(6107)'],
    products = ['S(252)(251)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CC(O)=[C]C[C]=O(6402)'],
    products = ['S(252)(251)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CC([O])=CC=C[O](5355)'],
    products = ['S(252)(251)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(O)=C[CH]C=O(6403)'],
    products = ['S(252)(251)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C[C](O)C=C[C]=O(6109)'],
    products = ['S(252)(251)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['CC([O])=[C]CC=O(6404)'],
    products = ['S(252)(251)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C(O)=[C]CC=O(6405)'],
    products = ['S(252)(251)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['CC(O)[C]=C[C]=O(6406)'],
    products = ['S(252)(251)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['CC(=O)[C]=CC[O](6105)'],
    products = ['S(252)(251)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(O)=C=CC[O](6407)'],
    products = ['S(252)(251)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['CC([O])=[C]C=CO(5357)'],
    products = ['S(252)(251)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C(O)=C=C[CH]O(6408)'],
    products = ['S(252)(251)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['CH2(S)(21)(22)', 'O=CC=C=CO(6409)'],
    products = ['S(252)(251)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(884144,'m^3/(mol*s)'), n=0.112, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cd_sec] for rate rule [carbene;Cd/H/NonDeO]
Euclidian distance = 1.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction54',
    reactants = ['CC(O)=C[C]C=O(6412)'],
    products = ['S(252)(251)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.14647e+14,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CsJ2-C;CsJ2C;CH=C]
Euclidian distance = 0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction52',
    reactants = ['CO(10)(11)', 'C=C=C(C)O(6410)'],
    products = ['S(252)(251)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction27',
    reactants = ['S(252)(251)'],
    products = ['S(247)(246)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(104,'s^-1'), n=3.21, Ea=(82.0482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction53',
    reactants = ['CC(O)=C=C=CO(6411)'],
    products = ['S(252)(251)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

network(
    label = '388',
    isomers = [
        'S(252)(251)',
    ],
    reactants = [
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '388',
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

