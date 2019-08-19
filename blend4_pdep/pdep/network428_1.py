species(
    label = 'S(573)(572)',
    structure = SMILES('O=CC1=C(O)C2CC(=O)C12'),
    E0 = (-350.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(5120.61,'J/mol'), sigma=(7.38277,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=799.83 K, Pc=28.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.639949,0.078475,-6.94413e-05,3.18768e-08,-6.124e-12,-42057.3,23.6853], Tmin=(100,'K'), Tmax=(1201.89,'K')), NASAPolynomial(coeffs=[12.8628,0.0377962,-1.86729e-05,3.71657e-09,-2.66536e-13,-44995.4,-37.5309], Tmin=(1201.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-350.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1)"""),
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
    label = 'O=CC1=[C]C2CC(=O)C12(7389)',
    structure = SMILES('O=CC1=[C]C2CC(=O)C12'),
    E0 = (116.388,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16437,0.0698016,-9.25007e-05,8.76215e-08,-3.60381e-11,14093.2,22.2096], Tmin=(100,'K'), Tmax=(717.751,'K')), NASAPolynomial(coeffs=[3.30304,0.04767,-2.49054e-05,5.01291e-09,-3.59676e-13,14049.2,14.4336], Tmin=(717.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(cyclobutene-vinyl)"""),
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
    label = 'O=C[C]1C(=O)C2CC(=O)C12(7336)',
    structure = SMILES('[O]C1=C(C=O)C2C(=O)CC12'),
    E0 = (-212.852,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96656,0.0746056,-7.58449e-05,4.59121e-08,-1.23495e-11,-25497.7,24.0404], Tmin=(100,'K'), Tmax=(860.901,'K')), NASAPolynomial(coeffs=[7.57091,0.0439195,-2.23779e-05,4.50767e-09,-3.25724e-13,-26634.8,-6.83245], Tmin=(860.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-212.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=CC1=C(O)[C]2CC(=O)C21(7390)',
    structure = SMILES('O=CC1=C(O)[C]2CC(=O)C21'),
    E0 = (-198.083,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.869878,0.0723997,-6.03803e-05,2.53479e-08,-4.42206e-12,-23714.5,23.7415], Tmin=(100,'K'), Tmax=(1311.82,'K')), NASAPolynomial(coeffs=[13.2385,0.0346851,-1.72555e-05,3.43179e-09,-2.45379e-13,-26959.6,-39.287], Tmin=(1311.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(CCJ(C)CO)"""),
)

species(
    label = 'O=CC1[C]2C(=O)CC2C=1O(7391)',
    structure = SMILES('O=CC1[C]2C(=O)CC2C=1O'),
    E0 = (-230.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.211578,0.0720795,-5.4846e-05,1.97095e-08,-2.78582e-12,-27631.9,25.4036], Tmin=(100,'K'), Tmax=(1676.7,'K')), NASAPolynomial(coeffs=[21.1955,0.0220196,-1.00617e-05,1.90303e-09,-1.30838e-13,-34668.6,-86.6766], Tmin=(1676.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-230.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(C=CCJ(C)C=O)"""),
)

species(
    label = 'O=CC1=C(O)C2[CH]C(=O)C12(7392)',
    structure = SMILES('O=CC1=C(O)C2[CH]C(=O)C12'),
    E0 = (-188.528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.178798,0.0778256,-6.65935e-05,2.74822e-08,-4.50251e-12,-22532.2,25.7296], Tmin=(100,'K'), Tmax=(1449.91,'K')), NASAPolynomial(coeffs=[19.3524,0.0249309,-1.18727e-05,2.32221e-09,-1.64395e-13,-28092.3,-73.8951], Tmin=(1449.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(CCJC=O)"""),
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
    label = 'O=C1CC2C(O)=[C]C12(7393)',
    structure = SMILES('O=C1CC2C(O)=[C]C12'),
    E0 = (32.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21066,0.0602552,-5.31942e-05,2.46494e-08,-4.67595e-12,4026.33,19.482], Tmin=(100,'K'), Tmax=(1244.82,'K')), NASAPolynomial(coeffs=[12.0867,0.0253069,-1.10814e-05,2.09565e-09,-1.46394e-13,1318.62,-35.3701], Tmin=(1244.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'O=[C]C1=C(O)C2CC(=O)C12(7394)',
    structure = SMILES('O=[C]C1=C(O)C2CC(=O)C12'),
    E0 = (-187.053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582938,0.0814812,-8.45115e-05,4.68579e-08,-1.08529e-11,-22379.4,24.4664], Tmin=(100,'K'), Tmax=(1018.11,'K')), NASAPolynomial(coeffs=[11.9498,0.0368219,-1.87131e-05,3.77185e-09,-2.72784e-13,-24693.9,-30.5758], Tmin=(1018.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-187.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(C=C(C)CJ=O)"""),
)

species(
    label = '[O]C1[C](C=O)C2C(=O)CC12(7341)',
    structure = SMILES('[O]C=C1C([O])C2CC(=O)C12'),
    E0 = (94.0538,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.336155,0.0641445,-1.36917e-05,-3.11486e-08,1.67129e-11,11458.4,28.5539], Tmin=(100,'K'), Tmax=(1015.85,'K')), NASAPolynomial(coeffs=[19.4204,0.0254637,-1.04194e-05,2.03931e-09,-1.50536e-13,5699.51,-73.0774], Tmin=(1015.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.0538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsOsH) + polycyclic(s2_4_4_ane) + radical(CC(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O=C[C]1C(O)[C]2CC(=O)C12(7395)',
    structure = SMILES('[O]C=C1C(O)[C]2CC(=O)C21'),
    E0 = (16.2669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.45173,0.0614163,-8.03242e-06,-3.60576e-08,1.82795e-11,2098.8,29.219], Tmin=(100,'K'), Tmax=(1012.12,'K')), NASAPolynomial(coeffs=[18.9728,0.0256271,-1.0431e-05,2.03936e-09,-1.50622e-13,-3566.31,-69.8228], Tmin=(1012.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.2669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsOsH) + polycyclic(s2_4_4_ane) + radical(C=COJ) + radical(CCJ(C)CO)"""),
)

species(
    label = 'O=CC1[C](O)C2CC(=O)[C]12(7352)',
    structure = SMILES('O=CC1[C](O)C2CC(=O)[C]12'),
    E0 = (63.7475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.19753,0.0905621,-0.000122555,1.01549e-07,-3.4821e-11,7797.34,27.4927], Tmin=(100,'K'), Tmax=(790.859,'K')), NASAPolynomial(coeffs=[7.94815,0.0440707,-2.05494e-05,3.90528e-09,-2.70195e-13,6799.4,-6.63963], Tmin=(790.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.7475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + polycyclic(s2_4_4_ane) + radical(C2CsJOH) + radical(C2CJCHO)"""),
)

species(
    label = '[O]C1CC2[C]1C(C=O)=C2O(7396)',
    structure = SMILES('[O]C1CC2[C]1C(C=O)=C2O'),
    E0 = (-25.7423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.103776,0.0698338,-2.84462e-05,-1.67163e-08,1.15732e-11,-2942.09,28.0668], Tmin=(100,'K'), Tmax=(1044.54,'K')), NASAPolynomial(coeffs=[20.5037,0.0246144,-1.07565e-05,2.14817e-09,-1.59007e-13,-8998.66,-79.8318], Tmin=(1044.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.7423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(CCJ(C)CO) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1[CH]C2C(O)=C(C=O)C12(7397)',
    structure = SMILES('[O]C1[CH]C2C(O)=C(C=O)C12'),
    E0 = (21.5858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0260061,0.0712999,-3.13951e-05,-1.59904e-08,1.2003e-11,2753.3,29.9526], Tmin=(100,'K'), Tmax=(1028.71,'K')), NASAPolynomial(coeffs=[21.3567,0.0223148,-9.48068e-06,1.8947e-09,-1.41323e-13,-3432.03,-82.2921], Tmin=(1028.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.5858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = 'O=[C]C1[C](O)C2CC(=O)C12(7355)',
    structure = SMILES('O=[C]C1[C](O)C2CC(=O)C12'),
    E0 = (65.6084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.894172,0.0747296,-6.97411e-05,3.8041e-08,-9.13503e-12,7997.33,28.0727], Tmin=(100,'K'), Tmax=(959.288,'K')), NASAPolynomial(coeffs=[8.38146,0.0435095,-2.09237e-05,4.11493e-09,-2.93588e-13,6560.83,-7.73794], Tmin=(959.288,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.6084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + polycyclic(s2_4_4_ane) + radical(CC(C)CJ=O) + radical(C2CsJOH)"""),
)

species(
    label = 'O=C[C]1[C]2C(=O)CC2C1O(7398)',
    structure = SMILES('[O]C=C1[C]2C(=O)CC2C1O'),
    E0 = (-16.6065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.661667,0.0514112,2.89453e-05,-7.93178e-08,3.47362e-11,-1857.4,27.7273], Tmin=(100,'K'), Tmax=(979.238,'K')), NASAPolynomial(coeffs=[20.7962,0.0224642,-8.35672e-06,1.6602e-09,-1.27824e-13,-8356.13,-82.0358], Tmin=(979.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.6065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsOsH) + polycyclic(s2_4_4_ane) + radical(C=CCJ(C)C=O) + radical(C=COJ)"""),
)

species(
    label = 'O=C[C]1C(O)C2[CH]C(=O)C12(7399)',
    structure = SMILES('[O]C=C1C(O)C2C=C([O])C12'),
    E0 = (-118.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.142293,0.0733981,-2.7667e-05,-3.00992e-08,2.0374e-11,-14077.3,28.8157], Tmin=(100,'K'), Tmax=(954.331,'K')), NASAPolynomial(coeffs=[23.8527,0.0157303,-4.46306e-06,8.00874e-10,-6.16404e-14,-20610.9,-96.0608], Tmin=(954.331,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + polycyclic(s2_4_4_ene_1) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=CC1[C](O)[C]2CC(=O)C21(7359)',
    structure = SMILES('O=CC1[C](O)[C]2CC(=O)C21'),
    E0 = (59.4936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23826,0.0680258,-5.22423e-05,2.24754e-08,-4.39117e-12,7248.7,26.5591], Tmin=(100,'K'), Tmax=(1111.82,'K')), NASAPolynomial(coeffs=[7.47987,0.0455709,-2.19482e-05,4.31099e-09,-3.06888e-13,5860.76,-4.21477], Tmin=(1111.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.4936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + polycyclic(s2_4_4_ane) + radical(C2CsJOH) + radical(CCJ(C)CO)"""),
)

species(
    label = '[O]C1C[C]2C(O)=C(C=O)C21(7400)',
    structure = SMILES('[O]C1C[C]2C(O)=C(C=O)C21'),
    E0 = (-25.7423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.103776,0.0698338,-2.84462e-05,-1.67163e-08,1.15732e-11,-2942.09,28.0668], Tmin=(100,'K'), Tmax=(1044.54,'K')), NASAPolynomial(coeffs=[20.5037,0.0246144,-1.07565e-05,2.14817e-09,-1.59007e-13,-8998.66,-79.8318], Tmin=(1044.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.7423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(CC(C)OJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[O]CC1[C]2C(=O)CC2C=1O(7401)',
    structure = SMILES('[O]CC1[C]2C(=O)CC2C=1O'),
    E0 = (-66.6359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.731865,0.0703543,-5.17437e-05,1.93248e-08,-2.98445e-12,-7896.08,26.2475], Tmin=(100,'K'), Tmax=(1475,'K')), NASAPolynomial(coeffs=[13.9093,0.0346185,-1.54019e-05,2.89898e-09,-2.00389e-13,-11783.4,-42.4473], Tmin=(1475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.6359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + polycyclic(s2_4_4_ene_1) + radical(C=CCJ(C)C=O) + radical(CCOJ)"""),
)

species(
    label = 'O=CC1[C](O)C2[CH]C(=O)C21(7353)',
    structure = SMILES('[O]C1=CC2[C](O)C(C=O)C12'),
    E0 = (-75.7787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.156161,0.0800633,-7.35385e-05,3.50866e-08,-6.73255e-12,-8971.99,29.618], Tmin=(100,'K'), Tmax=(1249.56,'K')), NASAPolynomial(coeffs=[16.2016,0.0287,-1.18813e-05,2.19135e-09,-1.51233e-13,-12982,-51.367], Tmin=(1249.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.7787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + polycyclic(s2_4_4_ene_1) + radical(C=C(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'O=CC1=C(O)[C]2C[C](O)C21(7402)',
    structure = SMILES('O=CC1=C(O)[C]2C[C](O)C21'),
    E0 = (-79.4753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2782.5,750,1395,475,1775,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.279081,0.083163,-7.26218e-05,3.07672e-08,-5.12309e-12,-9395.58,29.7612], Tmin=(100,'K'), Tmax=(1446.41,'K')), NASAPolynomial(coeffs=[21.5673,0.0227479,-9.96877e-06,1.88988e-09,-1.3193e-13,-15715.4,-83.6985], Tmin=(1446.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.4753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(C2CsJOH) + radical(CCJ(C)CO)"""),
)

species(
    label = 'O=C1CC2[C]1C([CH]O)=C2O(7403)',
    structure = SMILES('O=C1CC2[C]1C([CH]O)=C2O'),
    E0 = (-175.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.191239,0.0666076,-1.57415e-05,-3.42803e-08,1.94172e-11,-20900.7,27.9271], Tmin=(100,'K'), Tmax=(983.556,'K')), NASAPolynomial(coeffs=[20.9444,0.0218912,-8.06685e-06,1.53986e-09,-1.145e-13,-26902.6,-81.6087], Tmin=(983.556,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-175.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + polycyclic(s2_4_4_ene_1) + radical(C=CCJO) + radical(C=CCJ(C)C=O)"""),
)

species(
    label = '[O]CC1=C(O)[C]2CC(=O)C21(7404)',
    structure = SMILES('[O]CC1=C(O)[C]2CC(=O)C21'),
    E0 = (-33.7625,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803989,0.0767317,-7.49754e-05,4.3894e-08,-1.132e-11,-3951.08,26.7518], Tmin=(100,'K'), Tmax=(900.826,'K')), NASAPolynomial(coeffs=[8.06199,0.0445029,-2.13089e-05,4.17669e-09,-2.97323e-13,-5258.69,-7.50566], Tmin=(900.826,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.7625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + polycyclic(s2_4_4_ene_1) + radical(CCJ(C)CO) + radical(CCOJ)"""),
)

species(
    label = '[O]C1CC2C(O)=C([C]=O)C12(7405)',
    structure = SMILES('[O]C1CC2C(O)=C([C]=O)C12'),
    E0 = (-14.7125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.255538,0.0800164,-5.73085e-05,1.16759e-08,2.06467e-12,-1604.7,29.0313], Tmin=(100,'K'), Tmax=(1086.41,'K')), NASAPolynomial(coeffs=[21.4084,0.0230263,-1.00654e-05,1.98029e-09,-1.44247e-13,-7655.83,-83.4655], Tmin=(1086.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.7125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(C=C(C)CJ=O) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C[C]1C(=O)C2CC(=O)C12(7342)',
    structure = SMILES('[O]CC1=C([O])C2CC(=O)C12'),
    E0 = (-48.5316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.519427,0.0842594,-0.000112818,9.87463e-08,-3.62448e-11,-5719.1,28.3616], Tmin=(100,'K'), Tmax=(777.407,'K')), NASAPolynomial(coeffs=[5.42681,0.0484408,-2.33145e-05,4.50504e-09,-3.14938e-13,-6162.75,7.97592], Tmin=(777.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.5316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + polycyclic(s2_4_4_ene_1) + radical(CCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=CC1=C(O)C2[CH][C](O)C12(7406)',
    structure = SMILES('O=CC1=C(O)C2[CH][C](O)C12'),
    E0 = (-32.1472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2782.5,750,1395,475,1775,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.41243,0.0852829,-7.78555e-05,3.44597e-08,-5.95864e-12,-3697.76,31.8466], Tmin=(100,'K'), Tmax=(1403.46,'K')), NASAPolynomial(coeffs=[22.2568,0.0206743,-8.80357e-06,1.65932e-09,-1.15953e-13,-10060.9,-85.2029], Tmin=(1403.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.1472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(C2CsJOH) + radical(CCJCO)"""),
)

species(
    label = 'O=CC1[C]2[C](O)CC2C=1O(7407)',
    structure = SMILES('O=CC1[C]2[C](O)CC2C=1O'),
    E0 = (-79.4753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2782.5,750,1395,475,1775,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.279081,0.083163,-7.26218e-05,3.07672e-08,-5.12309e-12,-9395.58,29.7612], Tmin=(100,'K'), Tmax=(1446.41,'K')), NASAPolynomial(coeffs=[21.5673,0.0227479,-9.96877e-06,1.88988e-09,-1.3193e-13,-15715.4,-83.6985], Tmin=(1446.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.4753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(CCJ(C)CO) + radical(C2CsJOH)"""),
)

species(
    label = 'O=C1C[C]2C(O)=C([CH]O)C21(7408)',
    structure = SMILES('O=C1C[C]2C(O)=C([CH]O)C21'),
    E0 = (-142.173,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0492218,0.0769572,-5.38377e-05,1.02797e-08,2.4822e-12,-16943.2,29.5292], Tmin=(100,'K'), Tmax=(1052.04,'K')), NASAPolynomial(coeffs=[19.3564,0.0246653,-9.92165e-06,1.86797e-09,-1.33114e-13,-22215.6,-70.7287], Tmin=(1052.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-142.173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + polycyclic(s2_4_4_ene_1) + radical(CCJ(C)CO) + radical(C=CCJO)"""),
)

species(
    label = 'O=[C]C1=C(O)C2C[C](O)C12(7409)',
    structure = SMILES('O=[C]C1=C(O)C2C[C](O)C12'),
    E0 = (-68.4454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1855,455,950,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.252407,0.0890189,-8.74202e-05,4.2454e-08,-8.15909e-12,-8075.27,29.3262], Tmin=(100,'K'), Tmax=(1255.77,'K')), NASAPolynomial(coeffs=[19.6336,0.0256768,-1.17598e-05,2.28766e-09,-1.62824e-13,-13069.8,-71.1415], Tmin=(1255.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.4454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(C2CsJOH) + radical(C=C(C)CJ=O)"""),
)

species(
    label = '[O]C1=C([CH]O)C2C(=O)CC12(7410)',
    structure = SMILES('[O]C1=C([CH]O)C2C(=O)CC12'),
    E0 = (-156.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.168791,0.0781484,-6.73004e-05,2.96589e-08,-5.25624e-12,-18732.8,29.3605], Tmin=(100,'K'), Tmax=(1344.67,'K')), NASAPolynomial(coeffs=[16.8562,0.028508,-1.19256e-05,2.20485e-09,-1.51981e-13,-23220.6,-56.0885], Tmin=(1344.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-156.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + polycyclic(s2_4_4_ene_1) + radical(C=C(C)OJ) + radical(C=CCJO)"""),
)

species(
    label = '[O]C1CC2C(=O)[C](C=O)C12(7345)',
    structure = SMILES('[O]C1=C(C=O)C2C([O])CC12'),
    E0 = (-40.5115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.121544,0.0732616,-4.91838e-05,1.13748e-08,4.11445e-13,-4722.76,28.6253], Tmin=(100,'K'), Tmax=(1186.14,'K')), NASAPolynomial(coeffs=[18.6049,0.0274798,-1.22161e-05,2.35986e-09,-1.67706e-13,-10271.7,-68.6086], Tmin=(1186.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.5115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(CC(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]CC1=C(O)C2[CH]C(=O)C12(7411)',
    structure = SMILES('[O]CC1=C(O)C2[CH]C(=O)C12'),
    E0 = (-24.2077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.518634,0.0778409,-6.8081e-05,3.14236e-08,-6.02341e-12,-2787.37,27.2511], Tmin=(100,'K'), Tmax=(1217.7,'K')), NASAPolynomial(coeffs=[13.2693,0.0359558,-1.64851e-05,3.17559e-09,-2.23873e-13,-5892.64,-36.775], Tmin=(1217.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.2077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + polycyclic(s2_4_4_ene_1) + radical(CCOJ) + radical(CCJC=O)"""),
)

species(
    label = 'O=C[C]1C(=O)C2C[C](O)C12(7351)',
    structure = SMILES('[O]C1=C(C=O)C2[C](O)CC12'),
    E0 = (-94.2444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.367505,0.0795603,-7.05812e-05,3.18943e-08,-5.87694e-12,-11204.1,28.0383], Tmin=(100,'K'), Tmax=(1276.59,'K')), NASAPolynomial(coeffs=[15.4196,0.032397,-1.51644e-05,2.95434e-09,-2.09537e-13,-15047.2,-48.2551], Tmin=(1276.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.2444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s2_4_4_ene_1) + radical(C=C(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'O=C1[CH]C2C(O)=C([CH]O)C12(7412)',
    structure = SMILES('O=C1[CH]C2C(O)=C([CH]O)C12'),
    E0 = (-132.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.176349,0.0759464,-3.85488e-05,-1.40353e-08,1.30763e-11,-15785.7,29.4801], Tmin=(100,'K'), Tmax=(986.291,'K')), NASAPolynomial(coeffs=[22.4506,0.0197657,-7.22632e-06,1.37426e-09,-1.02047e-13,-21979.9,-88.144], Tmin=(986.291,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + polycyclic(s2_4_4_ene_1) + radical(CCJC=O) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]C1C(O)=C(C=O)C1[C]=O(7414)',
    structure = SMILES('[CH2]C1C(O)=C(C=O)C1[C]=O'),
    E0 = (17.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,2750,3150,900,1100,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.566323,0.088287,-8.5155e-05,4.04593e-08,-7.48287e-12,2297.07,35.4034], Tmin=(100,'K'), Tmax=(1323.01,'K')), NASAPolynomial(coeffs=[22.0779,0.0198242,-7.53317e-06,1.34548e-09,-9.17858e-14,-3694.62,-80.18], Tmin=(1323.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(Isobutyl) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'O=[C]CC1[CH]C(C=O)=C1O(7415)',
    structure = SMILES('O=[C]CC1[CH]C(C=O)=C1O'),
    E0 = (20.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,2750,3150,900,1100,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.29662,0.0821242,-6.67147e-05,2.07375e-08,-6.28182e-13,2647,36.3062], Tmin=(100,'K'), Tmax=(1075.56,'K')), NASAPolynomial(coeffs=[21.3735,0.0210511,-8.76067e-06,1.6874e-09,-1.21839e-13,-3143.43,-75.066], Tmin=(1075.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C[C]=C(O)C1[CH]C(=O)C1(7416)',
    structure = SMILES('[O]C=C=C(O)C1[CH]C(=O)C1'),
    E0 = (-22.8608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.171844,0.0649271,1.78256e-05,-9.03389e-08,4.48487e-11,-2574.62,33.028], Tmin=(100,'K'), Tmax=(936.172,'K')), NASAPolynomial(coeffs=[29.1229,0.0078534,1.6772e-07,-6.82609e-11,-5.99986e-15,-11043.6,-122.307], Tmin=(936.172,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.8608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclobutanone) + radical(C=COJ) + radical(CCJC=O)"""),
)

species(
    label = 'O=CC(=[C]O)C1[CH]CC1=O(7417)',
    structure = SMILES('O=CC(=[C]O)C1[CH]CC1=O'),
    E0 = (57.4447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,350,440,435,1725,2782.5,750,1395,475,1775,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.150705,0.078144,-6.4244e-05,2.57282e-08,-4.05163e-12,7069.44,37.9544], Tmin=(100,'K'), Tmax=(1526.81,'K')), NASAPolynomial(coeffs=[21.2386,0.022107,-9.19056e-06,1.68948e-09,-1.1549e-13,538.023,-74.2878], Tmin=(1526.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.4447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-OdCsCs) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutanone) + radical(C=CJO) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]C(=O)C1[CH]C(O)=C1C=O(7418)',
    structure = SMILES('C=C([O])C1C=C(O)C1=C[O]'),
    E0 = (-62.0466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,3010,987.5,1337.5,450,1655,2750,3150,900,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.79246,0.108652,-0.000115244,5.61149e-08,-9.90773e-12,-7182.71,37.1273], Tmin=(100,'K'), Tmax=(1672.54,'K')), NASAPolynomial(coeffs=[30.0022,0.00314861,3.65525e-06,-9.55861e-10,6.93929e-14,-14366.1,-126.636], Tmin=(1672.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.0466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C1=C[C](C=O)C(O)=CC1(7220)',
    structure = SMILES('[O]C=C1C=C([O])CC=C1O'),
    E0 = (-204.926,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.399781,0.075623,-2.01323e-05,-4.8084e-08,2.95152e-11,-24469,27.9206], Tmin=(100,'K'), Tmax=(932.398,'K')), NASAPolynomial(coeffs=[27.5204,0.00972261,-7.90772e-07,6.02684e-11,-1.0247e-14,-32017.5,-117.383], Tmin=(932.398,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(14cyclohexadiene3methylene) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CH2CO(26)(27)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=CC1=CC=C1O(7419)',
    structure = SMILES('O=CC1=CC=C1O'),
    E0 = (95.6777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,3150,900,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.821255,0.0574881,-3.67553e-05,-6.17105e-09,9.54463e-12,11633,18.6539], Tmin=(100,'K'), Tmax=(969.725,'K')), NASAPolynomial(coeffs=[20.2627,0.00607611,-1.74996e-06,3.70325e-10,-3.23728e-14,6509.19,-81.5197], Tmin=(969.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.6777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + ring(cyclobutadiene_13)"""),
)

species(
    label = 'S(599)(598)',
    structure = SMILES('O=CC1=CC(=O)CC=C1O'),
    E0 = (-347.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(5238.5,'J/mol'), sigma=(7.47626,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=818.24 K, Pc=28.44 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.240181,0.0799814,-5.77922e-05,1.19872e-08,2.12357e-12,-41675.6,27.5419], Tmin=(100,'K'), Tmax=(1072.24,'K')), NASAPolynomial(coeffs=[21.1684,0.0228702,-9.72789e-06,1.89403e-09,-1.375e-13,-47574.6,-83.3346], Tmin=(1072.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-347.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(24cyclohexadienone)"""),
)

species(
    label = 'O=C[C]1C2[C]3CC2C1(O)O3(7426)',
    structure = SMILES('[O]C=C1C2[C]3CC2C1(O)O3'),
    E0 = (100.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.151956,0.0709992,-2.13607e-05,-2.89595e-08,1.59867e-11,12313.7,22.7739], Tmin=(100,'K'), Tmax=(1070.18,'K')), NASAPolynomial(coeffs=[25.0277,0.0202407,-1.09835e-05,2.43079e-09,-1.89192e-13,4441.62,-112.011], Tmin=(1070.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + polycyclic(s3_4_5_ane) + polycyclic(s2_4_5_ane) + polycyclic(s2_4_4_ane) - ring(Cyclobutane) - ring(Cyclobutane) - ring(Tetrahydrofuran) + radical(C=COJ) + radical(C2CsJOCs)"""),
)

species(
    label = 'O=CC12O[C]3CC([C]1O)C32(7427)',
    structure = SMILES('O=CC12O[C]3CC([C]1O)C32'),
    E0 = (248.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.12721,0.0977992,-0.000158365,1.61682e-07,-6.74901e-11,30038.9,22.8211], Tmin=(100,'K'), Tmax=(729.823,'K')), NASAPolynomial(coeffs=[3.2513,0.0589188,-3.37357e-05,6.98846e-09,-5.07229e-13,30162.3,12.703], Tmin=(729.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + polycyclic(s2_4_4_ane) + polycyclic(s2_4_4_ane) + polycyclic(s2_4_4_ane) - ring(Oxetane) - ring(Cyclobutane) - ring(Cyclobutane) + radical(C2CsJOCs) + radical(C2CsJOH)"""),
)

species(
    label = 'OC1=C2[CH]OO[C]3CC1C32(7428)',
    structure = SMILES('OC1=C2[CH]OO[C]3CC1C32'),
    E0 = (169.859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33701,0.0579066,-1.72051e-05,-7.87222e-09,3.30743e-12,20525.3,15.9567], Tmin=(100,'K'), Tmax=(1573.2,'K')), NASAPolynomial(coeffs=[20.5616,0.0352731,-2.06503e-05,4.19278e-09,-2.95108e-13,11228.5,-95.8246], Tmin=(1573.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + polycyclic(s2_4_6_ane) + Estimated bicyclic component: polycyclic(s2_4_6_ane) - ring(12dioxane) - ring(Cyclobutane) + ring(12dioxane) + ring(Cyclobutene) + polycyclic(s2_4_4_ene_1) - ring(Cyclobutene) - ring(12dioxane) - ring(Cyclobutane) + radical(C2CsJOOC) + radical(C=CCJO)"""),
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
    label = 'O=C1CC2C(O)=CC12(7413)',
    structure = SMILES('O=C1CC2C(O)=CC12'),
    E0 = (-219.851,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.887068,0.0597536,-4.22545e-05,1.26822e-08,-8.22011e-13,-26322.6,20.1764], Tmin=(100,'K'), Tmax=(1172.93,'K')), NASAPolynomial(coeffs=[14.5053,0.0237956,-9.67669e-06,1.78611e-09,-1.23813e-13,-30238.4,-50.7701], Tmin=(1172.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-219.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'S(571)(570)',
    structure = SMILES('O=CC1C(=O)C2CC(=O)C12'),
    E0 = (-209.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(5053.51,'J/mol'), sigma=(7.09538,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=789.35 K, Pc=32.1 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08501,0.0601037,-3.32837e-05,7.858e-09,-6.59487e-13,-25050.5,26.695], Tmin=(100,'K'), Tmax=(1928.53,'K')), NASAPolynomial(coeffs=[20.7416,0.0257765,-1.15953e-05,2.09294e-09,-1.36712e-13,-33830.3,-84.1526], Tmin=(1928.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-209.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-OdCsCs) + group(Cds-OdCsH) + polycyclic(s2_4_4_ane)"""),
)

species(
    label = 'O=CC1C=C2CC(O2)C=1O(7420)',
    structure = SMILES('O=CC1C=C2CC(O2)C=1O'),
    E0 = (-13.2047,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.217441,0.0611799,3.44456e-05,-1.10013e-07,5.14106e-11,-1407.54,24.6425], Tmin=(100,'K'), Tmax=(954.393,'K')), NASAPolynomial(coeffs=[32.0955,0.00432375,3.14919e-07,8.86639e-11,-2.60223e-14,-11153.8,-148.488], Tmin=(954.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.2047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s3_4_6_diene_0_2)"""),
)

species(
    label = 'O=CC1OC2=CC(C2)C=1O(7421)',
    structure = SMILES('O=CC1OC2=CC(C2)C=1O'),
    E0 = (-13.0898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.118902,0.065346,-7.83756e-06,-4.46013e-08,2.30064e-11,-1416.87,25.1994], Tmin=(100,'K'), Tmax=(998.893,'K')), NASAPolynomial(coeffs=[23.138,0.019268,-7.8711e-06,1.62364e-09,-1.26058e-13,-8315.5,-97.3412], Tmin=(998.893,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.0898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)H) + polycyclic(s3_4_6_diene_1_4)"""),
)

species(
    label = 'O=CC1=C(O)C2CC(O)=C12(7422)',
    structure = SMILES('O=CC1=C(O)C2CC(O)=C12'),
    E0 = (-174.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14368,0.095489,-7.36799e-05,5.21592e-09,1.04916e-11,-20844.7,23.6861], Tmin=(100,'K'), Tmax=(958.193,'K')), NASAPolynomial(coeffs=[29.7725,0.00794993,-1.64137e-06,3.18309e-10,-2.97533e-14,-28675.5,-134.093], Tmin=(958.193,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(Cyclobutane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(Cyclobutene)"""),
)

species(
    label = 'C=C1OC2C(O)=C(C=O)C12(7423)',
    structure = SMILES('C=C1OC2C(O)=C(C=O)C12'),
    E0 = (-299.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,3150,900,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.479928,0.0812604,-4.8217e-05,-6.48695e-09,1.03485e-11,-35886.1,25.0365], Tmin=(100,'K'), Tmax=(1019.74,'K')), NASAPolynomial(coeffs=[24.9625,0.0175098,-7.4686e-06,1.53989e-09,-1.18203e-13,-42949.4,-107.396], Tmin=(1019.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-299.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'O=CC1=C(O)C2C=C(O)C12(7424)',
    structure = SMILES('O=CC1=C(O)C2C=C(O)C12'),
    E0 = (-170.343,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.685892,0.0966156,-0.000102289,5.26991e-08,-1.05819e-11,-20313.3,23.6424], Tmin=(100,'K'), Tmax=(1216.98,'K')), NASAPolynomial(coeffs=[22.2414,0.021258,-9.40625e-06,1.81802e-09,-1.29625e-13,-25893.8,-91.4709], Tmin=(1216.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-170.343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(Cyclobutane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(Cyclobutene)"""),
)

species(
    label = 'O=C1CC2C(O)=C=COC12(7425)',
    structure = SMILES('O=C1CC2C(O)=C=COC12'),
    E0 = (-186.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.87695,0.124446,-0.000135037,6.81573e-08,-1.33786e-11,-22163.1,15.0649], Tmin=(100,'K'), Tmax=(1239.32,'K')), NASAPolynomial(coeffs=[29.1744,0.0242257,-1.37368e-05,2.90682e-09,-2.16079e-13,-29859.6,-141.402], Tmin=(1239.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-186.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(s2_4_6_ane) - ring(Oxane) - ring(Cyclobutane) + ring(six-inringtwodouble-12) + ring(Cyclobutane)"""),
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
    E0 = (144.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-1.05965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (14.2283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-19.1639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (26.2263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (65.1133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (24.7393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (116.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (39.4128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (86.8934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-3.44145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (43.8867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (88.4697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (76.4875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-18.7287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (122.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (37.6578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-3.23575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-15.1107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-71.1073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-166.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-25.3945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (10.2608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-23.5584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-7.1739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-40.2503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-133.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-33.9274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-131.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-15.5382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (15.0173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-59.7264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-93.3927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (25.9313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (28.5418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-15.3296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (65.3525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-54.1389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (-197.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (156.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (-225.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (251.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (248.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (169.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (-33.7633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (-207.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (300.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (215.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (-31.8718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (-40.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (-27.2292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (127.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction455',
    reactants = ['OH(5)(5)', 'O=CC1=[C]C2CC(=O)C12(7389)'],
    products = ['S(573)(572)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/NonDe] for rate rule [O_pri_rad;Cd_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction456',
    reactants = ['H(3)(3)', 'O=C[C]1C(=O)C2CC(=O)C12(7336)'],
    products = ['S(573)(572)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction457',
    reactants = ['H(3)(3)', 'O=CC1=C(O)[C]2CC(=O)C21(7390)'],
    products = ['S(573)(572)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction458',
    reactants = ['H(3)(3)', 'O=CC1[C]2C(=O)CC2C=1O(7391)'],
    products = ['S(573)(572)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.62e+13,'cm^3/(mol*s)'), n=0.228, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 13 used for C_rad/TwoDeCs;H_rad
Exact match found for rate rule [C_rad/TwoDeCs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction459',
    reactants = ['H(3)(3)', 'O=CC1=C(O)C2[CH]C(=O)C12(7392)'],
    products = ['S(573)(572)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H/OneDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction460',
    reactants = ['HCO(14)(15)', 'O=C1CC2C(O)=[C]C12(7393)'],
    products = ['S(573)(572)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34536e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_rad;CO_pri_rad] + [Cd_rad/NonDe;Y_rad] for rate rule [Cd_rad/NonDe;CO_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction461',
    reactants = ['H(3)(3)', 'O=[C]C1=C(O)C2CC(=O)C12(7394)'],
    products = ['S(573)(572)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_sec_rad;Y_rad] for rate rule [CO_rad/OneDe;H_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction462',
    reactants = ['[O]C1[C](C=O)C2C(=O)CC12(7341)'],
    products = ['S(573)(572)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction463',
    reactants = ['O=C[C]1C(O)[C]2CC(=O)C12(7395)'],
    products = ['S(573)(572)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.5515e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction464',
    reactants = ['O=CC1[C](O)C2CC(=O)[C]12(7352)'],
    products = ['S(573)(572)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.5515e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction465',
    reactants = ['[O]C1CC2[C]1C(C=O)=C2O(7396)'],
    products = ['S(573)(572)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction466',
    reactants = ['[O]C1[CH]C2C(O)=C(C=O)C12(7397)'],
    products = ['S(573)(572)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction467',
    reactants = ['O=[C]C1[C](O)C2CC(=O)C12(7355)'],
    products = ['S(573)(572)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction468',
    reactants = ['O=C[C]1[C]2C(=O)CC2C1O(7398)'],
    products = ['S(573)(572)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction469',
    reactants = ['O=C[C]1C(O)C2[CH]C(=O)C12(7399)'],
    products = ['S(573)(572)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.47101e+09,'s^-1'), n=0.37, Ea=(99.6901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad_NDe] + [R3radExo;Y_rad;XH_Rrad_NDe] for rate rule [R3radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction470',
    reactants = ['O=CC1[C](O)[C]2CC(=O)C21(7359)'],
    products = ['S(573)(572)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction471',
    reactants = ['[O]C1C[C]2C(O)=C(C=O)C21(7400)'],
    products = ['S(573)(572)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction472',
    reactants = ['[O]CC1[C]2C(=O)CC2C=1O(7401)'],
    products = ['S(573)(572)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction473',
    reactants = ['O=CC1[C](O)C2[CH]C(=O)C21(7353)'],
    products = ['S(573)(572)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(60.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction474',
    reactants = ['O=CC1=C(O)[C]2C[C](O)C21(7402)'],
    products = ['S(573)(572)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction475',
    reactants = ['O=C1CC2[C]1C([CH]O)=C2O(7403)'],
    products = ['S(573)(572)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction476',
    reactants = ['[O]CC1=C(O)[C]2CC(=O)C21(7404)'],
    products = ['S(573)(572)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction477',
    reactants = ['[O]C1CC2C(O)=C([C]=O)C12(7405)'],
    products = ['S(573)(572)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction478',
    reactants = ['[O]C[C]1C(=O)C2CC(=O)C12(7342)'],
    products = ['S(573)(572)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction479',
    reactants = ['O=CC1=C(O)C2[CH][C](O)C12(7406)'],
    products = ['S(573)(572)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction480',
    reactants = ['O=CC1[C]2[C](O)CC2C=1O(7407)'],
    products = ['S(573)(572)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction481',
    reactants = ['O=C1C[C]2C(O)=C([CH]O)C21(7408)'],
    products = ['S(573)(572)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction482',
    reactants = ['O=[C]C1=C(O)C2C[C](O)C12(7409)'],
    products = ['S(573)(572)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction483',
    reactants = ['[O]C1=C([CH]O)C2C(=O)CC12(7410)'],
    products = ['S(573)(572)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction484',
    reactants = ['[O]C1CC2C(=O)[C](C=O)C12(7345)'],
    products = ['S(573)(572)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction485',
    reactants = ['[O]CC1=C(O)C2[CH]C(=O)C12(7411)'],
    products = ['S(573)(572)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction486',
    reactants = ['O=C[C]1C(=O)C2C[C](O)C12(7351)'],
    products = ['S(573)(572)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_NDe] for rate rule [R6radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction487',
    reactants = ['O=C1[CH]C2C(O)=C([CH]O)C12(7412)'],
    products = ['S(573)(572)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction489',
    reactants = ['[CH2]C1C(O)=C(C=O)C1[C]=O(7414)'],
    products = ['S(573)(572)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction490',
    reactants = ['O=[C]CC1[CH]C(C=O)=C1O(7415)'],
    products = ['S(573)(572)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction491',
    reactants = ['O=C[C]=C(O)C1[CH]C(=O)C1(7416)'],
    products = ['S(573)(572)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriDe_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction492',
    reactants = ['O=CC(=[C]O)C1[CH]CC1=O(7417)'],
    products = ['S(573)(572)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriND_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction493',
    reactants = ['[CH2]C(=O)C1[CH]C(O)=C1C=O(7418)'],
    products = ['S(573)(572)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction494',
    reactants = ['[O]C1=C[C](C=O)C(O)=CC1(7220)'],
    products = ['S(573)(572)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] + [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction495',
    reactants = ['CH2CO(26)(27)', 'O=CC1=CC=C1O(7419)'],
    products = ['S(573)(572)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.77e-07,'m^3/(mol*s)'), n=2.94, Ea=(121.336,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [diene_out;diene_in;allene_unsub] for rate rule [diene_out;diene_in_DeNd;allene_unsub]
Euclidian distance = 1.0
family: Diels_alder_addition"""),
)

reaction(
    label = 'reaction502',
    reactants = ['S(599)(598)'],
    products = ['S(573)(572)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction503',
    reactants = ['O=C[C]1C2[C]3CC2C1(O)O3(7426)'],
    products = ['S(573)(572)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction504',
    reactants = ['O=CC12O[C]3CC([C]1O)C32(7427)'],
    products = ['S(573)(572)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction505',
    reactants = ['OC1=C2[CH]OO[C]3CC1C32(7428)'],
    products = ['S(573)(572)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction488',
    reactants = ['CO(10)(11)', 'O=C1CC2C(O)=CC12(7413)'],
    products = ['S(573)(572)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(0.0591985,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd/H/NonDeC]
Euclidian distance = 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction441',
    reactants = ['S(573)(572)'],
    products = ['S(571)(570)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction496',
    reactants = ['O=CC1C=C2CC(O2)C=1O(7420)'],
    products = ['S(573)(572)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_Cs;R_O_C]
Euclidian distance = 1.41421356237
family: ketoenol"""),
)

reaction(
    label = 'reaction497',
    reactants = ['O=CC1OC2=CC(C2)C=1O(7421)'],
    products = ['S(573)(572)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(5460.43,'s^-1'), n=2.73, Ea=(228.457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.2360679775
family: ketoenol"""),
)

reaction(
    label = 'reaction498',
    reactants = ['O=CC1=C(O)C2CC(O)=C12(7422)'],
    products = ['S(573)(572)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction499',
    reactants = ['C=C1OC2C(O)=C(C=O)C12(7423)'],
    products = ['S(573)(572)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction500',
    reactants = ['O=CC1=C(O)C2C=C(O)C12(7424)'],
    products = ['S(573)(572)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CsC;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction501',
    reactants = ['O=C1CC2C(O)=C=COC12(7425)'],
    products = ['S(573)(572)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

network(
    label = '428',
    isomers = [
        'S(573)(572)',
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
    label = '428',
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

