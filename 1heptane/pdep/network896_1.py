species(
    label = '[O]C(C=O)C([O])C=O(1310)',
    structure = SMILES('[O]C(C=O)C([O])C=O'),
    E0 = (-162.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1019.1,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0281104,'amu*angstrom^2'), symmetry=1, barrier=(20.7159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900935,'amu*angstrom^2'), symmetry=1, barrier=(20.7143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.418182,'amu*angstrom^2'), symmetry=1, barrier=(20.7158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21395,0.0534636,-4.18756e-05,1.44531e-08,-1.45425e-12,-19401.3,34.1683], Tmin=(100,'K'), Tmax=(1160.74,'K')), NASAPolynomial(coeffs=[14.2824,0.0175632,-7.28639e-06,1.36663e-09,-9.58967e-14,-23050.4,-33.4779], Tmin=(1160.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=OCOJ)"""),
)

species(
    label = 'OCHCHO(48)',
    structure = SMILES('O=CC=O'),
    E0 = (-225.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,180,857.339,857.946,861.419,1717.87],'cm^-1')),
        HinderedRotor(inertia=(0.00964493,'amu*angstrom^2'), symmetry=1, barrier=(5.0669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.68412,0.000478013,4.26391e-05,-5.79018e-08,2.31669e-11,-27198.5,4.51187], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[8.72507,0.00633097,-2.35575e-06,3.89783e-10,-2.37487e-14,-29102.4,-20.3904], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-225.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""OCHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[O]C(C=O)C1OC1[O](5634)',
    structure = SMILES('[O]C(C=O)C1OC1[O]'),
    E0 = (-110.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06322,0.0574975,-5.41941e-05,2.80741e-08,-5.75063e-12,-13200.5,30.4939], Tmin=(100,'K'), Tmax=(1282.48,'K')), NASAPolynomial(coeffs=[12.3406,0.019081,-5.46875e-06,7.73648e-10,-4.44759e-14,-15826.4,-25.6787], Tmin=(1282.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C1OC(C=O)C1[O](4883)',
    structure = SMILES('[O]C1OC(C=O)C1[O]'),
    E0 = (-129.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71058,0.0395893,5.95998e-06,-4.02855e-08,2.01204e-11,-15493.5,26.8773], Tmin=(100,'K'), Tmax=(919.065,'K')), NASAPolynomial(coeffs=[13.0128,0.0183569,-5.01642e-06,7.75164e-10,-5.21025e-14,-18751.8,-33.1193], Tmin=(919.065,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + ring(Oxetane) + radical(CC(C)OJ) + radical(CCOJ)"""),
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
    label = '[O]C(C=O)C(=O)C=O(5635)',
    structure = SMILES('[O]C(C=O)C(=O)C=O'),
    E0 = (-319.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1390,370,380,2900,435,375,552.5,462.5,1710,414.727,2227.29],'cm^-1')),
        HinderedRotor(inertia=(0.0573439,'amu*angstrom^2'), symmetry=1, barrier=(6.99913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.056653,'amu*angstrom^2'), symmetry=1, barrier=(7.00509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319644,'amu*angstrom^2'), symmetry=1, barrier=(39.2824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56042,0.0587646,-7.18915e-05,5.2131e-08,-1.60772e-11,-38316.2,27.524], Tmin=(100,'K'), Tmax=(775.467,'K')), NASAPolynomial(coeffs=[7.29972,0.0291609,-1.46295e-05,2.90404e-09,-2.07413e-13,-39206.3,1.29457], Tmin=(775.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH) + group(Cds-O2d(Cds-O2d)H) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C=C[O](2537)',
    structure = SMILES('[O]C=C[O]'),
    E0 = (-18.8461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,714.947,715.113],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7018,-0.0027881,5.18222e-05,-5.96735e-08,2.03775e-11,-2247.83,13.3237], Tmin=(100,'K'), Tmax=(1035.3,'K')), NASAPolynomial(coeffs=[6.65087,0.0113883,-5.76543e-06,1.26602e-09,-9.87778e-14,-4228.83,-7.62441], Tmin=(1035.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.8461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'HCO(16)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[O]C(C=O)C=O(3641)',
    structure = SMILES('[O]C(C=O)C=O'),
    E0 = (-211.279,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1930.68],'cm^-1')),
        HinderedRotor(inertia=(0.0134357,'amu*angstrom^2'), symmetry=1, barrier=(6.41037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.474193,'amu*angstrom^2'), symmetry=1, barrier=(10.9026,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63661,0.0319192,-2.26797e-05,8.51925e-09,-1.38613e-12,-25363.8,20.5163], Tmin=(100,'K'), Tmax=(1344.19,'K')), NASAPolynomial(coeffs=[6.77299,0.0196105,-8.94442e-06,1.70721e-09,-1.19208e-13,-26475.8,-0.662934], Tmin=(1344.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-211.279,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C(C=O)[C](O)C=O(5636)',
    structure = SMILES('[O]C=C(O)C([O])C=O'),
    E0 = (-308.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.451657,0.0635926,-3.67959e-05,-1.66077e-08,1.6229e-11,-36962,31.0991], Tmin=(100,'K'), Tmax=(925.667,'K')), NASAPolynomial(coeffs=[23.2326,0.00262012,1.29108e-06,-3.10181e-10,1.76627e-14,-42784.8,-85.7168], Tmin=(925.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=COJ)"""),
)

species(
    label = '[O][C](C=O)C(O)C=O(5637)',
    structure = SMILES('[O]C=C([O])C(O)C=O'),
    E0 = (-414.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.592428,0.06237,-3.99609e-05,-8.25752e-09,1.19559e-11,-49709.1,31.314], Tmin=(100,'K'), Tmax=(934.48,'K')), NASAPolynomial(coeffs=[21.1811,0.0060351,-5.68719e-07,4.85755e-11,-6.59706e-15,-54945.3,-74.0471], Tmin=(934.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-414.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C(C=O)C(O)[C]=O(5638)',
    structure = SMILES('[O]C(C=O)C(O)[C]=O'),
    E0 = (-245.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.982542,0.058462,-5.04008e-05,1.76754e-08,-1.16368e-12,-29468.1,35.7254], Tmin=(100,'K'), Tmax=(1041.52,'K')), NASAPolynomial(coeffs=[15.9534,0.014107,-5.44622e-06,1.01449e-09,-7.23183e-14,-33299.4,-40.5308], Tmin=(1041.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-245.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C([C]=O)C(O)C=O(5639)',
    structure = SMILES('[O]C([C]=O)C(O)C=O'),
    E0 = (-245.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.982542,0.058462,-5.04008e-05,1.76754e-08,-1.16368e-12,-29468.1,35.7254], Tmin=(100,'K'), Tmax=(1041.52,'K')), NASAPolynomial(coeffs=[15.9534,0.014107,-5.44622e-06,1.01449e-09,-7.23183e-14,-33299.4,-40.5308], Tmin=(1041.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-245.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(C=O)C1[CH]OO1(5640)',
    structure = SMILES('[O]C(C=O)C1[CH]OO1'),
    E0 = (102.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63429,0.0546283,-4.70648e-05,2.18311e-08,-4.29285e-12,12421,25.5673], Tmin=(100,'K'), Tmax=(1171.52,'K')), NASAPolynomial(coeffs=[9.37817,0.0281876,-1.32103e-05,2.56575e-09,-1.81627e-13,10606.6,-13.0182], Tmin=(1171.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(12dioxetane) + radical(CCsJOO) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C1[CH]OOC1C=O(4837)',
    structure = SMILES('[O]C1[CH]OOC1C=O'),
    E0 = (3.2769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47472,0.0386033,2.76591e-05,-7.70304e-08,3.7127e-11,500.986,25.4803], Tmin=(100,'K'), Tmax=(897.943,'K')), NASAPolynomial(coeffs=[18.3154,0.00920124,5.72238e-07,-3.44542e-10,2.51187e-14,-4362.47,-64.1938], Tmin=(897.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.2769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(12dioxolane) + radical(CC(C)OJ) + radical(CCsJOOC)"""),
)

species(
    label = 'O=CC(=O)C(O)C=O(5641)',
    structure = SMILES('O=CC(=O)C(O)C=O'),
    E0 = (-563.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49923,0.0595752,-6.06742e-05,3.41797e-08,-8.15436e-12,-67627.5,26.963], Tmin=(100,'K'), Tmax=(984.985,'K')), NASAPolynomial(coeffs=[8.98183,0.0291887,-1.43999e-05,2.86006e-09,-2.05127e-13,-69101.6,-9.02304], Tmin=(984.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-563.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH) + group(Cds-O2d(Cds-O2d)H)"""),
)

species(
    label = 'CO(12)',
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
    label = '[O]CC([O])C=O(970)',
    structure = SMILES('[O]CC([O])C=O'),
    E0 = (-63.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,440.965,1697,2474.95],'cm^-1')),
        HinderedRotor(inertia=(0.0566694,'amu*angstrom^2'), symmetry=1, barrier=(6.95388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0548441,'amu*angstrom^2'), symmetry=1, barrier=(6.74531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4282.73,'J/mol'), sigma=(6.77277,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=668.95 K, Pc=31.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80966,0.0319521,-2.9017e-07,-6.00266e-08,5.97141e-11,-7611.74,22.9413], Tmin=(100,'K'), Tmax=(441.444,'K')), NASAPolynomial(coeffs=[4.24755,0.0273331,-1.31716e-05,2.58322e-09,-1.83678e-13,-7820.63,16.2519], Tmin=(441.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCOJ)"""),
)

species(
    label = 'O=CC1OOC1C=O(1313)',
    structure = SMILES('O=CC1OOC1C=O'),
    E0 = (-248.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68928,0.0376839,9.73746e-06,-3.91228e-08,1.68674e-11,-29850.6,24.7644], Tmin=(100,'K'), Tmax=(1026.29,'K')), NASAPolynomial(coeffs=[14.0811,0.0198348,-8.67732e-06,1.74773e-09,-1.3048e-13,-33997.6,-43.1526], Tmin=(1026.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + ring(12dioxetane)"""),
)

species(
    label = '[O]C(C=O)O[CH]C=O(1311)',
    structure = SMILES('[O]C=COC([O])C=O'),
    E0 = (-290.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,266.809,267.46,267.479,267.595,267.737],'cm^-1')),
        HinderedRotor(inertia=(0.43173,'amu*angstrom^2'), symmetry=1, barrier=(21.9157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.429668,'amu*angstrom^2'), symmetry=1, barrier=(21.9216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.432166,'amu*angstrom^2'), symmetry=1, barrier=(21.9199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4632.59,'J/mol'), sigma=(7.12501,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=723.60 K, Pc=29.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.386483,0.0624346,-2.59843e-05,-3.11712e-08,2.19752e-11,-34774.4,29.7735], Tmin=(100,'K'), Tmax=(928.8,'K')), NASAPolynomial(coeffs=[24.7494,0.00133247,1.92598e-06,-4.08465e-10,2.24894e-14,-41190.2,-96.1391], Tmin=(928.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C(=CO)C([O])C=O(5642)',
    structure = SMILES('[O]C(=CO)C([O])C=O'),
    E0 = (-312.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,346.712,347.451,348.133],'cm^-1')),
        HinderedRotor(inertia=(0.158227,'amu*angstrom^2'), symmetry=1, barrier=(13.4582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157154,'amu*angstrom^2'), symmetry=1, barrier=(13.4591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159362,'amu*angstrom^2'), symmetry=1, barrier=(13.4623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.551851,0.0777728,-8.64014e-05,4.40374e-08,-8.19617e-12,-37360,34.7994], Tmin=(100,'K'), Tmax=(1553.23,'K')), NASAPolynomial(coeffs=[22.9372,0.00151404,2.47149e-06,-6.43793e-10,4.7468e-14,-42754.8,-82.7422], Tmin=(1553.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=C(C)OJ)"""),
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
    label = '[O]C([CH]C=O)C=O(4744)',
    structure = SMILES('[O]C=CC([O])C=O'),
    E0 = (-93.9344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,461.745,461.829,461.843],'cm^-1')),
        HinderedRotor(inertia=(0.0970111,'amu*angstrom^2'), symmetry=1, barrier=(14.6777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0970098,'amu*angstrom^2'), symmetry=1, barrier=(14.6782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44256,0.0449914,-1.47293e-05,-2.06564e-08,1.30766e-11,-11195.4,27.4845], Tmin=(100,'K'), Tmax=(961.149,'K')), NASAPolynomial(coeffs=[16.1571,0.0107158,-3.31513e-06,6.12131e-10,-4.67717e-14,-15269.4,-49.4004], Tmin=(961.149,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-93.9344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(C=COJ)"""),
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
    E0 = (-162.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-59.9093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-40.1568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-75.0528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-162.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-84.8701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-8.29166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-86.9191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-86.9191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-34.8062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-37.6922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (103.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (3.2769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-98.7996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (122.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-153.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (23.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-162.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (149.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(C=O)C([O])C=O(1310)'],
    products = ['OCHCHO(48)', 'OCHCHO(48)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C(C=O)C([O])C=O(1310)'],
    products = ['[O]C(C=O)C1OC1[O](5634)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.557e+12,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(C=O)C([O])C=O(1310)'],
    products = ['[O]C1OC(C=O)C1[O](4883)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.448e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[O]C(C=O)C(=O)C=O(5635)'],
    products = ['[O]C(C=O)C([O])C=O(1310)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(25.1243,'m^3/(mol*s)'), n=1.86, Ea=(32.426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-DeNd_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OCHCHO(48)', '[O]C=C[O](2537)'],
    products = ['[O]C(C=O)C([O])C=O(1310)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(421.062,'m^3/(mol*s)'), n=1.35673, Ea=(81.9464,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-DeH_O;YJ] for rate rule [CO-DeH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond
Ea raised from 78.3 to 81.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['HCO(16)', '[O]C(C=O)C=O(3641)'],
    products = ['[O]C(C=O)C([O])C=O(1310)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.04e+12,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;CO_pri_rad] for rate rule [CO-CsH_O;CO_pri_rad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(C=O)C([O])C=O(1310)'],
    products = ['[O]C(C=O)[C](O)C=O(5636)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.40446e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(C=O)C([O])C=O(1310)'],
    products = ['[O][C](C=O)C(O)C=O(5637)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(C=O)C([O])C=O(1310)'],
    products = ['[O]C(C=O)C(O)[C]=O(5638)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(C=O)C([O])C=O(1310)'],
    products = ['[O]C([C]=O)C(O)C=O(5639)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.50344e+06,'s^-1'), n=1.80068, Ea=(127.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;CO_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C=C[O](2537)', '[O]C=C[O](2537)'],
    products = ['[O]C(C=O)C([O])C=O(1310)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(C=O)C([O])C=O(1310)'],
    products = ['[O]C(C=O)C1[CH]OO1(5640)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.006e+11,'s^-1'), n=0.221, Ea=(265.786,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(C=O)C([O])C=O(1310)'],
    products = ['[O]C1[CH]OOC1C=O(4837)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.62368e+09,'s^-1'), n=0.551229, Ea=(165.477,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 164.4 to 165.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(C=O)C([O])C=O(1310)'],
    products = ['O=CC(=O)C(O)C=O(5641)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CO(12)', '[O]CC([O])C=O(970)'],
    products = ['[O]C(C=O)C([O])C=O(1310)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(C=O)C([O])C=O(1310)'],
    products = ['O=CC1OOC1C=O(1313)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O]C(C=O)C([O])C=O(1310)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(=CO)C([O])C=O(5642)'],
    products = ['[O]C(C=O)C([O])C=O(1310)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(149.951,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol
Ea raised from 147.5 to 150.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['O(4)', '[O]C([CH]C=O)C=O(4744)'],
    products = ['[O]C(C=O)C([O])C=O(1310)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '896',
    isomers = [
        '[O]C(C=O)C([O])C=O(1310)',
    ],
    reactants = [
        ('OCHCHO(48)', 'OCHCHO(48)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '896',
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

