species(
    label = '[O]C(O)C([O])C=O(1315)',
    structure = SMILES('[O]C(O)C([O])C=O'),
    E0 = (-263.825,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,342.93,343.178,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00545562,'amu*angstrom^2'), symmetry=1, barrier=(8.94958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00545545,'amu*angstrom^2'), symmetry=1, barrier=(8.95022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265483,'amu*angstrom^2'), symmetry=1, barrier=(22.2133,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8123,0.0527175,-6.2343e-05,4.46689e-08,-1.37608e-11,-31656.2,29.5754], Tmin=(100,'K'), Tmax=(773.215,'K')), NASAPolynomial(coeffs=[6.66135,0.0276317,-1.36765e-05,2.70743e-09,-1.93253e-13,-32406,7.42883], Tmin=(773.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(C=OCOJ)"""),
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
    label = '[O]C(O)C1OC1[O](8455)',
    structure = SMILES('[O]C(O)C1OC1[O]'),
    E0 = (-212.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55204,0.0581035,-7.95598e-05,6.45995e-08,-2.06046e-11,-25450.8,25.5957], Tmin=(100,'K'), Tmax=(923.394,'K')), NASAPolynomial(coeffs=[6.34573,0.0263938,-1.02709e-05,1.73927e-09,-1.10693e-13,-25869.6,5.37743], Tmin=(923.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-212.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C1OC(O)C1[O](7934)',
    structure = SMILES('[O]C1OC(O)C1[O]'),
    E0 = (-213.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34415,0.0487002,-4.14789e-05,1.93693e-08,-3.55292e-12,-25536.4,24.9663], Tmin=(100,'K'), Tmax=(1488.49,'K')), NASAPolynomial(coeffs=[11.4754,0.016666,-4.3513e-06,5.70275e-10,-3.09998e-14,-28019.8,-26.1518], Tmin=(1488.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(Oxetane) + radical(CC(C)OJ) + radical(CCOJ)"""),
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
    label = '[O]C(O)C(=O)C=O(8456)',
    structure = SMILES('[O]C(O)C(=O)C=O'),
    E0 = (-402.07,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,375,552.5,462.5,1710,3615,1277.5,1000,1380,1390,370,380,2900,435,193.18,627.668],'cm^-1')),
        HinderedRotor(inertia=(0.00451831,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140663,'amu*angstrom^2'), symmetry=1, barrier=(3.72232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140512,'amu*angstrom^2'), symmetry=1, barrier=(3.72157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57306,0.0580324,-8.3496e-05,6.8954e-08,-2.34371e-11,-48274.8,25.8157], Tmin=(100,'K'), Tmax=(729.787,'K')), NASAPolynomial(coeffs=[7.64699,0.0240197,-1.21044e-05,2.38325e-09,-1.68395e-13,-49142.2,-1.44271], Tmin=(729.787,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-402.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C(C=O)C(=O)O(8457)',
    structure = SMILES('[O]C(C=O)C(=O)O'),
    E0 = (-476.966,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,388.732,388.916,388.951,388.994,389.369,1103.18],'cm^-1')),
        HinderedRotor(inertia=(0.0493264,'amu*angstrom^2'), symmetry=1, barrier=(42.81,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00103319,'amu*angstrom^2'), symmetry=1, barrier=(5.78282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0215462,'amu*angstrom^2'), symmetry=1, barrier=(42.8064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22608,0.0411567,-3.22677e-05,1.33031e-08,-2.33885e-12,-57303.8,23.9929], Tmin=(100,'K'), Tmax=(1280.46,'K')), NASAPolynomial(coeffs=[8.12031,0.0227436,-1.06973e-05,2.07238e-09,-1.46113e-13,-58813.3,-5.90047], Tmin=(1280.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-476.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(C=OCOJ)"""),
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
    label = '[O]C(O)C=O(4774)',
    structure = SMILES('[O]C(O)C=O'),
    E0 = (-294.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,3550.67],'cm^-1')),
        HinderedRotor(inertia=(0.0510999,'amu*angstrom^2'), symmetry=1, barrier=(2.09244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.625001,'amu*angstrom^2'), symmetry=1, barrier=(14.37,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.73016,0.0299168,-2.84183e-05,1.57731e-08,-3.78628e-12,-35325.3,19.2328], Tmin=(100,'K'), Tmax=(971.508,'K')), NASAPolynomial(coeffs=[6.01114,0.0164081,-7.56104e-06,1.46047e-09,-1.03199e-13,-35962.9,3.49881], Tmin=(971.508,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-294.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(C=OCOJ)"""),
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
    label = '[O]C(O)[C](O)C=O(8458)',
    structure = SMILES('[O]C=C(O)C([O])O'),
    E0 = (-410.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.620987,0.0680187,-7.57653e-05,3.79118e-08,-6.46141e-12,-49198.5,27.3435], Tmin=(100,'K'), Tmax=(981.131,'K')), NASAPolynomial(coeffs=[18.1293,0.00841079,-2.63148e-06,4.47423e-10,-3.12906e-14,-53200.7,-59.6778], Tmin=(981.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-410.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C(C=O)[C](O)O(8459)',
    structure = SMILES('[O]C(C=O)[C](O)O'),
    E0 = (-284.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44172,0.0568813,-6.35393e-05,3.65746e-08,-8.40099e-12,-34099.8,30.0299], Tmin=(100,'K'), Tmax=(1056.03,'K')), NASAPolynomial(coeffs=[11.7518,0.0178279,-8.06544e-06,1.55307e-09,-1.09885e-13,-36277.2,-20.272], Tmin=(1056.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-284.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O][C](O)C(O)C=O(8460)',
    structure = SMILES('[O][C](O)C(O)C=O'),
    E0 = (-302.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69762,0.0552921,-6.50212e-05,4.32965e-08,-1.21033e-11,-36280.8,29.8968], Tmin=(100,'K'), Tmax=(854.558,'K')), NASAPolynomial(coeffs=[8.0547,0.025536,-1.27906e-05,2.54983e-09,-1.82925e-13,-37367.3,0.226727], Tmin=(854.558,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-302.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[O]C(O)C(O)[C]=O(8461)',
    structure = SMILES('[O]C(O)C(O)[C]=O'),
    E0 = (-347.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52728,0.0582852,-7.25015e-05,4.93371e-08,-1.37178e-11,-41720.7,30.6353], Tmin=(100,'K'), Tmax=(869.714,'K')), NASAPolynomial(coeffs=[9.28835,0.02259,-1.09372e-05,2.14539e-09,-1.52412e-13,-43070.6,-5.72388], Tmin=(869.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-347.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C](C=O)C(O)O(8462)',
    structure = SMILES('[O]C=C([O])C(O)O'),
    E0 = (-498.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.330931,0.0702969,-8.34651e-05,4.65079e-08,-9.69062e-12,-59756.8,28.3293], Tmin=(100,'K'), Tmax=(1302.07,'K')), NASAPolynomial(coeffs=[19.6421,0.00438138,6.33995e-08,-1.46831e-10,1.36271e-14,-64226.9,-67.7875], Tmin=(1302.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-498.019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C([C]=O)C(O)O(8463)',
    structure = SMILES('[O]C([C]=O)C(O)O'),
    E0 = (-329.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25088,0.0600874,-7.16436e-05,4.32914e-08,-1.02641e-11,-39538.7,30.8439], Tmin=(100,'K'), Tmax=(1034.07,'K')), NASAPolynomial(coeffs=[12.8924,0.0150544,-6.31798e-06,1.17479e-09,-8.16174e-14,-41946.2,-25.7094], Tmin=(1034.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-329.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C([O])C(O)C=O(8464)',
    structure = SMILES('[O]C([O])C(O)C=O'),
    E0 = (-281.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78361,0.0549511,-7.95749e-05,7.56604e-08,-2.97811e-11,-33825.3,29.0477], Tmin=(100,'K'), Tmax=(761.865,'K')), NASAPolynomial(coeffs=[3.85677,0.0337268,-1.74301e-05,3.46735e-09,-2.46178e-13,-33841.1,21.5792], Tmin=(761.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-281.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C(O)C1[CH]OO1(8465)',
    structure = SMILES('[O]C(O)C1[CH]OO1'),
    E0 = (0.959555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52387,0.0624676,-9.8767e-05,9.46117e-08,-3.58514e-11,196.783,22.8101], Tmin=(100,'K'), Tmax=(812.337,'K')), NASAPolynomial(coeffs=[3.88415,0.0346607,-1.75352e-05,3.42011e-09,-2.38719e-13,347.325,15.2006], Tmin=(812.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.959555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(12dioxetane) + radical(CCsJOO) + radical(CCOJ)"""),
)

species(
    label = '[O]C1[CH]OOC1O(7968)',
    structure = SMILES('[O]C1[CH]OOC1O'),
    E0 = (-80.3204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56682,0.0422964,-7.09147e-07,-4.24442e-08,2.43741e-11,-9561.96,21.9243], Tmin=(100,'K'), Tmax=(873.458,'K')), NASAPolynomial(coeffs=[16.0317,0.00883726,4.52368e-07,-3.61073e-10,3.0437e-14,-13339.4,-53.062], Tmin=(873.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.3204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(12dioxolane) + radical(CCsJOOC) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=CC(=O)C(O)O(8466)',
    structure = SMILES('O=CC(=O)C(O)O'),
    E0 = (-645.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55816,0.0580969,-6.86685e-05,4.47684e-08,-1.20908e-11,-77588,24.408], Tmin=(100,'K'), Tmax=(888.06,'K')), NASAPolynomial(coeffs=[8.96902,0.0247163,-1.22852e-05,2.44069e-09,-1.74807e-13,-78904.2,-10.4651], Tmin=(888.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-645.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H)"""),
)

species(
    label = 'O=CC(O)C(=O)O(8467)',
    structure = SMILES('O=CC(O)C(=O)O'),
    E0 = (-720.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80941,0.0461047,-3.5258e-05,1.34958e-08,-2.10836e-12,-86599.7,24.7101], Tmin=(100,'K'), Tmax=(1479.11,'K')), NASAPolynomial(coeffs=[11.3258,0.0203691,-9.15885e-06,1.73236e-09,-1.20086e-13,-89414.9,-24.9263], Tmin=(1479.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-720.699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsH) + group(Cds-OdCsOs) + group(Cds-OdCsH)"""),
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
    label = '[O]CC([O])O(980)',
    structure = SMILES('[O]CC([O])O'),
    E0 = (-165.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1380,1390,370,380,2900,435,180,2144.82,2145.09],'cm^-1')),
        HinderedRotor(inertia=(0.1642,'amu*angstrom^2'), symmetry=1, barrier=(3.77529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16403,'amu*angstrom^2'), symmetry=1, barrier=(3.77137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4356.85,'J/mol'), sigma=(6.93803,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=680.53 K, Pc=29.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42536,0.0451945,-8.59012e-05,9.42384e-08,-3.78131e-11,-19824.8,21.0512], Tmin=(100,'K'), Tmax=(848.613,'K')), NASAPolynomial(coeffs=[-0.435924,0.032337,-1.66085e-05,3.22031e-09,-2.22264e-13,-18390.6,39.9745], Tmin=(848.613,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = 'O=CC1OOC1O(1319)',
    structure = SMILES('O=CC1OOC1O'),
    E0 = (-332.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73305,0.0419214,-2.04177e-05,-2.36495e-09,3.22942e-12,-39911.4,22.0766], Tmin=(100,'K'), Tmax=(1129.11,'K')), NASAPolynomial(coeffs=[12.3018,0.0186574,-8.34579e-06,1.62757e-09,-1.1675e-13,-43201.7,-34.1966], Tmin=(1129.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-332.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + ring(12dioxetane)"""),
)

species(
    label = '[O]C([CH]OO)C=O(8468)',
    structure = SMILES('[O]C([CH]OO)C=O'),
    E0 = (-29.2116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,464.704,2618.01],'cm^-1')),
        HinderedRotor(inertia=(0.0915826,'amu*angstrom^2'), symmetry=1, barrier=(13.1246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257086,'amu*angstrom^2'), symmetry=1, barrier=(38.5852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0435295,'amu*angstrom^2'), symmetry=1, barrier=(2.25105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67846,'amu*angstrom^2'), symmetry=1, barrier=(38.5922,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44118,0.0603413,-7.14048e-05,4.64068e-08,-1.24276e-11,-3424.68,27.7987], Tmin=(100,'K'), Tmax=(897.851,'K')), NASAPolynomial(coeffs=[9.40787,0.0248498,-1.21123e-05,2.38238e-09,-1.69651e-13,-4855.29,-9.77781], Tmin=(897.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.2116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCsJOOH)"""),
)

species(
    label = '[O][CH]C(C=O)OO(8469)',
    structure = SMILES('[O][CH]C(C=O)OO'),
    E0 = (-55.5242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,346.79,1788.45],'cm^-1')),
        HinderedRotor(inertia=(0.0769666,'amu*angstrom^2'), symmetry=1, barrier=(6.5684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.434893,'amu*angstrom^2'), symmetry=1, barrier=(37.1137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0769664,'amu*angstrom^2'), symmetry=1, barrier=(6.56842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.434883,'amu*angstrom^2'), symmetry=1, barrier=(37.1136,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58052,0.0617287,-6.84077e-05,2.08154e-08,1.47843e-11,-6599.3,26.8502], Tmin=(100,'K'), Tmax=(528.745,'K')), NASAPolynomial(coeffs=[7.73553,0.0290228,-1.49362e-05,2.96217e-09,-2.10081e-13,-7443.89,-0.753633], Tmin=(528.745,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.5242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O]C(O)O[CH]C=O(1317)',
    structure = SMILES('[O]C=COC([O])O'),
    E0 = (-369.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180,180,910.925,911.789],'cm^-1')),
        HinderedRotor(inertia=(0.230122,'amu*angstrom^2'), symmetry=1, barrier=(18.5603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.807226,'amu*angstrom^2'), symmetry=1, barrier=(18.5597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0314561,'amu*angstrom^2'), symmetry=1, barrier=(18.555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4708.06,'J/mol'), sigma=(7.29164,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=735.39 K, Pc=27.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41963,0.0577655,-5.76161e-05,2.89678e-08,-5.88519e-12,-44357.2,23.5691], Tmin=(100,'K'), Tmax=(1172.93,'K')), NASAPolynomial(coeffs=[12.2793,0.020731,-1.02544e-05,2.04837e-09,-1.47522e-13,-46904.7,-30.5547], Tmin=(1172.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-369.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(OCOJ)"""),
)

species(
    label = '[O]C(=CO)C([O])O(8470)',
    structure = SMILES('[O]C(=CO)C([O])O'),
    E0 = (-413.776,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,328.653,328.792,329.055],'cm^-1')),
        HinderedRotor(inertia=(0.158597,'amu*angstrom^2'), symmetry=1, barrier=(12.0004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153371,'amu*angstrom^2'), symmetry=1, barrier=(11.9887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151449,'amu*angstrom^2'), symmetry=1, barrier=(11.9814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.597275,0.0707567,-8.57782e-05,4.76335e-08,-9.27726e-12,-49639.5,27.5201], Tmin=(100,'K'), Tmax=(920.725,'K')), NASAPolynomial(coeffs=[17.581,0.00819188,-2.12837e-06,2.99768e-10,-1.84103e-14,-53242.5,-55.5962], Tmin=(920.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-413.776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(CCOJ)"""),
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
    label = '[O]C(O)[CH]C=O(7809)',
    structure = SMILES('[O]C=CC([O])O'),
    E0 = (-195.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,442.992,443.288,445.231],'cm^-1')),
        HinderedRotor(inertia=(0.0864937,'amu*angstrom^2'), symmetry=1, barrier=(12.2135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0874108,'amu*angstrom^2'), symmetry=1, barrier=(12.2373,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80761,0.0470676,-4.52296e-05,2.24065e-08,-4.46932e-12,-23440.5,23.029], Tmin=(100,'K'), Tmax=(1201.25,'K')), NASAPolynomial(coeffs=[10.9078,0.0167648,-7.39014e-06,1.40626e-09,-9.87829e-14,-25626.8,-22.5429], Tmin=(1201.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C([CH]O)C=O(3643)',
    structure = SMILES('[O]C([CH]O)C=O'),
    E0 = (-109.001,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.130976,'amu*angstrom^2'), symmetry=1, barrier=(3.01139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0583656,'amu*angstrom^2'), symmetry=1, barrier=(12.0817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0919641,'amu*angstrom^2'), symmetry=1, barrier=(19.9596,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76206,0.0505058,-5.90343e-05,3.6861e-08,-9.2422e-12,-13030.3,25.2761], Tmin=(100,'K'), Tmax=(969.377,'K')), NASAPolynomial(coeffs=[9.83735,0.0171844,-7.47374e-06,1.40174e-09,-9.74014e-14,-14595.9,-13.4314], Tmin=(969.377,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCsJOH) + radical(C=OCOJ)"""),
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
    E0 = (-263.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-161.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-141.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-157.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-218.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-180.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-167.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-263.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-88.1814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-109.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-105.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-188.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-188.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-188.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-136.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-158.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-13.4679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1.96138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-80.3204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-200.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-200.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (20.8682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-255.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (63.5335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (51.5863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-55.7718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-263.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (47.4451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (134.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['HOCHO(59)', 'OCHCHO(48)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['[O]C(O)C1OC1[O](8455)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['[O]C1OC(O)C1[O](7934)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[O]C(O)C(=O)C=O(8456)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(25.1243,'m^3/(mol*s)'), n=1.86, Ea=(32.426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-DeNd_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[O]C(C=O)C(=O)O(8457)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O][CH]O(171)', 'OCHCHO(48)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(210.531,'m^3/(mol*s)'), n=1.35673, Ea=(38.9443,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-DeH_O;YJ] for rate rule [CO-DeH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['HCO(16)', '[O]C(O)C=O(4774)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;CO_pri_rad] for rate rule [CO-CsH_O;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HOCHO(59)', '[O]C=C[O](2537)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(13.7458,'m^3/(mol*s)'), n=1.39198, Ea=(144.232,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-NdH_O;YJ] for rate rule [CO-NdH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 143.1 to 144.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)', '[O]C(C=O)C=O(3641)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0849,'cm^3/(mol*s)'), n=3.486, Ea=(94.7258,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [CO-CsH_O;OJ] for rate rule [CO-CsH_O;OJ_pri]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['[O]C(O)[C](O)C=O(8458)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['[O]C(C=O)[C](O)O(8459)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['[O][C](O)C(O)C=O(8460)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['[O]C(O)C(O)[C]=O(8461)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['[O][C](C=O)C(O)O(8462)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['[O]C([C]=O)C(O)O(8463)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.75172e+06,'s^-1'), n=1.80068, Ea=(127.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['[O]C([O])C(O)C=O(8464)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(94.0113,'s^-1'), n=2.81534, Ea=(105.641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][CH]O(171)', '[O]C=C[O](2537)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['[O]C(O)C1[CH]OO1(8465)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(265.786,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['[O]C1[CH]OOC1O(7968)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.81184e+09,'s^-1'), n=0.551229, Ea=(183.505,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 182.5 to 183.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['O=CC(=O)C(O)O(8466)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['O=CC(O)C(=O)O(8467)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CO(12)', '[O]CC([O])O(980)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C(O)C([O])C=O(1315)'],
    products = ['O=CC1OOC1O(1319)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C([CH]OO)C=O(8468)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.263e+10,'s^-1'), n=0, Ea=(92.7451,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnOOH;C_rad_out_H/NonDeC] for rate rule [ROOH;C_rad_out_H/NonDeC]
Euclidian distance = 1.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][CH]C(C=O)OO(8469)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OOH_S;Y_rad_out]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C(=CO)C([O])O(8470)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(149.951,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol
Ea raised from 147.5 to 150.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['O(4)', '[O]C(O)[CH]C=O(7809)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['O(4)', '[O]C([CH]O)C=O(3643)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '901',
    isomers = [
        '[O]C(O)C([O])C=O(1315)',
    ],
    reactants = [
        ('HOCHO(59)', 'OCHCHO(48)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '901',
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

