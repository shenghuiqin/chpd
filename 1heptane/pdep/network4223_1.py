species(
    label = 'C=C([CH]C)C([O])O(24164)',
    structure = SMILES('[CH2]C(=CC)C([O])O'),
    E0 = (-51.8111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,222.441,2441.78],'cm^-1')),
        HinderedRotor(inertia=(0.00340648,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231634,'amu*angstrom^2'), symmetry=1, barrier=(8.13497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231679,'amu*angstrom^2'), symmetry=1, barrier=(8.13502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.991753,'amu*angstrom^2'), symmetry=1, barrier=(34.8284,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47893,0.0634877,-4.95477e-05,5.76648e-10,2.02979e-11,-6148.44,25.87], Tmin=(100,'K'), Tmax=(533.4,'K')), NASAPolynomial(coeffs=[5.98219,0.0412814,-1.96197e-05,3.81552e-09,-2.69739e-13,-6793.36,5.43256], Tmin=(533.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.8111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCOJ)"""),
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
    label = '[CH2]C1([CH]C)OC1O(25335)',
    structure = SMILES('[CH2]C1([CH]C)OC1O'),
    E0 = (41.2064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.230231,0.0799189,-8.50696e-05,4.69036e-08,-9.73782e-12,5119.61,28.9007], Tmin=(100,'K'), Tmax=(1374.42,'K')), NASAPolynomial(coeffs=[17.1553,0.0166529,-2.19699e-06,-4.25118e-13,1.36329e-14,1537.17,-56.1506], Tmin=(1374.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.2064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCJCO) + radical(CJC(C)OC)"""),
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
    label = 'C=C([CH]C)C(=O)O(18229)',
    structure = SMILES('[CH2]C(=CC)C(=O)O'),
    E0 = (-176.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40434,0.0620366,-6.63328e-05,4.5094e-08,-1.35262e-11,-21182.6,23.9649], Tmin=(100,'K'), Tmax=(785.718,'K')), NASAPolynomial(coeffs=[6.46974,0.0362489,-1.71013e-05,3.32146e-09,-2.34831e-13,-21978.5,0.748917], Tmin=(785.718,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-176.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)O2s) + radical(C=C(C=O)CJ)"""),
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
    label = '[CH2]C(=CC)[C](O)O(25371)',
    structure = SMILES('[CH2]C([CH]C)=C(O)O'),
    E0 = (-128.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.778897,0.0839354,-8.37235e-05,4.03682e-08,-7.30482e-12,-15312.2,26.8218], Tmin=(100,'K'), Tmax=(1544.24,'K')), NASAPolynomial(coeffs=[22.8843,0.0110749,-1.71519e-06,1.13867e-10,-2.74316e-15,-21241.5,-93.1567], Tmin=(1544.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-128.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'CC=C(C)[C]([O])O(25372)',
    structure = SMILES('C[CH][C](C)C(=O)O'),
    E0 = (-169.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23808,0.0526879,-2.52971e-05,1.54006e-09,1.50176e-12,-20322.8,27.6693], Tmin=(100,'K'), Tmax=(1248.22,'K')), NASAPolynomial(coeffs=[11.7582,0.0306392,-1.28173e-05,2.36078e-09,-1.61968e-13,-23857.7,-29.0565], Tmin=(1248.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCJCC=O) + radical(CCJ(C)CO)"""),
)

species(
    label = 'C[C]=C(C)C([O])O(25373)',
    structure = SMILES('C[C]=C(C)C([O])O'),
    E0 = (34.5315,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1685,370,3615,1277.5,1000,1380,1390,370,380,2900,435,188.223,2913.08],'cm^-1')),
        HinderedRotor(inertia=(0.300285,'amu*angstrom^2'), symmetry=1, barrier=(7.54936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.300269,'amu*angstrom^2'), symmetry=1, barrier=(7.54921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.300401,'amu*angstrom^2'), symmetry=1, barrier=(7.54946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.300342,'amu*angstrom^2'), symmetry=1, barrier=(7.54951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.883825,0.0788644,-0.000122201,1.16582e-07,-4.37512e-11,4255.38,28.1075], Tmin=(100,'K'), Tmax=(829.863,'K')), NASAPolynomial(coeffs=[3.26062,0.045435,-2.20591e-05,4.22724e-09,-2.91896e-13,4617.51,21.6428], Tmin=(829.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.5315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'CC=C(C)C([O])[O](25374)',
    structure = SMILES('CC=C(C)C([O])[O]'),
    E0 = (22.3949,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,3082.91],'cm^-1')),
        HinderedRotor(inertia=(0.147579,'amu*angstrom^2'), symmetry=1, barrier=(3.39314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147687,'amu*angstrom^2'), symmetry=1, barrier=(3.39561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147549,'amu*angstrom^2'), symmetry=1, barrier=(3.39245,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9976,0.0423622,-1.77626e-05,2.18686e-09,4.75761e-15,2649.64,16.3843], Tmin=(100,'K'), Tmax=(2743.04,'K')), NASAPolynomial(coeffs=[49.8,-0.00785776,6.38113e-07,-8.29112e-11,1.09048e-14,-28712,-262.167], Tmin=(2743.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.3949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(=[C]C)C(O)O(25375)',
    structure = SMILES('[CH2]C(=[C]C)C(O)O'),
    E0 = (-39.6745,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1685,370,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,328.217],'cm^-1')),
        HinderedRotor(inertia=(0.00156487,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11696,'amu*angstrom^2'), symmetry=1, barrier=(8.941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11696,'amu*angstrom^2'), symmetry=1, barrier=(8.941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11696,'amu*angstrom^2'), symmetry=1, barrier=(8.941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201327,'amu*angstrom^2'), symmetry=1, barrier=(15.3905,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.732682,0.0763695,-9.04073e-05,6.17016e-08,-1.74638e-11,-4658.07,27.3928], Tmin=(100,'K'), Tmax=(851.436,'K')), NASAPolynomial(coeffs=[9.80393,0.0337536,-1.53302e-05,2.91736e-09,-2.03629e-13,-6202.8,-14.9118], Tmin=(851.436,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.6745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C[C](C)C([O])O(21234)',
    structure = SMILES('[CH2]C=C(C)C([O])O'),
    E0 = (-51.8111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,222.441,2441.78],'cm^-1')),
        HinderedRotor(inertia=(0.00340648,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231634,'amu*angstrom^2'), symmetry=1, barrier=(8.13497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231679,'amu*angstrom^2'), symmetry=1, barrier=(8.13502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.991753,'amu*angstrom^2'), symmetry=1, barrier=(34.8284,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47893,0.0634877,-4.95477e-05,5.76648e-10,2.02979e-11,-6148.44,25.87], Tmin=(100,'K'), Tmax=(533.4,'K')), NASAPolynomial(coeffs=[5.98219,0.0412814,-1.96197e-05,3.81552e-09,-2.69739e-13,-6793.36,5.43256], Tmin=(533.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.8111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C](C=C)C(O)O(21238)',
    structure = SMILES('[CH2]C=C([CH2])C(O)O'),
    E0 = (-126.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.550469,0.0719644,-6.70508e-05,3.2548e-08,-6.34851e-12,-15028.9,27.8414], Tmin=(100,'K'), Tmax=(1231.69,'K')), NASAPolynomial(coeffs=[14.8562,0.0255053,-1.04708e-05,1.92326e-09,-1.32476e-13,-18552.9,-44.1566], Tmin=(1231.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-126.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'CC1C[C]1C([O])O(25376)',
    structure = SMILES('CC1C[C]1C([O])O'),
    E0 = (-5.35034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55287,0.0421716,7.81743e-06,-3.49073e-08,1.46459e-11,-545.353,26.5474], Tmin=(100,'K'), Tmax=(1040.32,'K')), NASAPolynomial(coeffs=[11.6939,0.029603,-1.21595e-05,2.30973e-09,-1.65004e-13,-4085.18,-29.6502], Tmin=(1040.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.35034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CCOJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2][C]1C(C)OC1O(25377)',
    structure = SMILES('[CH2][C]1C(C)OC1O'),
    E0 = (-8.13915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59117,0.038547,3.38873e-05,-7.99115e-08,3.76508e-11,-878.444,24.4308], Tmin=(100,'K'), Tmax=(877.9,'K')), NASAPolynomial(coeffs=[14.5052,0.0198945,-2.90893e-06,1.75373e-10,-4.74739e-15,-4694.55,-45.0105], Tmin=(877.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.13915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(Isobutyl) + radical(CCJ(C)CO)"""),
)

species(
    label = 'CC=C(C)C(=O)O(25378)',
    structure = SMILES('CC=C(C)C(=O)O'),
    E0 = (-333.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29215,0.0660362,-8.41357e-05,7.55547e-08,-2.8889e-11,-39985.1,24.8501], Tmin=(100,'K'), Tmax=(774.12,'K')), NASAPolynomial(coeffs=[3.64612,0.0436436,-2.09249e-05,4.04822e-09,-2.83526e-13,-40043.1,16.0759], Tmin=(774.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-333.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)O2s)"""),
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
    label = '[CH2]C(=C)C([O])O(25379)',
    structure = SMILES('[CH2]C(=C)C([O])O'),
    E0 = (-15.7855,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,290.536,1979.31],'cm^-1')),
        HinderedRotor(inertia=(0.00203006,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130681,'amu*angstrom^2'), symmetry=1, barrier=(7.74076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.607176,'amu*angstrom^2'), symmetry=1, barrier=(35.4766,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76335,0.0535403,-5.74598e-05,3.88681e-08,-1.15784e-11,-1821.9,22.2762], Tmin=(100,'K'), Tmax=(791.134,'K')), NASAPolynomial(coeffs=[6.21877,0.0310127,-1.47457e-05,2.87284e-09,-2.03466e-13,-2526.84,1.82536], Tmin=(791.134,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-15.7855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(Allyl_P)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4287.55,'J/mol'), sigma=(7.04051,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=669.70 K, Pc=27.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10312,0.0704205,-9.28953e-05,8.23788e-08,-3.07053e-11,6224.63,29.478], Tmin=(100,'K'), Tmax=(781.226,'K')), NASAPolynomial(coeffs=[4.4893,0.0429091,-2.05377e-05,3.96225e-09,-2.76788e-13,6006.01,15.9647], Tmin=(781.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'CC=C1COC1O(25334)',
    structure = SMILES('CC=C1COC1O'),
    E0 = (-233.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50312,0.0441335,7.75272e-07,-2.65439e-08,1.14001e-11,-27971,23.8559], Tmin=(100,'K'), Tmax=(1070.3,'K')), NASAPolynomial(coeffs=[11.5536,0.0301506,-1.26728e-05,2.41542e-09,-1.71938e-13,-31473,-31.6238], Tmin=(1070.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C([CH]OO)=CC(25380)',
    structure = SMILES('[CH2]C([CH]OO)=CC'),
    E0 = (118.207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1000.1],'cm^-1')),
        HinderedRotor(inertia=(0.0331587,'amu*angstrom^2'), symmetry=1, barrier=(23.5354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02357,'amu*angstrom^2'), symmetry=1, barrier=(23.5339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02351,'amu*angstrom^2'), symmetry=1, barrier=(23.5326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.913508,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.609889,0.0639826,-3.79664e-05,2.82203e-09,3.43266e-12,14348.1,27.295], Tmin=(100,'K'), Tmax=(1073.94,'K')), NASAPolynomial(coeffs=[15.9141,0.025656,-1.0519e-05,1.97575e-09,-1.39657e-13,9983.93,-52.6447], Tmin=(1073.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(C=CCJO) + radical(Allyl_P)"""),
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
    label = '[CH2]C([CH]O)=CC(24557)',
    structure = SMILES('[CH2]C([CH]C)=CO'),
    E0 = (28.2424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,355.272],'cm^-1')),
        HinderedRotor(inertia=(0.263834,'amu*angstrom^2'), symmetry=1, barrier=(23.6945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256376,'amu*angstrom^2'), symmetry=1, barrier=(23.6399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266715,'amu*angstrom^2'), symmetry=1, barrier=(23.6255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261247,'amu*angstrom^2'), symmetry=1, barrier=(23.6444,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.980735,0.0481604,1.39064e-05,-6.56205e-08,3.25648e-11,3522.37,20.7053], Tmin=(100,'K'), Tmax=(925.665,'K')), NASAPolynomial(coeffs=[20.6065,0.0103929,-1.11864e-06,9.84218e-11,-1.08904e-14,-2126.31,-83.3472], Tmin=(925.665,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.2424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Allyl_S)"""),
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
    label = 'CC=[C]C([O])O(25381)',
    structure = SMILES('CC=[C]C([O])O'),
    E0 = (73.5865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,180,2732.86],'cm^-1')),
        HinderedRotor(inertia=(0.335777,'amu*angstrom^2'), symmetry=1, barrier=(7.72017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.335529,'amu*angstrom^2'), symmetry=1, barrier=(7.71448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.335586,'amu*angstrom^2'), symmetry=1, barrier=(7.71579,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68891,0.0593122,-9.31679e-05,9.13307e-08,-3.50199e-11,8925.39,24.4611], Tmin=(100,'K'), Tmax=(827.711,'K')), NASAPolynomial(coeffs=[2.546,0.0364827,-1.79296e-05,3.45422e-09,-2.39227e-13,9423.66,24.3552], Tmin=(827.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.5865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'C=CC(=C)C([O])O(21232)',
    structure = SMILES('C=CC(=C)C([O])O'),
    E0 = (-78.0227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,397.009,397.011,397.011],'cm^-1')),
        HinderedRotor(inertia=(0.0999754,'amu*angstrom^2'), symmetry=1, barrier=(11.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0999757,'amu*angstrom^2'), symmetry=1, barrier=(11.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0999754,'amu*angstrom^2'), symmetry=1, barrier=(11.1821,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19147,0.0641413,-6.094e-05,3.20755e-08,-7.07393e-12,-9284.8,25.4075], Tmin=(100,'K'), Tmax=(1068.6,'K')), NASAPolynomial(coeffs=[10.1746,0.0305157,-1.37397e-05,2.62886e-09,-1.84902e-13,-11204.7,-18.527], Tmin=(1068.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.0227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]CC(=C)C([O])O(25382)',
    structure = SMILES('[CH2]CC(=C)C([O])O'),
    E0 = (15.3162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,370.664,1482.84,1482.84],'cm^-1')),
        HinderedRotor(inertia=(0.0762109,'amu*angstrom^2'), symmetry=1, barrier=(7.43025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0762112,'amu*angstrom^2'), symmetry=1, barrier=(7.43025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0762108,'amu*angstrom^2'), symmetry=1, barrier=(7.43025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.076211,'amu*angstrom^2'), symmetry=1, barrier=(7.43025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63698,0.0608584,-3.78085e-05,-2.72259e-08,4.30286e-11,1918.32,27.9513], Tmin=(100,'K'), Tmax=(499.64,'K')), NASAPolynomial(coeffs=[5.86524,0.0408707,-1.94205e-05,3.77001e-09,-2.65915e-13,1322.76,8.75469], Tmin=(499.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.3162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C](CC)C(=O)O(2673)',
    structure = SMILES('[CH2][C](CC)C(=O)O'),
    E0 = (-159.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38598,0.0607249,-4.54771e-05,1.86185e-08,-3.32491e-12,-19061.7,25.6718], Tmin=(100,'K'), Tmax=(1244.6,'K')), NASAPolynomial(coeffs=[8.8328,0.0367917,-1.66327e-05,3.16809e-09,-2.21427e-13,-20915.3,-11.8843], Tmin=(1244.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-159.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CJC(C)C=O) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH]=C(CC)C([O])O(25383)',
    structure = SMILES('[CH]=C(CC)C([O])O'),
    E0 = (57.166,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,266.807,4000],'cm^-1')),
        HinderedRotor(inertia=(0.134469,'amu*angstrom^2'), symmetry=1, barrier=(6.79953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134685,'amu*angstrom^2'), symmetry=1, barrier=(6.79917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.585934,'amu*angstrom^2'), symmetry=1, barrier=(29.5689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134596,'amu*angstrom^2'), symmetry=1, barrier=(6.79993,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48783,0.0642604,-4.73486e-05,-1.46129e-08,3.55935e-11,6957,27.0948], Tmin=(100,'K'), Tmax=(512.186,'K')), NASAPolynomial(coeffs=[6.43148,0.0401215,-1.90303e-05,3.68388e-09,-2.59135e-13,6260.8,4.6995], Tmin=(512.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCOJ)"""),
)

species(
    label = 'C=C(CC)C([O])[O](25384)',
    structure = SMILES('C=C(CC)C([O])[O]'),
    E0 = (35.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,180,1254.83,1262.63,4000],'cm^-1')),
        HinderedRotor(inertia=(0.258642,'amu*angstrom^2'), symmetry=1, barrier=(5.94668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30483,'amu*angstrom^2'), symmetry=1, barrier=(30.0005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259687,'amu*angstrom^2'), symmetry=1, barrier=(5.97072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34031,0.0457218,-2.17575e-05,4.00188e-09,-2.52674e-13,4293.13,19.9002], Tmin=(100,'K'), Tmax=(2938.86,'K')), NASAPolynomial(coeffs=[49.3243,-0.00797921,1.11594e-06,-1.57955e-10,1.36646e-14,-26572.6,-258.047], Tmin=(2938.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C([CH]C)C(O)O(25385)',
    structure = SMILES('[CH]C(=CC)C(O)O'),
    E0 = (-58.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.858725,0.0730824,-6.98531e-05,3.99751e-08,-9.92792e-12,-6905.86,27.2959], Tmin=(100,'K'), Tmax=(942.756,'K')), NASAPolynomial(coeffs=[8.51141,0.0406128,-1.81909e-05,3.44202e-09,-2.40011e-13,-8348.77,-9.17274], Tmin=(942.756,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C[CH][C]1COC1O(25386)',
    structure = SMILES('C[CH][C]1COC1O'),
    E0 = (-6.03642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75185,0.038077,2.26771e-05,-5.68307e-08,2.5655e-11,-634.51,25.6621], Tmin=(100,'K'), Tmax=(911.453,'K')), NASAPolynomial(coeffs=[11.0834,0.0268488,-7.76225e-06,1.21391e-09,-7.99757e-14,-3570.24,-25.2653], Tmin=(911.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.03642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCJ(C)CO) + radical(Cs_S)"""),
)

species(
    label = 'C=C(CC)C(=O)O(883)',
    structure = SMILES('C=C(CC)C(=O)O'),
    E0 = (-319.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57897,0.05856,-5.08894e-05,2.86954e-08,-7.63048e-12,-38384.5,24.9497], Tmin=(100,'K'), Tmax=(846.972,'K')), NASAPolynomial(coeffs=[5.36873,0.0406621,-1.91922e-05,3.7461e-09,-2.66232e-13,-39026.5,7.29574], Tmin=(846.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(=C)C(O)O(21245)',
    structure = SMILES('C=CC(=C)C(O)O'),
    E0 = (-303.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.356792,0.0725237,-6.86352e-05,3.33023e-08,-6.3788e-12,-36392.5,26.8441], Tmin=(100,'K'), Tmax=(1269.9,'K')), NASAPolynomial(coeffs=[16.7067,0.0210242,-7.80454e-06,1.36792e-09,-9.20367e-14,-40545,-55.9414], Tmin=(1269.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-303.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=[C]C(C)C([O])O(21237)',
    structure = SMILES('C=[C]C(C)C([O])O'),
    E0 = (55.2145,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,417.085,417.141,2298.36],'cm^-1')),
        HinderedRotor(inertia=(0.061065,'amu*angstrom^2'), symmetry=1, barrier=(7.53853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.061025,'amu*angstrom^2'), symmetry=1, barrier=(7.53817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0610319,'amu*angstrom^2'), symmetry=1, barrier=(7.53858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125336,'amu*angstrom^2'), symmetry=1, barrier=(15.466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4232,'J/mol'), sigma=(7.00504,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=661.03 K, Pc=27.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31677,0.0641451,-6.35813e-05,4.03656e-08,-1.15084e-11,6732.88,28.9242], Tmin=(100,'K'), Tmax=(816.931,'K')), NASAPolynomial(coeffs=[6.32979,0.0395998,-1.85132e-05,3.58769e-09,-2.53614e-13,5913.81,5.75288], Tmin=(816.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.2145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'C=C1C(C)OC1O(25325)',
    structure = SMILES('C=C1C(C)OC1O'),
    E0 = (-239.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40964,0.0398403,2.85601e-05,-6.5621e-08,2.74821e-11,-28718.2,23.42], Tmin=(100,'K'), Tmax=(991.121,'K')), NASAPolynomial(coeffs=[15.7326,0.0233725,-9.07842e-06,1.77748e-09,-1.33084e-13,-33587.7,-55.795], Tmin=(991.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-239.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
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
    label = 'C=[C]C([O])O(7844)',
    structure = SMILES('C=[C]C([O])O'),
    E0 = (109.612,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,3615,1277.5,1000,1380,1390,370,380,2900,435,197.172,2159.89],'cm^-1')),
        HinderedRotor(inertia=(0.240952,'amu*angstrom^2'), symmetry=1, barrier=(6.66133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238298,'amu*angstrom^2'), symmetry=1, barrier=(6.65741,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32057,0.0438937,-7.30106e-05,7.30931e-08,-2.81388e-11,13237,20.3884], Tmin=(100,'K'), Tmax=(834.17,'K')), NASAPolynomial(coeffs=[2.64552,0.0264948,-1.32391e-05,2.55884e-09,-1.77136e-13,13733.9,22.1831], Tmin=(834.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S)"""),
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
    E0 = (-51.8111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (41.7447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (67.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (29.0007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (183.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (43.2621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (73.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (150.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (196.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (112.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (4.6341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (97.9761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (37.4127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (366.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (179.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (74.4244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (11.5891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (404.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (220.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-43.5268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (210.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (271.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (455.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (133.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (129.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (157.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (202.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (89.6649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (367.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (74.4244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (11.5891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-26.8378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (149.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-43.5268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (453.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['HOCHO(59)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['[CH2]C1([CH]C)OC1O(25335)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(93.5558,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C([CH]C)C(=O)O(18229)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(25.1243,'m^3/(mol*s)'), n=1.86, Ea=(32.426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-DeNd_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['HOCHO(59)', 'C=[C][CH]C(18176)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.87291,'m^3/(mol*s)'), n=1.39198, Ea=(57.1551,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-NdH_O;YJ] for rate rule [CO-NdH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O][CH]O(171)', 'CH3CHCCH2(18175)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OH(5)', 'C=C([CH]C)C=O(18210)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.5e+06,'cm^3/(mol*s)'), n=2.16, Ea=(19.2464,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdH_O;YJ] for rate rule [CO-CdH_O;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['[CH2]C(=CC)[C](O)O(25371)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['CC=C(C)[C]([O])O(25372)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.4947e+07,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[C]=C(C)C([O])O(25373)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC=C(C)C([O])[O](25374)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.7265e+07,'s^-1'), n=1.395, Ea=(89.7468,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;O_rad_out;Cs_H_out_2H] + [R4H_SS(Cd)S;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_SS(Cd)S;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(=[C]C)C(O)O(25375)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[C](C)C([O])O(21234)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['[CH2][C](C=C)C(O)O(21238)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.45388e+06,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_SSMS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][CH]O(171)', 'C=[C][CH]C(18176)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['CC1C[C]1C([O])O(25376)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['[CH2][C]1C(C)OC1O(25377)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['CC=C(C)C(=O)O(25378)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH2(S)(23)', '[CH2]C(=C)C([O])O(25379)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC=[C]CC([O])O(21299)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['CC=C1COC1O(25334)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([CH]OO)=CC(25380)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.263e+10,'s^-1'), n=0, Ea=(92.7451,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnOOH;C_rad_out_1H] for rate rule [ROOH;C_rad_out_H/OneDe]
Euclidian distance = 1.41421356237
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(4)', '[CH2]C([CH]O)=CC(24557)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(19)', 'CC=[C]C([O])O(25381)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', 'C=CC(=C)C([O])O(21232)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]CC(=C)C([O])O(25382)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['[CH2][C](CC)C(=O)O(2673)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.18083e+09,'s^-1'), n=1.04667, Ea=(209.2,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C(CC)C([O])O(25383)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C(CC)C([O])[O](25384)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(420000,'s^-1'), n=1.76, Ea=(53.8899,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R4H_SS(Cd)S;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['[CH]=C([CH]C)C(O)O(25385)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.4207e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_singleH] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['C[CH][C]1COC1O(25386)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['C=C(CC)C(=O)O(883)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['C=CC(=C)C(O)O(21245)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=[C]C(C)C([O])O(21237)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C([CH]C)C([O])O(24164)'],
    products = ['C=C1C(C)OC1O(25325)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CHCH3(T)(95)', 'C=[C]C([O])O(7844)'],
    products = ['C=C([CH]C)C([O])O(24164)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4223',
    isomers = [
        'C=C([CH]C)C([O])O(24164)',
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
    label = '4223',
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

