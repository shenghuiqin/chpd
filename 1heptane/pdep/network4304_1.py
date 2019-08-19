species(
    label = '[CH2]C(O)(C[C]=CC)C(=O)CC(19830)',
    structure = SMILES('[CH2]C(O)(C[C]=CC)C(=O)CC'),
    E0 = (6.24915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,513.982,636.666,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.4157,0.149235,-0.000195039,1.4321e-07,-4.26962e-11,975.177,44.2594], Tmin=(100,'K'), Tmax=(816.956,'K')), NASAPolynomial(coeffs=[16.461,0.0568149,-2.53555e-05,4.74827e-09,-3.26924e-13,-2109.24,-42.9941], Tmin=(816.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.24915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJC(C)(C=O)O) + radical(Cds_S)"""),
)

species(
    label = 'C=C(O)C(=O)CC(4626)',
    structure = SMILES('C=C(O)C(=O)CC'),
    E0 = (-358.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,332.629,332.633],'cm^-1')),
        HinderedRotor(inertia=(0.152178,'amu*angstrom^2'), symmetry=1, barrier=(11.9524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152152,'amu*angstrom^2'), symmetry=1, barrier=(11.9478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152157,'amu*angstrom^2'), symmetry=1, barrier=(11.9534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152209,'amu*angstrom^2'), symmetry=1, barrier=(11.954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4071.48,'J/mol'), sigma=(6.49965,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.96 K, Pc=33.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.785843,0.0702816,-6.54766e-05,3.23543e-08,-6.53762e-12,-42997.7,22.9463], Tmin=(100,'K'), Tmax=(1175.98,'K')), NASAPolynomial(coeffs=[12.9689,0.0288414,-1.26178e-05,2.38822e-09,-1.6711e-13,-45863,-37.8049], Tmin=(1175.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(O)=C([O])CC(4557)',
    structure = SMILES('[CH2]C(O)=C([O])CC'),
    E0 = (-184.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,402.631,402.631],'cm^-1')),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100034,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4303.36,'J/mol'), sigma=(7.05412,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=672.18 K, Pc=27.82 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.516209,0.0797578,-7.96062e-05,3.90546e-08,-7.19184e-12,-22064.1,30.4076], Tmin=(100,'K'), Tmax=(1528.49,'K')), NASAPolynomial(coeffs=[20.9225,0.0119894,-1.65417e-06,6.23996e-11,2.30461e-15,-27255.3,-77.6604], Tmin=(1528.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC=[C]CC1(O)CC1([O])CC(25697)',
    structure = SMILES('CC=[C]CC1(O)CC1([O])CC'),
    E0 = (100.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.05063,0.122802,-0.000109174,5.06503e-08,-9.43112e-12,12347.1,42.2292], Tmin=(100,'K'), Tmax=(1289.68,'K')), NASAPolynomial(coeffs=[23.462,0.0436732,-1.71402e-05,3.07589e-09,-2.08945e-13,5766.5,-87.3447], Tmin=(1289.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(Cds_S) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1(O)CC(=CC)C1([O])CC(25698)',
    structure = SMILES('[CH2]C1(O)CC(=CC)C1([O])CC'),
    E0 = (45.8055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[9.13945,0.0453452,6.46628e-05,-8.85869e-08,2.10516e-11,5184.8,-5.53165], Tmin=(100,'K'), Tmax=(1811.35,'K')), NASAPolynomial(coeffs=[103.555,0.0212895,-6.81555e-05,1.65126e-08,-1.21294e-12,-59276.7,-600.644], Tmin=(1811.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.8055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(methylenecyclobutane) + radical(CJC(C)2O) + radical(C=CC(C)2OJ)"""),
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
    label = '[CH2]C(O)(C=C=CC)C(=O)CC(25699)',
    structure = SMILES('[CH2]C(O)(C=C=CC)C(=O)CC'),
    E0 = (-66.7564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1705.35,2812.68,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154791,'amu*angstrom^2'), symmetry=1, barrier=(3.55895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154791,'amu*angstrom^2'), symmetry=1, barrier=(3.55895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154791,'amu*angstrom^2'), symmetry=1, barrier=(3.55895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154791,'amu*angstrom^2'), symmetry=1, barrier=(3.55895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154791,'amu*angstrom^2'), symmetry=1, barrier=(3.55895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154791,'amu*angstrom^2'), symmetry=1, barrier=(3.55895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154791,'amu*angstrom^2'), symmetry=1, barrier=(3.55895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56233,0.137768,-0.000156439,6.70696e-08,1.09436e-11,-7843.43,39.2286], Tmin=(100,'K'), Tmax=(574.9,'K')), NASAPolynomial(coeffs=[13.9078,0.0591872,-2.72203e-05,5.13731e-09,-3.53958e-13,-10102.3,-31.0181], Tmin=(574.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.7564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = '[CH2]C(O)(CC#CC)C(=O)CC(25700)',
    structure = SMILES('[CH2]C(O)(CC#CC)C(=O)CC'),
    E0 = (-71.6448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2100,2250,500,550,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,230.053,919.456,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152331,'amu*angstrom^2'), symmetry=1, barrier=(3.50239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152331,'amu*angstrom^2'), symmetry=1, barrier=(3.50239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152331,'amu*angstrom^2'), symmetry=1, barrier=(3.50239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152331,'amu*angstrom^2'), symmetry=1, barrier=(3.50239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152331,'amu*angstrom^2'), symmetry=1, barrier=(3.50239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152331,'amu*angstrom^2'), symmetry=1, barrier=(3.50239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152331,'amu*angstrom^2'), symmetry=1, barrier=(3.50239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152331,'amu*angstrom^2'), symmetry=1, barrier=(3.50239,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.18714,0.140731,-0.000172182,1.12719e-07,-2.82325e-11,-8398.15,42.4273], Tmin=(100,'K'), Tmax=(772.155,'K')), NASAPolynomial(coeffs=[18.2996,0.0482397,-1.8996e-05,3.33051e-09,-2.20652e-13,-11968.4,-53.7443], Tmin=(772.155,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.6448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CJC(C)(C=O)O)"""),
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
    label = 'C2H5CO(71)',
    structure = SMILES('CC[C]=O'),
    E0 = (-44.1874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.192008,'amu*angstrom^2'), symmetry=1, barrier=(4.41465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191563,'amu*angstrom^2'), symmetry=1, barrier=(4.40441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91313,0.0225012,-8.59328e-06,4.80319e-10,2.48743e-13,-5274.24,14.655], Tmin=(100,'K'), Tmax=(1673.18,'K')), NASAPolynomial(coeffs=[7.46306,0.0158512,-6.42141e-06,1.12499e-09,-7.32055e-14,-7388.54,-11.406], Tmin=(1673.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.1874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2CO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C(O)C[C]=CC(25072)',
    structure = SMILES('C=C(O)C[C]=CC'),
    E0 = (76.0074,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,243.467,243.872],'cm^-1')),
        HinderedRotor(inertia=(0.305623,'amu*angstrom^2'), symmetry=1, barrier=(12.8619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298837,'amu*angstrom^2'), symmetry=1, barrier=(12.8646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.29943,'amu*angstrom^2'), symmetry=1, barrier=(12.8639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296566,'amu*angstrom^2'), symmetry=1, barrier=(12.8554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.500333,0.067629,-4.63891e-05,8.81413e-09,2.3976e-12,9275.75,27.7218], Tmin=(100,'K'), Tmax=(1014.13,'K')), NASAPolynomial(coeffs=[15.9809,0.0244511,-8.97317e-06,1.60437e-09,-1.11148e-13,5216.35,-51.7138], Tmin=(1014.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.0074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = 'C#CCC([CH2])(O)C(=O)CC(25701)',
    structure = SMILES('C#CCC([CH2])(O)C(=O)CC'),
    E0 = (-29.4637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,216.738,1534.09,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149247,'amu*angstrom^2'), symmetry=1, barrier=(3.43149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149247,'amu*angstrom^2'), symmetry=1, barrier=(3.43149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149247,'amu*angstrom^2'), symmetry=1, barrier=(3.43149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149247,'amu*angstrom^2'), symmetry=1, barrier=(3.43149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149247,'amu*angstrom^2'), symmetry=1, barrier=(3.43149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149247,'amu*angstrom^2'), symmetry=1, barrier=(3.43149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149247,'amu*angstrom^2'), symmetry=1, barrier=(3.43149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.172,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80063,0.132931,-0.000174186,1.19002e-07,-3.10237e-11,-3339.54,37.6187], Tmin=(100,'K'), Tmax=(761.012,'K')), NASAPolynomial(coeffs=[18.3544,0.0398886,-1.62125e-05,2.88047e-09,-1.91886e-13,-6780.62,-56.5669], Tmin=(761.012,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.4637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)(C=O)O)"""),
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
    label = 'C=C(C[C]=CC)C(=O)CC(25702)',
    structure = SMILES('C=C(C[C]=CC)C(=O)CC'),
    E0 = (89.1883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.199,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.572228,0.0893009,-3.4008e-05,-7.92618e-08,8.99511e-11,10836.2,34.5588], Tmin=(100,'K'), Tmax=(477.497,'K')), NASAPolynomial(coeffs=[5.33674,0.0692098,-3.31609e-05,6.48988e-09,-4.61034e-13,10155.3,12.7285], Tmin=(477.497,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.1883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)([CH]C=CC)C(=O)CC(19833)',
    structure = SMILES('[CH2]C(O)([CH]C=CC)C(=O)CC'),
    E0 = (-114.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.44027,0.14664,-0.000170873,1.07951e-07,-2.75899e-11,-13564.6,41.0766], Tmin=(100,'K'), Tmax=(948.498,'K')), NASAPolynomial(coeffs=[19.6413,0.0535197,-2.36105e-05,4.44758e-09,-3.09483e-13,-17753.5,-64.2869], Tmin=(948.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJCO) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C(O)(CC=[C]C)C(=O)CC(19827)',
    structure = SMILES('[CH2]C(O)(CC=[C]C)C(=O)CC'),
    E0 = (6.24915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,513.982,636.666,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.4157,0.149235,-0.000195039,1.4321e-07,-4.26962e-11,975.177,44.2594], Tmin=(100,'K'), Tmax=(816.956,'K')), NASAPolynomial(coeffs=[16.461,0.0568149,-2.53555e-05,4.74827e-09,-3.26924e-13,-2109.24,-42.9941], Tmin=(816.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.24915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJC(C)(C=O)O) + radical(Cds_S)"""),
)

species(
    label = 'CC=[C]CC(C)([O])C(=O)CC(25703)',
    structure = SMILES('CC=[C]CC(C)([O])C(=O)CC'),
    E0 = (36.6205,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1685,370,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,375,552.5,462.5,1710,180,180,180,180,884.708,1341.81,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.168852,'amu*angstrom^2'), symmetry=1, barrier=(3.88225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168852,'amu*angstrom^2'), symmetry=1, barrier=(3.88225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168852,'amu*angstrom^2'), symmetry=1, barrier=(3.88225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168852,'amu*angstrom^2'), symmetry=1, barrier=(3.88225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168852,'amu*angstrom^2'), symmetry=1, barrier=(3.88225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168852,'amu*angstrom^2'), symmetry=1, barrier=(3.88225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168852,'amu*angstrom^2'), symmetry=1, barrier=(3.88225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.76341,0.137422,-0.000186362,1.55197e-07,-5.3188e-11,4601.73,44.703], Tmin=(100,'K'), Tmax=(805.601,'K')), NASAPolynomial(coeffs=[9.629,0.0673997,-3.09262e-05,5.83277e-09,-4.01498e-13,3202.85,-5.08589], Tmin=(805.601,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.6205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'CC=[C][CH]C(C)(O)C(=O)CC(25704)',
    structure = SMILES('CC=[C][CH]C(C)(O)C(=O)CC'),
    E0 = (-88.4455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00822,0.139873,-0.000164285,1.09904e-07,-3.03727e-11,-10427.9,40.3881], Tmin=(100,'K'), Tmax=(872.344,'K')), NASAPolynomial(coeffs=[15.3749,0.060169,-2.72417e-05,5.17685e-09,-3.61211e-13,-13460.9,-41.1023], Tmin=(872.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-88.4455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)',
    structure = SMILES('[CH2]C=CCC([CH2])(O)C(=O)CC'),
    E0 = (-80.0934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,431.394,712.751,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4713.9,'J/mol'), sigma=(7.86321,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=736.30 K, Pc=22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.28638,0.141016,-0.000157743,9.51904e-08,-2.31948e-11,-9408.87,43.6005], Tmin=(100,'K'), Tmax=(994.152,'K')), NASAPolynomial(coeffs=[20.0925,0.0509721,-2.18793e-05,4.08028e-09,-2.82779e-13,-13858.4,-64.2329], Tmin=(994.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.0934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C([O])(CC=CC)C(=O)CC(19836)',
    structure = SMILES('[CH2]C([O])(CC=CC)C(=O)CC'),
    E0 = (10.39,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,180,180,180,180,965.576,1262.15,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.64891,0.136361,-0.000157514,9.15868e-08,-1.42864e-11,1441.84,43.5103], Tmin=(100,'K'), Tmax=(626.327,'K')), NASAPolynomial(coeffs=[13.1203,0.062199,-2.81876e-05,5.32467e-09,-3.68776e-13,-803.683,-23.9894], Tmin=(626.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJC(C)(C=O)O) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'C[CH]C(=O)C(C)(O)C[C]=CC(25705)',
    structure = SMILES('CC=[C]CC(C)(O)C([O])=CC'),
    E0 = (-68.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.73662,0.12415,-0.000121014,6.42887e-08,-1.39034e-11,-8077.23,45.0739], Tmin=(100,'K'), Tmax=(1110.54,'K')), NASAPolynomial(coeffs=[19.0357,0.0493333,-1.99628e-05,3.6286e-09,-2.48243e-13,-12691.1,-57.3188], Tmin=(1110.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)(CC=CC)C(=O)[CH]C(19837)',
    structure = SMILES('[CH2]C(O)(CC=CC)C([O])=CC'),
    E0 = (-93.2884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.86455,0.125971,-0.000123382,6.52485e-08,-1.39761e-11,-11006.2,45.934], Tmin=(100,'K'), Tmax=(1123.85,'K')), NASAPolynomial(coeffs=[20.049,0.0479755,-1.92803e-05,3.49494e-09,-2.38838e-13,-15931.7,-62.3444], Tmin=(1123.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-93.2884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'C[C]=[C]CC(C)(O)C(=O)CC(25706)',
    structure = SMILES('C[C]=[C]CC(C)(O)C(=O)CC'),
    E0 = (32.4796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,1670,1700,300,440,3615,1277.5,1000,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,180,180,180,437.765,714.554,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161267,'amu*angstrom^2'), symmetry=1, barrier=(3.70784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161267,'amu*angstrom^2'), symmetry=1, barrier=(3.70784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161267,'amu*angstrom^2'), symmetry=1, barrier=(3.70784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161267,'amu*angstrom^2'), symmetry=1, barrier=(3.70784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161267,'amu*angstrom^2'), symmetry=1, barrier=(3.70784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161267,'amu*angstrom^2'), symmetry=1, barrier=(3.70784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161267,'amu*angstrom^2'), symmetry=1, barrier=(3.70784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161267,'amu*angstrom^2'), symmetry=1, barrier=(3.70784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12647,0.144408,-0.000196528,1.57661e-07,-5.18107e-11,4117.81,44.067], Tmin=(100,'K'), Tmax=(798.088,'K')), NASAPolynomial(coeffs=[12.6911,0.0625584,-2.84378e-05,5.34308e-09,-3.67196e-13,1994.2,-22.5642], Tmin=(798.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.4796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC(=O)C(C)(O)C[C]=CC(25707)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)C[C]=CC'),
    E0 = (6.22689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.14241,0.143623,-0.000185862,1.38886e-07,-4.2593e-11,962.222,44.2941], Tmin=(100,'K'), Tmax=(792.011,'K')), NASAPolynomial(coeffs=[14.4808,0.059669,-2.68635e-05,5.05242e-09,-3.48698e-13,-1670.97,-32.0275], Tmin=(792.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.22689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJCC=O) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)CC=CC(19838)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)CC=CC'),
    E0 = (-20.0035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,553.142,596.05,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.4762,0.149084,-0.000187158,1.28954e-07,-3.5865e-11,-2178.54,44.6396], Tmin=(100,'K'), Tmax=(875.607,'K')), NASAPolynomial(coeffs=[18.4246,0.0536016,-2.35821e-05,4.40834e-09,-3.04197e-13,-5838.62,-53.4177], Tmin=(875.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.0035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJCC=O) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C=[C]CC(C)(O)C(=O)CC(19831)',
    structure = SMILES('[CH2]C=[C]CC(C)(O)C(=O)CC'),
    E0 = (-53.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,488.561,660.355,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7944,0.133513,-0.000148452,9.34466e-08,-2.43116e-11,-6274.75,42.6991], Tmin=(100,'K'), Tmax=(923.946,'K')), NASAPolynomial(coeffs=[15.6761,0.0578777,-2.56587e-05,4.84455e-09,-3.37411e-13,-9503.05,-40.2036], Tmin=(923.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'CC=[C]CC1(O)CO[C]1CC(23228)',
    structure = SMILES('CC=[C]CC1(O)CO[C]1CC'),
    E0 = (91.9846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,3150,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67108,0.123883,-0.000128707,7.82667e-08,-1.95111e-11,11268.1,41.5313], Tmin=(100,'K'), Tmax=(971.297,'K')), NASAPolynomial(coeffs=[15.5037,0.0531531,-1.94771e-05,3.29476e-09,-2.14184e-13,7931.71,-40.8269], Tmin=(971.297,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.9846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Oxetane) + radical(C2CsJOCs) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1(O)CC(=CC)O[C]1CC(25708)',
    structure = SMILES('[CH2]C1(O)CC(=CC)O[C]1CC'),
    E0 = (-47.9848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98108,0.113356,-6.9775e-05,1.50198e-09,9.72194e-12,-5539.82,41.5467], Tmin=(100,'K'), Tmax=(996.319,'K')), NASAPolynomial(coeffs=[26.2959,0.038,-1.37909e-05,2.49426e-09,-1.75786e-13,-13068.9,-104.277], Tmin=(996.319,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.9848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopentane) + radical(C2CsJOC(O)) + radical(CJC(C)2O)"""),
)

species(
    label = 'CH3CH2OH(54)',
    structure = SMILES('CCO'),
    E0 = (-248.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00248522,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168417,'amu*angstrom^2'), symmetry=1, barrier=(8.09857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0684,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3014.83,'J/mol'), sigma=(4.53,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.8587,-0.00374017,6.95554e-05,-8.86548e-08,3.51688e-11,-29996.1,4.80185], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.56244,0.0152042,-5.38968e-06,8.6225e-10,-5.12898e-14,-31525.6,-9.47302], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-248.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2OH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]C(=C=O)C[C]=CC(25709)',
    structure = SMILES('C=C([C]=O)C[C]=CC'),
    E0 = (323.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,1855,455,950,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,284.694,284.917],'cm^-1')),
        HinderedRotor(inertia=(0.190734,'amu*angstrom^2'), symmetry=1, barrier=(10.9863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189501,'amu*angstrom^2'), symmetry=1, barrier=(10.9857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190783,'amu*angstrom^2'), symmetry=1, barrier=(10.9862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190418,'amu*angstrom^2'), symmetry=1, barrier=(10.9883,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0574,0.072013,-6.96937e-05,4.11883e-08,-1.099e-11,38991.4,28.8802], Tmin=(100,'K'), Tmax=(861.162,'K')), NASAPolynomial(coeffs=[6.92652,0.0447516,-2.22087e-05,4.42792e-09,-3.18277e-13,37980.6,1.44237], Tmin=(861.162,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(C=C(C)CJ=O)"""),
)

species(
    label = 'CC=C=CC(C)(O)C(=O)CC(25710)',
    structure = SMILES('CC=C=CC(C)(O)C(=O)CC'),
    E0 = (-279.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.79289,0.135616,-0.000162438,1.14089e-07,-3.33595e-11,-33472.6,39.8871], Tmin=(100,'K'), Tmax=(824.452,'K')), NASAPolynomial(coeffs=[13.4529,0.0616513,-2.78741e-05,5.28343e-09,-3.67542e-13,-35986.6,-30.7225], Tmin=(824.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds)"""),
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
    label = '[CH2]C(O)(C[C]=CC)C(C)=O(25711)',
    structure = SMILES('[CH2]C(O)(C[C]=CC)C(C)=O'),
    E0 = (31.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1765.48,2745.71,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156603,'amu*angstrom^2'), symmetry=1, barrier=(3.60062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156603,'amu*angstrom^2'), symmetry=1, barrier=(3.60062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156603,'amu*angstrom^2'), symmetry=1, barrier=(3.60062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156603,'amu*angstrom^2'), symmetry=1, barrier=(3.60062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156603,'amu*angstrom^2'), symmetry=1, barrier=(3.60062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156603,'amu*angstrom^2'), symmetry=1, barrier=(3.60062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156603,'amu*angstrom^2'), symmetry=1, barrier=(3.60062,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49963,0.125421,-0.000150032,9.80065e-08,-2.5857e-11,4002.75,39.6408], Tmin=(100,'K'), Tmax=(921.124,'K')), NASAPolynomial(coeffs=[16.9301,0.0453878,-1.96977e-05,3.67443e-09,-2.53917e-13,607.628,-47.7574], Tmin=(921.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJC(C)(C=O)O) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)(C[C]=C)C(=O)CC(25688)',
    structure = SMILES('[CH2]C(O)(C[C]=C)C(=O)CC'),
    E0 = (42.2748,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,555.264,594.351,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80138,0.134055,-0.000175888,1.26558e-07,-3.66325e-11,5287.53,40.2469], Tmin=(100,'K'), Tmax=(843.587,'K')), NASAPolynomial(coeffs=[16.6162,0.0467241,-2.06022e-05,3.83739e-09,-2.63506e-13,2180.2,-45.4743], Tmin=(843.587,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.2748,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)(C=O)O) + radical(Cds_S)"""),
)

species(
    label = 'CC=[C]CCC(O)=C([O])CC(20035)',
    structure = SMILES('CC=[C]CCC(O)=C([O])CC'),
    E0 = (-84.4672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.19373,0.126097,-0.000118494,5.81967e-08,-1.13953e-11,-9927.56,45.8874], Tmin=(100,'K'), Tmax=(1236.02,'K')), NASAPolynomial(coeffs=[23.9646,0.0414437,-1.57619e-05,2.78675e-09,-1.88069e-13,-16394,-85.8546], Tmin=(1236.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.4672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC=[C]C[C](O)CC(=O)CC(25712)',
    structure = SMILES('CC=[C]C[C](O)CC(=O)CC'),
    E0 = (-21.1957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00146,0.146725,-0.000218865,1.94482e-07,-6.89663e-11,-2347.36,45.1227], Tmin=(100,'K'), Tmax=(829.28,'K')), NASAPolynomial(coeffs=[8.0455,0.0712974,-3.36557e-05,6.37835e-09,-4.38116e-13,-3086.46,4.12336], Tmin=(829.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.1957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C2CsJOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)',
    structure = SMILES('[CH2]C(=CC)C([CH2])(O)C(=O)CC'),
    E0 = (-94.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1589.73,1600,2920.46,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4732.62,'J/mol'), sigma=(7.88591,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=739.22 K, Pc=21.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39668,0.149333,-0.000195449,1.45429e-07,-4.41529e-11,-11190.4,42.2789], Tmin=(100,'K'), Tmax=(801.412,'K')), NASAPolynomial(coeffs=[15.6661,0.0591689,-2.66711e-05,5.01446e-09,-3.45918e-13,-14085.2,-40.863], Tmin=(801.412,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = 'CC=C1CC(O)(C1)C(=O)CC(25713)',
    structure = SMILES('CC=C1CC(O)(C1)C(=O)CC'),
    E0 = (-302.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[9.68117,0.0384734,7.46587e-05,-9.40387e-08,2.20453e-11,-36736.4,-9.25191], Tmin=(100,'K'), Tmax=(1809.89,'K')), NASAPolynomial(coeffs=[101.685,0.0239961,-6.98659e-05,1.68512e-08,-1.23581e-12,-100972,-593.156], Tmin=(1809.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-302.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + ring(methylenecyclobutane)"""),
)

species(
    label = '[CH2]C(C[C]=CC)=C(CC)OO(25714)',
    structure = SMILES('[CH2]C(C[C]=CC)=C(CC)OO'),
    E0 = (254.275,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1685,370,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19811,0.119094,-0.000108859,5.54499e-08,-1.19474e-11,30765.5,43.6301], Tmin=(100,'K'), Tmax=(1085.14,'K')), NASAPolynomial(coeffs=[14.8659,0.0598807,-2.70105e-05,5.1664e-09,-3.63148e-13,27279,-35.1822], Tmin=(1085.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(O)=C(CC)OC[C]=CC(20006)',
    structure = SMILES('[CH2]C(O)=C(CC)OC[C]=CC'),
    E0 = (-30.2364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.70021,0.133275,-0.000129935,6.48846e-08,-1.27388e-11,-3383.46,45.5872], Tmin=(100,'K'), Tmax=(1243.98,'K')), NASAPolynomial(coeffs=[27.4462,0.0363376,-1.30452e-05,2.24048e-09,-1.491e-13,-10883.6,-106.433], Tmin=(1243.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.2364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C(O)(C[C]=CC)C(=C)OC(25715)',
    structure = SMILES('[CH2]C(O)(C[C]=CC)C(=C)OC'),
    E0 = (61.4583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1087.58,1130.78,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161683,'amu*angstrom^2'), symmetry=1, barrier=(3.71741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161683,'amu*angstrom^2'), symmetry=1, barrier=(3.71741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161683,'amu*angstrom^2'), symmetry=1, barrier=(3.71741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161683,'amu*angstrom^2'), symmetry=1, barrier=(3.71741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161683,'amu*angstrom^2'), symmetry=1, barrier=(3.71741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161683,'amu*angstrom^2'), symmetry=1, barrier=(3.71741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161683,'amu*angstrom^2'), symmetry=1, barrier=(3.71741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161683,'amu*angstrom^2'), symmetry=1, barrier=(3.71741,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16907,0.131115,-0.000131184,6.92592e-08,-1.46562e-11,7617.77,45.828], Tmin=(100,'K'), Tmax=(1142.87,'K')), NASAPolynomial(coeffs=[22.4486,0.0449526,-1.80966e-05,3.29116e-09,-2.25626e-13,1990.89,-76.2252], Tmin=(1142.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.4583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(O)(C[C]=CC)C(O)=CC(25716)',
    structure = SMILES('[CH2]C(O)(C[C]=CC)C(O)=CC'),
    E0 = (6.74859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1600,1726.44,2779.54,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155349,'amu*angstrom^2'), symmetry=1, barrier=(3.57179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155349,'amu*angstrom^2'), symmetry=1, barrier=(3.57179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155349,'amu*angstrom^2'), symmetry=1, barrier=(3.57179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155349,'amu*angstrom^2'), symmetry=1, barrier=(3.57179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155349,'amu*angstrom^2'), symmetry=1, barrier=(3.57179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155349,'amu*angstrom^2'), symmetry=1, barrier=(3.57179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155349,'amu*angstrom^2'), symmetry=1, barrier=(3.57179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155349,'amu*angstrom^2'), symmetry=1, barrier=(3.57179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.28746,0.134818,-0.000141596,7.86528e-08,-1.74452e-11,1041.11,46.425], Tmin=(100,'K'), Tmax=(1095.28,'K')), NASAPolynomial(coeffs=[22.505,0.0442735,-1.75916e-05,3.17333e-09,-2.16483e-13,-4389.73,-75.4403], Tmin=(1095.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.74859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(Cds_S)"""),
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
    label = 'CC=[C]C[C](O)C(=O)CC(25717)',
    structure = SMILES('CC=[C]CC(O)=C([O])CC'),
    E0 = (-59.4175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.61916,0.112019,-0.000106859,5.27128e-08,-1.02711e-11,-6934.21,41.6086], Tmin=(100,'K'), Tmax=(1249.91,'K')), NASAPolynomial(coeffs=[23.1183,0.0328532,-1.18532e-05,2.03919e-09,-1.35667e-13,-13118.1,-83.2534], Tmin=(1249.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.4175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([CH2])(O)C(=O)CC(5988)',
    structure = SMILES('[CH2]C([CH2])(O)C(=O)CC'),
    E0 = (-64.1774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24882,0.119041,-0.00016428,1.17739e-07,-3.31244e-11,-7533.04,31.7431], Tmin=(100,'K'), Tmax=(875.483,'K')), NASAPolynomial(coeffs=[18.1268,0.0305211,-1.26247e-05,2.26293e-09,-1.51404e-13,-10925.9,-59.1577], Tmin=(875.483,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.1774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CJC(C)(C=O)O)"""),
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
    E0 = (6.24915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (100.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (86.7265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (153.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (158.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (18.1955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (125.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (6.24915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (136.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (131.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (163.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (240.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (175.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (124.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (168.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (155.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (92.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (42.3362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (194.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (58.9676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (72.3215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (58.9676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (176.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (192.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (57.1894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (403.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (84.4962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (451.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (462.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (163.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (163.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (176.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (14.5335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (397.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (112.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (320.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (149.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (322.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (499.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['C=C(O)C(=O)CC(4626)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['CC=[C]CC1(O)CC1([O])CC(25697)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(94.5285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 90.6 to 94.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['[CH2]C1(O)CC(=CC)C1([O])CC(25698)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.98674e+07,'s^-1'), n=1.31443, Ea=(80.4773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_cddouble]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(O)(C=C=CC)C(=O)CC(25699)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.23e+08,'cm^3/(mol*s)'), n=1.64, Ea=(8.49352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2566 used for Cds-CsH_Ca;HJ
Exact match found for rate rule [Cds-CsH_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C(O)(CC#CC)C(=O)CC(25700)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.02e+09,'cm^3/(mol*s)'), n=1.64, Ea=(18.4933,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2702 used for Ct-Cs_Ct-Cs;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C(O)C(=O)CC(4626)', 'C=[C][CH]C(18176)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H5CO(71)', 'C=C(O)C[C]=CC(25072)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cds-OsCs_Cds;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'CH3CHCCH2(18175)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00429749,'m^3/(mol*s)'), n=2.45395, Ea=(45.5771,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 44.4 to 45.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(17)', 'C#CCC([CH2])(O)C(=O)CC(25701)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(178000,'cm^3/(mol*s)'), n=2.41, Ea=(30.1248,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2271 used for Ct-H_Ct-Cs;CsJ-HHH
Exact match found for rate rule [Ct-H_Ct-Cs;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OH(5)', 'C=C(C[C]=CC)C(=O)CC(25702)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.52832,'m^3/(mol*s)'), n=2.02802, Ea=(14.4047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;YJ] for rate rule [Cds-COCs_Cds;OJ_pri]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['[CH2]C(O)([CH]C=CC)C(=O)CC(19833)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/O)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)(CC=[C]C)C(=O)CC(19827)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CC=[C]CC(C)([O])C(=O)CC(25703)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['CC=[C][CH]C(C)(O)C(=O)CC(25704)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([O])(CC=CC)C(=O)CC(19836)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['C[CH]C(=O)C(C)(O)C[C]=CC(25705)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['[CH2]C(O)(CC=CC)C(=O)[CH]C(19837)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.44776e+08,'s^-1'), n=0.81, Ea=(36.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_CC(O2d)CC;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C[C]=[C]CC(C)(O)C(=O)CC(25706)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['[CH2]CC(=O)C(C)(O)C[C]=CC(25707)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['[CH2]CC(=O)C([CH2])(O)CC=CC(19838)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(60905.8,'s^-1'), n=1.5925, Ea=(66.0723,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_SSSSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['[CH2]C=[C]CC(C)(O)C(=O)CC(19831)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5436.63,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R6HJ_3;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'C=[C][CH]C(18176)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['CC=[C]CC1(O)CO[C]1CC(23228)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['[CH2]C1(O)CC(=CC)O[C]1CC(25708)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.67521e+10,'s^-1'), n=0.355, Ea=(50.9402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cddouble] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_cddouble]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH3CH2OH(54)', '[CH2]C(=C=O)C[C]=CC(25709)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['CC=C=CC(C)(O)C(=O)CC(25710)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(S)(23)', '[CH2]C(O)(C[C]=CC)C(C)=O(25711)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(S)(23)', '[CH2]C(O)(C[C]=C)C(=O)CC(25688)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['CC=[C]CCC(O)=C([O])CC(20035)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['CC=[C]C[C](O)CC(=O)CC(25712)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['CC=C1CC(O)(C1)C(=O)CC(25713)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(C[C]=CC)=C(CC)OO(25714)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(O)=C(CC)OC[C]=CC(20006)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=C)OC(25715)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(O)(C[C]=CC)C(O)=CC(25716)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CH2(19)', 'CC=[C]C[C](O)C(=O)CC(25717)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C([CH2])(O)C(=O)CC(5988)', '[C]=CC(24199)'],
    products = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4304',
    isomers = [
        '[CH2]C(O)(C[C]=CC)C(=O)CC(19830)',
    ],
    reactants = [
        ('C=C(O)C(=O)CC(4626)', 'CH3CHCCH2(18175)'),
        ('[CH2]C(O)=C([O])CC(4557)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4304',
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

