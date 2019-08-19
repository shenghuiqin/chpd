species(
    label = 'C=[C]C(C)C([O])(CC)C(=C)O(20152)',
    structure = SMILES('C=[C]C(C)C([O])(CC)C(=C)O'),
    E0 = (38.9266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1004.78,1208.81,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39537,0.125449,-0.000112778,5.22387e-08,-9.59442e-12,4924.78,47.788], Tmin=(100,'K'), Tmax=(1319.05,'K')), NASAPolynomial(coeffs=[26.1448,0.0389017,-1.43578e-05,2.49583e-09,-1.66652e-13,-2604.41,-97.805], Tmin=(1319.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.9266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)2OJ)"""),
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
    label = '[CH2]C1(O)OC1(CC)C(C)[C]=C(26107)',
    structure = SMILES('[CH2]C1(O)OC1(CC)C(C)[C]=C'),
    E0 = (85.1289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.3293,0.146384,-0.000147165,7.33518e-08,-1.36727e-11,10570.1,47.658], Tmin=(100,'K'), Tmax=(1531.16,'K')), NASAPolynomial(coeffs=[33.5177,0.0230701,-2.41479e-06,-9.83378e-11,2.21307e-14,1845.38,-141.699], Tmin=(1531.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.1289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(CJC(O)2C)"""),
)

species(
    label = '[CH2]C1(O)C(=C)C(C)C1([O])CC(25664)',
    structure = SMILES('[CH2]C1(O)C(=C)C(C)C1([O])CC'),
    E0 = (51.3196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[8.30878,0.0492268,6.07581e-05,-8.75673e-08,2.1155e-11,5890.71,-3.82484], Tmin=(100,'K'), Tmax=(1789.8,'K')), NASAPolynomial(coeffs=[100.402,0.0231125,-6.7963e-05,1.64774e-08,-1.21394e-12,-55857.9,-582.138], Tmin=(1789.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.3196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(C=CC(C)(O)CJ) + radical(CC(C)2OJ)"""),
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
    label = 'C=C=C(C)C([O])(CC)C(=C)O(26108)',
    structure = SMILES('C=C=C(C)C([O])(CC)C(=C)O'),
    E0 = (-48.0652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,180,180,180,413.339,724.225,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15526,'amu*angstrom^2'), symmetry=1, barrier=(3.56974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15526,'amu*angstrom^2'), symmetry=1, barrier=(3.56974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15526,'amu*angstrom^2'), symmetry=1, barrier=(3.56974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15526,'amu*angstrom^2'), symmetry=1, barrier=(3.56974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15526,'amu*angstrom^2'), symmetry=1, barrier=(3.56974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15526,'amu*angstrom^2'), symmetry=1, barrier=(3.56974,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.81622,0.12002,-0.000105917,4.80402e-08,-8.76723e-12,-5565.09,41.8378], Tmin=(100,'K'), Tmax=(1307.82,'K')), NASAPolynomial(coeffs=[23.1415,0.0436864,-1.83659e-05,3.41072e-09,-2.35966e-13,-12093.1,-85.2666], Tmin=(1307.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.0652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C#CC(C)C([O])(CC)C(=C)O(26109)',
    structure = SMILES('C#CC(C)C([O])(CC)C(=C)O'),
    E0 = (-33.9652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,180,180,180,191.58,935.208,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150437,'amu*angstrom^2'), symmetry=1, barrier=(3.45885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150437,'amu*angstrom^2'), symmetry=1, barrier=(3.45885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150437,'amu*angstrom^2'), symmetry=1, barrier=(3.45885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150437,'amu*angstrom^2'), symmetry=1, barrier=(3.45885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150437,'amu*angstrom^2'), symmetry=1, barrier=(3.45885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150437,'amu*angstrom^2'), symmetry=1, barrier=(3.45885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150437,'amu*angstrom^2'), symmetry=1, barrier=(3.45885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10112,0.121428,-0.000103073,3.75166e-08,-2.58513e-12,-3854.08,44.0563], Tmin=(100,'K'), Tmax=(984.963,'K')), NASAPolynomial(coeffs=[24.7494,0.0373198,-1.29565e-05,2.22219e-09,-1.50093e-13,-10352.9,-91.2147], Tmin=(984.963,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.9652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CC(C)2OJ)"""),
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
    label = 'C2H5(29)',
    structure = SMILES('C[CH2]'),
    E0 = (107.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1190.6,1642.82,1642.96,3622.23,3622.39],'cm^-1')),
        HinderedRotor(inertia=(0.866817,'amu*angstrom^2'), symmetry=1, barrier=(19.9298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24186,-0.00356905,4.82667e-05,-5.85401e-08,2.25805e-11,12969,4.44704], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.32196,0.0123931,-4.39681e-06,7.0352e-10,-4.18435e-14,12175.9,0.171104], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(107.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H5""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=[C]C(C)C(=O)C(=C)O(26110)',
    structure = SMILES('C=[C]C(C)C(=O)C(=C)O'),
    E0 = (-47.2994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.235266,0.0975538,-0.000109349,6.70158e-08,-1.68132e-11,-5539.97,31.2561], Tmin=(100,'K'), Tmax=(959.233,'K')), NASAPolynomial(coeffs=[13.7943,0.0390501,-1.78629e-05,3.43273e-09,-2.41792e-13,-8231.48,-35.8446], Tmin=(959.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.2994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'CH2COH(99)',
    structure = SMILES('C=[C]O'),
    E0 = (103.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.989114,'amu*angstrom^2'), symmetry=1, barrier=(22.7417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1624,0.0134245,5.56346e-06,-1.95511e-08,9.36369e-12,12455.2,10.1544], Tmin=(100,'K'), Tmax=(925.618,'K')), NASAPolynomial(coeffs=[8.19875,0.00453462,-8.93448e-07,1.26083e-10,-9.46513e-15,10971.3,-16.733], Tmin=(925.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2COH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=[C]C(C)C(=O)CC(24949)',
    structure = SMILES('C=[C]C(C)C(=O)CC'),
    E0 = (27.5574,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.214577,0.0905039,-0.000112418,9.35951e-08,-3.3297e-11,3443.76,30.5337], Tmin=(100,'K'), Tmax=(777.38,'K')), NASAPolynomial(coeffs=[5.85739,0.0532382,-2.46298e-05,4.68995e-09,-3.25708e-13,2815.14,6.33103], Tmin=(777.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.5574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = 'C=C=CC([O])(CC)C(=C)O(26111)',
    structure = SMILES('C=C=CC([O])(CC)C(=C)O'),
    E0 = (-9.0102,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,379.07,757.399,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154065,'amu*angstrom^2'), symmetry=1, barrier=(3.54226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154065,'amu*angstrom^2'), symmetry=1, barrier=(3.54226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154065,'amu*angstrom^2'), symmetry=1, barrier=(3.54226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154065,'amu*angstrom^2'), symmetry=1, barrier=(3.54226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154065,'amu*angstrom^2'), symmetry=1, barrier=(3.54226,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.172,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35328,0.104348,-8.98181e-05,3.88403e-08,-6.65797e-12,-879.839,39.4301], Tmin=(100,'K'), Tmax=(1403.81,'K')), NASAPolynomial(coeffs=[23.8792,0.0324512,-1.29947e-05,2.35706e-09,-1.60798e-13,-7964.18,-90.8608], Tmin=(1403.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.0102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C[C](C)C([O])(CC)C(=C)O(20149)',
    structure = SMILES('[CH2]C=C(C)C([O])(CC)C(=C)O'),
    E0 = (-73.1703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,262.04,866.396,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.24133,0.122242,-0.00010271,4.37536e-08,-7.43058e-12,-8563.24,44.6026], Tmin=(100,'K'), Tmax=(1410.73,'K')), NASAPolynomial(coeffs=[26.3926,0.0410524,-1.63822e-05,2.95775e-09,-2.00987e-13,-16642.1,-103.392], Tmin=(1410.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.1703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH]=CC(C)C([O])(CC)C(=C)O(20158)',
    structure = SMILES('[CH]=CC(C)C([O])(CC)C(=C)O'),
    E0 = (48.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,454.304,681.32,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25135,0.122328,-9.8184e-05,3.30406e-08,-1.91601e-12,6033.04,47.1925], Tmin=(100,'K'), Tmax=(1037.77,'K')), NASAPolynomial(coeffs=[25.6991,0.0395897,-1.47213e-05,2.63245e-09,-1.81647e-13,-1114.11,-95.1734], Tmin=(1037.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=[C][C](C)C(O)(CC)C(=C)O(26112)',
    structure = SMILES('[CH2][C]=C(C)C(O)(CC)C(=C)O'),
    E0 = (-64.4312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45192,0.13093,-0.000123957,6.01856e-08,-1.1597e-11,-7507.8,44.6059], Tmin=(100,'K'), Tmax=(1256.85,'K')), NASAPolynomial(coeffs=[25.9976,0.0403878,-1.58975e-05,2.86796e-09,-1.95914e-13,-14659.1,-99.1504], Tmin=(1256.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.4312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C(C)C(O)([CH]C)C(=C)O(26113)',
    structure = SMILES('C=[C]C(C)C(O)([CH]C)C(=C)O'),
    E0 = (9.726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,180,180,180,523.054,610.775,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.31615,0.123569,-9.79371e-05,2.65283e-08,2.40075e-12,1410.89,48.7761], Tmin=(100,'K'), Tmax=(974.711,'K')), NASAPolynomial(coeffs=[27.2311,0.0349815,-1.18811e-05,2.05402e-09,-1.41163e-13,-5900.93,-100.976], Tmin=(974.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C(C=C)C([O])(CC)C(=C)O(19804)',
    structure = SMILES('[CH2]C(C=C)C([O])(CC)C(=C)O'),
    E0 = (6.16714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,891.493,1315.34,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4819.15,'J/mol'), sigma=(8.10355,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=752.74 K, Pc=20.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10588,0.119884,-9.2049e-05,2.58179e-08,1.27997e-12,974.314,48.0474], Tmin=(100,'K'), Tmax=(991.066,'K')), NASAPolynomial(coeffs=[24.8373,0.0397338,-1.4018e-05,2.44045e-09,-1.66611e-13,-5770.46,-88.7809], Tmin=(991.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.16714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=CC(C)C([O])([CH]C)C(=C)O(20153)',
    structure = SMILES('C=CC(C)C([O])([CH]C)C(=C)O'),
    E0 = (0.986866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.9345,0.11304,-7.10988e-05,4.06188e-09,8.6232e-12,347.804,48.1475], Tmin=(100,'K'), Tmax=(996.605,'K')), NASAPolynomial(coeffs=[25.6259,0.0387884,-1.40757e-05,2.53058e-09,-1.77273e-13,-6951.51,-93.7823], Tmin=(996.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.986866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C([C]=C)C(O)(CC)C(=C)O(20154)',
    structure = SMILES('[CH2]C([C]=C)C(O)(CC)C(=C)O'),
    E0 = (14.9063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,411.055,724.254,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.46676,0.130155,-0.000117921,4.69382e-08,-4.32561e-12,2036.52,48.6023], Tmin=(100,'K'), Tmax=(957.092,'K')), NASAPolynomial(coeffs=[26.3948,0.0360113,-1.18734e-05,1.9759e-09,-1.31513e-13,-4700.86,-95.7079], Tmin=(957.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.9063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC(O)(C(=C)O)C(C)[C]=C(26114)',
    structure = SMILES('[CH2]CC(O)(C(=C)O)C(C)[C]=C'),
    E0 = (15.0703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,469.957,666.409,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15644,'amu*angstrom^2'), symmetry=1, barrier=(3.59687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15644,'amu*angstrom^2'), symmetry=1, barrier=(3.59687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15644,'amu*angstrom^2'), symmetry=1, barrier=(3.59687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15644,'amu*angstrom^2'), symmetry=1, barrier=(3.59687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15644,'amu*angstrom^2'), symmetry=1, barrier=(3.59687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15644,'amu*angstrom^2'), symmetry=1, barrier=(3.59687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15644,'amu*angstrom^2'), symmetry=1, barrier=(3.59687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15644,'amu*angstrom^2'), symmetry=1, barrier=(3.59687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.83704,0.133697,-0.000131072,6.55848e-08,-1.28314e-11,2072.82,49.9414], Tmin=(100,'K'), Tmax=(1253.47,'K')), NASAPolynomial(coeffs=[28.5429,0.0335605,-1.12426e-05,1.85318e-09,-1.20477e-13,-5794.04,-108.538], Tmin=(1253.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.0703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C)C(O)(CC)C(=C)[O](26115)',
    structure = SMILES('C=[C]C(C)C(O)(CC)C(=C)[O]'),
    E0 = (-52.3713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.27306,0.126265,-0.000118816,5.86105e-08,-1.14951e-11,-6063.09,47.629], Tmin=(100,'K'), Tmax=(1237.66,'K')), NASAPolynomial(coeffs=[24.3027,0.0403742,-1.47187e-05,2.53837e-09,-1.68836e-13,-12641.4,-86.2501], Tmin=(1237.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.3713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C(C)[C]=C(26116)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C(C)[C]=C'),
    E0 = (56.9201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1600,1663.96,2832.55,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.63281,0.132873,-0.000125169,5.58442e-08,-8.35221e-12,7096.11,47.8192], Tmin=(100,'K'), Tmax=(1012.68,'K')), NASAPolynomial(coeffs=[27.2175,0.0359216,-1.26032e-05,2.17335e-09,-1.46948e-13,-24.1565,-101.873], Tmin=(1012.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.9201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC([O])(C(=C)O)C(C)C=C(20159)',
    structure = SMILES('[CH2]CC([O])(C(=C)O)C(C)C=C'),
    E0 = (6.33117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,961.091,1246.59,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.13233,0.119442,-9.16307e-05,2.74485e-08,-2.34878e-13,995.628,48.1488], Tmin=(100,'K'), Tmax=(1035.85,'K')), NASAPolynomial(coeffs=[25.2102,0.0402001,-1.50282e-05,2.69832e-09,-1.86709e-13,-6082.2,-91.5478], Tmin=(1035.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.33117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=CC(C)C([O])(CC)C(=C)[O](20160)',
    structure = SMILES('C=CC(C)C([O])(CC)C(=C)[O]'),
    E0 = (-61.1104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.0554,0.117492,-9.72816e-05,4.18401e-08,-7.20398e-12,-7118.8,47.6006], Tmin=(100,'K'), Tmax=(1393.41,'K')), NASAPolynomial(coeffs=[24.4317,0.041457,-1.54313e-05,2.67984e-09,-1.78063e-13,-14500.3,-88.9716], Tmin=(1393.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.1104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)C(C)C=C(20161)',
    structure = SMILES('[CH]=C(O)C([O])(CC)C(C)C=C'),
    E0 = (48.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,454.305,681.32,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25134,0.122328,-9.81835e-05,3.30399e-08,-1.91574e-12,6033.04,47.1924], Tmin=(100,'K'), Tmax=(1037.76,'K')), NASAPolynomial(coeffs=[25.6991,0.0395898,-1.47213e-05,2.63246e-09,-1.81648e-13,-1114.08,-95.173], Tmin=(1037.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C(C)C(O)(CC)C(=C)O(26117)',
    structure = SMILES('[CH]=[C]C(C)C(O)(CC)C(=C)O'),
    E0 = (56.9201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1600,1663.96,2832.55,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153414,'amu*angstrom^2'), symmetry=1, barrier=(3.52729,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.63281,0.132873,-0.000125169,5.58442e-08,-8.35221e-12,7096.11,47.8192], Tmin=(100,'K'), Tmax=(1012.68,'K')), NASAPolynomial(coeffs=[27.2175,0.0359216,-1.26032e-05,2.17335e-09,-1.46948e-13,-24.1565,-101.873], Tmin=(1012.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.9201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C(C)C1(CC)OC[C]1O(23189)',
    structure = SMILES('C=[C]C(C)C1(CC)OC[C]1O'),
    E0 = (92.1821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,3150,900,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77085,0.11899,-9.91238e-05,3.70844e-08,-2.30393e-12,11302,41.5609], Tmin=(100,'K'), Tmax=(911.667,'K')), NASAPolynomial(coeffs=[20.0204,0.0456562,-1.51193e-05,2.4586e-09,-1.58837e-13,6403,-66.6315], Tmin=(911.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.1821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(Cds_S) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C1C[C](O)C([O])(CC)C1C(25812)',
    structure = SMILES('C=C1C[C](O)C([O])(CC)C1C'),
    E0 = (-43.5578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5167,0.1037,-5.2294e-05,-5.60425e-09,9.20097e-12,-5024.99,38.5515], Tmin=(100,'K'), Tmax=(1049.1,'K')), NASAPolynomial(coeffs=[22.4157,0.0461603,-1.8223e-05,3.3733e-09,-2.37138e-13,-11901.6,-86.8968], Tmin=(1049.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.5578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclopentane) + radical(CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C=C(C)C(O)(CC)C(=C)O(26118)',
    structure = SMILES('C=C=C(C)C(O)(CC)C(=C)O'),
    E0 = (-277.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.47376,0.129545,-0.000120355,5.71062e-08,-1.07333e-11,-33091.8,42.9859], Tmin=(100,'K'), Tmax=(1288.79,'K')), NASAPolynomial(coeffs=[26.6245,0.0392325,-1.52427e-05,2.73327e-09,-1.86054e-13,-40592.1,-104.779], Tmin=(1288.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-277.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = 'C=[C]CC([O])(CC)C(=C)O(26119)',
    structure = SMILES('C=[C]CC([O])(CC)C(=C)O'),
    E0 = (70.6787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1036.49,1179.57,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162038,'amu*angstrom^2'), symmetry=1, barrier=(3.72556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162038,'amu*angstrom^2'), symmetry=1, barrier=(3.72556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162038,'amu*angstrom^2'), symmetry=1, barrier=(3.72556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162038,'amu*angstrom^2'), symmetry=1, barrier=(3.72556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162038,'amu*angstrom^2'), symmetry=1, barrier=(3.72556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162038,'amu*angstrom^2'), symmetry=1, barrier=(3.72556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5578,0.111128,-0.000102613,4.89134e-08,-9.25811e-12,8710.09,42.7781], Tmin=(100,'K'), Tmax=(1279.66,'K')), NASAPolynomial(coeffs=[22.8455,0.0348461,-1.31952e-05,2.32871e-09,-1.56992e-13,2464.6,-80.9709], Tmin=(1279.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(70.6787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C)C(C)([O])C(=C)O(26120)',
    structure = SMILES('C=[C]C(C)C(C)([O])C(=C)O'),
    E0 = (62.7068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,180,180,180,463.016,673.86,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156422,'amu*angstrom^2'), symmetry=1, barrier=(3.59645,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156422,'amu*angstrom^2'), symmetry=1, barrier=(3.59645,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156422,'amu*angstrom^2'), symmetry=1, barrier=(3.59645,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156422,'amu*angstrom^2'), symmetry=1, barrier=(3.59645,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156422,'amu*angstrom^2'), symmetry=1, barrier=(3.59645,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156422,'amu*angstrom^2'), symmetry=1, barrier=(3.59645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47491,0.107561,-8.8712e-05,3.19145e-08,-2.67736e-12,7750.15,42.2314], Tmin=(100,'K'), Tmax=(1037.92,'K')), NASAPolynomial(coeffs=[22.9108,0.0342559,-1.26494e-05,2.24978e-09,-1.54665e-13,1574.51,-81.6873], Tmin=(1037.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.7068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C(=C)[CH]C(25632)',
    structure = SMILES('[CH2]C(=CC)C([O])(CC)C(=C)O'),
    E0 = (-73.1703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,262.04,866.396,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4889.28,'J/mol'), sigma=(8.15989,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=763.69 K, Pc=20.42 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.24133,0.122242,-0.00010271,4.37536e-08,-7.43059e-12,-8563.24,44.6026], Tmin=(100,'K'), Tmax=(1410.71,'K')), NASAPolynomial(coeffs=[26.3925,0.0410524,-1.63822e-05,2.95775e-09,-2.00988e-13,-16642.1,-103.392], Tmin=(1410.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.1703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(C)[CH]C([O])(CC)C(=C)O(26121)',
    structure = SMILES('C=C(C)[CH]C([O])(CC)C(=C)O'),
    E0 = (-136.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.01065,0.111231,-5.49639e-05,-1.71996e-08,1.6872e-11,-16130.9,45.4805], Tmin=(100,'K'), Tmax=(986.88,'K')), NASAPolynomial(coeffs=[27.5833,0.0376565,-1.36216e-05,2.48847e-09,-1.77735e-13,-24230.4,-108.344], Tmin=(986.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C=C1OC(CC)(C(=C)O)C1C(26058)',
    structure = SMILES('C=C1OC(CC)(C(=C)O)C1C'),
    E0 = (-306.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8767,0.0969757,1.14247e-05,-1.11918e-07,5.88713e-11,-36670.5,35.6738], Tmin=(100,'K'), Tmax=(911.619,'K')), NASAPolynomial(coeffs=[35.4603,0.0186955,-5.32207e-07,-2.35935e-10,1.4264e-14,-47032.7,-160.497], Tmin=(911.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-306.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(2methyleneoxetane)"""),
)

species(
    label = 'C=[C]C(C)C([C]=C)(CC)OO(26122)',
    structure = SMILES('C=[C]C(C)C([C]=C)(CC)OO'),
    E0 = (333.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,180,180,180,235.506,900.155,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152144,'amu*angstrom^2'), symmetry=1, barrier=(3.49809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152144,'amu*angstrom^2'), symmetry=1, barrier=(3.49809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152144,'amu*angstrom^2'), symmetry=1, barrier=(3.49809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152144,'amu*angstrom^2'), symmetry=1, barrier=(3.49809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152144,'amu*angstrom^2'), symmetry=1, barrier=(3.49809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152144,'amu*angstrom^2'), symmetry=1, barrier=(3.49809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152144,'amu*angstrom^2'), symmetry=1, barrier=(3.49809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152144,'amu*angstrom^2'), symmetry=1, barrier=(3.49809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.73812,0.125812,-0.000119196,6.10106e-08,-1.28104e-11,40345.3,45.3903], Tmin=(100,'K'), Tmax=(1134.32,'K')), NASAPolynomial(coeffs=[18.9861,0.0527321,-2.25568e-05,4.21352e-09,-2.92605e-13,35643.7,-57.2041], Tmin=(1134.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C)C([O])(CC)C(C)=O(26123)',
    structure = SMILES('C=[C]C(C)C([O])(CC)C(C)=O'),
    E0 = (42.5306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25901,0.121885,-0.000125224,7.55283e-08,-1.92938e-11,5299.33,44.0201], Tmin=(100,'K'), Tmax=(929.599,'K')), NASAPolynomial(coeffs=[12.8843,0.061029,-2.7029e-05,5.10937e-09,-3.56319e-13,2669.74,-23.1811], Tmin=(929.599,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.5306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CC(C)(C=O)OJ)"""),
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
    label = 'C=[C]C(C)[C](CC)C(=C)O(26124)',
    structure = SMILES('[CH2]C(O)=C(CC)C(C)[C]=C'),
    E0 = (139.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75648,0.121974,-0.000119808,6.30129e-08,-1.3339e-11,17007.7,39.7686], Tmin=(100,'K'), Tmax=(1141.07,'K')), NASAPolynomial(coeffs=[20.5016,0.0439484,-1.72388e-05,3.08673e-09,-2.09565e-13,11928.1,-70.5511], Tmin=(1141.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(Cds_S)"""),
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
    label = 'C=C(O)C([O])([CH]C)CC(6518)',
    structure = SMILES('C=C(O)C([O])([CH]C)CC'),
    E0 = (-71.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,524.731,610.588,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62983,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.840076,0.0926017,-5.93304e-05,5.06435e-09,6.37595e-12,-8384.46,39.5015], Tmin=(100,'K'), Tmax=(995.569,'K')), NASAPolynomial(coeffs=[21.324,0.0322481,-1.16347e-05,2.07897e-09,-1.44926e-13,-14219.8,-74.4714], Tmin=(995.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
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
    E0 = (38.9266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (86.6097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (101.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (175.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (190.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (38.9266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (88.6074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (178.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (38.9266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (159.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (221.508,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (153.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (114.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (115.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (190.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (83.2352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (98.5863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (98.7503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (140.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (101.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (74.7291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (81.6034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (81.2212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (89.9603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (176.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (163.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (82.0218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (117.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (490.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (482.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (133.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (183.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (47.2109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (440.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (243.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (382.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (329.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=C(O)C(=O)CC(4626)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['[CH2]C1(O)OC1(CC)C(C)[C]=C(26107)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(47.6831,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['[CH2]C1(O)C(=C)C(C)C1([O])CC(25664)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(62.76,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C=C(C)C([O])(CC)C(=C)O(26108)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.51e+07,'cm^3/(mol*s)'), n=1.64, Ea=(11.8407,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2579 used for Cds-CsCs_Ca;HJ
Exact match found for rate rule [Cds-CsCs_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#CC(C)C([O])(CC)C(=C)O(26109)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C(O)C(=O)CC(4626)', 'C=[C][CH]C(18176)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(36.3378,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 32.2 to 36.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H5(29)', 'C=[C]C(C)C(=O)C(=C)O(26110)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2COH(99)', 'C=[C]C(C)C(=O)CC(24949)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.16e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using template [CO-CsCs_O;CJ] for rate rule [CO-CsCs_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'CH3CHCCH2(18175)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00472174,'m^3/(mol*s)'), n=2.41, Ea=(78.2545,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 73.2 to 78.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH3(17)', 'C=C=CC([O])(CC)C(=C)O(26111)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(10800,'cm^3/(mol*s)'), n=2.41, Ea=(32.1331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 597 used for Cds-CsH_Ca;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Ca;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=C[C](C)C([O])(CC)C(=C)O(20149)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_Cs2]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC(C)C([O])(CC)C(=C)O(20158)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=[C][C](C)C(O)(CC)C(=C)O(26112)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]C(C)C(O)([CH]C)C(=C)O(26113)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=CC(C)C([O])([CH]C)C(=C)O(20153)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([C]=C)C(O)(CC)C(=C)O(20154)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC(O)(C(=C)O)C(C)[C]=C(26114)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=[C]C(C)C(O)(CC)C(=C)[O](26115)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(O)C(O)(CC)C(C)[C]=C(26116)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['[CH2]CC([O])(C(=C)O)C(C)C=C(20159)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(561575,'s^-1'), n=1.6076, Ea=(35.8025,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_CCC;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=CC(C)C([O])(CC)C(=C)[O](20160)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);Y_rad_out;XH_out] for rate rule [R5H_CCC(Cd);Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(O)C([O])(CC)C(C)C=C(20161)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 3.60555127546
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]C(C)C(O)(CC)C(=C)O(26117)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'C=[C][CH]C(18176)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=[C]C(C)C1(CC)OC[C]1O(23189)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=C1C[C](O)C([O])(CC)C1C(25812)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.6e+11,'s^-1'), n=0.27, Ea=(43.0952,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 8 used for R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cddouble
Exact match found for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=C=C(C)C(O)(CC)C(=C)O(26118)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(S)(23)', 'C=[C]CC([O])(CC)C(=C)O(26119)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(S)(23)', 'C=[C]C(C)C(C)([O])C(=C)O(26120)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=C(C)[CH]C([O])(CC)C(=C)O(26121)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.95888e+10,'s^-1'), n=0.7315, Ea=(144.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HC)CJ;CJ;CH3] + [cCs(-HC)CJ;CdsJ;C] for rate rule [cCs(-HC)CJ;CdsJ;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=C1OC(CC)(C(=C)O)C1C(26058)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C]C(C)C([C]=C)(CC)OO(26122)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=[C]C(C)C([O])(CC)C(C)=O(26123)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(4)', 'C=[C]C(C)[C](CC)C(=C)O(26124)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction37',
    reactants = ['H2CC(41)', 'C=C(O)C([O])([CH]C)CC(6518)'],
    products = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4294',
    isomers = [
        'C=[C]C(C)C([O])(CC)C(=C)O(20152)',
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
    label = '4294',
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

