species(
    label = 'O=C[CH]OO[CH]O(1318)',
    structure = SMILES('O=C[CH]OO[CH]O'),
    E0 = (-116.242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,270.967],'cm^-1')),
        HinderedRotor(inertia=(0.283772,'amu*angstrom^2'), symmetry=1, barrier=(14.7415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.65145,'amu*angstrom^2'), symmetry=1, barrier=(33.8502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00230112,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.908989,'amu*angstrom^2'), symmetry=1, barrier=(47.2209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.912758,'amu*angstrom^2'), symmetry=1, barrier=(47.2172,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.802195,0.0751136,-0.000104076,7.53955e-08,-2.18262e-11,-13869.8,25.9926], Tmin=(100,'K'), Tmax=(843.879,'K')), NASAPolynomial(coeffs=[11.791,0.023028,-1.14962e-05,2.25873e-09,-1.59894e-13,-15724.5,-25.1568], Tmin=(843.879,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-116.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-OsOsHH) + group(Cds-OdCsH) + radical(OCJO) + radical(OCJC=O)"""),
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
    label = 'O=C[CH]OOC=O(7969)',
    structure = SMILES('O=C[CH]OOC=O'),
    E0 = (-355.612,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(1.16884,'amu*angstrom^2'), symmetry=1, barrier=(44.3568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1687,'amu*angstrom^2'), symmetry=1, barrier=(44.3859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15573,'amu*angstrom^2'), symmetry=1, barrier=(44.3551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17829,'amu*angstrom^2'), symmetry=1, barrier=(44.3463,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47849,0.043211,-9.12559e-06,-2.25736e-08,1.19202e-11,-42668.8,24.3935], Tmin=(100,'K'), Tmax=(1028.76,'K')), NASAPolynomial(coeffs=[15.7984,0.0145457,-6.71614e-06,1.38848e-09,-1.05324e-13,-47044.6,-52.0451], Tmin=(1028.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-355.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(OCJC=O)"""),
)

species(
    label = 'O=C=COO[CH]O(7970)',
    structure = SMILES('O=C=COO[CH]O'),
    E0 = (-72.5377,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2120,512.5,787.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.677652,'amu*angstrom^2'), symmetry=1, barrier=(15.5805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.678156,'amu*angstrom^2'), symmetry=1, barrier=(15.5921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26781,'amu*angstrom^2'), symmetry=1, barrier=(29.1496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.678711,'amu*angstrom^2'), symmetry=1, barrier=(15.6049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.754039,0.0733151,-0.000104697,7.3117e-08,-1.97621e-11,-8609.05,25.5373], Tmin=(100,'K'), Tmax=(914.425,'K')), NASAPolynomial(coeffs=[14.5174,0.0131092,-5.93682e-06,1.11464e-09,-7.68986e-14,-11126.1,-39.6317], Tmin=(914.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-72.5377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cds-(Cdd-O2d)OsH) + radical(OCJO)"""),
)

species(
    label = '[O]COO[CH]C=O(7971)',
    structure = SMILES('[O]COO[CH]C=O'),
    E0 = (-77.8047,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,340.912,340.916],'cm^-1')),
        HinderedRotor(inertia=(0.537753,'amu*angstrom^2'), symmetry=1, barrier=(44.3502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.53772,'amu*angstrom^2'), symmetry=1, barrier=(44.3504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0014504,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.537773,'amu*angstrom^2'), symmetry=1, barrier=(44.3498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44355,0.0555879,-4.70652e-05,1.9451e-08,-3.25623e-12,-9265.52,24.363], Tmin=(100,'K'), Tmax=(1393.08,'K')), NASAPolynomial(coeffs=[13.127,0.0220407,-1.09432e-05,2.16454e-09,-1.54039e-13,-12520.7,-35.8762], Tmin=(1393.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.8047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-OsOsHH) + group(Cds-OdCsH) + radical(OCOJ) + radical(OCJC=O)"""),
)

species(
    label = 'O=[C]COO[CH]O(7972)',
    structure = SMILES('O=[C]COO[CH]O'),
    E0 = (-100.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,3615,1277.5,1000,1855,455,950,2750,2850,1437.5,1250,1305,750,350,805.493],'cm^-1')),
        HinderedRotor(inertia=(0.548407,'amu*angstrom^2'), symmetry=1, barrier=(12.609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0583633,'amu*angstrom^2'), symmetry=1, barrier=(12.6087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.548404,'amu*angstrom^2'), symmetry=1, barrier=(12.6089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55329,'amu*angstrom^2'), symmetry=1, barrier=(35.7132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55328,'amu*angstrom^2'), symmetry=1, barrier=(35.713,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.908452,0.0721047,-9.57563e-05,6.5483e-08,-1.7848e-11,-11926.2,26.5337], Tmin=(100,'K'), Tmax=(895.242,'K')), NASAPolynomial(coeffs=[12.2765,0.0213116,-1.06517e-05,2.10793e-09,-1.50321e-13,-13961.7,-27.0528], Tmin=(895.242,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-OsOsHH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(OCJO)"""),
)

species(
    label = 'O=[C][CH]OOCO(7973)',
    structure = SMILES('O=[C][CH]OOCO'),
    E0 = (-150.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.610127,0.0627557,-4.66228e-05,1.85289e-09,7.23336e-12,-17972,27.5574], Tmin=(100,'K'), Tmax=(963.465,'K')), NASAPolynomial(coeffs=[20.7213,0.00663028,-1.85392e-06,3.60278e-10,-3.01304e-14,-23117.7,-75.3116], Tmin=(963.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-150.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-OsOsHH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[O][CH]OOCC=O(7974)',
    structure = SMILES('[O][CH]OOCC=O'),
    E0 = (-27.3294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,180,1822.69],'cm^-1')),
        HinderedRotor(inertia=(0.00404478,'amu*angstrom^2'), symmetry=1, barrier=(9.52792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78158,'amu*angstrom^2'), symmetry=1, barrier=(40.9621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.7813,'amu*angstrom^2'), symmetry=1, barrier=(40.9556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78078,'amu*angstrom^2'), symmetry=1, barrier=(40.9436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20006,0.0716758,-0.00012187,1.20385e-07,-4.65046e-11,-3195.99,25.2645], Tmin=(100,'K'), Tmax=(803.057,'K')), NASAPolynomial(coeffs=[3.76603,0.038155,-2.05189e-05,4.08791e-09,-2.8832e-13,-2939.36,17.6118], Tmin=(803.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.3294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-OsOsHH) + group(Cds-OdCsH) + radical(OCJO) + radical(OCOJ)"""),
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
    label = 'O[CH]OOC1[CH]O1(7975)',
    structure = SMILES('O[CH]OOC1[CH]O1'),
    E0 = (-15.0833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.520063,0.073072,-8.07503e-05,3.40514e-08,-9.26299e-13,-1685.2,25.4398], Tmin=(100,'K'), Tmax=(819.054,'K')), NASAPolynomial(coeffs=[18.0554,0.00832303,-4.23875e-07,-1.93601e-10,2.25223e-14,-5258.3,-59.9349], Tmin=(819.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-15.0833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Ethylene_oxide) + radical(OCJO) + radical(CCsJO)"""),
)

species(
    label = 'OC1O[CH][CH]OO1(7966)',
    structure = SMILES('OC1O[CH][CH]OO1'),
    E0 = (-104.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.330749,0.0617989,-6.48273e-05,3.46096e-08,-6.71715e-12,-12416,24.1324], Tmin=(100,'K'), Tmax=(1565.3,'K')), NASAPolynomial(coeffs=[13.4943,0.0105958,1.07164e-06,-6.25764e-10,5.53795e-14,-14385.2,-38.3989], Tmin=(1565.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-104.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-OsOsOsH) + ring(124trioxane) + radical(CCsJOCs) + radical(CCsJOOC)"""),
)

species(
    label = 'O=C=COOCO(7976)',
    structure = SMILES('O=C=COOCO'),
    E0 = (-261.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.660599,0.0587245,-3.06832e-05,-2.16895e-08,1.79662e-11,-31317,26.077], Tmin=(100,'K'), Tmax=(923.247,'K')), NASAPolynomial(coeffs=[23.0808,-1.31684e-05,2.36115e-06,-5.01705e-10,3.04301e-14,-37093.4,-89.1602], Tmin=(923.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-261.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cds-(Cdd-O2d)OsH)"""),
)

species(
    label = 'O=CCOOC=O(7977)',
    structure = SMILES('O=CCOOC=O'),
    E0 = (-494.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49997,0.0402662,7.09241e-06,-4.16151e-08,1.89915e-11,-59323.1,24.5585], Tmin=(100,'K'), Tmax=(1006.64,'K')), NASAPolynomial(coeffs=[16.4731,0.0151237,-6.63545e-06,1.38007e-09,-1.06435e-13,-64078.3,-56.4233], Tmin=(1006.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-494.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdOsH)"""),
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
    label = 'CHCHO(47)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (245.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06864,0.0187233,-1.21319e-05,-3.33727e-10,2.32882e-12,29739.4,14.7866], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.96288,0.00799899,-4.30606e-06,1.11076e-09,-1.11415e-13,28725.6,-5.17392], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(245.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[O]O[CH]O(4538)',
    structure = SMILES('[O]O[CH]O'),
    E0 = (9.28844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,492.5,1135,1000],'cm^-1')),
        HinderedRotor(inertia=(0.156407,'amu*angstrom^2'), symmetry=1, barrier=(3.59612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156804,'amu*angstrom^2'), symmetry=1, barrier=(3.60523,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58857,0.0452524,-0.00011197,1.22556e-07,-4.65581e-11,1154.51,14.2974], Tmin=(100,'K'), Tmax=(891.504,'K')), NASAPolynomial(coeffs=[-0.0459,0.0212294,-1.12418e-05,2.13324e-09,-1.41938e-13,3048.61,34.6933], Tmin=(891.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.28844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = 'HCOH(T)(285)',
    structure = SMILES('[CH]O'),
    E0 = (205.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,403.876,3308.82],'cm^-1')),
        HinderedRotor(inertia=(0.0103144,'amu*angstrom^2'), symmetry=1, barrier=(22.7121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75938,0.0029613,8.90411e-06,-1.35016e-08,5.39816e-12,24775.6,6.76286], Tmin=(100,'K'), Tmax=(940.429,'K')), NASAPolynomial(coeffs=[5.09112,0.00321239,-9.31686e-07,1.59615e-10,-1.15729e-14,24263.5,-0.971], Tmin=(940.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HCOH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]O[CH]C=O(3604)',
    structure = SMILES('[O]O[CH]C=O'),
    E0 = (38.8345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.304071,'amu*angstrom^2'), symmetry=1, barrier=(6.9912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17371,'amu*angstrom^2'), symmetry=1, barrier=(56.7408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0355,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48924,0.0359331,-5.02567e-05,4.15061e-08,-1.40856e-11,4722.53,16.7287], Tmin=(100,'K'), Tmax=(782.809,'K')), NASAPolynomial(coeffs=[5.94539,0.0159797,-7.6283e-06,1.46013e-09,-1.01349e-13,4251.69,1.34978], Tmin=(782.809,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.8345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(ROOJ)"""),
)

species(
    label = '[O][CH]C1OOC1O(7978)',
    structure = SMILES('[O][CH]C1OOC1O'),
    E0 = (-6.18977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26898,0.0686498,-0.000113436,1.0666e-07,-3.91421e-11,-654.43,22.3588], Tmin=(100,'K'), Tmax=(826.759,'K')), NASAPolynomial(coeffs=[5.19125,0.0324224,-1.64101e-05,3.18443e-09,-2.20941e-13,-713.412,7.74792], Tmin=(826.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.18977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(12dioxetane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'O[C]=COO[CH]O(7979)',
    structure = SMILES('O[C]=COO[CH]O'),
    E0 = (-14.2693,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,1685,370,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.835041,'amu*angstrom^2'), symmetry=1, barrier=(19.1992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834589,'amu*angstrom^2'), symmetry=1, barrier=(19.1889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834932,'amu*angstrom^2'), symmetry=1, barrier=(19.1967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835152,'amu*angstrom^2'), symmetry=1, barrier=(19.2018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835286,'amu*angstrom^2'), symmetry=1, barrier=(19.2049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.153909,0.0797073,-9.99704e-05,5.49871e-08,-1.01712e-11,-1573.13,30.0108], Tmin=(100,'K'), Tmax=(907.565,'K')), NASAPolynomial(coeffs=[20.5527,0.00467329,-5.35442e-07,6.77714e-13,1.97587e-15,-5888.24,-69.7973], Tmin=(907.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.2693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(OCJO) + radical(C=CJO)"""),
)

species(
    label = 'O[CH]OO[C]=CO(7980)',
    structure = SMILES('O[CH]OO[C]=CO'),
    E0 = (-14.2693,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,1685,370,3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.836087,'amu*angstrom^2'), symmetry=1, barrier=(19.2233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835877,'amu*angstrom^2'), symmetry=1, barrier=(19.2185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.83495,'amu*angstrom^2'), symmetry=1, barrier=(19.1971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834152,'amu*angstrom^2'), symmetry=1, barrier=(19.1788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834062,'amu*angstrom^2'), symmetry=1, barrier=(19.1767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.153909,0.0797073,-9.99704e-05,5.49871e-08,-1.01712e-11,-1573.13,30.0108], Tmin=(100,'K'), Tmax=(907.565,'K')), NASAPolynomial(coeffs=[20.5527,0.00467329,-5.35442e-07,6.77714e-13,1.97587e-15,-5888.24,-69.7973], Tmin=(907.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.2693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(OCJO)"""),
)

species(
    label = '[O]C=[C]OOCO(7981)',
    structure = SMILES('[O]C=[C]OOCO'),
    E0 = (-61.7668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,1685,370,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01487,'amu*angstrom^2'), symmetry=1, barrier=(23.3338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01484,'amu*angstrom^2'), symmetry=1, barrier=(23.3332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01482,'amu*angstrom^2'), symmetry=1, barrier=(23.3327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01871,'amu*angstrom^2'), symmetry=1, barrier=(23.4222,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61541,0.0563764,-1.60062e-05,-4.2575e-08,2.68427e-11,-7290.13,29.9744], Tmin=(100,'K'), Tmax=(917.28,'K')), NASAPolynomial(coeffs=[25.3655,-0.00389843,4.63355e-06,-9.40181e-10,5.96118e-14,-13835.5,-98.2212], Tmin=(917.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.7668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[O][CH]OOC=CO(7982)',
    structure = SMILES('[O][CH]OOC=CO'),
    E0 = (-26.6156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.10571,'amu*angstrom^2'), symmetry=1, barrier=(25.4225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10589,'amu*angstrom^2'), symmetry=1, barrier=(25.4265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10575,'amu*angstrom^2'), symmetry=1, barrier=(25.4234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10621,'amu*angstrom^2'), symmetry=1, barrier=(25.434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676034,0.072912,-8.93233e-05,5.39062e-08,-1.26999e-11,-3081.14,25.4879], Tmin=(100,'K'), Tmax=(1041.63,'K')), NASAPolynomial(coeffs=[15.508,0.0159544,-7.30034e-06,1.40902e-09,-9.99211e-14,-6170.98,-46.6726], Tmin=(1041.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.6156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(OCOJ) + radical(OCJO)"""),
)

species(
    label = 'O=COOC=CO(7983)',
    structure = SMILES('O=COOC=CO'),
    E0 = (-493.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.906536,0.0425544,3.45574e-05,-9.87523e-08,4.72035e-11,-59205.2,25.0181], Tmin=(100,'K'), Tmax=(925.543,'K')), NASAPolynomial(coeffs=[27.2389,-0.00542298,5.63107e-06,-1.07416e-09,6.33398e-14,-66899,-115.214], Tmin=(925.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-493.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsH)"""),
)

species(
    label = 'C1=COO1(3609)',
    structure = SMILES('C1=COO1'),
    E0 = (143.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31094,-0.00195743,7.27275e-05,-1.01582e-07,4.10226e-11,17349.4,10.0685], Tmin=(100,'K'), Tmax=(927.669,'K')), NASAPolynomial(coeffs=[13.7806,-0.00309387,3.40676e-06,-6.26613e-10,3.47414e-14,13513.3,-49.8617], Tmin=(927.669,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclobutene)"""),
)

species(
    label = 'OC1OC=COO1(7961)',
    structure = SMILES('OC1OC=COO1'),
    E0 = (-365.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66914,0.0338819,2.69995e-05,-7.35317e-08,3.503e-11,-43902.1,19.0065], Tmin=(100,'K'), Tmax=(913.087,'K')), NASAPolynomial(coeffs=[19.0592,0.00338907,2.03653e-06,-5.05497e-10,3.15999e-14,-48982.4,-73.7385], Tmin=(913.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(124trioxene)"""),
)

species(
    label = '[O]C(C=O)O[CH]O(1316)',
    structure = SMILES('[O]C(C=O)O[CH]O'),
    E0 = (-276.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,363.765,363.86,363.977,4000],'cm^-1')),
        HinderedRotor(inertia=(0.101988,'amu*angstrom^2'), symmetry=1, barrier=(9.5885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101978,'amu*angstrom^2'), symmetry=1, barrier=(9.58886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00626724,'amu*angstrom^2'), symmetry=1, barrier=(22.6483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240708,'amu*angstrom^2'), symmetry=1, barrier=(22.6444,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4424.28,'J/mol'), sigma=(6.96605,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.06 K, Pc=29.7 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998809,0.0706265,-0.000106738,8.61771e-08,-2.76818e-11,-33131,28.4801], Tmin=(100,'K'), Tmax=(794.665,'K')), NASAPolynomial(coeffs=[10.297,0.0217016,-1.03825e-05,1.98122e-09,-1.36815e-13,-34541.7,-13.8195], Tmin=(794.665,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-OsOsHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(OCJO)"""),
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
    label = '[CH]=COO[CH]O(7817)',
    structure = SMILES('[CH]=COO[CH]O'),
    E0 = (201.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(0.854399,'amu*angstrom^2'), symmetry=1, barrier=(19.6443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.854351,'amu*angstrom^2'), symmetry=1, barrier=(19.6432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.854325,'amu*angstrom^2'), symmetry=1, barrier=(19.6426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.854364,'amu*angstrom^2'), symmetry=1, barrier=(19.6435,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.983753,0.064797,-8.19779e-05,5.04179e-08,-1.19582e-11,24390,24.4484], Tmin=(100,'K'), Tmax=(1041.89,'K')), NASAPolynomial(coeffs=[15.009,0.010952,-4.45796e-06,8.16091e-10,-5.63824e-14,21467.4,-43.7908], Tmin=(1041.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCJO) + radical(Cds_P)"""),
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
    E0 = (-116.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-80.3204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-107.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (164.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (36.1673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (52.4373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-89.3603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (33.8465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-13.4679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (67.2128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-89.8833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-107.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-77.0175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-108.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (255.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (244.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-6.18977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-113.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (148.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (136.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-17.4583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (49.7732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-91.2692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (150.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-108.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (197.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (444.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['HOCHO(59)', 'OCHCHO(48)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['[O]C1[CH]OOC1O(7968)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(13013.2,'s^-1'), n=1.81618, Ea=(35.9221,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_csHNd] for rate rule [R6;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 31.4 to 35.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'O=C[CH]OOC=O(7969)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Od_CO-NdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'O=C=COO[CH]O(7970)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.3e-15,'cm^3/(molecule*s)'), n=1.43, Ea=(25.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ck_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]COO[CH]C=O(7971)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.3e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using template [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=[C]COO[CH]O(7972)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.67823e+08,'s^-1'), n=1.48018, Ea=(152.494,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out] for rate rule [R2H_S;CO_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['O=[C][CH]OOCO(7973)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_1H;XH_out] for rate rule [R5HJ_3;C_rad_out_H/NonDeO;CO_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][CH]OOCC=O(7974)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(60863,'s^-1'), n=1.86284, Ea=(61.1759,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][CH]O(171)', '[O]C=C[O](2537)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['O[CH]OOC1[CH]O1(7975)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra_H;radadd_intra_csHO]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['OC1O[CH][CH]OO1(7966)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(487000,'s^-1'), n=1.17, Ea=(26.3592,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHNd] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_csHO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['O=C=COOCO(7976)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['O=CCOOC=O(7977)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['O=CC1OOC1O(1319)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] + [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CHCHO(47)', '[O]O[CH]O(4538)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['HCOH(T)(285)', '[O]O[CH]C=O(3604)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['[O][CH]C1OOC1O(7978)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.73e+06,'s^-1'), n=1.31, Ea=(110.053,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 108.4 to 110.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['HOCHO(59)', '[O]C=C[O](2537)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.6e+11,'cm^3/(mol*s)'), n=0, Ea=(294.2,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_R;O_rad/OneDe] for rate rule [Od_CO-NdH;O_rad/OneDe]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O[C]=COO[CH]O(7979)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O[CH]OO[C]=CO(7980)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C=[C]OOCO(7981)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O][CH]OOC=CO(7982)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(970995,'s^-1'), n=0.905106, Ea=(76.3888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;O_rad_out;XH_out] for rate rule [R7HJ_1;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['O=COOC=CO(7983)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['[O][CH]O(171)', 'C1=COO1(3609)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(266.548,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO_intra] for rate rule [R3OO_SD;Y_rad_intra;OO_intra]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['OC1OC=COO1(7961)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_1H;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_H/NonDeO;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=C[CH]OO[CH]O(1318)'],
    products = ['[O]C(C=O)O[CH]O(1316)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_R]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O(4)', '[CH]=COO[CH]O(7817)'],
    products = ['O=C[CH]OO[CH]O(1318)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '904',
    isomers = [
        'O=C[CH]OO[CH]O(1318)',
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
    label = '904',
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

