species(
    label = '[CH2]C(O)(O[C]=O)C(O)=CC(9923)',
    structure = SMILES('[CH2]C(O)(O[C]=O)C(O)=CC'),
    E0 = (-396.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3010,987.5,1337.5,450,1655,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,251.719,1470.1,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.146812,'amu*angstrom^2'), symmetry=1, barrier=(3.37549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146812,'amu*angstrom^2'), symmetry=1, barrier=(3.37549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146812,'amu*angstrom^2'), symmetry=1, barrier=(3.37549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146812,'amu*angstrom^2'), symmetry=1, barrier=(3.37549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146812,'amu*angstrom^2'), symmetry=1, barrier=(3.37549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146812,'amu*angstrom^2'), symmetry=1, barrier=(3.37549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146812,'amu*angstrom^2'), symmetry=1, barrier=(3.37549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.75858,0.127704,-0.00014879,8.09369e-08,-1.64559e-11,-47482.7,39.8047], Tmin=(100,'K'), Tmax=(1335.93,'K')), NASAPolynomial(coeffs=[33.3255,0.00738682,8.70932e-08,-2.34831e-10,2.13629e-14,-56028.5,-140.632], Tmin=(1335.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-396.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(CJC(C)OC) + radical((O)CJOC)"""),
)

species(
    label = 'CO2(13)',
    structure = SMILES('O=C=O'),
    E0 = (-403.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([459.166,1086.67,1086.68,2300.05],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2779,0.00275783,7.12787e-06,-1.07855e-08,4.14228e-12,-48475.6,5.97856], Tmin=(100,'K'), Tmax=(988.185,'K')), NASAPolynomial(coeffs=[4.55071,0.00290728,-1.14643e-06,2.25798e-10,-1.69526e-14,-48986,-1.45662], Tmin=(988.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.131,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: BurkeH2O2"""),
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
    label = 'C[CH]C1(O)CC1(O)O[C]=O(31633)',
    structure = SMILES('C[CH]C1(O)CC1(O)O[C]=O'),
    E0 = (-335.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17072,0.0878829,-2.34309e-05,-6.00825e-08,3.71163e-11,-40114.3,36.4121], Tmin=(100,'K'), Tmax=(926.303,'K')), NASAPolynomial(coeffs=[33.3572,0.00553958,1.80941e-06,-4.46763e-10,2.34189e-14,-49374.9,-142.981], Tmin=(926.303,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdOsH) + ring(Cyclopropane) + radical(CCJCO) + radical((O)CJOC)"""),
)

species(
    label = '[CH2]C1(O)OC(=O)C1(O)[CH]C(31634)',
    structure = SMILES('[CH2]C1(O)OC(=O)C1(O)[CH]C'),
    E0 = (-362.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.08578,0.0844707,-7.04129e-06,-8.75666e-08,5.1455e-11,-43355,36.6425], Tmin=(100,'K'), Tmax=(889.647,'K')), NASAPolynomial(coeffs=[35.0286,-0.00117322,7.98444e-06,-1.87758e-09,1.31899e-13,-52817.4,-150.432], Tmin=(889.647,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-362.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(CCJCO) + radical(CJC(C)OC)"""),
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
    label = 'C=C(O)O[C]=O(5061)',
    structure = SMILES('C=C(O)O[C]=O'),
    E0 = (-254.888,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1855,455,950,3615,1277.5,1000,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06293,'amu*angstrom^2'), symmetry=1, barrier=(24.4389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06302,'amu*angstrom^2'), symmetry=1, barrier=(24.4408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06295,'amu*angstrom^2'), symmetry=1, barrier=(24.4394,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29067,0.0594252,-7.71778e-05,4.88549e-08,-1.1969e-11,-30558.1,17.2691], Tmin=(100,'K'), Tmax=(1007.08,'K')), NASAPolynomial(coeffs=[13.5381,0.0107798,-4.72263e-06,8.90961e-10,-6.23478e-14,-33024.9,-41.9042], Tmin=(1007.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-254.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = 'CC=[C]O(31159)',
    structure = SMILES('CC=[C]O'),
    E0 = (72.7982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.344811,'amu*angstrom^2'), symmetry=1, barrier=(7.92788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346062,'amu*angstrom^2'), symmetry=1, barrier=(7.95665,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.692,0.0334424,-4.69125e-05,4.57787e-08,-1.77014e-11,8798.16,15.2755], Tmin=(100,'K'), Tmax=(837.57,'K')), NASAPolynomial(coeffs=[2.44779,0.0239997,-1.10022e-05,2.07309e-09,-1.42159e-13,9211.19,18.6318], Tmin=(837.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.7982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=CJO)"""),
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
    label = 'C=C(O[C]=O)C(O)=CC(31635)',
    structure = SMILES('C=C(O[C]=O)C(O)=CC'),
    E0 = (-300.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,1855,455,950,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.782532,0.103128,-0.000122814,7.3154e-08,-1.70087e-11,-36027.6,27.992], Tmin=(100,'K'), Tmax=(1056.62,'K')), NASAPolynomial(coeffs=[20.2646,0.0234516,-9.70494e-06,1.78987e-09,-1.23918e-13,-40475.4,-74.7076], Tmin=(1056.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = 'CC=C(O)C(C)([O])O[C]=O(31636)',
    structure = SMILES('CC=C(O)C(C)([O])O[C]=O'),
    E0 = (-374.43,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,303.869,828.562,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152416,'amu*angstrom^2'), symmetry=1, barrier=(3.50435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152416,'amu*angstrom^2'), symmetry=1, barrier=(3.50435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152416,'amu*angstrom^2'), symmetry=1, barrier=(3.50435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152416,'amu*angstrom^2'), symmetry=1, barrier=(3.50435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152416,'amu*angstrom^2'), symmetry=1, barrier=(3.50435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152416,'amu*angstrom^2'), symmetry=1, barrier=(3.50435,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.883,0.113858,-0.00012529,6.7085e-08,-1.37061e-11,-44808.5,37.7907], Tmin=(100,'K'), Tmax=(1269.12,'K')), NASAPolynomial(coeffs=[27.6369,0.0161334,-4.25017e-06,5.94625e-10,-3.54631e-14,-51924.1,-110.175], Tmin=(1269.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-374.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(C=CC(C)(O)OJ) + radical((O)CJOC)"""),
)

species(
    label = '[CH2]C([O])(OC=O)C(O)=CC(31637)',
    structure = SMILES('[CH2]C([O])(OC=O)C(O)=CC'),
    E0 = (-360.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2782.5,750,1395,475,1775,1000,180,180,180,375.33,758.876,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154215,'amu*angstrom^2'), symmetry=1, barrier=(3.5457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154215,'amu*angstrom^2'), symmetry=1, barrier=(3.5457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154215,'amu*angstrom^2'), symmetry=1, barrier=(3.5457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154215,'amu*angstrom^2'), symmetry=1, barrier=(3.5457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154215,'amu*angstrom^2'), symmetry=1, barrier=(3.5457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154215,'amu*angstrom^2'), symmetry=1, barrier=(3.5457,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.64148,0.11438,-0.000125276,6.39717e-08,-1.15716e-11,-43130.6,38.0465], Tmin=(100,'K'), Tmax=(985.383,'K')), NASAPolynomial(coeffs=[26.0263,0.0191424,-6.29362e-06,1.05877e-09,-7.1611e-14,-49412.3,-99.2334], Tmin=(985.383,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-360.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(CJC(C)OC) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = 'C[CH]C(=O)C(C)(O)O[C]=O(9911)',
    structure = SMILES('CC=C([O])C(C)(O)O[C]=O'),
    E0 = (-469.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.19957,0.11266,-0.000117932,5.88747e-08,-1.10608e-11,-56245,39.5595], Tmin=(100,'K'), Tmax=(1453.38,'K')), NASAPolynomial(coeffs=[30.0031,0.0122589,-2.15914e-06,1.96252e-10,-8.66808e-15,-64362.1,-123.562], Tmin=(1453.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-469.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical((O)CJOC) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[C]=C(O)C(C)(O)O[C]=O(31638)',
    structure = SMILES('C[C]=C(O)C(C)(O)O[C]=O'),
    E0 = (-369.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,3580,3650,1210,1345,900,1100,1855,455,950,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.54394,0.120672,-0.000133612,6.94822e-08,-1.35229e-11,-44201.3,39.7622], Tmin=(100,'K'), Tmax=(1403.85,'K')), NASAPolynomial(coeffs=[32.2254,0.00892915,-6.73663e-07,-7.94412e-11,1.0013e-14,-52714.6,-135.326], Tmin=(1403.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-369.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(Cds_S) + radical((O)CJOC)"""),
)

species(
    label = '[CH2]C(O)(OC=O)C(=O)[CH]C(9912)',
    structure = SMILES('[CH2]C(O)(OC=O)C([O])=CC'),
    E0 = (-455.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.11555,0.11508,-0.000124675,6.45716e-08,-1.26637e-11,-54560.5,40.3766], Tmin=(100,'K'), Tmax=(1352.94,'K')), NASAPolynomial(coeffs=[29.7984,0.0130056,-2.94646e-06,3.71569e-10,-2.13293e-14,-62489.4,-120.625], Tmin=(1352.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-455.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(C=C(C)OJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C(O)(OC=O)C(O)=[C]C(31639)',
    structure = SMILES('[CH2]C(O)(OC=O)C(O)=[C]C'),
    E0 = (-355.573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2782.5,750,1395,475,1775,1000,180,180,180,214.767,1511.85,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.148018,'amu*angstrom^2'), symmetry=1, barrier=(3.40323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148018,'amu*angstrom^2'), symmetry=1, barrier=(3.40323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148018,'amu*angstrom^2'), symmetry=1, barrier=(3.40323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148018,'amu*angstrom^2'), symmetry=1, barrier=(3.40323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148018,'amu*angstrom^2'), symmetry=1, barrier=(3.40323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148018,'amu*angstrom^2'), symmetry=1, barrier=(3.40323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148018,'amu*angstrom^2'), symmetry=1, barrier=(3.40323,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.48086,0.12331,-0.000140996,7.58606e-08,-1.53612e-11,-42515.8,40.6566], Tmin=(100,'K'), Tmax=(1315.19,'K')), NASAPolynomial(coeffs=[32.0273,0.00965503,-1.44607e-06,9.19615e-11,-2.30426e-15,-50840.1,-132.419], Tmin=(1315.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-355.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(Cds_S) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C=C(O)C(C)(O)O[C]=O(31640)',
    structure = SMILES('[CH2]C=C(O)C(C)(O)O[C]=O'),
    E0 = (-455.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.29938,0.122497,-0.000129874,6.32543e-08,-1.13085e-11,-54545.9,42.3065], Tmin=(100,'K'), Tmax=(1607.23,'K')), NASAPolynomial(coeffs=[34.6322,0.00455109,2.17526e-06,-6.32492e-10,4.65058e-14,-63698,-149.23], Tmin=(1607.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-455.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(Allyl_P) + radical((O)CJOC)"""),
)

species(
    label = '[CH2]C=C(O)C([CH2])(O)OC=O(31641)',
    structure = SMILES('[CH2]C=C(O)C([CH2])(O)OC=O'),
    E0 = (-441.916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.16805,0.124434,-0.000135233,6.75159e-08,-1.24287e-11,-52863.6,42.9481], Tmin=(100,'K'), Tmax=(1538,'K')), NASAPolynomial(coeffs=[35.1726,0.00427945,1.88592e-06,-5.60912e-10,4.16489e-14,-62239.9,-150.669], Tmin=(1538,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-441.916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(Allyl_P) + radical(CJC(C)OC)"""),
)

species(
    label = '[O][C]=O(669)',
    structure = SMILES('[O][C]=O'),
    E0 = (31.9507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90478,-0.000175995,8.15126e-06,-1.13656e-08,4.4768e-12,3848.25,8.04855], Tmin=(100,'K'), Tmax=(975.388,'K')), NASAPolynomial(coeffs=[5.59398,-0.00122084,7.11747e-07,-9.7712e-11,3.97995e-15,3238.91,-1.49318], Tmin=(975.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.9507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(OJC=O) + radical((O)CJOH)"""),
)

species(
    label = 'CC1CC(O)(O[C]=O)[C]1O(31642)',
    structure = SMILES('CC1CC(O)(O[C]=O)[C]1O'),
    E0 = (-354.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.937598,0.0852109,-2.47385e-05,-4.9549e-08,3.06851e-11,-42380.4,34.4693], Tmin=(100,'K'), Tmax=(941.894,'K')), NASAPolynomial(coeffs=[29.8785,0.0127291,-2.29134e-06,3.75612e-10,-3.41757e-14,-50775.5,-126.104], Tmin=(941.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-354.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdOsH) + ring(Cyclobutane) + radical((O)CJOC) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1(O)OC(=O)C(C)[C]1O(31643)',
    structure = SMILES('[CH2]C1(O)OC(=O)C(C)[C]1O'),
    E0 = (-445.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.848359,0.0849591,-2.3967e-05,-5.6665e-08,3.67442e-11,-53394.4,34.7883], Tmin=(100,'K'), Tmax=(895.081,'K')), NASAPolynomial(coeffs=[29.7032,0.00887189,2.24903e-06,-7.46974e-10,5.42562e-14,-61284.9,-122.744], Tmin=(895.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-445.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(butyrolactone) + radical(C2CsJOH) + radical(CJC(C)OC)"""),
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
    label = '[CH2]C(O)(O[C]=O)C(=C)O(31644)',
    structure = SMILES('[CH2]C(O)(O[C]=O)C(=C)O'),
    E0 = (-360.941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.32394,0.114625,-0.00013698,7.3949e-08,-1.4641e-11,-43162.4,36.4386], Tmin=(100,'K'), Tmax=(1425.88,'K')), NASAPolynomial(coeffs=[32.3953,-0.00106871,3.97568e-06,-9.53938e-10,6.96093e-14,-51203.5,-136.857], Tmin=(1425.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-360.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(CJC(C)OC) + radical((O)CJOC)"""),
)

species(
    label = 'CC=C(O)C[C](O)O[C]=O(31645)',
    structure = SMILES('CC=C(O)C[C](O)O[C]=O'),
    E0 = (-376.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.784973,0.103382,-0.000110473,6.02774e-08,-1.30674e-11,-45110.5,38.0064], Tmin=(100,'K'), Tmax=(1119.44,'K')), NASAPolynomial(coeffs=[19.4515,0.0310723,-1.35804e-05,2.57425e-09,-1.80747e-13,-49641.2,-61.906], Tmin=(1119.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-376.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical((O)CJOC) + radical(Cs_P)"""),
)

species(
    label = 'CC=C(O)C1(O)CC(=O)O1(31646)',
    structure = SMILES('CC=C(O)C1(O)CC(=O)O1'),
    E0 = (-692.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.590279,0.0714728,2.38799e-05,-1.12539e-07,5.79808e-11,-83128.2,32.9329], Tmin=(100,'K'), Tmax=(901.976,'K')), NASAPolynomial(coeffs=[32.7664,0.00300728,5.5932e-06,-1.35123e-09,9.14293e-14,-92378,-142.471], Tmin=(901.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-692.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone)"""),
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
    label = '[CH2]C([O])(O)C(O)=CC(31647)',
    structure = SMILES('[CH2]C([O])(O)C(O)=CC'),
    E0 = (-206.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57438,0.103411,-0.000118869,6.43734e-08,-1.29463e-11,-24595.2,35.268], Tmin=(100,'K'), Tmax=(1387.58,'K')), NASAPolynomial(coeffs=[27.1597,0.00660338,8.89175e-07,-4.23158e-10,3.58034e-14,-31224,-107.922], Tmin=(1387.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-206.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)OJ) + radical(C=CC(O)2CJ)"""),
)

species(
    label = '[CH2]C(O)(O[C]=O)C(=O)CC(8749)',
    structure = SMILES('[CH2]C(O)(O[C]=O)C(=O)CC'),
    E0 = (-410.546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
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
    molecularWeight = (144.125,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4712.61,'J/mol'), sigma=(7.47664,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=736.10 K, Pc=25.59 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.63452,0.116926,-0.000134099,7.56879e-08,-1.65435e-11,-49168.4,37.6907], Tmin=(100,'K'), Tmax=(1127.37,'K')), NASAPolynomial(coeffs=[24.9136,0.0227301,-8.767e-06,1.57248e-09,-1.078e-13,-55154.2,-93.5711], Tmin=(1127.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-410.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdOsH) + radical(CJC(C)OC) + radical((O)CJOC)"""),
)

species(
    label = '[C]=O(361)',
    structure = SMILES('[C]=O'),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08916,0.00200416,-1.61661e-05,2.55058e-08,-1.16424e-11,52802.7,4.52505], Tmin=(100,'K'), Tmax=(856.11,'K')), NASAPolynomial(coeffs=[0.961625,0.00569045,-3.48044e-06,7.19202e-10,-5.08041e-14,53738.7,21.4663], Tmin=(856.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CdCdJ2_triplet)"""),
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
    label = 'CC=C(O)[C](O)O[C]=O(31648)',
    structure = SMILES('CC=C(O)[C](O)O[C]=O'),
    E0 = (-353.871,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.091724,0.0909408,-0.000105889,6.29203e-08,-1.4825e-11,-42414.2,32.6094], Tmin=(100,'K'), Tmax=(1033.93,'K')), NASAPolynomial(coeffs=[16.6679,0.026103,-1.18249e-05,2.26987e-09,-1.60221e-13,-45879.9,-48.8058], Tmin=(1033.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-353.871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical((O)CJOC) + radical(Cs_P)"""),
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
    E0 = (-396.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-335.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-339.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-396.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-174.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-260.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-235.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-232.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-313.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-325.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-355.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-322.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-341.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-316.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-90.8951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-271.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-337.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (58.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-130.916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-388.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-312.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-269.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (232.784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (27.6919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    products = ['CO2(13)', 'C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    products = ['C[CH]C1(O)CC1(O)O[C]=O(31633)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(61.6907,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra_cs2H] for rate rule [R4_S_D;doublebond_intra_HNd_secNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 57.6 to 61.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    products = ['[CH2]C1(O)OC(=O)C1(O)[CH]C(31634)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.65e+06,'s^-1'), n=1.3, Ea=(57.7392,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_HNd;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_HNd_secNd;radadd_intra_CO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CO2(13)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(23.3993,'m^3/(mol*s)'), n=2.021, Ea=(129.01,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd-O2d;CJ]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 127.8 to 129.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C(O)O[C]=O(5061)', 'CC=[C]O(31159)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00712612,'m^3/(mol*s)'), n=2.40979, Ea=(7.81798,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CdsJ] for rate rule [Cds-OsOs_Cds;CdsJ-O2s]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OH(5)', 'C=C(O[C]=O)C(O)=CC(31635)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.812049,'m^3/(mol*s)'), n=2.01336, Ea=(12.3541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;YJ] for rate rule [Cds-CdOs_Cds;OJ_pri]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC=C(O)C(C)([O])O[C]=O(31636)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C([O])(OC=O)C(O)=CC(31637)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.75172e+06,'s^-1'), n=1.80068, Ea=(127.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    products = ['C[CH]C(=O)C(C)(O)O[C]=O(9911)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00963743,'s^-1'), n=3.795, Ea=(83.0524,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_2H;O_H_out] + [R4H_SS(Cd)S;C_rad_out_2H;XH_out] for rate rule [R4H_SS(Cd)S;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[C]=C(O)C(C)(O)O[C]=O(31638)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    products = ['[CH2]C(O)(OC=O)C(=O)[CH]C(9912)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(67170.6,'s^-1'), n=1.77845, Ea=(41.7861,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;Y_rad_out;XH_out] for rate rule [R5H_SSSS;CO_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)(OC=O)C(O)=[C]C(31639)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_Cs;CO_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    products = ['[CH2]C=C(O)C(C)(O)O[C]=O(31640)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(121000,'s^-1'), n=1.9, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 92 used for R5H_SSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    products = ['[CH2]C=C(O)C([CH2])(O)OC=O(31641)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(37753.8,'s^-1'), n=1.925, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_RSSMS;CO_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][C]=O(669)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    products = ['CC1CC(O)(O[C]=O)[C]1O(31642)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.82e+08,'s^-1'), n=0.91, Ea=(125.52,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra_secNd;radadd_intra_cs2H] for rate rule [R4_Cs_RR_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    products = ['[CH2]C1(O)OC(=O)C(C)[C]1O(31643)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.31572e+09,'s^-1'), n=0.656667, Ea=(59.3431,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_secNd;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_secNd_HNd;radadd_intra_CO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH2(S)(23)', '[CH2]C(O)(O[C]=O)C(=C)O(31644)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.97e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC=C(O)C[C](O)O[C]=O(31645)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    products = ['CC=C(O)C1(O)CC(=O)O1(31646)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CO(12)', '[CH2]C([O])(O)C(O)=CC(31647)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.41e+07,'cm^3/(mol*s)'), n=0, Ea=(12.552,'kJ/mol'), T0=(1,'K'), Tmin=(250,'K'), Tmax=(2500,'K'), comment="""From training reaction 9 used for COm;O_rad/NonDe
Exact match found for rate rule [COm;O_rad/NonDe]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    products = ['[CH2]C(O)(O[C]=O)C(=O)CC(8749)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1713.16,'s^-1'), n=2.9, Ea=(127.847,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_Cs;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_Cs;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[C]=O(361)', '[CH2]C([O])(O)C(O)=CC(31647)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(19)', 'CC=C(O)[C](O)O[C]=O(31648)'],
    products = ['[CH2]C(O)(O[C]=O)C(O)=CC(9923)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '5110',
    isomers = [
        '[CH2]C(O)(O[C]=O)C(O)=CC(9923)',
    ],
    reactants = [
        ('CO2(13)', 'C=C(O)C(O)=CC(5562)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '5110',
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

