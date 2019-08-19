species(
    label = '[CH2]C(O)(C(=O)CC)C(C)[O](6328)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C(C)[O]'),
    E0 = (-241.679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,566.073,582.003,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.90639,0.134962,-0.000171496,1.17148e-07,-3.20109e-11,-28859,39.5985], Tmin=(100,'K'), Tmax=(893.955,'K')), NASAPolynomial(coeffs=[18.3767,0.0442012,-1.91969e-05,3.56532e-09,-2.45143e-13,-32485.2,-55.9808], Tmin=(893.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CC(C)OJ) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'CH3CHO(37)',
    structure = SMILES('CC=O'),
    E0 = (-178.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1305.64,1305.66,1305.67,3976.84],'cm^-1')),
        HinderedRotor(inertia=(0.136163,'amu*angstrom^2'), symmetry=1, barrier=(3.13064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-178.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'CCC1([O])CC1(O)C(C)[O](6862)',
    structure = SMILES('CCC1([O])CC1(O)C(C)[O]'),
    E0 = (-147.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58198,0.109075,-8.80036e-05,2.86511e-08,-1.04433e-12,-17485.2,37.7114], Tmin=(100,'K'), Tmax=(1024.62,'K')), NASAPolynomial(coeffs=[23.9158,0.0334606,-1.23339e-05,2.20712e-09,-1.52909e-13,-23966.2,-92.0498], Tmin=(1024.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-147.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1(O)C(C)OC1([O])CC(6863)',
    structure = SMILES('[CH2]C1(O)C(C)OC1([O])CC'),
    E0 = (-175.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.51198,0.118053,-0.00011214,5.47891e-08,-1.02295e-11,-20832.7,38.8145], Tmin=(100,'K'), Tmax=(1485.28,'K')), NASAPolynomial(coeffs=[26.007,0.0270284,-5.85255e-06,6.36355e-10,-2.95861e-14,-27735.9,-104.775], Tmin=(1485.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-175.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CJC(C)2O) + radical(CC(C)(O)OJ)"""),
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
    label = '[CH2]C(O)(C(C)=O)C(=O)CC(6864)',
    structure = SMILES('[CH2]C(O)(C(C)=O)C(=O)CC'),
    E0 = (-420.323,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,365,385,505,600,445,480,1700,1720,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180,1600,1695.41,2815.45,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154336,'amu*angstrom^2'), symmetry=1, barrier=(3.54849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154336,'amu*angstrom^2'), symmetry=1, barrier=(3.54849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154336,'amu*angstrom^2'), symmetry=1, barrier=(3.54849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154336,'amu*angstrom^2'), symmetry=1, barrier=(3.54849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154336,'amu*angstrom^2'), symmetry=1, barrier=(3.54849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154336,'amu*angstrom^2'), symmetry=1, barrier=(3.54849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154336,'amu*angstrom^2'), symmetry=1, barrier=(3.54849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.16,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13792,0.11908,-0.000143804,9.66892e-08,-2.65453e-11,-50373.5,38.142], Tmin=(100,'K'), Tmax=(882.092,'K')), NASAPolynomial(coeffs=[14.76,0.046984,-2.11983e-05,4.02145e-09,-2.80281e-13,-53178,-36.5611], Tmin=(882.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-420.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'C[CH][O](176)',
    structure = SMILES('C[CH][O]'),
    E0 = (157.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1642.51],'cm^-1')),
        HinderedRotor(inertia=(0.123965,'amu*angstrom^2'), symmetry=1, barrier=(2.85019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65562,0.0114444,2.34936e-06,-4.83164e-09,1.17966e-12,18963.9,10.3625], Tmin=(100,'K'), Tmax=(1718.65,'K')), NASAPolynomial(coeffs=[6.06294,0.0136322,-6.35953e-06,1.18407e-09,-7.90642e-14,16985.9,-5.90233], Tmin=(1718.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CCsJOH) + radical(CCOJ)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91313,0.0225012,-8.59328e-06,4.80319e-10,2.48743e-13,-5274.24,14.655], Tmin=(100,'K'), Tmax=(1673.18,'K')), NASAPolynomial(coeffs=[7.46306,0.0158512,-6.42141e-06,1.12499e-09,-7.32055e-14,-7388.54,-11.406], Tmin=(1673.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.1874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2CO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C(O)C(C)[O](4716)',
    structure = SMILES('C=C(O)C(C)[O]'),
    E0 = (-173.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,331.439,332.514],'cm^-1')),
        HinderedRotor(inertia=(0.198362,'amu*angstrom^2'), symmetry=1, barrier=(15.5316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199032,'amu*angstrom^2'), symmetry=1, barrier=(15.5563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198135,'amu*angstrom^2'), symmetry=1, barrier=(15.5413,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11637,0.0495766,-8.68826e-06,-3.4174e-08,1.94555e-11,-20753.9,23.5605], Tmin=(100,'K'), Tmax=(942.703,'K')), NASAPolynomial(coeffs=[18.2893,0.0112553,-2.68025e-06,4.49556e-10,-3.52443e-14,-25526.7,-66.4172], Tmin=(942.703,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ)"""),
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
    label = '[CH2]C(O)(C=O)C(=O)CC(5970)',
    structure = SMILES('[CH2]C(O)(C=O)C(=O)CC'),
    E0 = (-365.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3615,1277.5,1000,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.500114,0.10405,-0.000128411,8.70522e-08,-2.39327e-11,-43796.2,34.2575], Tmin=(100,'K'), Tmax=(883.269,'K')), NASAPolynomial(coeffs=[13.9292,0.0387082,-1.74495e-05,3.30532e-09,-2.30069e-13,-46345.3,-33.5653], Tmin=(883.269,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CJC(C)(C=O)O)"""),
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
    label = 'C=C(C(=O)CC)C(C)[O](6865)',
    structure = SMILES('C=C(C(=O)CC)C(C)[O]'),
    E0 = (-160.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.161,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.587209,0.0823862,-6.19134e-05,2.53366e-08,-4.58649e-12,-19168.6,32.3339], Tmin=(100,'K'), Tmax=(1209.11,'K')), NASAPolynomial(coeffs=[9.66521,0.0523542,-2.46563e-05,4.79418e-09,-3.3907e-13,-21363.9,-13.186], Tmin=(1209.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-160.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(O)([C](C)O)C(=O)CC(6866)',
    structure = SMILES('[CH2]C(O)([C](C)O)C(=O)CC'),
    E0 = (-295.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43779,0.150406,-0.000224877,1.79374e-07,-5.63098e-11,-35306.4,41.1097], Tmin=(100,'K'), Tmax=(851.937,'K')), NASAPolynomial(coeffs=[17.6282,0.0450777,-1.98575e-05,3.62676e-09,-2.43221e-13,-38322,-50.1146], Tmin=(851.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-295.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(C2CsJOH)"""),
)

species(
    label = 'CCC(=O)C(C)([O])C(C)[O](6867)',
    structure = SMILES('CCC(=O)C(C)([O])C(C)[O]'),
    E0 = (-211.307,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,180,180,180,180,914.171,1312.9,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.170168,'amu*angstrom^2'), symmetry=1, barrier=(3.9125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170168,'amu*angstrom^2'), symmetry=1, barrier=(3.9125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170168,'amu*angstrom^2'), symmetry=1, barrier=(3.9125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170168,'amu*angstrom^2'), symmetry=1, barrier=(3.9125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170168,'amu*angstrom^2'), symmetry=1, barrier=(3.9125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170168,'amu*angstrom^2'), symmetry=1, barrier=(3.9125,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.685726,0.114779,-0.0001234,5.70221e-08,2.19933e-12,-25256.7,38.0961], Tmin=(100,'K'), Tmax=(591.923,'K')), NASAPolynomial(coeffs=[11.0791,0.0556773,-2.53256e-05,4.78961e-09,-3.31825e-13,-27006.9,-15.5122], Tmin=(591.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-211.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CC(C)OJ) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'CCC(=O)C(C)(O)[C](C)[O](6868)',
    structure = SMILES('CCC(=O)C(C)(O)[C](C)[O]'),
    E0 = (-276.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75472,0.134427,-0.000185601,1.44689e-07,-4.56524e-11,-33075,38.8527], Tmin=(100,'K'), Tmax=(798.139,'K')), NASAPolynomial(coeffs=[14.4358,0.0505065,-2.2661e-05,4.22686e-09,-2.8917e-13,-35570.9,-35.052], Tmin=(798.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(O)C([CH2])(O)C(=O)CC(6869)',
    structure = SMILES('[CH2]C(O)C([CH2])(O)C(=O)CC'),
    E0 = (-260.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35196,0.14577,-0.00020325,1.49875e-07,-4.37062e-11,-31101.8,42.0175], Tmin=(100,'K'), Tmax=(842.814,'K')), NASAPolynomial(coeffs=[19.5929,0.0416236,-1.79014e-05,3.26972e-09,-2.21019e-13,-34801,-60.1018], Tmin=(842.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCO) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C([O])(C(=O)CC)C(C)O(6870)',
    structure = SMILES('[CH2]C([O])(C(=O)CC)C(C)O'),
    E0 = (-230.057,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,379.326,771.708,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.90963,0.138724,-0.000200266,1.60796e-07,-5.15458e-11,-27465,42.2054], Tmin=(100,'K'), Tmax=(833.394,'K')), NASAPolynomial(coeffs=[14.6024,0.0496085,-2.21188e-05,4.08804e-09,-2.77059e-13,-29874.7,-32.3916], Tmin=(833.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-230.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH2]C([O])C(C)(O)C(=O)CC(6871)',
    structure = SMILES('[CH2]C([O])C(C)(O)C(=O)CC'),
    E0 = (-241.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59665,0.128874,-0.000160422,1.09978e-07,-3.04938e-11,-28873.5,39.5053], Tmin=(100,'K'), Tmax=(877.496,'K')), NASAPolynomial(coeffs=[16.3233,0.047187,-2.07838e-05,3.88866e-09,-2.68543e-13,-32018.4,-44.6063], Tmin=(877.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C[CH]C(=O)C(C)(O)C(C)[O](6872)',
    structure = SMILES('CC=C([O])C(C)(O)C(C)[O]'),
    E0 = (-316.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58314,0.114106,-0.000112433,5.79739e-08,-1.18537e-11,-37896,41.6873], Tmin=(100,'K'), Tmax=(1190.59,'K')), NASAPolynomial(coeffs=[21.9353,0.0350912,-1.28843e-05,2.23179e-09,-1.48942e-13,-43496.1,-75.8781], Tmin=(1190.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-316.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(O)(C(=O)[CH]C)C(C)O(6873)',
    structure = SMILES('[CH2]C(O)(C([O])=CC)C(C)O'),
    E0 = (-333.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.63024,0.121226,-0.000133891,7.82001e-08,-1.81314e-11,-39934,42.9264], Tmin=(100,'K'), Tmax=(1053.16,'K')), NASAPolynomial(coeffs=[20.4633,0.0373109,-1.43708e-05,2.54099e-09,-1.71135e-13,-44587.6,-64.8063], Tmin=(1053.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-333.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]CC(=O)C(C)(O)C(C)[O](6874)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)C(C)[O]'),
    E0 = (-241.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59665,0.128874,-0.000160422,1.09978e-07,-3.04938e-11,-28873.5,39.5053], Tmin=(100,'K'), Tmax=(877.496,'K')), NASAPolynomial(coeffs=[16.3233,0.047187,-2.07838e-05,3.88866e-09,-2.68543e-13,-32018.4,-44.6063], Tmin=(877.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCC=O) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)C(C)O(6875)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)C(C)O'),
    E0 = (-260.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35196,0.14577,-0.00020325,1.49875e-07,-4.37062e-11,-31101.8,42.0175], Tmin=(100,'K'), Tmax=(842.814,'K')), NASAPolynomial(coeffs=[19.5929,0.0416236,-1.79014e-05,3.26972e-09,-2.21019e-13,-34801,-60.1018], Tmin=(842.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CJCC=O)"""),
)

species(
    label = 'CC[C]1OCC1(O)C(C)[O](5574)',
    structure = SMILES('CC[C]1OCC1(O)C(C)[O]'),
    E0 = (-155.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09107,0.108669,-0.000101358,4.65145e-08,-6.05645e-12,-18569,36.6237], Tmin=(100,'K'), Tmax=(854.569,'K')), NASAPolynomial(coeffs=[17.1005,0.0411144,-1.36644e-05,2.19607e-09,-1.39555e-13,-22320.7,-52.0401], Tmin=(854.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-155.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1(O)[C](CC)OOC1C(6464)',
    structure = SMILES('[CH2]C1(O)[C](CC)OOC1C'),
    E0 = (-37.2451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.81348,0.113766,-0.000109072,5.61347e-08,-1.13831e-11,-4258.27,36.609], Tmin=(100,'K'), Tmax=(1268.75,'K')), NASAPolynomial(coeffs=[21.9732,0.0344559,-1.02029e-05,1.50176e-09,-8.95014e-14,-9946.67,-82.4407], Tmin=(1268.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-37.2451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(CJC(C)2O) + radical(C2CsJOOC)"""),
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
    label = '[CH2]C(=C=O)C(C)[O](6876)',
    structure = SMILES('C=C([C]=O)C(C)[O]'),
    E0 = (73.8302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,410.947,416.093],'cm^-1')),
        HinderedRotor(inertia=(0.116213,'amu*angstrom^2'), symmetry=1, barrier=(13.8558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111697,'amu*angstrom^2'), symmetry=1, barrier=(13.9163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11518,'amu*angstrom^2'), symmetry=1, barrier=(13.9338,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39823,0.0576084,-4.63023e-05,1.83257e-08,-2.96323e-12,8972.61,25.6754], Tmin=(100,'K'), Tmax=(1423.86,'K')), NASAPolynomial(coeffs=[12.8693,0.0253831,-1.23538e-05,2.43071e-09,-1.72414e-13,5705.94,-33.7197], Tmin=(1423.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.8302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)CJ=O) + radical(CC(C)OJ)"""),
)

species(
    label = 'CCC(=O)C(C)(O)C(C)=O(6877)',
    structure = SMILES('CCC(=O)C(C)(O)C(C)=O'),
    E0 = (-631.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.566479,0.106503,-0.000108382,6.42718e-08,-1.61954e-11,-75844.9,36.4765], Tmin=(100,'K'), Tmax=(939.085,'K')), NASAPolynomial(coeffs=[11.773,0.0539428,-2.44252e-05,4.66868e-09,-3.27741e-13,-78162.4,-22.2783], Tmin=(939.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-631.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsCs)"""),
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
    label = '[CH2]C(O)(C[O])C(=O)CC(5530)',
    structure = SMILES('[CH2]C(O)(C[O])C(=O)CC'),
    E0 = (-209.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,544.926,612.52,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160633,'amu*angstrom^2'), symmetry=1, barrier=(3.69326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160633,'amu*angstrom^2'), symmetry=1, barrier=(3.69326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160633,'amu*angstrom^2'), symmetry=1, barrier=(3.69326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160633,'amu*angstrom^2'), symmetry=1, barrier=(3.69326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160633,'amu*angstrom^2'), symmetry=1, barrier=(3.69326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160633,'amu*angstrom^2'), symmetry=1, barrier=(3.69326,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4774.5,'J/mol'), sigma=(7.75786,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=745.77 K, Pc=23.2 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.23669,0.1259,-0.000195729,1.6687e-07,-5.58239e-11,-25068.2,36.334], Tmin=(100,'K'), Tmax=(850.788,'K')), NASAPolynomial(coeffs=[11.7735,0.0458102,-2.11646e-05,3.94237e-09,-2.66959e-13,-26597.1,-20.3056], Tmin=(850.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-209.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(O)(C(C)=O)C(C)[O](6878)',
    structure = SMILES('[CH2]C(O)(C(C)=O)C(C)[O]'),
    E0 = (-216.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180,180,180,180,1600,1718.1,2788.04,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154798,'amu*angstrom^2'), symmetry=1, barrier=(3.55912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154798,'amu*angstrom^2'), symmetry=1, barrier=(3.55912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154798,'amu*angstrom^2'), symmetry=1, barrier=(3.55912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154798,'amu*angstrom^2'), symmetry=1, barrier=(3.55912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154798,'amu*angstrom^2'), symmetry=1, barrier=(3.55912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154798,'amu*angstrom^2'), symmetry=1, barrier=(3.55912,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1414,0.113051,-0.000133722,8.22071e-08,-1.99723e-11,-25825,35.5139], Tmin=(100,'K'), Tmax=(1007.36,'K')), NASAPolynomial(coeffs=[19.2558,0.0320573,-1.31175e-05,2.39051e-09,-1.63675e-13,-29934.4,-63.04], Tmin=(1007.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-216.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CC(C)OJ) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'CCC(=O)[C](O)CC(C)[O](6327)',
    structure = SMILES('CCC([O])=C(O)CC(C)[O]'),
    E0 = (-332.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.66644,0.111665,-9.47744e-05,3.27473e-08,-1.46842e-12,-39762.4,41.1599], Tmin=(100,'K'), Tmax=(992.265,'K')), NASAPolynomial(coeffs=[24.5877,0.0309804,-1.08254e-05,1.88962e-09,-1.29825e-13,-46210.8,-91.5364], Tmin=(992.265,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-332.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'CCC(=O)C[C](O)C(C)[O](6879)',
    structure = SMILES('CCC(=O)C[C](O)C(C)[O]'),
    E0 = (-269.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49461,0.132445,-0.000195121,1.67918e-07,-5.79658e-11,-32181.4,40.4729], Tmin=(100,'K'), Tmax=(823.813,'K')), NASAPolynomial(coeffs=[9.86786,0.0588591,-2.76056e-05,5.22239e-09,-3.58661e-13,-33428.6,-8.34954], Tmin=(823.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-269.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'CCC(=O)C1(O)COC1C(6334)',
    structure = SMILES('CCC(=O)C1(O)COC1C'),
    E0 = (-499.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39112,0.107923,-0.000102609,5.39769e-08,-1.13317e-11,-59887.1,35.0209], Tmin=(100,'K'), Tmax=(1196.77,'K')), NASAPolynomial(coeffs=[19.0392,0.0374689,-1.15838e-05,1.75598e-09,-1.06506e-13,-64621.8,-66.5641], Tmin=(1196.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-499.619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Oxetane)"""),
)

species(
    label = '[CH2][C](C(=O)CC)C(C)OO(6880)',
    structure = SMILES('[CH2]C(=C([O])CC)C(C)OO'),
    E0 = (-140.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.64177,0.113093,-0.000103226,4.8146e-08,-8.92445e-12,-16722.6,41.6765], Tmin=(100,'K'), Tmax=(1302.65,'K')), NASAPolynomial(coeffs=[23.5558,0.0357207,-1.41329e-05,2.55041e-09,-1.73979e-13,-23287.3,-86.5496], Tmin=(1302.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-140.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(=C(CC)OO)C(C)[O](6881)',
    structure = SMILES('[CH2]C(=C(CC)OO)C(C)[O]'),
    E0 = (4.74164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.12328,0.107546,-9.43126e-05,4.27877e-08,-7.86484e-12,758.843,41.4], Tmin=(100,'K'), Tmax=(1291.89,'K')), NASAPolynomial(coeffs=[20.1358,0.0417218,-1.78843e-05,3.34721e-09,-2.32439e-13,-4733.98,-66.6072], Tmin=(1291.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.74164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(O)[C](CC)OC(C)[O](6326)',
    structure = SMILES('[CH2]C(O)=C(CC)OC(C)[O]'),
    E0 = (-316.137,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.0006,0.118447,-0.000115096,5.64866e-08,-1.08645e-11,-37795.1,41.7168], Tmin=(100,'K'), Tmax=(1270.15,'K')), NASAPolynomial(coeffs=[25.9337,0.0304752,-1.12049e-05,1.95701e-09,-1.31655e-13,-44891.3,-99.7303], Tmin=(1270.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-316.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(O)(C(=C)OC)C(C)[O](6882)',
    structure = SMILES('[CH2]C(O)(C(=C)OC)C(C)[O]'),
    E0 = (-186.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180,1025.13,1188.45,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.16085,'amu*angstrom^2'), symmetry=1, barrier=(3.69825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16085,'amu*angstrom^2'), symmetry=1, barrier=(3.69825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16085,'amu*angstrom^2'), symmetry=1, barrier=(3.69825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16085,'amu*angstrom^2'), symmetry=1, barrier=(3.69825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16085,'amu*angstrom^2'), symmetry=1, barrier=(3.69825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16085,'amu*angstrom^2'), symmetry=1, barrier=(3.69825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16085,'amu*angstrom^2'), symmetry=1, barrier=(3.69825,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.04267,0.121354,-0.000123446,6.38575e-08,-1.29316e-11,-22199.7,42.5411], Tmin=(100,'K'), Tmax=(1209.76,'K')), NASAPolynomial(coeffs=[25.3852,0.0306654,-1.09996e-05,1.89124e-09,-1.26149e-13,-28836,-95.0052], Tmin=(1209.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-186.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(O)(C(O)=CC)C(C)[O](6883)',
    structure = SMILES('[CH2]C(O)(C(O)=CC)C(C)[O]'),
    E0 = (-241.179,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180,1600,1673.81,2826.18,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153429,'amu*angstrom^2'), symmetry=1, barrier=(3.52763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153429,'amu*angstrom^2'), symmetry=1, barrier=(3.52763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153429,'amu*angstrom^2'), symmetry=1, barrier=(3.52763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153429,'amu*angstrom^2'), symmetry=1, barrier=(3.52763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153429,'amu*angstrom^2'), symmetry=1, barrier=(3.52763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153429,'amu*angstrom^2'), symmetry=1, barrier=(3.52763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153429,'amu*angstrom^2'), symmetry=1, barrier=(3.52763,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12141,0.124631,-0.000132546,7.17787e-08,-1.51793e-11,-28778.2,42.993], Tmin=(100,'K'), Tmax=(1161.68,'K')), NASAPolynomial(coeffs=[25.3059,0.030191,-1.06019e-05,1.79695e-09,-1.18843e-13,-35150.6,-93.4386], Tmin=(1161.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C=CC(C)(O)CJ)"""),
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
    label = 'CCC(=O)[C](O)C(C)[O](6884)',
    structure = SMILES('CCC([O])=C(O)C(C)[O]'),
    E0 = (-308.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.86287,0.092376,-6.39541e-05,3.59795e-09,9.0903e-12,-36969.9,36.9406], Tmin=(100,'K'), Tmax=(955.525,'K')), NASAPolynomial(coeffs=[24.2153,0.0216407,-6.67343e-06,1.14214e-09,-8.08252e-14,-43325.9,-91.0871], Tmin=(955.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CC(C)OJ) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C(O)([CH]C)C(=O)CC(6379)',
    structure = SMILES('[CH2]C(O)([CH]C)C(=O)CC'),
    E0 = (-99.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1754.16,2755.36,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16691,0.116456,-0.000135482,8.5667e-08,-2.18244e-11,-11803.5,37.2714], Tmin=(100,'K'), Tmax=(953.997,'K')), NASAPolynomial(coeffs=[16.727,0.0414301,-1.75175e-05,3.23263e-09,-2.22295e-13,-15217.7,-48.2139], Tmin=(953.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CCJCO) + radical(CJC(C)(C=O)O)"""),
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
    E0 = (-241.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-147.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-119.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-137.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-185.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-123.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-241.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-197.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-117.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-76.0083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-72.6293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-123.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-103.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-124.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-158.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-155.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-215.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-188.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-148.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-27.3426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-55.2535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-37.2451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (154.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-178.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (209.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (203.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-84.3603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-84.3603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-233.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-33.6937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (147.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-173.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (72.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-98.0656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (72.6119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (143.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['CH3CHO(37)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['CCC1([O])CC1(O)C(C)[O](6862)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(94.5285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 90.6 to 94.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['[CH2]C1(O)C(C)OC1([O])CC(6863)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(O)(C(C)=O)C(=O)CC(6864)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0366254,'m^3/(mol*s)'), n=1.743, Ea=(71.4418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsCs_O;YJ] for rate rule [CO-CsCs_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C[CH][O](176)', 'C=C(O)C(=O)CC(4626)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5CO(71)', 'C=C(O)C(C)[O](4716)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cds-OsCs_Cds;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3CHO(37)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(122.03,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 121.1 to 122.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3(17)', '[CH2]C(O)(C=O)C(=O)CC(5970)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.61258,'m^3/(mol*s)'), n=1.485, Ea=(32.0285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-HHH] for rate rule [CO-CsH_O;CsJ-HHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)', 'C=C(C(=O)CC)C(C)[O](6865)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.52832,'m^3/(mol*s)'), n=2.02802, Ea=(14.4047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;YJ] for rate rule [Cds-COCs_Cds;OJ_pri]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['[CH2]C(O)([C](C)O)C(=O)CC(6866)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CCC(=O)C(C)([O])C(C)[O](6867)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['CCC(=O)C(C)(O)[C](C)[O](6868)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])(C(=O)CC)C(C)O(6870)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(94.0113,'s^-1'), n=2.81534, Ea=(105.641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['[CH2]C([O])C(C)(O)C(=O)CC(6871)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['C[CH]C(=O)C(C)(O)C(C)[O](6872)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['[CH2]C(O)(C(=O)[CH]C)C(C)O(6873)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8e+10,'s^-1'), n=0, Ea=(25.7316,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""Estimated using template [R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_CC(O2d)CC;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['[CH2]CC(=O)C(C)(O)C(C)[O](6874)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['[CH2]CC(=O)C([CH2])(O)C(C)O(6875)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 140 used for R6H_SSSSS;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R6H_SSSSS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C[CH][O](176)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['CC[C]1OCC1(O)C(C)[O](5574)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['[CH2]C1(O)[C](CC)OOC1C(6464)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.81184e+09,'s^-1'), n=0.551229, Ea=(204.434,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 200.6 to 204.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH3CH2OH(54)', '[CH2]C(=C=O)C(C)[O](6876)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['CCC(=O)C(C)(O)C(C)=O(6877)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(S)(23)', '[CH2]C(O)(C[O])C(=O)CC(5530)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(S)(23)', '[CH2]C(O)(C(C)=O)C(C)[O](6878)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['CCC(=O)[C](O)CC(C)[O](6327)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['CCC(=O)C[C](O)C(C)[O](6879)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['CCC(=O)C1(O)COC1C(6334)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][C](C(=O)CC)C(C)OO(6880)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OOH_S;Y_rad_out]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(=C(CC)OO)C(C)[O](6881)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C(O)[C](CC)OC(C)[O](6326)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(O)(C(=C)OC)C(C)[O](6882)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(O)(C(O)=CC)C(C)[O](6883)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CH2(19)', 'CCC(=O)[C](O)C(C)[O](6884)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(4)', '[CH2]C(O)([CH]C)C(=O)CC(6379)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '1365',
    isomers = [
        '[CH2]C(O)(C(=O)CC)C(C)[O](6328)',
    ],
    reactants = [
        ('CH3CHO(37)', 'C=C(O)C(=O)CC(4626)'),
        ('CH3CHO(37)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1365',
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

