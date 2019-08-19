species(
    label = '[CH2]C(O)(C=[C]O)C(=O)CC(24011)',
    structure = SMILES('[CH2]C(O)(C=[C]O)C(=O)CC'),
    E0 = (-140.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1699.48,2811.34,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154466,'amu*angstrom^2'), symmetry=1, barrier=(3.55148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154466,'amu*angstrom^2'), symmetry=1, barrier=(3.55148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154466,'amu*angstrom^2'), symmetry=1, barrier=(3.55148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154466,'amu*angstrom^2'), symmetry=1, barrier=(3.55148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154466,'amu*angstrom^2'), symmetry=1, barrier=(3.55148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154466,'amu*angstrom^2'), symmetry=1, barrier=(3.55148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154466,'amu*angstrom^2'), symmetry=1, barrier=(3.55148,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.81638,0.134619,-0.000186387,1.31754e-07,-3.5377e-11,-16677.7,40.3469], Tmin=(100,'K'), Tmax=(742.807,'K')), NASAPolynomial(coeffs=[18.8896,0.0365712,-1.55609e-05,2.82154e-09,-1.8979e-13,-20125,-55.8897], Tmin=(742.807,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-140.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = 'HCCOH(50)',
    structure = SMILES('C#CO'),
    E0 = (80.0402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3615,1277.5,1000,287.694,287.7,1588.42],'cm^-1')),
        HinderedRotor(inertia=(0.002037,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05541,0.0252003,-3.80822e-05,3.09891e-08,-9.898e-12,9768.72,12.2272], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.3751,0.00549429,-1.88137e-06,2.93804e-10,-1.71772e-14,8932.78,-8.24498], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(80.0402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'CCC1([O])CC1(O)C=[C]O(25498)',
    structure = SMILES('CCC1([O])CC1(O)C=[C]O'),
    E0 = (-53.6993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30371,0.0991636,-6.17605e-05,-9.77275e-09,1.65994e-11,-6251.46,39.4337], Tmin=(100,'K'), Tmax=(930.416,'K')), NASAPolynomial(coeffs=[27.7211,0.0189841,-4.40378e-06,6.53127e-10,-4.61843e-14,-13583,-108.875], Tmin=(930.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.6993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(C=CJO) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1(O)C=C(O)C1([O])CC(25488)',
    structure = SMILES('[CH2]C1(O)C=C(O)C1([O])CC'),
    E0 = (-94.5789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.2017,0.0968013,-5.69743e-05,-1.19694e-08,1.64501e-11,-11171.8,38.5652], Tmin=(100,'K'), Tmax=(940.661,'K')), NASAPolynomial(coeffs=[26.8625,0.0210443,-5.66623e-06,9.20741e-10,-6.58119e-14,-18379.7,-105.361], Tmin=(940.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.5789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=CC(C)2OJ) + radical(C=CC(C)(O)CJ)"""),
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
    label = '[CH2]C(O)(C=C=O)C(=O)CC(9038)',
    structure = SMILES('[CH2]C(O)(C=C=O)C(=O)CC'),
    E0 = (-283.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1614.55,2905.61,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151214,'amu*angstrom^2'), symmetry=1, barrier=(3.47671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151214,'amu*angstrom^2'), symmetry=1, barrier=(3.47671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151214,'amu*angstrom^2'), symmetry=1, barrier=(3.47671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151214,'amu*angstrom^2'), symmetry=1, barrier=(3.47671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151214,'amu*angstrom^2'), symmetry=1, barrier=(3.47671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151214,'amu*angstrom^2'), symmetry=1, barrier=(3.47671,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77139,0.140835,-0.000233616,2.04318e-07,-6.8831e-11,-33891.8,36.4951], Tmin=(100,'K'), Tmax=(866.601,'K')), NASAPolynomial(coeffs=[12.5489,0.0466439,-2.19532e-05,4.08104e-09,-2.74321e-13,-35318.9,-24.4552], Tmin=(866.601,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-283.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = '[CH2]C(O)(C#CO)C(=O)CC(25499)',
    structure = SMILES('[CH2]C(O)(C#CO)C(=O)CC'),
    E0 = (-148.437,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2100,2250,500,550,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,223.63,1526.13,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149013,'amu*angstrom^2'), symmetry=1, barrier=(3.4261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149013,'amu*angstrom^2'), symmetry=1, barrier=(3.4261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149013,'amu*angstrom^2'), symmetry=1, barrier=(3.4261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149013,'amu*angstrom^2'), symmetry=1, barrier=(3.4261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149013,'amu*angstrom^2'), symmetry=1, barrier=(3.4261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149013,'amu*angstrom^2'), symmetry=1, barrier=(3.4261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149013,'amu*angstrom^2'), symmetry=1, barrier=(3.4261,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.73313,0.13239,-0.000194733,1.50739e-07,-4.58449e-11,-17652.4,37.5587], Tmin=(100,'K'), Tmax=(855.778,'K')), NASAPolynomial(coeffs=[17.5171,0.0362913,-1.55624e-05,2.80394e-09,-1.86626e-13,-20723,-51.0047], Tmin=(855.778,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CtH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtOs) + radical(CJC(C)(C=O)O)"""),
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
    label = 'C=C(O)C=[C]O(8492)',
    structure = SMILES('C=C(O)C=[C]O'),
    E0 = (-89.0201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3580,3650,1210,1345,900,1100,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.999644,'amu*angstrom^2'), symmetry=1, barrier=(22.9838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.998892,'amu*angstrom^2'), symmetry=1, barrier=(22.9665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.999021,'amu*angstrom^2'), symmetry=1, barrier=(22.9695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.383059,0.0711318,-8.04555e-05,4.12386e-08,-7.59281e-12,-10527,25.9678], Tmin=(100,'K'), Tmax=(1617.64,'K')), NASAPolynomial(coeffs=[20.9015,-0.000781229,4.10752e-06,-9.80442e-10,7.0779e-14,-14890.4,-79.1573], Tmin=(1617.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.0201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]O(172)',
    structure = SMILES('[CH]=[C]O'),
    E0 = (348.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,1058.91,1059.8],'cm^-1')),
        HinderedRotor(inertia=(0.315045,'amu*angstrom^2'), symmetry=1, barrier=(7.24351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22339,0.0198737,-3.18466e-05,3.10436e-08,-1.17427e-11,41917.3,13.7606], Tmin=(100,'K'), Tmax=(829.31,'K')), NASAPolynomial(coeffs=[3.78101,0.0111997,-5.33323e-06,1.02849e-09,-7.13893e-14,42030.6,12.4155], Tmin=(829.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(Cds_P) + radical(C=CJO)"""),
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
    label = 'C=C(C=[C]O)C(=O)CC(25500)',
    structure = SMILES('C=C(C=[C]O)C(=O)CC'),
    E0 = (-59.6526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.761577,0.097184,-9.77803e-05,4.99937e-08,-1.00626e-11,-6996.72,33.6165], Tmin=(100,'K'), Tmax=(1211.25,'K')), NASAPolynomial(coeffs=[20.6251,0.0265574,-1.0317e-05,1.85425e-09,-1.26718e-13,-12177.6,-73.6604], Tmin=(1211.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.6526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(O)([CH]C=O)C(=O)CC(9040)',
    structure = SMILES('[CH2]C(O)(C=C[O])C(=O)CC'),
    E0 = (-238.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,506.646,639.981,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60403,0.125948,-0.000157352,1.02784e-07,-2.65938e-11,-28502.2,37.8558], Tmin=(100,'K'), Tmax=(946.303,'K')), NASAPolynomial(coeffs=[19.5637,0.0364722,-1.55234e-05,2.86671e-09,-1.97038e-13,-32508.4,-63.0979], Tmin=(946.303,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CC(O)(C=O)CJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(O)([C]=CO)C(=O)CC(25501)',
    structure = SMILES('[CH2]C(O)([C]=CO)C(=O)CC'),
    E0 = (-142.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1670.2,2836.16,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153444,'amu*angstrom^2'), symmetry=1, barrier=(3.52798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153444,'amu*angstrom^2'), symmetry=1, barrier=(3.52798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153444,'amu*angstrom^2'), symmetry=1, barrier=(3.52798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153444,'amu*angstrom^2'), symmetry=1, barrier=(3.52798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153444,'amu*angstrom^2'), symmetry=1, barrier=(3.52798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153444,'amu*angstrom^2'), symmetry=1, barrier=(3.52798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153444,'amu*angstrom^2'), symmetry=1, barrier=(3.52798,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.20475,0.139432,-0.000192581,1.3546e-07,-3.71456e-11,-16889.2,39.0704], Tmin=(100,'K'), Tmax=(900.186,'K')), NASAPolynomial(coeffs=[22.0856,0.0314945,-1.27168e-05,2.2513e-09,-1.49624e-13,-21262.2,-75.5614], Tmin=(900.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-142.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CC(O)(C=O)CJ) + radical(Cds_S)"""),
)

species(
    label = 'CCC(=O)C(C)([O])C=[C]O(25502)',
    structure = SMILES('CCC(=O)C(C)([O])C=[C]O'),
    E0 = (-109.295,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,397.953,753.713,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.163937,'amu*angstrom^2'), symmetry=1, barrier=(3.76923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163937,'amu*angstrom^2'), symmetry=1, barrier=(3.76923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163937,'amu*angstrom^2'), symmetry=1, barrier=(3.76923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163937,'amu*angstrom^2'), symmetry=1, barrier=(3.76923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163937,'amu*angstrom^2'), symmetry=1, barrier=(3.76923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163937,'amu*angstrom^2'), symmetry=1, barrier=(3.76923,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49487,0.129497,-0.000191743,1.55946e-07,-5.02257e-11,-12955.5,40.1105], Tmin=(100,'K'), Tmax=(846.913,'K')), NASAPolynomial(coeffs=[13.8721,0.0447036,-1.99264e-05,3.66724e-09,-2.47307e-13,-15120.3,-28.8867], Tmin=(846.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = 'CCC(=O)C(C)(O)[C]=[C]O(25503)',
    structure = SMILES('CCC(=O)C(C)(O)[C]=[C]O'),
    E0 = (-115.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1670,1700,300,440,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180,1600,1764.23,2748.21,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15658,'amu*angstrom^2'), symmetry=1, barrier=(3.60009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15658,'amu*angstrom^2'), symmetry=1, barrier=(3.60009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15658,'amu*angstrom^2'), symmetry=1, barrier=(3.60009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15658,'amu*angstrom^2'), symmetry=1, barrier=(3.60009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15658,'amu*angstrom^2'), symmetry=1, barrier=(3.60009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15658,'amu*angstrom^2'), symmetry=1, barrier=(3.60009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15658,'amu*angstrom^2'), symmetry=1, barrier=(3.60009,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7055,0.130777,-0.000182889,1.35793e-07,-3.98697e-11,-13719.7,40.3534], Tmin=(100,'K'), Tmax=(837.534,'K')), NASAPolynomial(coeffs=[17.812,0.0375569,-1.59228e-05,2.88135e-09,-1.93433e-13,-16988.8,-50.3454], Tmin=(837.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-115.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C([O])(C=CO)C(=O)CC(25504)',
    structure = SMILES('[CH2]C([O])(C=CO)C(=O)CC'),
    E0 = (-135.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,544.188,606.99,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160348,'amu*angstrom^2'), symmetry=1, barrier=(3.68672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160348,'amu*angstrom^2'), symmetry=1, barrier=(3.68672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160348,'amu*angstrom^2'), symmetry=1, barrier=(3.68672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160348,'amu*angstrom^2'), symmetry=1, barrier=(3.68672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160348,'amu*angstrom^2'), symmetry=1, barrier=(3.68672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160348,'amu*angstrom^2'), symmetry=1, barrier=(3.68672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.97421,0.137887,-0.000200353,1.53977e-07,-4.66932e-11,-16125.8,38.758], Tmin=(100,'K'), Tmax=(837.092,'K')), NASAPolynomial(coeffs=[18.0726,0.0387739,-1.68005e-05,3.05676e-09,-2.05162e-13,-19365.7,-53.6968], Tmin=(837.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-135.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CC(O)(C=O)CJ) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = 'C[CH]C(=O)C(C)(O)C=[C]O(25505)',
    structure = SMILES('CC=C([O])C(C)(O)C=[C]O'),
    E0 = (-222.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.63771,0.110522,-0.000112632,5.76909e-08,-1.14792e-11,-26507.7,42.4159], Tmin=(100,'K'), Tmax=(1236.77,'K')), NASAPolynomial(coeffs=[25.0181,0.0243122,-8.07532e-06,1.33163e-09,-8.69319e-14,-33101.3,-91.8478], Tmin=(1236.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(O)(C=CO)C(=O)[CH]C(25506)',
    structure = SMILES('[CH2]C(O)(C=CO)C([O])=CC'),
    E0 = (-248.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.28623,0.118408,-0.000123991,6.33514e-08,-1.2333e-11,-29642.1,42.4932], Tmin=(100,'K'), Tmax=(1358.41,'K')), NASAPolynomial(coeffs=[29.3799,0.0176494,-4.43246e-06,6.03202e-10,-3.55241e-14,-37552,-117.425], Tmin=(1358.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]CC(=O)C(C)(O)C=[C]O(25507)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)C=[C]O'),
    E0 = (-141.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80185,0.131069,-0.000176656,1.23859e-07,-3.41312e-11,-16871.9,40.8609], Tmin=(100,'K'), Tmax=(893.174,'K')), NASAPolynomial(coeffs=[19.8059,0.0342938,-1.41216e-05,2.53502e-09,-1.70188e-13,-20731.5,-60.9413], Tmin=(893.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-141.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CJCC=O)"""),
)

species(
    label = 'CCC(=O)C(C)(O)[CH][C]=O(9042)',
    structure = SMILES('CCC(=O)C(C)(O)[CH][C]=O'),
    E0 = (-234.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20694,0.121887,-0.000171853,1.3609e-07,-4.33802e-11,-27977.3,40.641], Tmin=(100,'K'), Tmax=(815.054,'K')), NASAPolynomial(coeffs=[13.3573,0.0449303,-2.01367e-05,3.74398e-09,-2.55211e-13,-30169.3,-25.5277], Tmin=(815.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-234.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)C=CO(25508)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)C=CO'),
    E0 = (-168.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.36193,0.140495,-0.000189312,1.2778e-07,-3.34191e-11,-20038.8,39.7927], Tmin=(100,'K'), Tmax=(945.419,'K')), NASAPolynomial(coeffs=[24.2385,0.0279505,-1.07492e-05,1.86488e-09,-1.23005e-13,-25068.5,-87.0461], Tmin=(945.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CJCC=O) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = 'CC[C]1OCC1(O)C=[C]O(25509)',
    structure = SMILES('CC[C]1OCC1(O)C=[C]O'),
    E0 = (-62.4924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75532,0.110174,-0.000116519,6.40645e-08,-1.33743e-11,-7294.76,41.7075], Tmin=(100,'K'), Tmax=(1343.08,'K')), NASAPolynomial(coeffs=[22.2537,0.0241466,-4.22144e-06,2.72294e-10,-1.41684e-15,-12434.1,-76.3277], Tmin=(1343.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.4924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Oxetane) + radical(C2CsJOCs) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1(O)C=C(O)O[C]1CC(25435)',
    structure = SMILES('[CH2]C1(O)C=C(O)O[C]1CC'),
    E0 = (-153.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00006,0.113053,-0.000118698,6.24252e-08,-1.24681e-11,-18269.5,36.5844], Tmin=(100,'K'), Tmax=(1369.02,'K')), NASAPolynomial(coeffs=[26.0155,0.0197044,-3.82483e-06,3.53825e-10,-1.32769e-14,-24863.2,-103.44], Tmin=(1369.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-153.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + ring(2,3-Dihydrofuran) + radical(C2CsJOC(O)) + radical(C=CC(C)(O)CJ)"""),
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
    label = '[CH2]C(=C=O)C=[C]O(25510)',
    structure = SMILES('C=C([C]=O)C=[C]O'),
    E0 = (174.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,1855,455,950,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.02126,'amu*angstrom^2'), symmetry=1, barrier=(23.4808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02139,'amu*angstrom^2'), symmetry=1, barrier=(23.4838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02112,'amu*angstrom^2'), symmetry=1, barrier=(23.4775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0358532,0.0724548,-8.19803e-05,4.24781e-08,-8.18795e-12,21145.6,27.0165], Tmin=(100,'K'), Tmax=(1418.32,'K')), NASAPolynomial(coeffs=[22.1712,0.00223701,5.20982e-07,-1.73757e-10,1.27392e-14,15650.2,-84.7468], Tmin=(1418.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)CJ=O) + radical(C=CJO)"""),
)

species(
    label = 'CCC(=O)C(C)(O)C=C=O(9051)',
    structure = SMILES('CCC(=O)C(C)(O)C=C=O'),
    E0 = (-496.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.45282,0.130061,-0.000195697,1.64007e-07,-5.44902e-11,-59544.6,35.2992], Tmin=(100,'K'), Tmax=(839.713,'K')), NASAPolynomial(coeffs=[12.3276,0.0487334,-2.24019e-05,4.18078e-09,-2.84198e-13,-61305.9,-25.4836], Tmin=(839.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-496.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH)"""),
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
    label = '[CH2]C(O)(C=[C]O)C(C)=O(25511)',
    structure = SMILES('[CH2]C(O)(C=[C]O)C(C)=O'),
    E0 = (-114.941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11512,0.113617,-0.000152748,1.04064e-07,-2.75786e-11,-13641,36.4822], Tmin=(100,'K'), Tmax=(931.832,'K')), NASAPolynomial(coeffs=[19.5936,0.0247223,-9.65036e-06,1.68638e-09,-1.11594e-13,-17500.4,-61.963], Tmin=(931.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CC(O)(C=O)CJ) + radical(C=CJO)"""),
)

species(
    label = 'CCC(=O)C[C](O)C=[C]O(25512)',
    structure = SMILES('CCC(=O)C[C](O)C=[C]O'),
    E0 = (-168.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20382,0.122615,-0.000163702,1.16562e-07,-3.19266e-11,-20099.9,41.6614], Tmin=(100,'K'), Tmax=(702.38,'K')), NASAPolynomial(coeffs=[14.7846,0.0433507,-1.95996e-05,3.68154e-09,-2.53517e-13,-22636.7,-31.8949], Tmin=(702.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C2CsJOH)"""),
)

species(
    label = 'CCC(=O)[C](O)CC=[C]O(24010)',
    structure = SMILES('CCC([O])=C(O)CC=[C]O'),
    E0 = (-230.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52724,0.118259,-0.000123794,6.26135e-08,-1.1886e-11,-27439.9,45.5806], Tmin=(100,'K'), Tmax=(1460.41,'K')), NASAPolynomial(coeffs=[30.0899,0.0142206,-1.83575e-06,4.72527e-11,4.40283e-15,-35399,-118.764], Tmin=(1460.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-230.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC(=O)C1(O)C=C(O)C1(25485)',
    structure = SMILES('CCC(=O)C1(O)C=C(O)C1'),
    E0 = (-438.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59348,0.109488,-0.000109401,5.50283e-08,-1.07792e-11,-52478.7,35.0739], Tmin=(100,'K'), Tmax=(1253.75,'K')), NASAPolynomial(coeffs=[24.8014,0.0252749,-8.64543e-06,1.45183e-09,-9.57165e-14,-59097,-98.2344], Tmin=(1253.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-438.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH2]C(C=[C]O)=C(CC)OO(25513)',
    structure = SMILES('[CH2]C(C=[C]O)=C(CC)OO'),
    E0 = (90.0465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.89489,0.117175,-0.000123514,6.49216e-08,-1.3233e-11,11053,41.8427], Tmin=(100,'K'), Tmax=(1208.28,'K')), NASAPolynomial(coeffs=[26.1095,0.0244661,-8.42199e-06,1.41922e-09,-9.39395e-14,4285.57,-98.5605], Tmin=(1208.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.0465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = 'C=C(O)[C](CC)OC=[C]O(24009)',
    structure = SMILES('[CH2]C(O)=C(CC)OC=[C]O'),
    E0 = (-155.339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.69996,0.123534,-0.000133837,6.97343e-08,-1.36285e-11,-18421.4,45.5466], Tmin=(100,'K'), Tmax=(1411.84,'K')), NASAPolynomial(coeffs=[31.0645,0.0134068,-1.46251e-06,-3.03937e-11,1.01785e-14,-26513.5,-123.886], Tmin=(1411.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-155.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=C(O)CJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(O)(C=[C]O)C(=C)OC(25514)',
    structure = SMILES('[CH2]C(O)(C=[C]O)C(=C)OC'),
    E0 = (-91.8322,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,428.706,704.993,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155574,'amu*angstrom^2'), symmetry=1, barrier=(3.57695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155574,'amu*angstrom^2'), symmetry=1, barrier=(3.57695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155574,'amu*angstrom^2'), symmetry=1, barrier=(3.57695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155574,'amu*angstrom^2'), symmetry=1, barrier=(3.57695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155574,'amu*angstrom^2'), symmetry=1, barrier=(3.57695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155574,'amu*angstrom^2'), symmetry=1, barrier=(3.57695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155574,'amu*angstrom^2'), symmetry=1, barrier=(3.57695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10525,0.117851,-0.000123871,6.38049e-08,-1.26343e-11,-10811.1,43.2995], Tmin=(100,'K'), Tmax=(1268.58,'K')), NASAPolynomial(coeffs=[28.4624,0.0199013,-6.20129e-06,9.93995e-10,-6.44032e-14,-18440.6,-110.947], Tmin=(1268.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.8322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(O)(C=[C]O)C(O)=CC(25515)',
    structure = SMILES('[CH2]C(O)(C=[C]O)C(O)=CC'),
    E0 = (-146.542,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3615,3650,1210,1277.5,1345,900,1000,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.84318,0.117187,-0.00011958,5.49651e-08,-8.00909e-12,-17404.4,42.5238], Tmin=(100,'K'), Tmax=(980.45,'K')), NASAPolynomial(coeffs=[26.8287,0.0219876,-7.24654e-06,1.23473e-09,-8.45245e-14,-24073.3,-100.572], Tmin=(980.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CC(C)(O)CJ)"""),
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
    label = 'CCC(=O)[C](O)C=[C]O(25516)',
    structure = SMILES('CCC([O])=C(O)C=[C]O'),
    E0 = (-224.445,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.1003,0.111016,-0.000126232,6.75244e-08,-1.33508e-11,-26755,38.3947], Tmin=(100,'K'), Tmax=(1430.55,'K')), NASAPolynomial(coeffs=[28.9952,0.00657607,1.62103e-06,-6.06034e-10,4.93974e-14,-33861.7,-116.501], Tmin=(1430.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-224.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=CJO)"""),
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
    E0 = (-140.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-53.6993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-86.8028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-40.9116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (78.4648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-117.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (5.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-82.2184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-31.2806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (22.4212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (91.7968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (29.3832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (36.1394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (283.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-53.7492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-96.9527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-87.6396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-66.7389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-68.8135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (163.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (46.0672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-94.7524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (255.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-115.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (304.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (16.9604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (16.9604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-132.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (233.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-12.2256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (167.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-3.42831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (157.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['HCCOH(50)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['CCC1([O])CC1(O)C=[C]O(25498)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(86.6587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 81.9 to 86.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['[CH2]C1(O)C=C(O)C1([O])CC(25488)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_CO;carbonylbond_intra_Nd;radadd_intra_cdsingleNd]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(O)(C=C=O)C(=O)CC(9038)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.185e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""Estimated using template [Od_R;HJ] for rate rule [Od_Cdd;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C(O)(C#CO)C(=O)CC(25499)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5CO(71)', 'C=C(O)C=[C]O(8492)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-OneDeOs_Cds;CO_rad/NonDe]
Euclidian distance = 3.16227766017
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]O(172)', 'C=C(O)C(=O)CC(4626)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HCCOH(50)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)', 'C=C(C=[C]O)C(=O)CC(25500)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(931.236,'m^3/(mol*s)'), n=1.015, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;OJ_pri] for rate rule [Cds-TwoDe_Cds;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -7.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['[CH2]C(O)([CH]C=O)C(=O)CC(9040)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(O)([C]=CO)C(=O)CC(25501)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CCC(=O)C(C)([O])C=[C]O(25502)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CCC(=O)C(C)(O)[C]=[C]O(25503)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])(C=CO)C(=O)CC(25504)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.71035e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_singleNd]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['C[CH]C(=O)C(C)(O)C=[C]O(25505)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['[CH2]C(O)(C=CO)C(=O)[CH]C(25506)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(366176,'s^-1'), n=1.54456, Ea=(43.4053,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_DSSS;Cd_rad_out_singleNd;Cs_H_out_H/NonDeC]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['[CH2]CC(=O)C(C)(O)C=[C]O(25507)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['CCC(=O)C(C)(O)[CH][C]=O(9042)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_3;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['[CH2]CC(=O)C([CH2])(O)C=CO(25508)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.17493e+06,'s^-1'), n=1.01509, Ea=(71.5445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_DSSSS;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]O(172)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['CC[C]1OCC1(O)C=[C]O(25509)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['[CH2]C1(O)C=C(O)O[C]1CC(25435)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.18021e+10,'s^-1'), n=0.3725, Ea=(45.6056,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingleNd] for rate rule [R5_DS_CO;carbonyl_intra_Nd;radadd_intra_cdsingleNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH3CH2OH(54)', '[CH2]C(=C=O)C=[C]O(25510)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['CCC(=O)C(C)(O)C=C=O(9051)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(S)(23)', '[CH2]C(O)(C=[C]O)C(C)=O(25511)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['CCC(=O)C[C](O)C=[C]O(25512)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['CCC(=O)[C](O)CC=[C]O(24010)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    products = ['CCC(=O)C1(O)C=C(O)C1(25485)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriND_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C=[C]O)=C(CC)OO(25513)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C(O)[C](CC)OC=[C]O(24009)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(O)(C=[C]O)C(=C)OC(25514)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(O)(C=[C]O)C(O)=CC(25515)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2(19)', 'CCC(=O)[C](O)C=[C]O(25516)'],
    products = ['[CH2]C(O)(C=[C]O)C(=O)CC(24011)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/TDMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4145',
    isomers = [
        '[CH2]C(O)(C=[C]O)C(=O)CC(24011)',
    ],
    reactants = [
        ('HCCOH(50)', 'C=C(O)C(=O)CC(4626)'),
        ('HCCOH(50)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4145',
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

