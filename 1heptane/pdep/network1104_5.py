species(
    label = 'CCC1OCC=1O(3474)',
    structure = SMILES('CCC1OCC=1O'),
    E0 = (-274.176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2750,3150,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4129.74,'J/mol'), sigma=(6.71856,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=645.06 K, Pc=30.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.829143,0.0463633,3.51866e-05,-9.43311e-08,4.38647e-11,-32840,22.4006], Tmin=(100,'K'), Tmax=(930.832,'K')), NASAPolynomial(coeffs=[23.6091,0.00832506,3.35504e-08,-7.60831e-11,-2.88049e-15,-39673.9,-99.7951], Tmin=(930.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-274.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene)"""),
)

species(
    label = 'CCC1OCC1=O(3198)',
    structure = SMILES('CCC1OCC1=O'),
    E0 = (-248.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4038.67,'J/mol'), sigma=(6.4294,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=630.83 K, Pc=34.48 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71322,0.0296461,5.9882e-05,-9.84158e-08,3.90783e-11,-29768.4,23.7596], Tmin=(100,'K'), Tmax=(979.251,'K')), NASAPolynomial(coeffs=[15.843,0.0231785,-8.71379e-06,1.72785e-09,-1.32224e-13,-34992.9,-56.6588], Tmin=(979.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane)"""),
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
    label = 'C2H5CHO(70)',
    structure = SMILES('CCC=O'),
    E0 = (-204.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.207559,'amu*angstrom^2'), symmetry=1, barrier=(4.77219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208362,'amu*angstrom^2'), symmetry=1, barrier=(4.79065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90578,0.0240644,-7.06356e-06,-9.81837e-10,5.55825e-13,-24535.9,13.5806], Tmin=(100,'K'), Tmax=(1712.49,'K')), NASAPolynomial(coeffs=[7.69109,0.0189242,-7.84934e-06,1.38273e-09,-8.99057e-14,-27060.1,-14.6647], Tmin=(1712.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""propanal""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2CO(28)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]C(O)=C=O(4607)',
    structure = SMILES('C=C(O)[C]=O'),
    E0 = (-127.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1855,455,950,3615,1277.5,1000,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(1.03686,'amu*angstrom^2'), symmetry=1, barrier=(23.8395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03759,'amu*angstrom^2'), symmetry=1, barrier=(23.8562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43271,0.0523177,-6.55172e-05,3.78375e-08,-8.27375e-12,-15211.5,15.8439], Tmin=(100,'K'), Tmax=(1138,'K')), NASAPolynomial(coeffs=[15.2759,0.00365942,-1.38033e-06,2.64519e-10,-1.95358e-14,-18362.2,-52.7311], Tmin=(1138,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
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
    label = 'C=C=C([O])CC(4608)',
    structure = SMILES('C=C=C([O])CC'),
    E0 = (47.3032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,220.708,222.018],'cm^-1')),
        HinderedRotor(inertia=(0.34305,'amu*angstrom^2'), symmetry=1, barrier=(11.7864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339408,'amu*angstrom^2'), symmetry=1, barrier=(11.791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36132,0.0508126,-3.08268e-05,2.54214e-09,3.1323e-12,5790.59,21.9667], Tmin=(100,'K'), Tmax=(1007,'K')), NASAPolynomial(coeffs=[12.5297,0.020834,-7.59769e-06,1.34845e-09,-9.28702e-14,2811.96,-35.6132], Tmin=(1007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.3032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C(=O)C(=O)CC(4620)',
    structure = SMILES('[CH2]C(=O)C(=O)CC'),
    E0 = (-162.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,365,385,505,600,445,480,1700,1720,352.088,2238.1],'cm^-1')),
        HinderedRotor(inertia=(0.0859868,'amu*angstrom^2'), symmetry=1, barrier=(7.56458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859642,'amu*angstrom^2'), symmetry=1, barrier=(7.56459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859615,'amu*angstrom^2'), symmetry=1, barrier=(7.56392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.37635,'amu*angstrom^2'), symmetry=1, barrier=(33.1076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15233,0.0680038,-7.99912e-05,5.95539e-08,-1.92307e-11,-19492.6,25.214], Tmin=(100,'K'), Tmax=(738.439,'K')), NASAPolynomial(coeffs=[6.798,0.0374228,-1.78728e-05,3.47419e-09,-2.45165e-13,-20326.4,-0.311287], Tmin=(738.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + radical(CJCC=O)"""),
)

species(
    label = 'C=C(O)C([O])=CC(4627)',
    structure = SMILES('C=C(O)C([O])=CC'),
    E0 = (-232.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.757823,'amu*angstrom^2'), symmetry=1, barrier=(17.4238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756881,'amu*angstrom^2'), symmetry=1, barrier=(17.4022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.759059,'amu*angstrom^2'), symmetry=1, barrier=(17.4522,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.145471,0.0794715,-8.63325e-05,4.63576e-08,-9.46463e-12,-27754.2,24.6653], Tmin=(100,'K'), Tmax=(1313.65,'K')), NASAPolynomial(coeffs=[19.8698,0.0123363,-2.60601e-06,2.80197e-10,-1.3067e-14,-32478.7,-75.3243], Tmin=(1313.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = 'C=C([O])C(=C)O(4628)',
    structure = SMILES('C=C([O])C(=C)O'),
    E0 = (-196.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.916777,'amu*angstrom^2'), symmetry=1, barrier=(21.0785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.920477,'amu*angstrom^2'), symmetry=1, barrier=(21.1636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.284274,0.0664439,-7.46756e-05,3.95333e-08,-7.70645e-12,-23433.6,21.3173], Tmin=(100,'K'), Tmax=(1467.22,'K')), NASAPolynomial(coeffs=[18.8952,0.00393949,1.25465e-06,-4.33245e-10,3.4763e-14,-27628.4,-71.2884], Tmin=(1467.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = 'C=C(O)[C](C)[O](3238)',
    structure = SMILES('[CH2]C(O)=C(C)[O]'),
    E0 = (-162.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,282.345],'cm^-1')),
        HinderedRotor(inertia=(0.210861,'amu*angstrom^2'), symmetry=1, barrier=(11.9288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210833,'amu*angstrom^2'), symmetry=1, barrier=(11.9289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210854,'amu*angstrom^2'), symmetry=1, barrier=(11.9288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.326794,0.0657206,-7.01616e-05,3.62865e-08,-6.97475e-12,-19374.7,24.4612], Tmin=(100,'K'), Tmax=(1479.78,'K')), NASAPolynomial(coeffs=[17.9574,0.00731207,-5.68459e-08,-2.06761e-10,2.00592e-14,-23415.5,-63.528], Tmin=(1479.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
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
    label = 'CCC([O])=C(C)[O](4610)',
    structure = SMILES('CCC([O])=C(C)[O]'),
    E0 = (-206.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676626,0.0660722,-5.25149e-05,1.70923e-08,-3.95645e-13,-24656.8,26.9281], Tmin=(100,'K'), Tmax=(982.295,'K')), NASAPolynomial(coeffs=[14.7481,0.022808,-7.88241e-06,1.34753e-09,-9.06905e-14,-28098.4,-44.1539], Tmin=(982.295,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-206.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC(O)=C([CH2])O(4611)',
    structure = SMILES('[CH2]CC(O)=C([CH2])O'),
    E0 = (-117.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,325,375,415,465,420,450,1700,1750,424.855],'cm^-1')),
        HinderedRotor(inertia=(0.112797,'amu*angstrom^2'), symmetry=1, barrier=(14.3925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112269,'amu*angstrom^2'), symmetry=1, barrier=(14.3897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112576,'amu*angstrom^2'), symmetry=1, barrier=(14.3902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112523,'amu*angstrom^2'), symmetry=1, barrier=(14.3941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112345,'amu*angstrom^2'), symmetry=1, barrier=(14.394,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06791,0.0870529,-9.14386e-05,4.55824e-08,-8.38467e-12,-13928.7,32.6758], Tmin=(100,'K'), Tmax=(1562.56,'K')), NASAPolynomial(coeffs=[23.7823,0.00723112,7.46054e-07,-3.86281e-10,3.20893e-14,-19716.1,-91.9721], Tmin=(1562.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]C([O])=C(C)O(4612)',
    structure = SMILES('C[CH]C([O])=C(C)O'),
    E0 = (-143.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,363.145,364.572],'cm^-1')),
        HinderedRotor(inertia=(0.00128459,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110513,'amu*angstrom^2'), symmetry=1, barrier=(10.4043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111849,'amu*angstrom^2'), symmetry=1, barrier=(10.4087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111367,'amu*angstrom^2'), symmetry=1, barrier=(10.4225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.539627,0.0644069,-3.47983e-05,-1.16762e-08,1.24888e-11,-17178.7,28.4164], Tmin=(100,'K'), Tmax=(930.903,'K')), NASAPolynomial(coeffs=[18.8584,0.0154804,-3.95891e-06,6.11722e-10,-4.24718e-14,-21880,-65.5818], Tmin=(930.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C([O])=C(O)CC(4613)',
    structure = SMILES('[CH2]C([O])=C(O)CC'),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.516209,0.0797578,-7.96062e-05,3.90546e-08,-7.19184e-12,-22064.1,30.4076], Tmin=(100,'K'), Tmax=(1528.49,'K')), NASAPolynomial(coeffs=[20.9225,0.0119894,-1.65417e-06,6.23996e-11,2.30461e-15,-27255.3,-77.6604], Tmin=(1528.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]CC([O])=C(C)O(4614)',
    structure = SMILES('[CH2]CC([O])=C(C)O'),
    E0 = (-138.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,322.053,322.255],'cm^-1')),
        HinderedRotor(inertia=(0.00162288,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143441,'amu*angstrom^2'), symmetry=1, barrier=(10.5683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143425,'amu*angstrom^2'), symmetry=1, barrier=(10.5664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143458,'amu*angstrom^2'), symmetry=1, barrier=(10.5666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.149023,0.0765103,-7.49082e-05,3.66676e-08,-6.8673e-12,-16509.4,30.1846], Tmin=(100,'K'), Tmax=(1435.78,'K')), NASAPolynomial(coeffs=[19.7458,0.0146968,-3.65663e-06,4.85292e-10,-2.76616e-14,-21563.9,-70.6997], Tmin=(1435.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(RCCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1(O)O[C]1CC(4615)',
    structure = SMILES('[CH2]C1(O)O[C]1CC'),
    E0 = (20.6862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20042,0.0878293,-9.26709e-05,4.788e-08,-9.04116e-12,2698.08,30.0842], Tmin=(100,'K'), Tmax=(1572.68,'K')), NASAPolynomial(coeffs=[20.7922,0.0104727,1.54077e-06,-7.17204e-10,6.07465e-14,-1570.47,-77.5535], Tmin=(1572.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.6862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(CJC(O)2C)"""),
)

species(
    label = 'CCC1([O])C[C]1O(4616)',
    structure = SMILES('CCC1([O])C[C]1O'),
    E0 = (37.1314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.530362,0.0680393,-5.50591e-05,2.00541e-08,-1.98015e-12,4597.89,25.8256], Tmin=(100,'K'), Tmax=(1065.38,'K')), NASAPolynomial(coeffs=[15.6667,0.023151,-8.67151e-06,1.54756e-09,-1.06188e-13,695.012,-51.3375], Tmin=(1065.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.1314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(C2CsJOH) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2][C]=C(CC)OO(4617)',
    structure = SMILES('[CH2][C]=C(CC)OO'),
    E0 = (246.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,3615,1310,387.5,850,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.506098,'amu*angstrom^2'), symmetry=1, barrier=(11.6362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.019633,'amu*angstrom^2'), symmetry=1, barrier=(11.674,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0557866,'amu*angstrom^2'), symmetry=1, barrier=(32.801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215355,'amu*angstrom^2'), symmetry=1, barrier=(4.95144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.504366,'amu*angstrom^2'), symmetry=1, barrier=(11.5964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04933,0.0675382,-6.15487e-05,3.10919e-08,-6.63485e-12,29778.8,27.3366], Tmin=(100,'K'), Tmax=(1095.28,'K')), NASAPolynomial(coeffs=[10.3014,0.0337491,-1.52738e-05,2.92545e-09,-2.05772e-13,27752.1,-18.1414], Tmin=(1095.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=O)C([O])CC(3472)',
    structure = SMILES('C=C([O])C([O])CC'),
    E0 = (-59.5017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,406.649,412.992,413.371,415.884],'cm^-1')),
        HinderedRotor(inertia=(0.00102253,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0832567,'amu*angstrom^2'), symmetry=1, barrier=(10.1116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0835688,'amu*angstrom^2'), symmetry=1, barrier=(10.0603,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4230.88,'J/mol'), sigma=(6.99701,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=660.85 K, Pc=28.02 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728515,0.0616743,-3.46694e-05,-2.82684e-09,6.53814e-12,-7029.53,28.6979], Tmin=(100,'K'), Tmax=(994.485,'K')), NASAPolynomial(coeffs=[16.0073,0.0221145,-8.02439e-06,1.44923e-09,-1.01988e-13,-11151.1,-50.3725], Tmin=(994.485,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.5017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C(O)=[C]CC(4618)',
    structure = SMILES('[CH2]C(O)=[C]CC'),
    E0 = (129.652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,363.856],'cm^-1')),
        HinderedRotor(inertia=(0.133029,'amu*angstrom^2'), symmetry=1, barrier=(12.4987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133026,'amu*angstrom^2'), symmetry=1, barrier=(12.4985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133031,'amu*angstrom^2'), symmetry=1, barrier=(12.4984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133053,'amu*angstrom^2'), symmetry=1, barrier=(12.4985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.906588,0.0600581,-4.15031e-05,4.57085e-09,4.70951e-12,15712.2,24.2074], Tmin=(100,'K'), Tmax=(943.802,'K')), NASAPolynomial(coeffs=[15.3144,0.0182743,-5.7356e-06,9.4915e-10,-6.41253e-14,12134,-49.0174], Tmin=(943.802,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ) + radical(Cds_S)"""),
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
    label = 'CCC([O])=[C]O(4619)',
    structure = SMILES('CCC([O])=[C]O'),
    E0 = (-62.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,297.522,297.522],'cm^-1')),
        HinderedRotor(inertia=(0.185706,'amu*angstrom^2'), symmetry=1, barrier=(11.6652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185707,'amu*angstrom^2'), symmetry=1, barrier=(11.6652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185706,'amu*angstrom^2'), symmetry=1, barrier=(11.6652,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.423469,0.0614354,-6.11557e-05,2.97291e-08,-5.38234e-12,-7352.28,27.4708], Tmin=(100,'K'), Tmax=(1581.22,'K')), NASAPolynomial(coeffs=[17.1254,0.00842786,-6.66531e-07,-7.63516e-11,1.03047e-14,-11289.4,-56.5069], Tmin=(1581.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C([O])C(=O)CC(4621)',
    structure = SMILES('[CH2]C([O])C(=O)CC'),
    E0 = (1.82978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,745.949,752.438,4000],'cm^-1')),
        HinderedRotor(inertia=(0.12237,'amu*angstrom^2'), symmetry=1, barrier=(2.81352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123979,'amu*angstrom^2'), symmetry=1, barrier=(2.85052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.568819,'amu*angstrom^2'), symmetry=1, barrier=(13.0783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.568584,'amu*angstrom^2'), symmetry=1, barrier=(13.0729,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15704,0.0648368,-5.84199e-05,2.96528e-08,-6.3845e-12,320.515,28.752], Tmin=(100,'K'), Tmax=(1085.18,'K')), NASAPolynomial(coeffs=[9.74089,0.0331967,-1.46854e-05,2.78524e-09,-1.94876e-13,-1542.5,-13.3619], Tmin=(1085.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.82978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(C=OCOJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)C(=O)[CH]C(4622)',
    structure = SMILES('[CH2]C(O)C([O])=CC'),
    E0 = (-90.5188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,307.36,309.933],'cm^-1')),
        HinderedRotor(inertia=(0.00176934,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131214,'amu*angstrom^2'), symmetry=1, barrier=(8.83614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131255,'amu*angstrom^2'), symmetry=1, barrier=(8.8339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130366,'amu*angstrom^2'), symmetry=1, barrier=(8.83892,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.345407,0.0723962,-6.93553e-05,3.42305e-08,-6.65407e-12,-10748.6,30.4624], Tmin=(100,'K'), Tmax=(1255.15,'K')), NASAPolynomial(coeffs=[16.6574,0.0204121,-7.2304e-06,1.23314e-09,-8.16962e-14,-14843.4,-51.9406], Tmin=(1255.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.5188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])O(4623)',
    structure = SMILES('[CH2]CC(=O)C([CH2])O'),
    E0 = (-30.3142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,340.872,340.947],'cm^-1')),
        HinderedRotor(inertia=(0.00145124,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00145011,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00145026,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115055,'amu*angstrom^2'), symmetry=1, barrier=(9.49054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115066,'amu*angstrom^2'), symmetry=1, barrier=(9.49109,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.735022,0.0732832,-7.55104e-05,4.2154e-08,-9.62104e-12,-3529.57,30.1206], Tmin=(100,'K'), Tmax=(1050.32,'K')), NASAPolynomial(coeffs=[12.2345,0.029489,-1.29662e-05,2.45533e-09,-1.71838e-13,-5945.18,-25.9221], Tmin=(1050.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.3142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCC=O) + radical(CJCO)"""),
)

species(
    label = 'CC[C]1OC[C]1O(3468)',
    structure = SMILES('CC[C]1OC[C]1O'),
    E0 = (28.3384,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2750,3150,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.796463,0.0704813,-7.93574e-05,5.36594e-08,-1.4466e-11,3523.56,25.5318], Tmin=(100,'K'), Tmax=(1032.13,'K')), NASAPolynomial(coeffs=[9.25965,0.0300685,-9.56005e-06,1.42913e-09,-8.37405e-14,2182.09,-13.6012], Tmin=(1032.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(O)([C]=O)CC(4624)',
    structure = SMILES('[CH2]C(O)([C]=O)CC'),
    E0 = (-62.9497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.046423,0.0909749,-0.000117436,7.92038e-08,-2.10729e-11,-7427.11,27.9044], Tmin=(100,'K'), Tmax=(922.848,'K')), NASAPolynomial(coeffs=[15.1456,0.0251272,-1.04076e-05,1.88718e-09,-1.27989e-13,-10231.1,-44.1685], Tmin=(922.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.9497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)(O)CJ=O) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2][C](O)C(=C)OC(4625)',
    structure = SMILES('[CH2]C(O)=C([CH2])OC'),
    E0 = (-122.502,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,424.68,424.68],'cm^-1')),
        HinderedRotor(inertia=(0.1154,'amu*angstrom^2'), symmetry=1, barrier=(14.7692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1154,'amu*angstrom^2'), symmetry=1, barrier=(14.7692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1154,'amu*angstrom^2'), symmetry=1, barrier=(14.7692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1154,'amu*angstrom^2'), symmetry=1, barrier=(14.7692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1154,'amu*angstrom^2'), symmetry=1, barrier=(14.7692,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.970434,0.0863008,-9.01089e-05,4.50471e-08,-8.34194e-12,-14534.8,29.9491], Tmin=(100,'K'), Tmax=(1544.82,'K')), NASAPolynomial(coeffs=[23.1991,0.00851401,1.84605e-07,-2.90257e-10,2.61109e-14,-20188,-91.2941], Tmin=(1544.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C(O)C([O])[CH]C(4629)',
    structure = SMILES('C=C(O)C([O])[CH]C'),
    E0 = (2.59554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,334.962,334.964,334.965],'cm^-1')),
        HinderedRotor(inertia=(0.00150265,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00150216,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212724,'amu*angstrom^2'), symmetry=1, barrier=(16.9365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212695,'amu*angstrom^2'), symmetry=1, barrier=(16.9363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582894,0.0601292,-1.74707e-05,-3.07624e-08,1.89857e-11,448.888,30.2159], Tmin=(100,'K'), Tmax=(949.697,'K')), NASAPolynomial(coeffs=[20.0817,0.0148432,-4.13141e-06,7.20283e-10,-5.43182e-14,-4916.08,-71.5952], Tmin=(949.697,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.59554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC([O])C(=C)O(4630)',
    structure = SMILES('[CH2]CC([O])C(=C)O'),
    E0 = (7.93984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,390.278,390.485,4000],'cm^-1')),
        HinderedRotor(inertia=(0.174445,'amu*angstrom^2'), symmetry=1, barrier=(18.8536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174546,'amu*angstrom^2'), symmetry=1, barrier=(18.8536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174581,'amu*angstrom^2'), symmetry=1, barrier=(18.8514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174468,'amu*angstrom^2'), symmetry=1, barrier=(18.8538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.425828,0.066073,-3.65313e-05,-9.05001e-09,1.07201e-11,1094.94,30.0697], Tmin=(100,'K'), Tmax=(963.282,'K')), NASAPolynomial(coeffs=[19.3512,0.0167769,-5.37961e-06,9.56946e-10,-6.94138e-14,-3910.11,-67.5799], Tmin=(963.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.93984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C(O)C([O])CC(4631)',
    structure = SMILES('[CH]=C(O)C([O])CC'),
    E0 = (49.7896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,366.065,366.946],'cm^-1')),
        HinderedRotor(inertia=(0.147874,'amu*angstrom^2'), symmetry=1, barrier=(14.1392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148321,'amu*angstrom^2'), symmetry=1, barrier=(14.1305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148582,'amu*angstrom^2'), symmetry=1, barrier=(14.1331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148361,'amu*angstrom^2'), symmetry=1, barrier=(14.1336,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4230.88,'J/mol'), sigma=(6.99701,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=660.85 K, Pc=28.02 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.312838,0.068886,-4.28217e-05,-3.80964e-09,9.19362e-12,6132.1,29.0919], Tmin=(100,'K'), Tmax=(960.921,'K')), NASAPolynomial(coeffs=[19.819,0.0162025,-5.09347e-06,8.95974e-10,-6.47584e-14,1066.86,-71.0865], Tmin=(960.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.7896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)[C](O)CC(4632)',
    structure = SMILES('[CH]C(O)=C(O)CC'),
    E0 = (-110.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00203323,0.0732248,-3.43625e-05,-2.16115e-08,1.78208e-11,-13190.3,27.458], Tmin=(100,'K'), Tmax=(926.829,'K')), NASAPolynomial(coeffs=[21.8543,0.0167143,-4.08005e-06,6.09833e-10,-4.25184e-14,-18864.5,-85.0647], Tmin=(926.829,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(C)([O])C(=C)O(4633)',
    structure = SMILES('[CH2]C(C)([O])C(=C)O'),
    E0 = (-10.1611,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.315962,0.0720543,-5.80802e-05,1.62467e-08,1.20349e-12,-1081.43,29.2422], Tmin=(100,'K'), Tmax=(974.616,'K')), NASAPolynomial(coeffs=[17.7672,0.0195244,-6.61846e-06,1.14579e-09,-7.89563e-14,-5389.88,-59.1532], Tmin=(974.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.1611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'CCC1=[C]CO1(4546)',
    structure = SMILES('CCC1=[C]CO1'),
    E0 = (178.224,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,3150,900,1100,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79869,0.0322241,3.29422e-05,-6.97051e-08,2.99827e-11,21529.4,19.3365], Tmin=(100,'K'), Tmax=(955.106,'K')), NASAPolynomial(coeffs=[15.3026,0.0159918,-4.89132e-06,9.0493e-10,-6.94473e-14,17110.7,-54.8197], Tmin=(955.106,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Cds_S)"""),
)

species(
    label = 'CC[C]1OCC1=O(3457)',
    structure = SMILES('CCC1OCC=1[O]'),
    E0 = (-136.371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,3150,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19029,0.0424584,2.73692e-05,-7.67279e-08,3.55461e-11,-16282.7,22.6051], Tmin=(100,'K'), Tmax=(935.736,'K')), NASAPolynomial(coeffs=[19.8913,0.0116493,-2.00315e-06,3.11153e-10,-2.79423e-14,-21933.6,-77.868], Tmin=(935.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1OCC=1O(4547)',
    structure = SMILES('[CH2]C1OCC=1O'),
    E0 = (-92.6139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,3150,900,1100,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44312,0.034996,3.46899e-05,-9.02202e-08,4.3116e-11,-11027,17.9755], Tmin=(100,'K'), Tmax=(910.411,'K')), NASAPolynomial(coeffs=[22.9347,-0.00334417,5.45201e-06,-1.14291e-09,7.37539e-14,-17264.5,-96.4572], Tmin=(910.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-92.6139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(C=C(O)CJ)"""),
)

species(
    label = 'OC1=[C]OC1(4548)',
    structure = SMILES('OC1=[C]OC1'),
    E0 = (30.0051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3615,1277.5,1000,589.169,589.17,589.17,589.17,589.17,589.17,589.171,589.171,589.171,589.171],'cm^-1')),
        HinderedRotor(inertia=(0.000485642,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31646,0.0174648,5.03519e-05,-9.59339e-08,4.33992e-11,3687.74,15.2757], Tmin=(100,'K'), Tmax=(907.409,'K')), NASAPolynomial(coeffs=[19.7699,-0.00788073,6.96492e-06,-1.39985e-09,9.12099e-14,-1603.77,-78.9351], Tmin=(907.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.0051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=CJO)"""),
)

species(
    label = 'C[CH]C1OCC=1O(4549)',
    structure = SMILES('CC=C1OC[C]1O'),
    E0 = (-87.4419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.974495,0.0523775,-4.48829e-06,-4.46483e-08,2.53992e-11,-10394.6,21.6603], Tmin=(100,'K'), Tmax=(903.188,'K')), NASAPolynomial(coeffs=[18.7744,0.0119669,-1.18354e-06,1.11399e-11,7.44309e-16,-15177,-71.0768], Tmin=(903.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.4419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(2methyleneoxetane) + radical(C2CsJOH)"""),
)

species(
    label = 'CCC1O[CH]C=1O(4550)',
    structure = SMILES('CCC1O[CH]C=1O'),
    E0 = (-80.2487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.464354,0.0557213,8.49605e-06,-7.20383e-08,3.76449e-11,-9503.92,23.1841], Tmin=(100,'K'), Tmax=(922.416,'K')), NASAPolynomial(coeffs=[26.0982,0.00186806,2.88089e-06,-6.28278e-10,3.6679e-14,-16670.9,-111.629], Tmin=(922.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.2487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(CCsJOC(O))"""),
)

species(
    label = '[CH2]CC1OCC=1O(4551)',
    structure = SMILES('[CH2]CC1OCC=1O'),
    E0 = (-68.9297,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,3150,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.870589,0.0470461,2.49167e-05,-8.2319e-08,3.95324e-11,-8157.53,24.0385], Tmin=(100,'K'), Tmax=(929.099,'K')), NASAPolynomial(coeffs=[23.377,0.00607434,7.77029e-07,-2.12843e-10,7.246e-15,-14753.4,-95.8765], Tmin=(929.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.9297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(RCCJ)"""),
)

species(
    label = 'CC[C]1OCC1[O](3462)',
    structure = SMILES('CC[C]1OCC1[O]'),
    E0 = (82.0713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.780465,0.0616481,-4.98802e-05,2.36947e-08,-4.56034e-12,9994.59,25.2813], Tmin=(100,'K'), Tmax=(1382.55,'K')), NASAPolynomial(coeffs=[11.2631,0.0270095,-7.62259e-06,1.06313e-09,-6.02284e-14,7507.98,-27.1972], Tmin=(1382.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.0713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CC(C)OJ) + radical(C2CsJOCs)"""),
)

species(
    label = 'C[CH]C1OC[C]1O(3465)',
    structure = SMILES('C[CH]C1OC[C]1O'),
    E0 = (47.5366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.266935,0.0640504,-5.23731e-05,2.36401e-08,-4.12509e-12,5867.27,29.304], Tmin=(100,'K'), Tmax=(1641.76,'K')), NASAPolynomial(coeffs=[13.2057,0.0222547,-4.80173e-06,5.1209e-10,-2.29604e-14,3003.07,-35.317], Tmin=(1641.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.5366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOH) + radical(CCJCO)"""),
)

species(
    label = 'CC[C]1O[CH]C1O(4552)',
    structure = SMILES('CC[C]1O[CH]C1O'),
    E0 = (32.1666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.364379,0.0722166,-7.39902e-05,4.20324e-08,-9.24211e-12,4006.32,26.4043], Tmin=(100,'K'), Tmax=(1269.23,'K')), NASAPolynomial(coeffs=[13.3925,0.0225748,-5.17846e-06,5.60644e-10,-2.39277e-14,1390.56,-36.8313], Tmin=(1269.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.1666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC1O[CH][C]1O(3467)',
    structure = SMILES('CCC1O[CH][C]1O'),
    E0 = (28.0906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0216247,0.0735913,-7.18966e-05,3.75778e-08,-7.476e-12,3533.69,27.1152], Tmin=(100,'K'), Tmax=(1433.3,'K')), NASAPolynomial(coeffs=[15.3621,0.0194771,-3.43539e-06,2.32777e-10,-2.53102e-15,297.182,-48.3657], Tmin=(1433.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.0906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOH) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]CC1OC[C]1O(3469)',
    structure = SMILES('[CH2]CC1OC[C]1O'),
    E0 = (52.8809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.482484,0.0657927,-5.75334e-05,2.82806e-08,-5.45921e-12,6496.62,27.8067], Tmin=(100,'K'), Tmax=(1431.34,'K')), NASAPolynomial(coeffs=[13.0339,0.0235833,-5.8237e-06,7.14277e-10,-3.62802e-14,3634.25,-34.6953], Tmin=(1431.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.8809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(RCCJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C[C]1OCC1O(4553)',
    structure = SMILES('[CH2]C[C]1OCC1O'),
    E0 = (56.9569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24311,0.0591239,-3.91518e-05,2.78762e-09,7.38113e-12,6951.37,25.6191], Tmin=(100,'K'), Tmax=(778.896,'K')), NASAPolynomial(coeffs=[10.2079,0.0281176,-8.38924e-06,1.23564e-09,-7.36881e-14,5098.84,-18.3181], Tmin=(778.896,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.9569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(RCCJ)"""),
)

species(
    label = 'CC1OCC=1O(4554)',
    structure = SMILES('CC1OCC=1O'),
    E0 = (-251.531,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2750,3150,900,1100,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49086,0.0343344,3.81925e-05,-8.95545e-08,4.12024e-11,-30142.5,17.1135], Tmin=(100,'K'), Tmax=(924.073,'K')), NASAPolynomial(coeffs=[21.2828,0.00257053,2.24547e-06,-4.88959e-10,2.6689e-14,-36102,-89.2625], Tmin=(924.073,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene)"""),
)

species(
    label = '[CH2]OC(=[C]O)CC(4555)',
    structure = SMILES('[CH2]OC(=[C]O)CC'),
    E0 = (10.5629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.390481,0.0665368,-3.58177e-05,-1.33251e-08,1.35292e-11,1412.18,30.1989], Tmin=(100,'K'), Tmax=(935.674,'K')), NASAPolynomial(coeffs=[20.2282,0.0140471,-3.47703e-06,5.44573e-10,-3.91697e-14,-3714.76,-71.7469], Tmin=(935.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.5629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=COCJ)"""),
)

species(
    label = 'CC[C]=C(O)C[O](4556)',
    structure = SMILES('CC[C]=C(O)C[O]'),
    E0 = (43.3052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,923.78,3585.25],'cm^-1')),
        HinderedRotor(inertia=(0.539279,'amu*angstrom^2'), symmetry=1, barrier=(12.3991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0569386,'amu*angstrom^2'), symmetry=1, barrier=(12.423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.540944,'amu*angstrom^2'), symmetry=1, barrier=(12.4374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.539344,'amu*angstrom^2'), symmetry=1, barrier=(12.4006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07904,0.0673401,-6.65958e-05,3.83338e-08,-9.35376e-12,5311,27.5262], Tmin=(100,'K'), Tmax=(969.514,'K')), NASAPolynomial(coeffs=[9.16894,0.0339626,-1.49547e-05,2.82338e-09,-1.969e-13,3742.36,-11.2523], Tmin=(969.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.3052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'O=C1[CH]OC1(3456)',
    structure = SMILES('[O]C1=COC1'),
    E0 = (-71.9343,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,655.682,655.763,655.969,656.277,656.288,656.368,656.427,656.458,656.465,656.505,656.709,1600.65],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.46246,0.0118602,6.92436e-05,-1.16956e-07,5.12739e-11,-8575.74,13.124], Tmin=(100,'K'), Tmax=(906.709,'K')), NASAPolynomial(coeffs=[20.2081,-0.00870958,7.79101e-06,-1.56811e-09,1.02188e-13,-14166.3,-83.8336], Tmin=(906.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.9343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1OCC1=O(3458)',
    structure = SMILES('[CH2]C1OCC1=O'),
    E0 = (-14.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07531,0.026826,3.30733e-05,-6.5517e-08,2.74969e-11,-1607.29,19.9672], Tmin=(100,'K'), Tmax=(970.74,'K')), NASAPolynomial(coeffs=[14.7552,0.0128592,-4.49779e-06,9.09177e-10,-7.22864e-14,-5872.75,-50.1195], Tmin=(970.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(CJC(C)OC)"""),
)

species(
    label = 'C[CH]C1OCC1=O(3459)',
    structure = SMILES('C[CH]C1OCC1=O'),
    E0 = (-48.4443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90613,0.0244605,6.8361e-05,-1.07633e-07,4.27657e-11,-5731.67,25.5632], Tmin=(100,'K'), Tmax=(970.177,'K')), NASAPolynomial(coeffs=[16.3285,0.0190137,-6.73251e-06,1.35672e-09,-1.07181e-13,-11072.2,-56.6811], Tmin=(970.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.4443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(CCJCO)"""),
)

species(
    label = 'CCC1O[CH]C1=O(3460)',
    structure = SMILES('CCC1OC=C1[O]'),
    E0 = (-138.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.961359,0.0392489,5.81125e-05,-1.23844e-07,5.60869e-11,-16466.3,22.4869], Tmin=(100,'K'), Tmax=(922.86,'K')), NASAPolynomial(coeffs=[26.1999,0.00125128,3.82899e-06,-8.00752e-10,4.56637e-14,-24164.9,-113.72], Tmin=(922.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC1OCC1=O(3461)',
    structure = SMILES('[CH2]CC1OCC1=O'),
    E0 = (-43.1,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2950,3150,900,1000,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75375,0.0303401,4.95709e-05,-8.63469e-08,3.47201e-11,-5085.82,25.4007], Tmin=(100,'K'), Tmax=(981.527,'K')), NASAPolynomial(coeffs=[15.6147,0.0209214,-7.96665e-06,1.59023e-09,-1.22026e-13,-10074.1,-52.7616], Tmin=(981.527,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(RCCJ)"""),
)

species(
    label = 'CCC1O[CH]C1[O](3463)',
    structure = SMILES('CCC1O[CH]C1[O]'),
    E0 = (81.8235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03925,0.0523869,1.85101e-06,-5.15712e-08,2.86092e-11,9959.86,23.1675], Tmin=(100,'K'), Tmax=(874.64,'K')), NASAPolynomial(coeffs=[16.6925,0.0178384,-2.41961e-06,1.00877e-10,6.02748e-16,5804.96,-58.3526], Tmin=(874.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.8235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CC(C)OJ) + radical(CCsJOCs)"""),
)

species(
    label = 'C[CH]C1OCC1[O](3464)',
    structure = SMILES('C[CH]C1OCC1[O]'),
    E0 = (101.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51096,0.0404778,2.83056e-05,-7.26402e-08,3.42077e-11,12283,24.5229], Tmin=(100,'K'), Tmax=(890.958,'K')), NASAPolynomial(coeffs=[14.5431,0.0208075,-3.96462e-06,4.3242e-10,-2.47064e-14,8419.3,-45.4955], Tmin=(890.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC1OCC1[O](3466)',
    structure = SMILES('[CH2]CC1OCC1[O]'),
    E0 = (106.614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36846,0.0462487,9.84825e-06,-5.16944e-08,2.62559e-11,12928.4,24.3246], Tmin=(100,'K'), Tmax=(887.903,'K')), NASAPolynomial(coeffs=[13.735,0.0228731,-5.28897e-06,6.87091e-10,-4.1297e-14,9457.73,-41.0445], Tmin=(887.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CC(C)OJ) + radical(RCCJ)"""),
)

species(
    label = 'CC1OCC1=O(1089)',
    structure = SMILES('CC1OCC1=O'),
    E0 = (-224.566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43871,0.0139223,7.70714e-05,-1.10914e-07,4.3161e-11,-26934.4,18.9107], Tmin=(100,'K'), Tmax=(966.418,'K')), NASAPolynomial(coeffs=[14.5226,0.0154476,-5.29354e-06,1.0888e-09,-8.83359e-14,-31676.9,-51.4269], Tmin=(966.418,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-224.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane)"""),
)

species(
    label = 'CH2O(15)',
    structure = SMILES('C=O'),
    E0 = (-119.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CCC=C=O(2396)',
    structure = SMILES('CCC=C=O'),
    E0 = (-101.683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,311.417],'cm^-1')),
        HinderedRotor(inertia=(0.144733,'amu*angstrom^2'), symmetry=1, barrier=(9.96061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143318,'amu*angstrom^2'), symmetry=1, barrier=(9.96905,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24367,0.0334186,-1.32661e-05,-1.87751e-09,1.89793e-12,-12162.2,17.1242], Tmin=(100,'K'), Tmax=(1185.57,'K')), NASAPolynomial(coeffs=[8.36996,0.0212602,-8.65138e-06,1.58277e-09,-1.08606e-13,-14213,-15.9971], Tmin=(1185.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.683,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = '[CH2]OC([C]=O)CC(3470)',
    structure = SMILES('[CH2]OC([C]=O)CC'),
    E0 = (-32.6937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.571335,0.0680979,-5.89244e-05,2.61449e-08,-4.63233e-12,-3802.65,30.3345], Tmin=(100,'K'), Tmax=(1355.1,'K')), NASAPolynomial(coeffs=[15.8154,0.0231002,-9.11486e-06,1.63998e-09,-1.11448e-13,-7934.07,-47.8414], Tmin=(1355.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.6937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-OdCsH) + radical(CsJOCC2) + radical(CCCJ=O)"""),
)

species(
    label = 'CC[CH]OC[C]=O(3471)',
    structure = SMILES('CC[CH]OC[C]=O'),
    E0 = (-28.6152,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3576.1,'J/mol'), sigma=(6.18358,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=558.58 K, Pc=34.32 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.837597,0.071101,-6.86937e-05,3.61878e-08,-7.88774e-12,-3329,26.9209], Tmin=(100,'K'), Tmax=(1088.74,'K')), NASAPolynomial(coeffs=[11.6602,0.0313394,-1.39128e-05,2.64399e-09,-1.85317e-13,-5685.6,-26.2118], Tmin=(1088.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.6152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC=C([O])C[O](3185)',
    structure = SMILES('CCC=C([O])C[O]'),
    E0 = (-56.7318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,306.647,306.815,310.866,3900.7],'cm^-1')),
        HinderedRotor(inertia=(0.00178854,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.095552,'amu*angstrom^2'), symmetry=1, barrier=(6.46406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.095449,'amu*angstrom^2'), symmetry=1, barrier=(6.46141,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5192,0.0582795,-4.76049e-05,2.39043e-08,-5.44871e-12,-6737.06,26.9742], Tmin=(100,'K'), Tmax=(991.402,'K')), NASAPolynomial(coeffs=[6.56758,0.03791,-1.67842e-05,3.17815e-09,-2.21998e-13,-7738.01,2.66249], Tmin=(991.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.7318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ)"""),
)

species(
    label = 'CCC=C1COO1(3194)',
    structure = SMILES('CCC=C1COO1'),
    E0 = (32.1945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75612,0.0347556,3.23581e-05,-6.02328e-08,2.32253e-11,3965.58,24.9512], Tmin=(100,'K'), Tmax=(1026.74,'K')), NASAPolynomial(coeffs=[11.905,0.0304365,-1.2785e-05,2.48762e-09,-1.81144e-13,25.1527,-33.3192], Tmin=(1026.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.1945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutane)"""),
)

species(
    label = 'CCOC1=COC1(3473)',
    structure = SMILES('CCOC1=COC1'),
    E0 = (-228.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723596,0.045054,4.95413e-05,-1.15312e-07,5.27018e-11,-27351.6,21.3281], Tmin=(100,'K'), Tmax=(927.925,'K')), NASAPolynomial(coeffs=[26.1043,0.00509,1.88569e-06,-4.22614e-10,1.93451e-14,-35051.6,-115.331], Tmin=(927.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene)"""),
)

species(
    label = 'C=C1OOC1CC(3475)',
    structure = SMILES('C=C1OOC1CC'),
    E0 = (24.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58461,0.0343002,4.55884e-05,-8.17567e-08,3.25374e-11,3082.66,24.1296], Tmin=(100,'K'), Tmax=(995.949,'K')), NASAPolynomial(coeffs=[15.3038,0.0252698,-1.01961e-05,2.02922e-09,-1.5282e-13,-1934.9,-53.4724], Tmin=(995.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'CCC1OC=C1O(3476)',
    structure = SMILES('CCC1OC=C1O'),
    E0 = (-275.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.598555,0.0431727,6.58694e-05,-1.4138e-07,6.43837e-11,-33023.5,22.2884], Tmin=(100,'K'), Tmax=(921.123,'K')), NASAPolynomial(coeffs=[29.93,-0.00209382,5.87769e-06,-1.19082e-09,7.09603e-14,-41910.3,-135.717], Tmin=(921.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-275.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH2]CC(=O)C(=C)O(5558)',
    structure = SMILES('[CH2]CC(=O)C(=C)O'),
    E0 = (-146.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,251.471,251.487],'cm^-1')),
        HinderedRotor(inertia=(0.00266583,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304411,'amu*angstrom^2'), symmetry=1, barrier=(13.6601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304398,'amu*angstrom^2'), symmetry=1, barrier=(13.6601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304419,'amu*angstrom^2'), symmetry=1, barrier=(13.6601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.661955,0.0750645,-8.3545e-05,4.91082e-08,-1.16253e-11,-17546.5,24.0304], Tmin=(100,'K'), Tmax=(1021.82,'K')), NASAPolynomial(coeffs=[13.1639,0.0261243,-1.1702e-05,2.23537e-09,-1.57235e-13,-20101.5,-36.554], Tmin=(1021.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(CJCC=O)"""),
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
    label = '[CH]=C(O)C(=O)CC(5559)',
    structure = SMILES('[CH]=C(O)C(=O)CC'),
    E0 = (-111.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,253.138],'cm^-1')),
        HinderedRotor(inertia=(0.269783,'amu*angstrom^2'), symmetry=1, barrier=(12.268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269834,'amu*angstrom^2'), symmetry=1, barrier=(12.2682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269735,'amu*angstrom^2'), symmetry=1, barrier=(12.268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269833,'amu*angstrom^2'), symmetry=1, barrier=(12.2682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.7613,0.0731453,-7.95332e-05,4.60811e-08,-1.08237e-11,-13279.8,23.4437], Tmin=(100,'K'), Tmax=(1025.9,'K')), NASAPolynomial(coeffs=[12.568,0.027111,-1.22252e-05,2.34201e-09,-1.64964e-13,-15702.3,-33.8185], Tmin=(1025.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C(C)=O(3244)',
    structure = SMILES('C=C(O)C(C)=O'),
    E0 = (-333.051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,336.43],'cm^-1')),
        HinderedRotor(inertia=(0.168438,'amu*angstrom^2'), symmetry=1, barrier=(13.528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168465,'amu*angstrom^2'), symmetry=1, barrier=(13.5276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16843,'amu*angstrom^2'), symmetry=1, barrier=(13.5277,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34047,0.0507425,-3.5558e-05,7.1326e-09,1.49075e-12,-39954.3,19.6241], Tmin=(100,'K'), Tmax=(1048.53,'K')), NASAPolynomial(coeffs=[13.8723,0.0167055,-6.56506e-06,1.22369e-09,-8.67526e-14,-43339.3,-45.0383], Tmin=(1048.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-333.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = 'CCC(=O)C(C)=O(5560)',
    structure = SMILES('CCC(=O)C(C)=O'),
    E0 = (-374.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51646,0.0601905,-5.04133e-05,2.65316e-08,-6.55815e-12,-44953.8,22.5887], Tmin=(100,'K'), Tmax=(902.454,'K')), NASAPolynomial(coeffs=[5.68542,0.0417124,-1.97003e-05,3.84325e-09,-2.73025e-13,-45706.3,2.90371], Tmin=(902.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-374.471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs)"""),
)

species(
    label = 'C=C(O)C(=C)OC(5561)',
    structure = SMILES('C=C(O)C(=C)OC'),
    E0 = (-315.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.853151,'amu*angstrom^2'), symmetry=1, barrier=(19.6156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853387,'amu*angstrom^2'), symmetry=1, barrier=(19.621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853436,'amu*angstrom^2'), symmetry=1, barrier=(19.6222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.852549,'amu*angstrom^2'), symmetry=1, barrier=(19.6018,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.879363,0.085231,-8.65755e-05,4.21603e-08,-7.66487e-12,-37712.8,25.643], Tmin=(100,'K'), Tmax=(1548.07,'K')), NASAPolynomial(coeffs=[23.4666,0.00958852,-9.41414e-07,-3.17674e-11,7.00573e-15,-43724.6,-97.5224], Tmin=(1548.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-315.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C(CC)OO(5563)',
    structure = SMILES('C=C=C(CC)OO'),
    E0 = (33.9932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,300.836],'cm^-1')),
        HinderedRotor(inertia=(0.232368,'amu*angstrom^2'), symmetry=1, barrier=(5.3426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691206,'amu*angstrom^2'), symmetry=1, barrier=(15.8922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691038,'amu*angstrom^2'), symmetry=1, barrier=(15.8883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691275,'amu*angstrom^2'), symmetry=1, barrier=(15.8938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09237,0.0654842,-5.60116e-05,2.59736e-08,-5.06682e-12,4191.76,25.4768], Tmin=(100,'K'), Tmax=(1189.76,'K')), NASAPolynomial(coeffs=[10.852,0.0326716,-1.46425e-05,2.79265e-09,-1.95839e-13,1869.45,-23.3036], Tmin=(1189.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.9932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH]C(O)C(=O)CC(5564)',
    structure = SMILES('[CH]C(O)C(=O)CC'),
    E0 = (-3.19103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3615,1277.5,1000,1380,1390,370,380,2900,435,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.725856,0.0785881,-0.000108626,9.24897e-08,-3.25033e-11,-272.212,27.444], Tmin=(100,'K'), Tmax=(789.086,'K')), NASAPolynomial(coeffs=[6.89798,0.0392746,-1.86372e-05,3.57145e-09,-2.48192e-13,-996.407,0.71246], Tmin=(789.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.19103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=[C]C(O)=CC(31155)',
    structure = SMILES('C=[C]C(O)=CC'),
    E0 = (42.9982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.02093,'amu*angstrom^2'), symmetry=1, barrier=(23.4731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01944,'amu*angstrom^2'), symmetry=1, barrier=(23.4389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0196,'amu*angstrom^2'), symmetry=1, barrier=(23.4427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.502882,0.0666463,-6.75e-05,3.49508e-08,-6.96103e-12,5306.16,21.1922], Tmin=(100,'K'), Tmax=(1342.26,'K')), NASAPolynomial(coeffs=[16.3641,0.0144229,-3.60075e-06,4.62785e-10,-2.5161e-14,1494.65,-58.3348], Tmin=(1342.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.9982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C([O])C(O)=CC(31156)',
    structure = SMILES('C=C([O])C(O)=CC'),
    E0 = (-232.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.759458,'amu*angstrom^2'), symmetry=1, barrier=(17.4614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756943,'amu*angstrom^2'), symmetry=1, barrier=(17.4036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756852,'amu*angstrom^2'), symmetry=1, barrier=(17.4015,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.145462,0.0794714,-8.63323e-05,4.63573e-08,-9.46453e-12,-27754.2,24.6653], Tmin=(100,'K'), Tmax=(1313.67,'K')), NASAPolynomial(coeffs=[19.8697,0.0123364,-2.60611e-06,2.80219e-10,-1.30687e-14,-32478.6,-75.3236], Tmin=(1313.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)[C]=CC(31157)',
    structure = SMILES('C=C(O)[C]=CC'),
    E0 = (42.9982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.02093,'amu*angstrom^2'), symmetry=1, barrier=(23.4731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01944,'amu*angstrom^2'), symmetry=1, barrier=(23.4389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0196,'amu*angstrom^2'), symmetry=1, barrier=(23.4427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.502882,0.0666463,-6.75e-05,3.49508e-08,-6.96103e-12,5306.16,21.1922], Tmin=(100,'K'), Tmax=(1342.26,'K')), NASAPolynomial(coeffs=[16.3641,0.0144229,-3.60075e-06,4.62785e-10,-2.5161e-14,1494.65,-58.3348], Tmin=(1342.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.9982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C(O)C(=C)O(8478)',
    structure = SMILES('[CH]=C(O)C(=C)O'),
    E0 = (-86.7684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3120,650,792.5,1650,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(1.02953,'amu*angstrom^2'), symmetry=1, barrier=(23.671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02396,'amu*angstrom^2'), symmetry=1, barrier=(23.543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02862,'amu*angstrom^2'), symmetry=1, barrier=(23.6499,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.322663,0.0759296,-9.08749e-05,4.91324e-08,-9.62861e-12,-10263.7,22.3962], Tmin=(100,'K'), Tmax=(1482.05,'K')), NASAPolynomial(coeffs=[22.3441,-0.00152328,3.98935e-06,-9.50257e-10,6.95844e-14,-15194.9,-89.8453], Tmin=(1482.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.7684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C(O)C(=C)O(31158)',
    structure = SMILES('[CH2]C=C(O)C(=C)O'),
    E0 = (-251.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.12181,'amu*angstrom^2'), symmetry=1, barrier=(25.7926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12332,'amu*angstrom^2'), symmetry=1, barrier=(25.8274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1216,'amu*angstrom^2'), symmetry=1, barrier=(25.7877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12186,'amu*angstrom^2'), symmetry=1, barrier=(25.7938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04714,0.0838822,-8.70686e-05,4.2758e-08,-7.70112e-12,-30083.6,27.9774], Tmin=(100,'K'), Tmax=(1619.92,'K')), NASAPolynomial(coeffs=[22.8798,0.00666956,1.21675e-06,-4.84385e-10,3.87218e-14,-35456.7,-91.6557], Tmin=(1619.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
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
    label = 'C=C(O)C(O)=[C]C(31160)',
    structure = SMILES('C=C(O)C(O)=[C]C'),
    E0 = (-132.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,1685,370,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.846386,'amu*angstrom^2'), symmetry=1, barrier=(19.4601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.846124,'amu*angstrom^2'), symmetry=1, barrier=(19.4541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.846658,'amu*angstrom^2'), symmetry=1, barrier=(19.4663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.846444,'amu*angstrom^2'), symmetry=1, barrier=(19.4614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.524415,0.0878455,-0.000103087,5.81171e-08,-1.23276e-11,-15708.8,24.9954], Tmin=(100,'K'), Tmax=(1274.34,'K')), NASAPolynomial(coeffs=[22.1925,0.00883947,-1.02661e-06,-1.71752e-11,7.3761e-15,-20873.4,-87.6545], Tmin=(1274.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C(O)C(O)=CC(31161)',
    structure = SMILES('[CH]=C(O)C(O)=CC'),
    E0 = (-122.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.892967,'amu*angstrom^2'), symmetry=1, barrier=(20.5311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.892677,'amu*angstrom^2'), symmetry=1, barrier=(20.5244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.892566,'amu*angstrom^2'), symmetry=1, barrier=(20.5219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.893372,'amu*angstrom^2'), symmetry=1, barrier=(20.5404,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.750107,0.088929,-0.000102429,5.58194e-08,-1.13288e-11,-14584.3,25.736], Tmin=(100,'K'), Tmax=(1366.7,'K')), NASAPolynomial(coeffs=[23.6362,0.00641342,3.64202e-07,-2.87664e-10,2.56898e-14,-20209.4,-95.7252], Tmin=(1366.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([O])C(O)=CC(31162)',
    structure = SMILES('[CH2]C([O])C(O)=CC'),
    E0 = (2.03722,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,366.473,366.476],'cm^-1')),
        HinderedRotor(inertia=(0.141294,'amu*angstrom^2'), symmetry=1, barrier=(13.4661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141293,'amu*angstrom^2'), symmetry=1, barrier=(13.466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141295,'amu*angstrom^2'), symmetry=1, barrier=(13.4661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141293,'amu*angstrom^2'), symmetry=1, barrier=(13.4661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283774,0.0707516,-5.05058e-05,5.49106e-09,5.59804e-12,388.701,28.9881], Tmin=(100,'K'), Tmax=(967.788,'K')), NASAPolynomial(coeffs=[19.2911,0.0169688,-5.54956e-06,9.76858e-10,-6.95203e-14,-4450.66,-68.0839], Tmin=(967.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.03722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(O)C(O)=[C]C(31163)',
    structure = SMILES('[CH2]C(O)C(O)=[C]C'),
    E0 = (9.51822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1685,370,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,200.168],'cm^-1')),
        HinderedRotor(inertia=(0.00386889,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.449267,'amu*angstrom^2'), symmetry=1, barrier=(12.9706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.426874,'amu*angstrom^2'), symmetry=1, barrier=(12.9173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.449783,'amu*angstrom^2'), symmetry=1, barrier=(12.9095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.445007,'amu*angstrom^2'), symmetry=1, barrier=(12.9281,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0392474,0.0808275,-8.62719e-05,4.61606e-08,-9.57735e-12,1297.04,30.8137], Tmin=(100,'K'), Tmax=(1187.87,'K')), NASAPolynomial(coeffs=[18.9048,0.0170356,-5.7176e-06,9.51058e-10,-6.2495e-14,-3203.56,-63.8418], Tmin=(1187.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.51822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CJCO)"""),
)

species(
    label = 'C=C([O])C(O)[CH]C(11048)',
    structure = SMILES('C=C([O])C(O)[CH]C'),
    E0 = (-89.9604,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,300.618,300.7,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00186589,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00186341,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279062,'amu*angstrom^2'), symmetry=1, barrier=(17.9225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279467,'amu*angstrom^2'), symmetry=1, barrier=(17.9245,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773408,0.0603072,-3.14832e-05,-7.79839e-09,8.95214e-12,-10694,31.2248], Tmin=(100,'K'), Tmax=(968.867,'K')), NASAPolynomial(coeffs=[16.492,0.0198534,-6.69236e-06,1.18047e-09,-8.316e-14,-14887,-50.0315], Tmin=(968.867,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.9604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)[CH]C(16341)',
    structure = SMILES('[CH]=C(O)C(O)[CH]C'),
    E0 = (19.3309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,350,440,435,1725,402.616],'cm^-1')),
        HinderedRotor(inertia=(0.00104002,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115404,'amu*angstrom^2'), symmetry=1, barrier=(13.2734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115396,'amu*angstrom^2'), symmetry=1, barrier=(13.2735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115397,'amu*angstrom^2'), symmetry=1, barrier=(13.2735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115395,'amu*angstrom^2'), symmetry=1, barrier=(13.2736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.35222,0.0675747,-3.97791e-05,-8.68595e-09,1.16147e-11,2467.87,31.639], Tmin=(100,'K'), Tmax=(944.712,'K')), NASAPolynomial(coeffs=[20.3737,0.0138257,-3.69596e-06,6.11973e-10,-4.46799e-14,-2699.43,-71.1412], Tmin=(944.712,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.3309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCJCO)"""),
)

species(
    label = 'C[C]=C(O)[C](C)O(31164)',
    structure = SMILES('C[C]C(O)=C(C)O'),
    E0 = (-33.4494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.743435,'amu*angstrom^2'), symmetry=1, barrier=(17.093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744187,'amu*angstrom^2'), symmetry=1, barrier=(17.1103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743187,'amu*angstrom^2'), symmetry=1, barrier=(17.0873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744501,'amu*angstrom^2'), symmetry=1, barrier=(17.1175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744324,'amu*angstrom^2'), symmetry=1, barrier=(17.1135,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.956428,0.0868802,-9.04217e-05,4.47944e-08,-8.25531e-12,-3825.61,29.7053], Tmin=(100,'K'), Tmax=(1531.6,'K')), NASAPolynomial(coeffs=[24.0031,0.00827556,-2.96485e-07,-1.55319e-10,1.54736e-14,-9897.26,-96.2123], Tmin=(1531.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.4494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C=C(O)C([CH2])O(31165)',
    structure = SMILES('[CH2]C=C(O)C([CH2])O'),
    E0 = (-76.8244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,311.91],'cm^-1')),
        HinderedRotor(inertia=(0.0017328,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25271,'amu*angstrom^2'), symmetry=1, barrier=(17.4465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.252707,'amu*angstrom^2'), symmetry=1, barrier=(17.4466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25271,'amu*angstrom^2'), symmetry=1, barrier=(17.4465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.252711,'amu*angstrom^2'), symmetry=1, barrier=(17.4465,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.14439,0.0719594,-4.67619e-05,-4.40244e-09,1.07502e-11,-9089.27,29.9614], Tmin=(100,'K'), Tmax=(941.8,'K')), NASAPolynomial(coeffs=[21.49,0.0129779,-3.27478e-06,5.2804e-10,-3.87317e-14,-14514.8,-79.1975], Tmin=(941.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.8244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CJCO)"""),
)

species(
    label = '[CH2]C=C(O)[C](C)O(31166)',
    structure = SMILES('[CH2]C=C(O)[C](C)O'),
    E0 = (-111.786,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,3000,3100,440,815,1455,1000,400.151],'cm^-1')),
        HinderedRotor(inertia=(0.124318,'amu*angstrom^2'), symmetry=1, barrier=(14.0911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124327,'amu*angstrom^2'), symmetry=1, barrier=(14.0898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124332,'amu*angstrom^2'), symmetry=1, barrier=(14.0949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124279,'amu*angstrom^2'), symmetry=1, barrier=(14.0919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124231,'amu*angstrom^2'), symmetry=1, barrier=(14.0903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.303811,0.0806589,-8.15741e-05,4.06288e-08,-7.74696e-12,-13277.8,30.367], Tmin=(100,'K'), Tmax=(1379.32,'K')), NASAPolynomial(coeffs=[21.1156,0.0138305,-3.77368e-06,5.48421e-10,-3.34848e-14,-18738.4,-78.2323], Tmin=(1379.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(O)C(=C)O(5772)',
    structure = SMILES('C=C(O)C(=C)O'),
    E0 = (-333.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.01466,'amu*angstrom^2'), symmetry=1, barrier=(23.3291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01413,'amu*angstrom^2'), symmetry=1, barrier=(23.3169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01166,'amu*angstrom^2'), symmetry=1, barrier=(23.2599,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.498033,0.0751665,-8.30988e-05,4.22746e-08,-7.8095e-12,-39972.2,21.9417], Tmin=(100,'K'), Tmax=(1586.91,'K')), NASAPolynomial(coeffs=[21.8153,0.0015529,2.90201e-06,-7.522e-10,5.55296e-14,-44866.9,-89.1201], Tmin=(1586.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-333.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(O)=C(C)O(31167)',
    structure = SMILES('C=CC(O)=C(C)O'),
    E0 = (-370.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1976,0.0877436,-9.10797e-05,4.46402e-08,-8.05724e-12,-44357.6,27.774], Tmin=(100,'K'), Tmax=(1598.72,'K')), NASAPolynomial(coeffs=[24.3035,0.00683599,8.79619e-07,-3.98758e-10,3.21892e-14,-50325.6,-100.383], Tmin=(1598.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-370.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC=C(O)C(C)=O(31168)',
    structure = SMILES('CC=C(O)C(C)=O'),
    E0 = (-369.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.706729,0.0660995,-5.49688e-05,2.3391e-08,-4.00086e-12,-44265.9,23.7086], Tmin=(100,'K'), Tmax=(1391.19,'K')), NASAPolynomial(coeffs=[15.2261,0.0243525,-9.95644e-06,1.82065e-09,-1.24604e-13,-48305.7,-51.1328], Tmin=(1391.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-369.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs)"""),
)

species(
    label = 'C=C(O)C(O)[C]C(31169)',
    structure = SMILES('C=C(O)C(O)[C]C'),
    E0 = (0.198124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.201235,0.0823602,-9.02139e-05,5.21228e-08,-1.20015e-11,161.704,37.0062], Tmin=(100,'K'), Tmax=(1057.37,'K')), NASAPolynomial(coeffs=[15.0406,0.0262241,-1.05795e-05,1.91443e-09,-1.3065e-13,-2976.48,-35.4132], Tmin=(1057.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.198124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]C(O)C(O)=CC(31170)',
    structure = SMILES('[CH]C(O)C(O)=CC'),
    E0 = (10.3887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.361348,0.0803021,-9.09026e-05,5.41511e-08,-1.28617e-11,1380.35,26.8385], Tmin=(100,'K'), Tmax=(1025.05,'K')), NASAPolynomial(coeffs=[14.4094,0.0254823,-1.06811e-05,1.97638e-09,-1.36615e-13,-1499.6,-41.2826], Tmin=(1025.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.3887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'CC1CC(O)=C1O(31171)',
    structure = SMILES('CC1CC(O)=C1O'),
    E0 = (-316.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716318,0.0465129,4.18254e-05,-1.08121e-07,5.08748e-11,-37912.2,22.5083], Tmin=(100,'K'), Tmax=(918.644,'K')), NASAPolynomial(coeffs=[26.071,0.00311907,3.26936e-06,-7.39926e-10,4.4025e-14,-45398,-113.051], Tmin=(918.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-316.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene)"""),
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
    E0 = (34.7824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (75.6751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (79.6163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-12.5194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-25.9799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (257.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-176.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (34.5352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (17.6323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-34.4491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (274.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-71.2221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-82.9661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (46.2732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (37.1314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (353.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-45.2695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (372.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (319.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (115.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (67.2152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (45.7103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (28.3384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (100.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (191.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-184.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (139.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (125.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (241.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (98.0345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (147.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (206.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (75.4209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (43.5736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (138.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (130.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (134.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (142.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (104.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (70.6825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (55.3124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (91.4907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (116.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (91.4749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (168.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (18.8472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (51.5895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (-146.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (35.9398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (75.9397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (122.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (163.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (76.7232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (168.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (104.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (104.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (164.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (55.9046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (131.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (53.0638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (36.7064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (87.3989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (195.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (57.4675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (101.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (-24.4094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (-20.3308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (-51.2174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (-48.4475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (253.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (85.1973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (228.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (-132.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (79.8717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (48.9099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (-59.8722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (-19.4004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (-17.3313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS79',
    E0 = (64.9136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS80',
    E0 = (59.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS81',
    E0 = (100.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS82',
    E0 = (24.6912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS83',
    E0 = (24.8964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS84',
    E0 = (71.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS85',
    E0 = (-27.1186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS86',
    E0 = (3.89844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS87',
    E0 = (113.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS88',
    E0 = (-99.7195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS89',
    E0 = (-135.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS90',
    E0 = (-150.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS91',
    E0 = (-76.4621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS92',
    E0 = (-5.34099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS93',
    E0 = (-104.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS94',
    E0 = (86.8112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS95',
    E0 = (-154.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS96',
    E0 = (-1.38046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS97',
    E0 = (-230.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS98',
    E0 = (177.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS99',
    E0 = (53.6573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS100',
    E0 = (75.5667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS101',
    E0 = (-20.2932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS102',
    E0 = (75.5667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS103',
    E0 = (-20.2932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS104',
    E0 = (49.4191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS105',
    E0 = (-40.0423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS106',
    E0 = (176.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS107',
    E0 = (79.7437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS108',
    E0 = (88.9981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS109',
    E0 = (24.8986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS110',
    E0 = (25.4569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS111',
    E0 = (-94.6401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS112',
    E0 = (-27.1186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS113',
    E0 = (72.9184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS114',
    E0 = (-26.5603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS115',
    E0 = (82.7311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS116',
    E0 = (-109.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS117',
    E0 = (-25.0814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS118',
    E0 = (-150.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS119',
    E0 = (-76.4621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS120',
    E0 = (-51.8511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS121',
    E0 = (-77.2676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS122',
    E0 = (85.9973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS123',
    E0 = (-206.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS124',
    E0 = (-165.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS125',
    E0 = (57.0465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS126',
    E0 = (67.2371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS127',
    E0 = (-247.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction44',
    reactants = ['C2H5(29)', '[CH2]C(O)=C=O(4607)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(8.04,'m^3/(mol*s)'), n=1.68, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Od;CsJ] for rate rule [Ck_O;CsJ-CsHH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction45',
    reactants = ['OH(5)', 'C=C=C([O])CC(4608)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(41610,'cm^3/(mol*s)'), n=2.487, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 214 used for Ca_Cds-HH;OJ_pri
Exact match found for rate rule [Ca_Cds-HH;OJ_pri]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -7.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction59',
    reactants = ['H(3)', '[CH2]C(=O)C(=O)CC(4620)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.185e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""From training reaction 184 used for Od_R;HJ
Exact match found for rate rule [Od_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction67',
    reactants = ['H(3)', 'C=C(O)C([O])=CC(4627)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction68',
    reactants = ['CH3(17)', 'C=C([O])C(=C)O(4628)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.125041,'m^3/(mol*s)'), n=2.08753, Ea=(33.8923,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(S)(23)', 'C=C(O)[C](C)[O](3238)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.87e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.146, Ea=(0.0118826,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 6 used for carbene;C_pri/Cd
Exact match found for rate rule [carbene;C_pri/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C(O)=C([O])CC(4557)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(187863,'s^-1'), n=2.03383, Ea=(157.381,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeC;O_H_out] + [R3H_SS_2Cd;C_rad_out_H/NonDeC;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C(O)=C([O])CC(4557)'],
    products = ['CCC([O])=C(C)[O](4610)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.4947e+07,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_2Cd;C_rad_out_2H;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]CC(O)=C([CH2])O(4611)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.00963743,'s^-1'), n=3.795, Ea=(83.0524,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_2H;O_H_out] + [R4H_SS(Cd)S;C_rad_out_2H;XH_out] for rate rule [R4H_SS(Cd)S;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C[CH]C([O])=C(C)O(4612)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C([O])=C(O)CC(4613)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.11e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4H_SDS;O_rad_out;XH_out] for rate rule [R4H_SDS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]CC([O])=C(C)O(4614)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(121000,'s^-1'), n=1.9, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 92 used for R5H_SSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C1(O)O[C]1CC(4615)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C(O)=C([O])CC(4557)'],
    products = ['CCC1([O])C[C]1O(4616)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.57691e+12,'s^-1'), n=0.0709135, Ea=(222.075,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3_D;doublebond_intra;radadd_intra_cs2H] + [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 221.5 to 222.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2][C]=C(CC)OO(4617)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH;Y_rad_out] for rate rule [R2OOH_D;Cd_rad_out]
Euclidian distance = 1.41421356237
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(=O)C([O])CC(3472)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction57',
    reactants = ['O(4)', '[CH2]C(O)=[C]CC(4618)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction58',
    reactants = ['CH2(19)', 'CCC([O])=[C]O(4619)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2]C([O])C(=O)CC(4621)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2]C(O)C(=O)[CH]C(4622)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.2544e+06,'s^-1'), n=1.86276, Ea=(157.734,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2]CC(=O)C([CH2])O(4623)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(463.959,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2]C(O)=C([O])CC(4557)'],
    products = ['CC[C]1OC[C]1O(3468)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.43115e+08,'s^-1'), n=0.895859, Ea=(213.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_linear;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 212.8 to 213.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2]C(O)([C]=O)CC(4624)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH2][C](O)C(=C)OC(4625)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.16e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Y_12_11] for rate rule [Y_12_11a]
Euclidian distance = 1.0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction69',
    reactants = ['C=C(O)C([O])[CH]C(4629)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(956916,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction70',
    reactants = ['[CH2]CC([O])C(=C)O(4630)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH]=C(O)C([O])CC(4631)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH]=C(O)[C](O)CC(4632)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(868165,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction73',
    reactants = ['[CH2]C(C)([O])C(=C)O(4633)'],
    products = ['[CH2]C(O)=C([O])CC(4557)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['OH(5)', 'CCC1=[C]CO1(4546)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/NonDe] for rate rule [O_pri_rad;Cd_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', 'CC[C]1OCC1=O(3457)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH3(17)', '[CH2]C1OCC=1O(4547)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.02e+14,'cm^3/(mol*s)'), n=-0.32, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 127 used for C_rad/H2/Cd;C_methyl
Exact match found for rate rule [C_rad/H2/Cd;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C2H5(29)', 'OC1=[C]OC1(4548)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;C_pri_rad] for rate rule [Cd_rad/NonDe;C_rad/H2/Cs]
Euclidian distance = 2.2360679775
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', 'C[CH]C1OCC=1O(4549)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.71464e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', 'CCC1O[CH]C=1O(4550)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;H_rad] for rate rule [C_rad/H/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH2]CC1OCC=1O(4551)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['CC[C]1OCC1[O](3462)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C[CH]C1OC[C]1O(3465)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.5515e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CC[C]1O[CH]C1O(4552)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.5515e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CCC1O[CH][C]1O(3467)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]CC1OC[C]1O(3469)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C[C]1OCC1O(4553)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['CH2(S)(23)', 'CC1OCC=1O(4554)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.87e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.146, Ea=(0.0118826,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 6 used for carbene;C_pri/Cd
Exact match found for rate rule [carbene;C_pri/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]OC(=[C]O)CC(4555)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriND_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CC[C]=C(O)C[O](4556)'],
    products = ['CCC1OCC=1O(3474)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriND_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CCC1OCC=1O(3474)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1713.16,'s^-1'), n=2.9, Ea=(127.847,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_Cs;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction1',
    reactants = ['C2H5(29)', 'O=C1[CH]OC1(3456)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(4.37973e+08,'m^3/(mol*s)'), n=-0.525, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_sec_rad;C_rad/H2/Cs] for rate rule [C_rad/H/OneDeO;C_rad/H2/Cs]
Euclidian distance = 2.0
family: R_Recombination
Ea raised from -1.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'CC[C]1OCC1=O(3457)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""Estimated using template [C_rad/OneDe;H_rad] for rate rule [C_rad/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH3(17)', '[CH2]C1OCC1=O(3458)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C[CH]C1OCC1=O(3459)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'CCC1O[CH]C1=O(3460)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;H_rad] for rate rule [C_rad/H/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH2]CC1OCC1=O(3461)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC[C]1OCC1[O](3462)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CCC1O[CH]C1[O](3463)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[CH]C1OCC1[O](3464)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[CH]C1OC[C]1O(3465)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]CC1OCC1[O](3466)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CCC1O[CH][C]1O(3467)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CC[C]1OC[C]1O(3468)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CC1OC[C]1O(3469)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction63',
    reactants = ['CH2(S)(23)', 'CC1OCC1=O(1089)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C2H5CHO(70)', 'CH2CO(28)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(2.319e-07,'m^3/(mol*s)'), n=3.416, Ea=(322.616,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;doublebond] for rate rule [CO_HNd;mb_COC_2H]
Euclidian distance = 2.2360679775
family: 2+2_cycloaddition_CO"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2O(15)', 'CCC=C=O(2396)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(2.319e-07,'m^3/(mol*s)'), n=3.416, Ea=(322.616,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;doublebond] for rate rule [CO_2H;mb_COC_HNd]
Euclidian distance = 2.2360679775
family: 2+2_cycloaddition_CO"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]OC([C]=O)CC(3470)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC[CH]OC[C]=O(3471)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=O)C([O])CC(3472)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CCC=C([O])C[O](3185)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CCC=C1COO1(3194)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(3472.84,'s^-1'), n=2.78, Ea=(220.824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond;R2_doublebond;R_O_R] + [R_ROR;R1_doublebond;R2_doublebond_Cs;R_O] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_Cs;R_O_R]
Euclidian distance = 1.41421356237
family: ketoenol"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CCOC1=COC1(3473)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_Cs;R_O_C]
Euclidian distance = 1.41421356237
family: ketoenol"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C1OOC1CC(3475)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CCC1OC=C1O(3476)'],
    products = ['CCC1OCC1=O(3198)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CsC;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction74',
    reactants = ['OH(5)', 'C=C=C([O])CC(4608)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(6.117e+08,'m^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/OneDe] for rate rule [O_pri_rad;Cd_rad/OneDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction75',
    reactants = ['H(3)', '[CH2]C(=O)C(=O)CC(4620)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction76',
    reactants = ['CH3(17)', 'C=C([O])C(=C)O(4628)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(3.54203e+08,'m^3/(mol*s)'), n=-0.441, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;C_methyl] for rate rule [C_rad/H2/CO;C_methyl]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction77',
    reactants = ['C2H5(29)', '[CH2]C(O)=C=O(4607)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(3.12e+14,'cm^3/(mol*s)','*|/',3), n=-0.5, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [CO_sec_rad;C_rad/H2/Cs] for rate rule [CO_rad/OneDe;C_rad/H2/Cs]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction78',
    reactants = ['H(3)', 'C=C(O)C([O])=CC(4627)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS78',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H/OneDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction79',
    reactants = ['H(3)', '[CH2]CC(=O)C(=C)O(5558)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS79',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction80',
    reactants = ['C2H5CO(71)', 'CH2COH(99)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS80',
    kinetics = Arrhenius(A=(1.81e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;CO_rad] for rate rule [Cd_rad/NonDe;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Recombination"""),
)

reaction(
    label = 'reaction81',
    reactants = ['H(3)', '[CH]=C(O)C(=O)CC(5559)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS81',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction82',
    reactants = ['[CH2]C([O])C(=O)CC(4621)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS82',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction83',
    reactants = ['C=C(O)C([O])[CH]C(4629)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS83',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction84',
    reactants = ['[CH2]CC([O])C(=C)O(4630)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS84',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction85',
    reactants = ['[CH2]C(O)C(=O)[CH]C(4622)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS85',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction86',
    reactants = ['[CH2]C(=O)C([O])CC(3472)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS86',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction87',
    reactants = ['[CH]=C(O)C([O])CC(4631)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS87',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction88',
    reactants = ['[CH2]CC(O)=C([CH2])O(4611)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS88',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radEndo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction89',
    reactants = ['C[CH]C([O])=C(C)O(4612)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS89',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction90',
    reactants = ['[CH2]C([O])=C(O)CC(4613)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS90',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction91',
    reactants = ['[CH]=C(O)[C](O)CC(4632)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS91',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction92',
    reactants = ['[CH2]CC(=O)C([CH2])O(4623)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS92',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction93',
    reactants = ['[CH2]CC([O])=C(C)O(4614)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS93',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction94',
    reactants = ['CH2(S)(23)', 'C=C(O)C(C)=O(3244)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS94',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction95',
    reactants = ['C=C(O)C(=O)CC(4626)'],
    products = ['CCC(=O)C(C)=O(5560)'],
    transitionState = 'TS95',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(204.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction96',
    reactants = ['C=C(O)C(=C)OC(5561)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS96',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction97',
    reactants = ['C=C(O)C(O)=CC(5562)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS97',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction98',
    reactants = ['C=C=C(CC)OO(5563)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS98',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction99',
    reactants = ['[CH]C(O)C(=O)CC(5564)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS99',
    kinetics = Arrhenius(A=(6.14647e+14,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH]
Euclidian distance = 1.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['OH(5)', 'C=[C]C(O)=CC(31155)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS100',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/Cd] for rate rule [O_pri_rad;Cd_rad/Cd]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', 'C=C([O])C(O)=CC(31156)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS101',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['OH(5)', 'C=C(O)[C]=CC(31157)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS102',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad/Cd;Y_rad] for rate rule [Cd_rad/Cd;O_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', 'C=C(O)C([O])=CC(4627)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS103',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [O_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH3(17)', '[CH]=C(O)C(=C)O(8478)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS104',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.81e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 69 used for C_methyl;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2]C=C(O)C(=C)O(31158)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS105',
    kinetics = Arrhenius(A=(3.156e+12,'cm^3/(mol*s)'), n=0.461, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 15 used for C_rad/H2/Cd;H_rad
Exact match found for rate rule [C_rad/H2/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2COH(99)', 'CC=[C]O(31159)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS106',
    kinetics = Arrhenius(A=(556926,'m^3/(mol*s)'), n=0.4, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;Cd_rad] for rate rule [Cd_rad/NonDe;Cd_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Recombination
Ea raised from -2.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', 'C=C(O)C(O)=[C]C(31160)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS107',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', '[CH]=C(O)C(O)=CC(31161)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS108',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([O])C(O)=CC(31162)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS109',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C(O)C([O])[CH]C(4629)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS110',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]CC(O)=C([CH2])O(4611)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS111',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(O)C(=O)[CH]C(4622)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS112',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(O)C(O)=[C]C(31163)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS113',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C([O])C(O)[CH]C(11048)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS114',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=C(O)C(O)[CH]C(16341)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS115',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C[CH]C([O])=C(C)O(4612)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS116',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C[C]=C(O)[C](C)O(31164)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS117',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([O])=C(O)CC(4613)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS118',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=C(O)[C](O)CC(4632)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS119',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C=C(O)C([CH2])O(31165)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS120',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C=C(O)[C](C)O(31166)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS121',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['CH2(S)(23)', 'C=C(O)C(=C)O(5772)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS122',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=C(O)C(O)=CC(5562)'],
    products = ['C=CC(O)=C(C)O(31167)'],
    transitionState = 'TS123',
    kinetics = Arrhenius(A=(1.02873e+09,'s^-1'), n=1.23767, Ea=(163.714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH3_1;CdH2_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=C(O)C(O)=CC(5562)'],
    products = ['CC=C(O)C(C)=O(31168)'],
    transitionState = 'TS124',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(204.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction52',
    reactants = ['C=C(O)C(O)[C]C(31169)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS125',
    kinetics = Arrhenius(A=(6.14647e+14,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CsJ2-C;CsJ2C;CH]
Euclidian distance = 0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]C(O)C(O)=CC(31170)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS126',
    kinetics = Arrhenius(A=(6.14647e+14,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH]
Euclidian distance = 1.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction54',
    reactants = ['C=C(O)C(O)=CC(5562)'],
    products = ['CC1CC(O)=C1O(31171)'],
    transitionState = 'TS127',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH2_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

network(
    label = '1104',
    isomers = [
        'CCC1OCC=1O(3474)',
        'CCC1OCC1=O(3198)',
        '[CH2]C(O)=C([O])CC(4557)',
        'C=C(O)C(=O)CC(4626)',
        'C=C(O)C(O)=CC(5562)',
    ],
    reactants = [
        ('C2H5CHO(70)', 'CH2CO(28)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1104',
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

