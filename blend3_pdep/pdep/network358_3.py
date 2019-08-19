species(
    label = 'S(542)(541)',
    structure = SMILES('C=CCC([O])=CC=CO'),
    E0 = (-145.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.850846,'amu*angstrom^2'), symmetry=1, barrier=(19.5626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853462,'amu*angstrom^2'), symmetry=1, barrier=(19.6228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.851928,'amu*angstrom^2'), symmetry=1, barrier=(19.5875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.847732,'amu*angstrom^2'), symmetry=1, barrier=(19.491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4631.72,'J/mol'), sigma=(7.29257,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=723.46 K, Pc=27.1 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.730394,0.0830463,-2.84649e-05,-4.53675e-08,3.03046e-11,-17252.6,33.8563], Tmin=(100,'K'), Tmax=(916.952,'K')), NASAPolynomial(coeffs=[28.6987,0.0101175,1.28803e-07,-2.08027e-10,1.19019e-14,-24980.7,-118.282], Tmin=(916.952,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'S(568)(567)',
    structure = SMILES('C=CCC(O)=CC=C[O]'),
    E0 = (-141.367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.905591,'amu*angstrom^2'), symmetry=1, barrier=(20.8213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.905117,'amu*angstrom^2'), symmetry=1, barrier=(20.8104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.905747,'amu*angstrom^2'), symmetry=1, barrier=(20.8249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4631.72,'J/mol'), sigma=(7.29257,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=723.46 K, Pc=27.1 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.67876,0.0799789,-1.73132e-05,-5.65189e-08,3.36965e-11,-16812.9,33.5796], Tmin=(100,'K'), Tmax=(927.692,'K')), NASAPolynomial(coeffs=[29.0993,0.0105876,-5.19111e-07,-2.61725e-11,-3.81472e-15,-24876.8,-121.532], Tmin=(927.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-141.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'S(781)(780)',
    structure = SMILES('[CH2]C=CC(O)=CC=CO'),
    E0 = (-194.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3000,3100,440,815,1455,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16475,'amu*angstrom^2'), symmetry=1, barrier=(26.7799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16456,'amu*angstrom^2'), symmetry=1, barrier=(26.7756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16271,'amu*angstrom^2'), symmetry=1, barrier=(26.7329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16313,'amu*angstrom^2'), symmetry=1, barrier=(26.7426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1629,'amu*angstrom^2'), symmetry=1, barrier=(26.7374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4704.45,'J/mol'), sigma=(7.35595,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=734.82 K, Pc=26.82 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09218,0.0869107,-2.13076e-05,-6.58298e-08,4.12007e-11,-23241.6,31.3771], Tmin=(100,'K'), Tmax=(903.05,'K')), NASAPolynomial(coeffs=[33.3096,0.00291046,4.63809e-06,-1.13343e-09,7.7052e-14,-32243.1,-146.521], Tmin=(903.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-194.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'H(3)(3)',
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
    label = 'S(778)(777)',
    structure = SMILES('C=CC=C(O)C=CC=O'),
    E0 = (-191.688,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.95864,'amu*angstrom^2'), symmetry=1, barrier=(22.041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.961854,'amu*angstrom^2'), symmetry=1, barrier=(22.1149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958529,'amu*angstrom^2'), symmetry=1, barrier=(22.0385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958395,'amu*angstrom^2'), symmetry=1, barrier=(22.0354,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4437.48,'J/mol'), sigma=(6.82012,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=693.12 K, Pc=31.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.785439,0.0927288,-8.19227e-05,2.91418e-08,-1.90073e-12,-22871.5,30.1277], Tmin=(100,'K'), Tmax=(1031.1,'K')), NASAPolynomial(coeffs=[23.4878,0.0201354,-7.69741e-06,1.43995e-09,-1.03445e-13,-29023.8,-93.2806], Tmin=(1031.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-191.688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'S(779)(778)',
    structure = SMILES('[CH2]C=CC(O)=CC=C[O]'),
    E0 = (-53.4996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3000,3100,440,815,1455,1000,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.25513,'amu*angstrom^2'), symmetry=1, barrier=(28.858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2548,'amu*angstrom^2'), symmetry=1, barrier=(28.8504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25508,'amu*angstrom^2'), symmetry=1, barrier=(28.8569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2556,'amu*angstrom^2'), symmetry=1, barrier=(28.8688,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4704.45,'J/mol'), sigma=(7.35595,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=734.82 K, Pc=26.82 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.677322,0.0799229,-1.79714e-05,-5.92811e-08,3.61766e-11,-6244.68,31.2969], Tmin=(100,'K'), Tmax=(912.898,'K')), NASAPolynomial(coeffs=[29.9465,0.0067817,1.90948e-06,-5.54001e-10,3.54203e-14,-14379.5,-127.586], Tmin=(912.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.4996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=COJ)"""),
)

species(
    label = 'S(459)(458)',
    structure = SMILES('C=CCC(=O)C=CC=O'),
    E0 = (-165.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,368.992,372.476,372.574],'cm^-1')),
        HinderedRotor(inertia=(0.0221847,'amu*angstrom^2'), symmetry=1, barrier=(2.09348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682691,'amu*angstrom^2'), symmetry=1, barrier=(15.6964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16382,'amu*angstrom^2'), symmetry=1, barrier=(15.6983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161514,'amu*angstrom^2'), symmetry=1, barrier=(15.6953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4218.9,'J/mol'), sigma=(6.48785,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=658.98 K, Pc=35.05 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39102,0.06549,-3.82654e-05,8.80535e-09,-7.19387e-13,-19784.4,27.477], Tmin=(100,'K'), Tmax=(2264.17,'K')), NASAPolynomial(coeffs=[35.3576,0.0135168,-9.15587e-06,1.8014e-09,-1.19076e-13,-37224.9,-168.697], Tmin=(2264.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=C=CC(O)=CC=CO(3156)',
    structure = SMILES('C=C=CC(O)=CC=CO'),
    E0 = (-136.414,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15902,'amu*angstrom^2'), symmetry=1, barrier=(26.6482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15901,'amu*angstrom^2'), symmetry=1, barrier=(26.648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.156,'amu*angstrom^2'), symmetry=1, barrier=(26.5787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1564,'amu*angstrom^2'), symmetry=1, barrier=(26.5879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.19335,0.117059,-0.000127922,6.37793e-08,-1.15109e-11,-16112.3,36.4375], Tmin=(100,'K'), Tmax=(1633.34,'K')), NASAPolynomial(coeffs=[32.3324,0.00127156,4.84942e-06,-1.20321e-09,8.6965e-14,-23877.6,-140.629], Tmin=(1633.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'OH(5)(5)',
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
    label = 'C=CC=C=CC=CO(5318)',
    structure = SMILES('C=CC=C=CC=CO'),
    E0 = (77.4792,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,540,610,2055,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.19048,'amu*angstrom^2'), symmetry=1, barrier=(27.3715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19144,'amu*angstrom^2'), symmetry=1, barrier=(27.3934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19031,'amu*angstrom^2'), symmetry=1, barrier=(27.3677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.27045,0.0725244,-1.4941e-05,-5.58318e-08,3.36067e-11,9492.31,26.6784], Tmin=(100,'K'), Tmax=(915.203,'K')), NASAPolynomial(coeffs=[27.9977,0.00559915,1.94256e-06,-5.27506e-10,3.29712e-14,1946.67,-120.15], Tmin=(915.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.4792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC=C(O)C=C=CO(3134)',
    structure = SMILES('C=CC=C(O)C=C=CO'),
    E0 = (-136.414,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15902,'amu*angstrom^2'), symmetry=1, barrier=(26.6482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15901,'amu*angstrom^2'), symmetry=1, barrier=(26.648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.156,'amu*angstrom^2'), symmetry=1, barrier=(26.5787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1564,'amu*angstrom^2'), symmetry=1, barrier=(26.5879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.19335,0.117059,-0.000127922,6.37793e-08,-1.15109e-11,-16112.3,36.4375], Tmin=(100,'K'), Tmax=(1633.34,'K')), NASAPolynomial(coeffs=[32.3324,0.00127156,4.84942e-06,-1.20321e-09,8.6965e-14,-23877.6,-140.629], Tmin=(1633.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]C=C[C]=CC=CO(5308)',
    structure = SMILES('[CH2]C=C[C]=CC=CO'),
    E0 = (217.926,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.41426,'amu*angstrom^2'), symmetry=1, barrier=(32.5167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40694,'amu*angstrom^2'), symmetry=1, barrier=(32.3483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41547,'amu*angstrom^2'), symmetry=1, barrier=(32.5445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40627,'amu*angstrom^2'), symmetry=1, barrier=(32.333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.107114,0.0704822,-1.1406e-05,-5.8112e-08,3.47059e-11,26377,28.1951], Tmin=(100,'K'), Tmax=(899.692,'K')), NASAPolynomial(coeffs=[26.1293,0.00824228,1.65446e-06,-5.75265e-10,4.08864e-14,19454.1,-107.844], Tmin=(899.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C[CH]C(=O)C=C[CH]O(2830)',
    structure = SMILES('[CH2]C=CC([O])=CC=CO'),
    E0 = (-57.1574,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3000,3100,440,815,1455,1000,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15639,'amu*angstrom^2'), symmetry=1, barrier=(26.5877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14897,'amu*angstrom^2'), symmetry=1, barrier=(26.417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15684,'amu*angstrom^2'), symmetry=1, barrier=(26.598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15,'amu*angstrom^2'), symmetry=1, barrier=(26.4407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.727654,0.0829658,-2.89863e-05,-4.84005e-08,3.29517e-11,-6684.47,31.5694], Tmin=(100,'K'), Tmax=(901.662,'K')), NASAPolynomial(coeffs=[29.5731,0.00626649,2.583e-06,-7.41832e-10,5.1628e-14,-14495.1,-124.489], Tmin=(901.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.1574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC=C(O)C=C[CH2](2092)',
    structure = SMILES('[CH]=CC=C(O)C=C[CH2]'),
    E0 = (260.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.14974,'amu*angstrom^2'), symmetry=1, barrier=(26.4348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15689,'amu*angstrom^2'), symmetry=1, barrier=(26.5992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15828,'amu*angstrom^2'), symmetry=1, barrier=(26.6312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15349,'amu*angstrom^2'), symmetry=1, barrier=(26.5209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.122173,0.0746292,-3.35115e-05,-2.68416e-08,2.09011e-11,31545.5,28.4372], Tmin=(100,'K'), Tmax=(919.333,'K')), NASAPolynomial(coeffs=[23.6744,0.0132278,-2.07945e-06,2.21195e-10,-1.59779e-14,25389.5,-94.0506], Tmin=(919.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[CH2](891)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,229.711,230.18,230.787],'cm^-1')),
        HinderedRotor(inertia=(1.33306,'amu*angstrom^2'), symmetry=1, barrier=(50.5153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,44095,-3.44606], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = 'O[C]=CC=CO(5309)',
    structure = SMILES('O[C]=CC=CO'),
    E0 = (-83.9199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.09535,'amu*angstrom^2'), symmetry=1, barrier=(25.1843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09523,'amu*angstrom^2'), symmetry=1, barrier=(25.1814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09648,'amu*angstrom^2'), symmetry=1, barrier=(25.2103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803519,0.0503097,1.529e-06,-6.70039e-08,3.84255e-11,-9959.14,22.3195], Tmin=(100,'K'), Tmax=(888.056,'K')), NASAPolynomial(coeffs=[26.7325,-0.0104688,9.5814e-06,-2.02664e-09,1.39785e-13,-16773.1,-112.13], Tmin=(888.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.9199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = 'C=C[C]=C(O)C=C[CH]O(3123)',
    structure = SMILES('[CH2]C=[C]C(O)=CC=CO'),
    E0 = (4.03321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.31565,'amu*angstrom^2'), symmetry=1, barrier=(30.2493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31479,'amu*angstrom^2'), symmetry=1, barrier=(30.2296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31591,'amu*angstrom^2'), symmetry=1, barrier=(30.2554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31725,'amu*angstrom^2'), symmetry=1, barrier=(30.2862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31664,'amu*angstrom^2'), symmetry=1, barrier=(30.2721,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.99783,0.114806,-0.000124446,6.27452e-08,-1.14742e-11,770.755,37.8274], Tmin=(100,'K'), Tmax=(1617.68,'K')), NASAPolynomial(coeffs=[29.9632,0.00474118,4.09789e-06,-1.14409e-09,8.61908e-14,-6156.03,-125.493], Tmin=(1617.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.03321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'CHCHOH(49)(49)',
    structure = SMILES('[CH]=CO'),
    E0 = (120.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,458.887,459.577,459.893,460.324],'cm^-1')),
        HinderedRotor(inertia=(0.000141718,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1868.27,'J/mol'), sigma=(4.162,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99181,0.0473233,-6.66059e-05,4.68997e-08,-1.30686e-11,15200.1,31.4259], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[8.78246,0.00524797,-1.71857e-06,2.59722e-10,-1.48227e-14,12883.6,-21.0851], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(298,'K'), Tmax=(6000,'K'), E0=(120.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCHOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH]=C(O)C=C[CH2](3149)',
    structure = SMILES('[CH]=C(O)C=C[CH2]'),
    E0 = (209.155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(1.22226,'amu*angstrom^2'), symmetry=1, barrier=(28.1022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22654,'amu*angstrom^2'), symmetry=1, barrier=(28.2007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23214,'amu*angstrom^2'), symmetry=1, barrier=(28.3293,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17422,0.0501348,-1.56149e-05,-2.77912e-08,1.81652e-11,25268.4,21.1489], Tmin=(100,'K'), Tmax=(912.551,'K')), NASAPolynomial(coeffs=[17.8093,0.00925475,-1.07841e-06,5.99114e-11,-4.11581e-15,20898.5,-64.892], Tmin=(912.551,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC=C(O)[C]=C[CH]O(3113)',
    structure = SMILES('[CH2]C=CC(O)=[C]C=CO'),
    E0 = (4.03321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.31565,'amu*angstrom^2'), symmetry=1, barrier=(30.2493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31479,'amu*angstrom^2'), symmetry=1, barrier=(30.2296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31591,'amu*angstrom^2'), symmetry=1, barrier=(30.2554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31725,'amu*angstrom^2'), symmetry=1, barrier=(30.2862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31664,'amu*angstrom^2'), symmetry=1, barrier=(30.2721,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.99783,0.114806,-0.000124446,6.27452e-08,-1.14742e-11,770.755,37.8274], Tmin=(100,'K'), Tmax=(1617.68,'K')), NASAPolynomial(coeffs=[29.9632,0.00474118,4.09789e-06,-1.14409e-09,8.61908e-14,-6156.03,-125.493], Tmin=(1617.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.03321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=CC(O)=C[C]=CO(3147)',
    structure = SMILES('[CH2]C=CC(O)=C[C]=CO'),
    E0 = (4.03321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.31565,'amu*angstrom^2'), symmetry=1, barrier=(30.2493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31479,'amu*angstrom^2'), symmetry=1, barrier=(30.2296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31591,'amu*angstrom^2'), symmetry=1, barrier=(30.2554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31725,'amu*angstrom^2'), symmetry=1, barrier=(30.2862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31664,'amu*angstrom^2'), symmetry=1, barrier=(30.2721,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.99783,0.114806,-0.000124446,6.27452e-08,-1.14742e-11,770.755,37.8274], Tmin=(100,'K'), Tmax=(1617.68,'K')), NASAPolynomial(coeffs=[29.9632,0.00474118,4.09789e-06,-1.14409e-09,8.61908e-14,-6156.03,-125.493], Tmin=(1617.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.03321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C]C=C(O)C=C[CH]O(3126)',
    structure = SMILES('C=[C]C=C(O)[CH]C=CO'),
    E0 = (33.0826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.15484,0.117963,-0.000127981,6.38524e-08,-1.15816e-11,4270.48,37.8751], Tmin=(100,'K'), Tmax=(1616.55,'K')), NASAPolynomial(coeffs=[32.173,0.00314763,3.98111e-06,-1.05405e-09,7.76993e-14,-3571.17,-138.456], Tmin=(1616.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.0826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=CC(O)=CC=[C]O(3145)',
    structure = SMILES('[CH2]C=CC(O)=CC=[C]O'),
    E0 = (44.7819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08975,'amu*angstrom^2'), symmetry=1, barrier=(25.0555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09009,'amu*angstrom^2'), symmetry=1, barrier=(25.0633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08867,'amu*angstrom^2'), symmetry=1, barrier=(25.0308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08854,'amu*angstrom^2'), symmetry=1, barrier=(25.0277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08878,'amu*angstrom^2'), symmetry=1, barrier=(25.0332,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.48897,0.107334,-0.00011244,5.53308e-08,-9.96048e-12,5650.15,39.5383], Tmin=(100,'K'), Tmax=(1624.19,'K')), NASAPolynomial(coeffs=[28.518,0.00715004,2.08383e-06,-7.06747e-10,5.49211e-14,-1280.1,-115.419], Tmin=(1624.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.7819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C[CH]C(O)=CC=CO(3163)',
    structure = SMILES('[CH]C=CC(O)=CC=CO'),
    E0 = (57.6665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,350,440,435,1725,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.35962,0.118672,-0.00012244,5.8992e-08,-1.03779e-11,7237.97,39.1012], Tmin=(100,'K'), Tmax=(1671.23,'K')), NASAPolynomial(coeffs=[31.2725,0.00779424,2.19709e-06,-7.47189e-10,5.75046e-14,-429.253,-134.071], Tmin=(1671.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.6665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C2H3(28)(29)',
    structure = SMILES('[CH]=C'),
    E0 = (286.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,677.08,1086.68,3788.01],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH]C(O)=CC=CO(5321)',
    structure = SMILES('[CH]C(O)=CC=CO'),
    E0 = (5.22887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.86599,'amu*angstrom^2'), symmetry=1, barrier=(42.9028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.86592,'amu*angstrom^2'), symmetry=1, barrier=(42.9011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.86574,'amu*angstrom^2'), symmetry=1, barrier=(42.8971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.86606,'amu*angstrom^2'), symmetry=1, barrier=(42.9044,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00887873,0.066945,-9.66871e-06,-6.19878e-08,3.68411e-11,793.217,25.2186], Tmin=(100,'K'), Tmax=(900.357,'K')), NASAPolynomial(coeffs=[28.0366,0.00033086,4.71042e-06,-1.1073e-09,7.56665e-14,-6607.16,-120.192], Tmin=(900.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.22887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC=[C]O(3098)',
    structure = SMILES('C=CC=[C]O'),
    E0 = (124.873,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.01736,'amu*angstrom^2'), symmetry=1, barrier=(23.3911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02157,'amu*angstrom^2'), symmetry=1, barrier=(23.4878,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85281,0.034453,6.0896e-06,-4.53851e-08,2.39475e-11,15108,18.6912], Tmin=(100,'K'), Tmax=(905.911,'K')), NASAPolynomial(coeffs=[16.7736,0.00282296,1.74834e-06,-4.54027e-10,3.0294e-14,10999.1,-59.5759], Tmin=(905.911,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C[CH]O(3569)',
    structure = SMILES('[CH]C=CO'),
    E0 = (203.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.16611,'amu*angstrom^2'), symmetry=1, barrier=(49.8031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.16414,'amu*angstrom^2'), symmetry=1, barrier=(49.7578,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43518,0.0321429,-2.11365e-05,7.43426e-09,-1.09824e-12,24561.7,14.1742], Tmin=(100,'K'), Tmax=(1529.28,'K')), NASAPolynomial(coeffs=[8.00322,0.0175792,-6.85177e-06,1.2071e-09,-8.02562e-14,22858.7,-15.0537], Tmin=(1529.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]CC=C(O)C=C=CO(5313)',
    structure = SMILES('[CH2]CC=C(O)C=C=CO'),
    E0 = (-41.6103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.60792,0.114483,-0.000122316,6.13507e-08,-1.1355e-11,-4740.17,38.1927], Tmin=(100,'K'), Tmax=(1543.31,'K')), NASAPolynomial(coeffs=[30.8587,0.00753571,1.27035e-06,-5.19136e-10,4.16316e-14,-12663.6,-129.99], Tmin=(1543.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.6103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C=C(O)C=CCO(5314)',
    structure = SMILES('C=C[C]=C(O)C=CCO'),
    E0 = (-58.3657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3389,0.105689,-0.000116258,6.43532e-08,-1.37036e-11,-6817.49,34.1901], Tmin=(100,'K'), Tmax=(1213.68,'K')), NASAPolynomial(coeffs=[23.4501,0.019993,-5.4056e-06,7.49414e-10,-4.32865e-14,-12540.3,-88.9903], Tmin=(1213.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.3657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C=CC(O)=CC[CH]O(5335)',
    structure = SMILES('C=C=CC(O)=CC[CH]O'),
    E0 = (-17.5937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25447,0.110775,-0.00012987,7.67191e-08,-1.76162e-11,-1922.63,33.9285], Tmin=(100,'K'), Tmax=(1073.45,'K')), NASAPolynomial(coeffs=[22.004,0.0241085,-8.76656e-06,1.50897e-09,-1.00497e-13,-6916.07,-79.9291], Tmin=(1073.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.5937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCsJOH)"""),
)

species(
    label = 'CC=CC(O)=[C]C=CO(5304)',
    structure = SMILES('CC=CC(O)=[C]C=CO'),
    E0 = (-114.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08872,'amu*angstrom^2'), symmetry=1, barrier=(25.0317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08895,'amu*angstrom^2'), symmetry=1, barrier=(25.0372,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08854,'amu*angstrom^2'), symmetry=1, barrier=(25.0276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08873,'amu*angstrom^2'), symmetry=1, barrier=(25.0322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08899,'amu*angstrom^2'), symmetry=1, barrier=(25.0381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.83028,0.118616,-0.000130647,6.74936e-08,-1.27875e-11,-13440.7,35.6551], Tmin=(100,'K'), Tmax=(1525.12,'K')), NASAPolynomial(coeffs=[31.0012,0.00646554,2.68841e-06,-8.58796e-10,6.72739e-14,-21036.5,-132.913], Tmin=(1525.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=CC1(O)C=CC1O(5315)',
    structure = SMILES('C=C[CH]C1(O)C=CC1O'),
    E0 = (-99.6008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.267443,0.0534453,5.10482e-05,-1.18914e-07,5.31829e-11,-11818.6,34.1344], Tmin=(100,'K'), Tmax=(942.783,'K')), NASAPolynomial(coeffs=[26.352,0.0150488,-2.8522e-06,5.13192e-10,-4.77606e-14,-19949,-107.206], Tmin=(942.783,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.6008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'OC=C[CH]C1(O)C=CC1(5299)',
    structure = SMILES('OC=C[CH]C1(O)C=CC1'),
    E0 = (-135.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0875306,0.0547951,6.89507e-05,-1.54476e-07,7.0557e-11,-16127.6,32.6948], Tmin=(100,'K'), Tmax=(922.868,'K')), NASAPolynomial(coeffs=[32.612,0.00391692,3.97888e-06,-8.67987e-10,4.83832e-14,-26031.9,-143.398], Tmin=(922.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-135.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'O[CH]C1C=C(O)C=CC1(5300)',
    structure = SMILES('O[CH]C1C=C(O)C=CC1'),
    E0 = (-136.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.226519,0.0767264,-3.09742e-05,-2.45617e-08,1.77316e-11,-16286.6,27.9345], Tmin=(100,'K'), Tmax=(957.049,'K')), NASAPolynomial(coeffs=[21.9967,0.0226876,-7.15828e-06,1.25666e-09,-9.05322e-14,-22319.3,-87.598], Tmin=(957.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(CCsJOH)"""),
)

species(
    label = 'C[C]=CC(O)=CC=CO(5301)',
    structure = SMILES('C[C]=CC(O)=CC=CO'),
    E0 = (-75.1762,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.978521,'amu*angstrom^2'), symmetry=1, barrier=(22.4981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.981693,'amu*angstrom^2'), symmetry=1, barrier=(22.571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.9843,'amu*angstrom^2'), symmetry=1, barrier=(22.631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.980911,'amu*angstrom^2'), symmetry=1, barrier=(22.5531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.984093,'amu*angstrom^2'), symmetry=1, barrier=(22.6262,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.83074,0.117097,-0.000127337,6.45172e-08,-1.1985e-11,-8767.19,36.5414], Tmin=(100,'K'), Tmax=(1553.91,'K')), NASAPolynomial(coeffs=[31.5379,0.00554594,2.62432e-06,-7.98329e-10,6.11594e-14,-16661.7,-135.45], Tmin=(1553.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.1762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S)"""),
)

species(
    label = 'CC=[C]C(O)=CC=CO(5302)',
    structure = SMILES('CC=[C]C(O)=CC=CO'),
    E0 = (-114.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08909,'amu*angstrom^2'), symmetry=1, barrier=(25.0403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08874,'amu*angstrom^2'), symmetry=1, barrier=(25.0323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08874,'amu*angstrom^2'), symmetry=1, barrier=(25.0322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08863,'amu*angstrom^2'), symmetry=1, barrier=(25.0298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08895,'amu*angstrom^2'), symmetry=1, barrier=(25.0372,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.83028,0.118616,-0.000130647,6.74936e-08,-1.27875e-11,-13440.7,35.6551], Tmin=(100,'K'), Tmax=(1525.12,'K')), NASAPolynomial(coeffs=[31.0012,0.00646554,2.68841e-06,-8.58796e-10,6.72739e-14,-21036.5,-132.913], Tmin=(1525.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C)"""),
)

species(
    label = 'CC=CC([O])=CC=CO(5303)',
    structure = SMILES('CC=CC([O])=CC=CO'),
    E0 = (-175.213,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.952452,'amu*angstrom^2'), symmetry=1, barrier=(21.8988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.955461,'amu*angstrom^2'), symmetry=1, barrier=(21.9679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.949807,'amu*angstrom^2'), symmetry=1, barrier=(21.8379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.95283,'amu*angstrom^2'), symmetry=1, barrier=(21.9074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52374,0.109463,-0.000112724,5.49949e-08,-9.87934e-12,-20809.1,36.4775], Tmin=(100,'K'), Tmax=(1610.97,'K')), NASAPolynomial(coeffs=[29.0115,0.00926663,9.57572e-07,-4.86553e-10,3.98668e-14,-28128.4,-121.881], Tmin=(1610.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-175.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC=CC(O)=C[C]=CO(5305)',
    structure = SMILES('CC=CC(O)=C[C]=CO'),
    E0 = (-114.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08909,'amu*angstrom^2'), symmetry=1, barrier=(25.0403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08874,'amu*angstrom^2'), symmetry=1, barrier=(25.0323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08874,'amu*angstrom^2'), symmetry=1, barrier=(25.0322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08863,'amu*angstrom^2'), symmetry=1, barrier=(25.0298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08895,'amu*angstrom^2'), symmetry=1, barrier=(25.0372,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.83028,0.118616,-0.000130647,6.74936e-08,-1.27875e-11,-13440.7,35.6551], Tmin=(100,'K'), Tmax=(1525.12,'K')), NASAPolynomial(coeffs=[31.0012,0.00646554,2.68841e-06,-8.58796e-10,6.72739e-14,-21036.5,-132.913], Tmin=(1525.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C)"""),
)

species(
    label = 'CC=CC(O)=CC=[C]O(5306)',
    structure = SMILES('CC=CC(O)=CC=[C]O'),
    E0 = (-73.2738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.905077,'amu*angstrom^2'), symmetry=1, barrier=(20.8095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.905303,'amu*angstrom^2'), symmetry=1, barrier=(20.8147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.904884,'amu*angstrom^2'), symmetry=1, barrier=(20.8051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.904836,'amu*angstrom^2'), symmetry=1, barrier=(20.804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.90483,'amu*angstrom^2'), symmetry=1, barrier=(20.8038,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.31655,0.111097,-0.00011852,5.99688e-08,-1.12419e-11,-8561.51,37.3477], Tmin=(100,'K'), Tmax=(1515.03,'K')), NASAPolynomial(coeffs=[29.3737,0.00913408,5.43367e-07,-3.93517e-10,3.38647e-14,-16064.3,-121.776], Tmin=(1515.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.2738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = 'CC=CC(O)=CC=C[O](5307)',
    structure = SMILES('CC=CC(O)=CC=C[O]'),
    E0 = (-171.555,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01814,'amu*angstrom^2'), symmetry=1, barrier=(23.409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02019,'amu*angstrom^2'), symmetry=1, barrier=(23.4561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01916,'amu*angstrom^2'), symmetry=1, barrier=(23.4325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01658,'amu*angstrom^2'), symmetry=1, barrier=(23.3731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.961686,0.0888998,-4.16337e-05,-3.25111e-08,2.56446e-11,-20435.9,30.7585], Tmin=(100,'K'), Tmax=(922.143,'K')), NASAPolynomial(coeffs=[29.21,0.0108513,-6.08858e-07,-4.49801e-11,1.59042e-16,-28246.6,-124.536], Tmin=(922.143,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'OC=CC=C(O)C1[CH]C1(5310)',
    structure = SMILES('OC=CC=C(O)C1[CH]C1'),
    E0 = (-31.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.556145,0.0709707,1.86112e-05,-1.03215e-07,5.30793e-11,-3578.5,32.8848], Tmin=(100,'K'), Tmax=(915.177,'K')), NASAPolynomial(coeffs=[32.6726,0.00268234,4.42138e-06,-1.00778e-09,6.3014e-14,-12882.8,-142.084], Tmin=(915.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(cyclopropane)"""),
)

species(
    label = 'OC=CC1CC=C[C]1O(5311)',
    structure = SMILES('OC=CC1C[CH]C=C1O'),
    E0 = (-202.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.220593,0.0645236,2.9363e-05,-1.05739e-07,5.13526e-11,-24158.6,26.2629], Tmin=(100,'K'), Tmax=(926.848,'K')), NASAPolynomial(coeffs=[29.4519,0.00936131,6.64054e-07,-2.40581e-10,8.04796e-15,-32790,-131.526], Tmin=(926.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-202.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopentene) + radical(cyclopentene-allyl)"""),
)

species(
    label = 'OC1C=CCC(O)[CH]C=1(5312)',
    structure = SMILES('OC1C=CCC(O)[CH]C=1'),
    E0 = (-196.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19575,0.098758,-7.44551e-05,1.02453e-08,7.4172e-12,-23444,13.6358], Tmin=(100,'K'), Tmax=(964.136,'K')), NASAPolynomial(coeffs=[26.4444,0.0194711,-6.1546e-06,1.0857e-09,-7.86687e-14,-30418.4,-127.232], Tmin=(964.136,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C[CH]C(=O)CC=CO(5175)',
    structure = SMILES('C=CC=C([O])CC=CO'),
    E0 = (-145.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.850846,'amu*angstrom^2'), symmetry=1, barrier=(19.5626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853462,'amu*angstrom^2'), symmetry=1, barrier=(19.6228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.851928,'amu*angstrom^2'), symmetry=1, barrier=(19.5875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.847732,'amu*angstrom^2'), symmetry=1, barrier=(19.491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.730394,0.0830463,-2.84649e-05,-4.53675e-08,3.03046e-11,-17252.6,33.8563], Tmin=(100,'K'), Tmax=(916.952,'K')), NASAPolynomial(coeffs=[28.6987,0.0101175,1.28803e-07,-2.08027e-10,1.19019e-14,-24980.7,-118.282], Tmin=(916.952,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C[CH]C(O)=CCC=O(5207)',
    structure = SMILES('[CH2]C=CC(O)=CCC=O'),
    E0 = (-173.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.206695,0.0759738,-3.0766e-05,-1.90534e-08,1.34326e-11,-20668.8,32.0085], Tmin=(100,'K'), Tmax=(1017.89,'K')), NASAPolynomial(coeffs=[21.2313,0.0271194,-1.09247e-05,2.10888e-09,-1.54219e-13,-26866.4,-80.8033], Tmin=(1017.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'CH2(17)(18)',
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
    label = '[CH]=CC(O)=CC=CO(5316)',
    structure = SMILES('[CH]=CC(O)=CC=CO'),
    E0 = (-29.8962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3120,650,792.5,1650,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.11825,'amu*angstrom^2'), symmetry=1, barrier=(25.7108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1187,'amu*angstrom^2'), symmetry=1, barrier=(25.7212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11965,'amu*angstrom^2'), symmetry=1, barrier=(25.743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1198,'amu*angstrom^2'), symmetry=1, barrier=(25.7464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.76881,0.10662,-0.000119297,5.99045e-08,-1.07704e-11,-3315.44,34.4596], Tmin=(100,'K'), Tmax=(1665.1,'K')), NASAPolynomial(coeffs=[29.863,-0.00238276,6.47588e-06,-1.49339e-09,1.05746e-13,-9938.68,-126.865], Tmin=(1665.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.8962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1C(O)=CC1[CH]O(5317)',
    structure = SMILES('C=CC1C(O)=CC1[CH]O'),
    E0 = (-9.67099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.590869,0.0880533,-6.39876e-05,8.80419e-09,6.20514e-12,-986.21,31.0322], Tmin=(100,'K'), Tmax=(962.527,'K')), NASAPolynomial(coeffs=[22.3576,0.0221383,-7.16479e-06,1.23815e-09,-8.67738e-14,-6768.26,-85.8918], Tmin=(962.527,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.67099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CCsJOH)"""),
)

species(
    label = 'C=[C]CC(O)=CC=CO(5164)',
    structure = SMILES('C=[C]CC(O)=CC=CO'),
    E0 = (-44.9877,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11773,0.0914691,-4.51463e-05,-3.40714e-08,2.78101e-11,-5206.81,34.2198], Tmin=(100,'K'), Tmax=(912.719,'K')), NASAPolynomial(coeffs=[31.2146,0.0063148,1.8754e-06,-5.43303e-10,3.53901e-14,-13464,-131.714], Tmin=(912.719,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.9877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CCC(O)=[C]C=CO(5163)',
    structure = SMILES('C=CCC(O)=[C]C=CO'),
    E0 = (-83.8341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.91979,0.113773,-0.000119256,5.85564e-08,-1.05019e-11,-9800.64,39.8351], Tmin=(100,'K'), Tmax=(1633.63,'K')), NASAPolynomial(coeffs=[30.1395,0.00701243,2.47436e-06,-7.93219e-10,6.08359e-14,-17157.4,-125.34], Tmin=(1633.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.8341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CCC(O)=CC=CO(5166)',
    structure = SMILES('[CH]=CCC(O)=CC=CO'),
    E0 = (-35.7334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3120,650,792.5,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.936941,'amu*angstrom^2'), symmetry=1, barrier=(21.5421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.93698,'amu*angstrom^2'), symmetry=1, barrier=(21.543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937213,'amu*angstrom^2'), symmetry=1, barrier=(21.5484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937011,'amu*angstrom^2'), symmetry=1, barrier=(21.5437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.936955,'amu*angstrom^2'), symmetry=1, barrier=(21.5424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16259,0.0904334,-3.71171e-05,-4.59137e-08,3.28902e-11,-4090.24,34.3105], Tmin=(100,'K'), Tmax=(911.534,'K')), NASAPolynomial(coeffs=[32.6806,0.00392106,3.22176e-06,-7.99202e-10,5.22552e-14,-12835.8,-139.957], Tmin=(911.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.7334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC(O)=C[C]=CO(5165)',
    structure = SMILES('C=CCC(O)=C[C]=CO'),
    E0 = (-83.8341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.91979,0.113773,-0.000119256,5.85564e-08,-1.05019e-11,-9800.64,39.8351], Tmin=(100,'K'), Tmax=(1633.63,'K')), NASAPolynomial(coeffs=[30.1395,0.00701243,2.47436e-06,-7.93219e-10,6.08359e-14,-17157.4,-125.34], Tmin=(1633.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.8341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CCC(O)=CC=[C]O(5167)',
    structure = SMILES('C=CCC(O)=CC=[C]O'),
    E0 = (-43.0853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.878689,0.0886789,-4.74582e-05,-2.42101e-08,2.2371e-11,-4989,36.0162], Tmin=(100,'K'), Tmax=(919.262,'K')), NASAPolynomial(coeffs=[28.2698,0.0109304,-6.87976e-07,-4.19822e-11,1.10907e-15,-12422,-113.435], Tmin=(919.262,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.0853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'C=CC1[C](O)C1C=CO(5319)',
    structure = SMILES('C=CC1[C](O)C1C=CO'),
    E0 = (-18.0629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.395087,0.0738719,-6.02134e-06,-6.49185e-08,3.60634e-11,-1993.13,32.3794], Tmin=(100,'K'), Tmax=(925.877,'K')), NASAPolynomial(coeffs=[27.5784,0.0118623,-8.8996e-07,2.75865e-11,-6.95038e-15,-9695.26,-114.042], Tmin=(925.877,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.0629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C2CsJOH)"""),
)

species(
    label = 'C=CC1C(O)=C[CH]C1O(5320)',
    structure = SMILES('C=CC1C(O)=C[CH]C1O'),
    E0 = (-176.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.335665,0.0756104,-1.77008e-05,-4.28742e-08,2.49318e-11,-21043.6,26.7868], Tmin=(100,'K'), Tmax=(958.559,'K')), NASAPolynomial(coeffs=[24.2455,0.0207149,-6.40942e-06,1.16446e-09,-8.75966e-14,-27946.6,-102.189], Tmin=(958.559,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-176.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopentene) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C1C=C(O)C1C=CO(5322)',
    structure = SMILES('[CH2]C1C=C(O)C1C=CO'),
    E0 = (-33.8517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.626487,0.0770909,-4.33206e-06,-7.6259e-08,4.31778e-11,-3881.72,30.552], Tmin=(100,'K'), Tmax=(906.151,'K')), NASAPolynomial(coeffs=[30.4678,0.00623225,3.04869e-06,-8.2286e-10,5.51804e-14,-12243,-131.438], Tmin=(906.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.8517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC=C(O)C[C]=CO(5323)',
    structure = SMILES('C=CC=C(O)C[C]=CO'),
    E0 = (-44.9877,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11773,0.0914691,-4.51463e-05,-3.40714e-08,2.78101e-11,-5206.81,34.2198], Tmin=(100,'K'), Tmax=(912.719,'K')), NASAPolynomial(coeffs=[31.2146,0.0063148,1.8754e-06,-5.43303e-10,3.53901e-14,-13464,-131.714], Tmin=(912.719,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.9877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=C[C]=C(O)CC=CO(5324)',
    structure = SMILES('C=C[C]=C(O)CC=CO'),
    E0 = (-83.8341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.91979,0.113773,-0.000119256,5.85564e-08,-1.05019e-11,-9800.64,39.8351], Tmin=(100,'K'), Tmax=(1633.63,'K')), NASAPolynomial(coeffs=[30.1395,0.00701243,2.47436e-06,-7.93219e-10,6.08359e-14,-17157.4,-125.34], Tmin=(1633.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.8341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC=C(O)CC=[C]O(5325)',
    structure = SMILES('C=CC=C(O)CC=[C]O'),
    E0 = (-43.0853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.878689,0.0886789,-4.74582e-05,-2.42101e-08,2.2371e-11,-4989,36.0162], Tmin=(100,'K'), Tmax=(919.262,'K')), NASAPolynomial(coeffs=[28.2698,0.0109304,-6.87976e-07,-4.19822e-11,1.10907e-15,-12422,-113.435], Tmin=(919.262,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.0853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]C=C(O)CC=CO(5326)',
    structure = SMILES('C=[C]C=C(O)CC=CO'),
    E0 = (-83.8341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.91979,0.113773,-0.000119256,5.85564e-08,-1.05019e-11,-9800.64,39.8351], Tmin=(100,'K'), Tmax=(1633.63,'K')), NASAPolynomial(coeffs=[30.1395,0.00701243,2.47436e-06,-7.93219e-10,6.08359e-14,-17157.4,-125.34], Tmin=(1633.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.8341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC=C(O)CC=C[O](5327)',
    structure = SMILES('C=CC=C(O)CC=C[O]'),
    E0 = (-141.367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.905591,'amu*angstrom^2'), symmetry=1, barrier=(20.8213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.905117,'amu*angstrom^2'), symmetry=1, barrier=(20.8104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.905747,'amu*angstrom^2'), symmetry=1, barrier=(20.8249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.67876,0.0799789,-1.73132e-05,-5.65189e-08,3.36965e-11,-16812.9,33.5796], Tmin=(100,'K'), Tmax=(927.692,'K')), NASAPolynomial(coeffs=[29.0993,0.0105876,-5.19111e-07,-2.61725e-11,-3.81472e-15,-24876.8,-121.532], Tmin=(927.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-141.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CC=C(O)CC=CO(5328)',
    structure = SMILES('[CH]=CC=C(O)CC=CO'),
    E0 = (-35.7334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3120,650,792.5,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.936941,'amu*angstrom^2'), symmetry=1, barrier=(21.5421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.93698,'amu*angstrom^2'), symmetry=1, barrier=(21.543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937213,'amu*angstrom^2'), symmetry=1, barrier=(21.5484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937011,'amu*angstrom^2'), symmetry=1, barrier=(21.5437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.936955,'amu*angstrom^2'), symmetry=1, barrier=(21.5424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16259,0.0904334,-3.71171e-05,-4.59137e-08,3.28902e-11,-4090.24,34.3105], Tmin=(100,'K'), Tmax=(911.534,'K')), NASAPolynomial(coeffs=[32.6806,0.00392106,3.22176e-06,-7.99202e-10,5.22552e-14,-12835.8,-139.957], Tmin=(911.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.7334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CC=C(O)C1[CH]C1O(5329)',
    structure = SMILES('C=CC=C(O)C1[CH]C1O'),
    E0 = (-21.0325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.223179,0.0670638,1.51475e-05,-8.84882e-08,4.47973e-11,-2353.67,33.6436], Tmin=(100,'K'), Tmax=(926.191,'K')), NASAPolynomial(coeffs=[28.5014,0.00961057,3.31662e-07,-1.84507e-10,5.56377e-15,-10531.2,-118.155], Tmin=(926.191,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.0325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1C=C(O)C=CC1O(5330)',
    structure = SMILES('[CH2]C1C=C(O)C=CC1O'),
    E0 = (-124.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0753033,0.0646483,9.81421e-06,-7.25632e-08,3.65598e-11,-14872.4,28.9551], Tmin=(100,'K'), Tmax=(932.688,'K')), NASAPolynomial(coeffs=[23.8317,0.0179351,-3.78666e-06,5.79068e-10,-4.48892e-14,-21703.5,-96.8642], Tmin=(932.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC=C(O)C=CC[O](5196)',
    structure = SMILES('C=CC=C(O)C=CC[O]'),
    E0 = (-31.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.723883,'amu*angstrom^2'), symmetry=1, barrier=(16.6435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.724829,'amu*angstrom^2'), symmetry=1, barrier=(16.6653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.72384,'amu*angstrom^2'), symmetry=1, barrier=(16.6425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.724985,'amu*angstrom^2'), symmetry=1, barrier=(16.6688,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.738134,0.0944042,-9.19376e-05,4.58744e-08,-9.00532e-12,-3628.46,33.3039], Tmin=(100,'K'), Tmax=(1243.5,'K')), NASAPolynomial(coeffs=[20.5358,0.025972,-9.38978e-06,1.61883e-09,-1.07952e-13,-8919.29,-73.9666], Tmin=(1243.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = 'C=CC=C(O)C=[C]CO(5331)',
    structure = SMILES('C=CC=C(O)C=[C]CO'),
    E0 = (-19.5194,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28969,0.103651,-0.000111412,5.97314e-08,-1.23289e-11,-2146.36,34.8934], Tmin=(100,'K'), Tmax=(1232.76,'K')), NASAPolynomial(coeffs=[24.0096,0.0190804,-5.4888e-06,8.16572e-10,-5.00732e-14,-8195.42,-91.6902], Tmin=(1232.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.5194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CC=C(O)[C]=CCO(5332)',
    structure = SMILES('[CH2]C=CC(O)=C=CCO'),
    E0 = (-86.5242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3580,3650,1210,1345,900,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.741509,'amu*angstrom^2'), symmetry=1, barrier=(17.0488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48907,'amu*angstrom^2'), symmetry=1, barrier=(34.2366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.741955,'amu*angstrom^2'), symmetry=1, barrier=(17.059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.741711,'amu*angstrom^2'), symmetry=1, barrier=(17.0534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.741761,'amu*angstrom^2'), symmetry=1, barrier=(17.0546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.651838,0.0931457,-8.42217e-05,3.42389e-08,-3.76125e-12,-10230.8,33.7179], Tmin=(100,'K'), Tmax=(982.666,'K')), NASAPolynomial(coeffs=[20.1484,0.0261734,-9.0036e-06,1.53492e-09,-1.03353e-13,-15173.2,-70.6147], Tmin=(982.666,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.5242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C[CH]C(=O)C=CCO(5187)',
    structure = SMILES('C=CC=C([O])C=CCO'),
    E0 = (-119.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.624113,'amu*angstrom^2'), symmetry=1, barrier=(14.3496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.62427,'amu*angstrom^2'), symmetry=1, barrier=(14.3532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.624294,'amu*angstrom^2'), symmetry=1, barrier=(14.3537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.625075,'amu*angstrom^2'), symmetry=1, barrier=(14.3717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.912223,0.0952927,-9.47025e-05,4.80158e-08,-9.47857e-12,-14191.6,34.5687], Tmin=(100,'K'), Tmax=(1268.54,'K')), NASAPolynomial(coeffs=[21.7634,0.0224595,-7.0051e-06,1.09982e-09,-6.93919e-14,-19837.5,-79.7994], Tmin=(1268.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C]C=C(O)C=CCO(5333)',
    structure = SMILES('C=C=CC(O)=C[CH]CO'),
    E0 = (-80.9749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,540,610,2055,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35607,0.107805,-0.000112821,5.66585e-08,-1.04282e-11,-9537.45,31.9418], Tmin=(100,'K'), Tmax=(1027.39,'K')), NASAPolynomial(coeffs=[24.077,0.0223304,-7.80287e-06,1.34528e-09,-9.12092e-14,-15478.2,-94.9241], Tmin=(1027.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.9749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=CC=C(O)C=CCO(5334)',
    structure = SMILES('[CH]=CC=C(O)C=CCO'),
    E0 = (-10.265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3120,650,792.5,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.787367,'amu*angstrom^2'), symmetry=1, barrier=(18.1031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.7861,'amu*angstrom^2'), symmetry=1, barrier=(18.074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788626,'amu*angstrom^2'), symmetry=1, barrier=(18.1321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.787812,'amu*angstrom^2'), symmetry=1, barrier=(18.1134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.787491,'amu*angstrom^2'), symmetry=1, barrier=(18.106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53009,0.104888,-0.000111203,5.78977e-08,-1.14819e-11,-1021.2,35.6881], Tmin=(100,'K'), Tmax=(1340.21,'K')), NASAPolynomial(coeffs=[25.8138,0.0161047,-3.80549e-06,4.81004e-10,-2.66031e-14,-7706.34,-101.834], Tmin=(1340.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'HCOH(T)(1710)',
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
    label = '[CH]=CC(O)=CC=C(3102)',
    structure = SMILES('[CH]=CC(O)=CC=C'),
    E0 = (178.896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.05134,'amu*angstrom^2'), symmetry=1, barrier=(24.1724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05475,'amu*angstrom^2'), symmetry=1, barrier=(24.2508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05767,'amu*angstrom^2'), symmetry=1, barrier=(24.3179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.259971,0.0677715,-3.5563e-05,-2.01267e-08,1.79606e-11,21664.4,23.7016], Tmin=(100,'K'), Tmax=(916.46,'K')), NASAPolynomial(coeffs=[22.8556,0.00762483,-9.18113e-08,-1.2112e-10,7.16957e-15,15907.1,-92.1531], Tmin=(916.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1CC(=CC=CO)O1(5161)',
    structure = SMILES('[CH2]C1CC(=CC=CO)O1'),
    E0 = (-38.1712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.946137,0.080054,4.28183e-06,-9.91269e-08,5.5528e-11,-4385.64,27.5135], Tmin=(100,'K'), Tmax=(891.504,'K')), NASAPolynomial(coeffs=[35.1894,-0.00229947,8.61252e-06,-1.98589e-09,1.383e-13,-13999,-160.451], Tmin=(891.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-38.1712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=CCC1=CC([CH]O)O1(5162)',
    structure = SMILES('C=CCC1=CC([CH]O)O1'),
    E0 = (23.8081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.335205,0.0753283,-1.88348e-05,-4.21243e-08,2.47714e-11,3037.66,33.2516], Tmin=(100,'K'), Tmax=(960.35,'K')), NASAPolynomial(coeffs=[24.8354,0.0187108,-5.72175e-06,1.05879e-09,-8.1327e-14,-4020.56,-98.7414], Tmin=(960.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.8081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CCsJOH)"""),
)

species(
    label = 'C3H5(89)(88)',
    structure = SMILES('[CH2]C=C'),
    E0 = (156.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.332071,'amu*angstrom^2'), symmetry=1, barrier=(25.4373,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.29606,0.00579326,4.33883e-05,-5.99838e-08,2.33791e-11,18908.2,9.02024], Tmin=(100,'K'), Tmax=(942.201,'K')), NASAPolynomial(coeffs=[8.0689,0.0101832,-2.84768e-06,5.00815e-10,-3.79574e-14,16914.6,-19.5287], Tmin=(942.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O=C=CC=CO(5155)',
    structure = SMILES('O=C=CC=CO'),
    E0 = (-180.379,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(1.07177,'amu*angstrom^2'), symmetry=1, barrier=(24.642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07268,'amu*angstrom^2'), symmetry=1, barrier=(24.6631,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.822541,0.0545621,-5.58541e-05,2.72142e-08,-4.91939e-12,-21567,21.5833], Tmin=(100,'K'), Tmax=(1576.29,'K')), NASAPolynomial(coeffs=[16.5566,0.00530553,-1.0316e-07,-1.19684e-10,1.12857e-14,-25368.2,-57.8079], Tmin=(1576.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-180.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CCO)H) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = '[CH]=CC=C([O])CC=C(2414)',
    structure = SMILES('[CH]C=CC(=O)CC=C'),
    E0 = (292.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.931118,0.0630806,-3.07415e-05,3.99648e-09,3.8936e-13,35345.7,28.997], Tmin=(100,'K'), Tmax=(1661.25,'K')), NASAPolynomial(coeffs=[18.9086,0.033344,-1.61263e-05,3.04128e-09,-2.05787e-13,27503,-72.4864], Tmin=(1661.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC[C]([O])C=CC=O(2769)',
    structure = SMILES('C=CCC([O])=CC=C[O]'),
    E0 = (-3.56209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.847619,'amu*angstrom^2'), symmetry=1, barrier=(19.4884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.846683,'amu*angstrom^2'), symmetry=1, barrier=(19.4669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.847147,'amu*angstrom^2'), symmetry=1, barrier=(19.4776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.318776,0.0760891,-2.51903e-05,-3.88268e-08,2.53343e-11,-255.512,33.7881], Tmin=(100,'K'), Tmax=(932.108,'K')), NASAPolynomial(coeffs=[25.3822,0.0139103,-2.55479e-06,3.60804e-10,-2.88539e-14,-7136.82,-99.6089], Tmin=(932.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.56209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[O][C]=CC=CO(5156)',
    structure = SMILES('O=[C]C=C[CH]O'),
    E0 = (7.41646,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.848579,'amu*angstrom^2'), symmetry=1, barrier=(19.5105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.851297,'amu*angstrom^2'), symmetry=1, barrier=(19.573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.848713,'amu*angstrom^2'), symmetry=1, barrier=(19.5136,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78346,0.048442,-4.55665e-05,2.11568e-08,-3.94802e-12,971.938,22.0766], Tmin=(100,'K'), Tmax=(1272.32,'K')), NASAPolynomial(coeffs=[11.8131,0.0169105,-8.39291e-06,1.67892e-09,-1.20833e-13,-1580.29,-28.7265], Tmin=(1272.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.41646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]C([O])=CC=CO(5157)',
    structure = SMILES('C=C([O])[CH]C=CO'),
    E0 = (-79.88,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23817,'amu*angstrom^2'), symmetry=1, barrier=(28.4679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23448,'amu*angstrom^2'), symmetry=1, barrier=(28.3831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23815,'amu*angstrom^2'), symmetry=1, barrier=(28.4676,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.347664,0.0625043,-1.82904e-05,-4.25449e-08,2.72266e-11,-9459.18,24.5509], Tmin=(100,'K'), Tmax=(911.865,'K')), NASAPolynomial(coeffs=[24.9673,0.00211818,2.72483e-06,-6.50021e-10,4.22106e-14,-15928.6,-102.807], Tmin=(911.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C=[C]CC(=O)C=C[CH]O(2833)',
    structure = SMILES('C=[C]CC([O])=CC=CO'),
    E0 = (92.8171,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.835352,'amu*angstrom^2'), symmetry=1, barrier=(19.2064,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835032,'amu*angstrom^2'), symmetry=1, barrier=(19.199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835386,'amu*angstrom^2'), symmetry=1, barrier=(19.2072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835631,'amu*angstrom^2'), symmetry=1, barrier=(19.2128,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.11863,0.103339,-0.000106984,5.23013e-08,-9.42349e-12,11410.5,39.332], Tmin=(100,'K'), Tmax=(1599.51,'K')), NASAPolynomial(coeffs=[27.9497,0.00823509,8.7432e-07,-4.35635e-10,3.55135e-14,4338.52,-111.891], Tmin=(1599.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.8171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([O])CC=C(4101)',
    structure = SMILES('[CH]=C([O])CC=C'),
    E0 = (259.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,350,440,435,1725,356.976,356.992],'cm^-1')),
        HinderedRotor(inertia=(0.121036,'amu*angstrom^2'), symmetry=1, barrier=(10.9441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121029,'amu*angstrom^2'), symmetry=1, barrier=(10.9437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51949,0.0464528,-2.33308e-05,-6.76514e-09,7.1228e-12,31258.2,23.688], Tmin=(100,'K'), Tmax=(965.886,'K')), NASAPolynomial(coeffs=[13.3363,0.0162286,-5.45376e-06,9.53765e-10,-6.66555e-14,28102.5,-37.4299], Tmin=(965.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC(=O)[C]=C[CH]O(2818)',
    structure = SMILES('C=CCC([O])=[C]C=CO'),
    E0 = (53.9707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.973677,'amu*angstrom^2'), symmetry=1, barrier=(22.3867,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97422,'amu*angstrom^2'), symmetry=1, barrier=(22.3992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973241,'amu*angstrom^2'), symmetry=1, barrier=(22.3767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.974311,'amu*angstrom^2'), symmetry=1, barrier=(22.4013,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10918,0.104768,-0.000110051,5.50396e-08,-1.01512e-11,6736.57,38.4121], Tmin=(100,'K'), Tmax=(1563.61,'K')), NASAPolynomial(coeffs=[27.4056,0.00917618,9.23397e-07,-4.92225e-10,4.12909e-14,-37.714,-109.32], Tmin=(1563.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.9707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CCC([O])=C[C]=CO(5158)',
    structure = SMILES('C=CCC([O])=C[C]=CO'),
    E0 = (53.9707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.973677,'amu*angstrom^2'), symmetry=1, barrier=(22.3867,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97422,'amu*angstrom^2'), symmetry=1, barrier=(22.3992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973241,'amu*angstrom^2'), symmetry=1, barrier=(22.3767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.974311,'amu*angstrom^2'), symmetry=1, barrier=(22.4013,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10918,0.104768,-0.000110051,5.50396e-08,-1.01512e-11,6736.57,38.4121], Tmin=(100,'K'), Tmax=(1563.61,'K')), NASAPolynomial(coeffs=[27.4056,0.00917618,9.23397e-07,-4.92225e-10,4.12909e-14,-37.714,-109.32], Tmin=(1563.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.9707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CCC([O])=CC=[C]O(5159)',
    structure = SMILES('C=CCC([O])=CC=[C]O'),
    E0 = (94.7195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.704664,'amu*angstrom^2'), symmetry=1, barrier=(16.2016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.705026,'amu*angstrom^2'), symmetry=1, barrier=(16.2099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.704571,'amu*angstrom^2'), symmetry=1, barrier=(16.1995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.70566,'amu*angstrom^2'), symmetry=1, barrier=(16.2245,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5966,0.0972605,-9.79509e-05,4.75378e-08,-8.61156e-12,11615.8,40.1091], Tmin=(100,'K'), Tmax=(1560.09,'K')), NASAPolynomial(coeffs=[25.8625,0.0117253,-1.16173e-06,-3.96746e-11,8.85396e-15,4889.44,-98.6758], Tmin=(1560.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.7195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=CCC([O])=CC=CO(5160)',
    structure = SMILES('[CH]=CCC([O])=CC=CO'),
    E0 = (102.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.884036,'amu*angstrom^2'), symmetry=1, barrier=(20.3257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.885501,'amu*angstrom^2'), symmetry=1, barrier=(20.3594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.882073,'amu*angstrom^2'), symmetry=1, barrier=(20.2806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.88397,'amu*angstrom^2'), symmetry=1, barrier=(20.3242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42509,0.105256,-0.000108751,5.25578e-08,-9.29792e-12,12538.8,40.3713], Tmin=(100,'K'), Tmax=(1649.68,'K')), NASAPolynomial(coeffs=[28.6485,0.00679402,1.79695e-06,-6.10922e-10,4.68092e-14,5432.05,-115.562], Tmin=(1649.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC[C]1OC1C=CO(5168)',
    structure = SMILES('C=CC[C]1OC1C=CO'),
    E0 = (15.3459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.24248,0.0993545,-9.48004e-05,4.48253e-08,-7.80125e-12,2103.3,39.9673], Tmin=(100,'K'), Tmax=(1708.29,'K')), NASAPolynomial(coeffs=[22.7817,0.0163433,-4.70876e-07,-3.54112e-10,3.50165e-14,-2883.7,-83.7322], Tmin=(1708.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.3459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C2CsJO)"""),
)

species(
    label = 'OC=CC=C1C[CH]CO1(5169)',
    structure = SMILES('OC=CC=C1C[CH]CO1'),
    E0 = (-106.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.110966,0.0462228,9.75768e-05,-1.89091e-07,8.4482e-11,-12677.1,30.7607], Tmin=(100,'K'), Tmax=(916.738,'K')), NASAPolynomial(coeffs=[34.4145,-0.00101905,7.26827e-06,-1.53056e-09,9.37957e-14,-23270.9,-155.228], Tmin=(916.738,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-106.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(CCJCO)"""),
)

species(
    label = 'C=CCC1=C[CH]C(O)O1(5170)',
    structure = SMILES('C=CCC1=C[CH]C(O)O1'),
    E0 = (-171.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.293212,0.0589982,2.48842e-05,-8.57798e-08,4.05699e-11,-20412.8,26.9711], Tmin=(100,'K'), Tmax=(934.949,'K')), NASAPolynomial(coeffs=[22.6924,0.0203296,-4.78594e-06,7.69374e-10,-5.8552e-14,-27099.6,-92.9462], Tmin=(934.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=CC=C(CC=C)OO(5171)',
    structure = SMILES('[CH]=CC=C(CC=C)OO'),
    E0 = (297.554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,235.981,235.981],'cm^-1')),
        HinderedRotor(inertia=(0.40083,'amu*angstrom^2'), symmetry=1, barrier=(15.8394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400829,'amu*angstrom^2'), symmetry=1, barrier=(15.8394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.40083,'amu*angstrom^2'), symmetry=1, barrier=(15.8394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.40083,'amu*angstrom^2'), symmetry=1, barrier=(15.8394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.40083,'amu*angstrom^2'), symmetry=1, barrier=(15.8394,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.465781,0.0905257,-8.33183e-05,3.91085e-08,-7.30349e-12,35954.7,35.9951], Tmin=(100,'K'), Tmax=(1291.33,'K')), NASAPolynomial(coeffs=[19.4326,0.028889,-1.17218e-05,2.14598e-09,-1.47623e-13,30815.6,-65.0909], Tmin=(1291.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC(=O)[CH]CC=O(4556)',
    structure = SMILES('C=CCC([O])=CCC=O'),
    E0 = (-123.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,317.438,317.52,317.524,317.585],'cm^-1')),
        HinderedRotor(inertia=(0.0016724,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00167241,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223221,'amu*angstrom^2'), symmetry=1, barrier=(15.9621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223145,'amu*angstrom^2'), symmetry=1, barrier=(15.9627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.114062,0.0723827,-3.7781e-05,-5.76727e-10,4.26865e-12,-14677.8,34.6478], Tmin=(100,'K'), Tmax=(1142.81,'K')), NASAPolynomial(coeffs=[17.994,0.0321103,-1.42041e-05,2.75192e-09,-1.96453e-13,-20221.3,-60.3732], Tmin=(1142.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-123.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O(4)(4)',
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
    label = 'C=CC[C]=CC=CO(5172)',
    structure = SMILES('C=CC[C]=CC=CO'),
    E0 = (169.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.912557,'amu*angstrom^2'), symmetry=1, barrier=(20.9815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.911731,'amu*angstrom^2'), symmetry=1, barrier=(20.9625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.912072,'amu*angstrom^2'), symmetry=1, barrier=(20.9703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.912437,'amu*angstrom^2'), symmetry=1, barrier=(20.9787,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.126168,0.0728766,-2.32027e-05,-3.77935e-08,2.44276e-11,20559.7,30.6018], Tmin=(100,'K'), Tmax=(930.179,'K')), NASAPolynomial(coeffs=[24.0689,0.0145272,-2.79698e-06,3.94439e-10,-3.03117e-14,14081.7,-95.0006], Tmin=(930.179,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CCC1([O])C=CC1O(5183)',
    structure = SMILES('C=CCC1([O])C=CC1O'),
    E0 = (59.3584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.363927,0.0616618,3.18831e-06,-5.17944e-08,2.49239e-11,7286.42,34.2043], Tmin=(100,'K'), Tmax=(982.077,'K')), NASAPolynomial(coeffs=[19.4779,0.0270104,-9.86852e-06,1.86018e-09,-1.36691e-13,1448.87,-68.2705], Tmin=(982.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.3584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1CC(=O)C=CC1O(5184)',
    structure = SMILES('[CH2]C1CC(=O)C=CC1O'),
    E0 = (-137.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.880153,0.0542882,8.72564e-06,-3.92093e-08,1.52726e-11,-16395.2,27.5827], Tmin=(100,'K'), Tmax=(1113.92,'K')), NASAPolynomial(coeffs=[13.3398,0.0422071,-1.89866e-05,3.69797e-09,-2.64664e-13,-21197.2,-42.9673], Tmin=(1113.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-137.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + ring(Cyclohexane) + radical(Isobutyl)"""),
)

species(
    label = 'C=CCC(=O)C=C=CO(2839)',
    structure = SMILES('C=CCC(=O)C=C=CO'),
    E0 = (-109.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.927847,'amu*angstrom^2'), symmetry=1, barrier=(21.333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.927603,'amu*angstrom^2'), symmetry=1, barrier=(21.3274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.928123,'amu*angstrom^2'), symmetry=1, barrier=(21.3394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.926834,'amu*angstrom^2'), symmetry=1, barrier=(21.3097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.392564,0.0829515,-6.24901e-05,1.82704e-08,-8.33914e-13,-13053.2,31.5134], Tmin=(100,'K'), Tmax=(1178.76,'K')), NASAPolynomial(coeffs=[21.9318,0.0250375,-1.14965e-05,2.27031e-09,-1.63739e-13,-19555.7,-85.1179], Tmin=(1178.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CCC(=O)C=CC[O](4558)',
    structure = SMILES('C=CCC(=O)C=CC[O]'),
    E0 = (-5.18239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97398,0.0614141,-3.04684e-05,5.22236e-09,-1.93412e-13,-566.311,28.6887], Tmin=(100,'K'), Tmax=(2252.36,'K')), NASAPolynomial(coeffs=[33.8429,0.0182164,-1.06233e-05,1.98964e-09,-1.27748e-13,-18321.1,-158.482], Tmin=(2252.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.18239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = 'C=CCC(=O)C=[C]CO(5185)',
    structure = SMILES('C=CCC(=O)C=[C]CO'),
    E0 = (6.95424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26565,0.071976,-5.27409e-05,2.11897e-08,-4.0623e-12,924.911,30.8862], Tmin=(100,'K'), Tmax=(1063.32,'K')), NASAPolynomial(coeffs=[5.93664,0.0544046,-2.79532e-05,5.64861e-09,-4.08384e-13,-68.4393,8.0646], Tmin=(1063.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.95424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CCC(=O)[C]=CCO(5186)',
    structure = SMILES('C=CCC([O])=C=CCO'),
    E0 = (-36.5867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,249.727,249.73,249.748,4000],'cm^-1')),
        HinderedRotor(inertia=(0.32813,'amu*angstrom^2'), symmetry=1, barrier=(14.5204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.32806,'amu*angstrom^2'), symmetry=1, barrier=(14.5204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.328127,'amu*angstrom^2'), symmetry=1, barrier=(14.5204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.328124,'amu*angstrom^2'), symmetry=1, barrier=(14.5204,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0701628,0.0866113,-8.15901e-05,4.11622e-08,-8.42038e-12,-4251.37,35.4125], Tmin=(100,'K'), Tmax=(1172.36,'K')), NASAPolynomial(coeffs=[15.5164,0.033431,-1.35472e-05,2.46923e-09,-1.69288e-13,-7905.98,-42.2624], Tmin=(1172.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-36.5867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C]CC(=O)C=CCO(5188)',
    structure = SMILES('C=[C]CC(=O)C=CCO'),
    E0 = (6.95424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26565,0.071976,-5.27409e-05,2.11897e-08,-4.0623e-12,924.911,30.8862], Tmin=(100,'K'), Tmax=(1063.32,'K')), NASAPolynomial(coeffs=[5.93664,0.0544046,-2.79532e-05,5.64861e-09,-4.08384e-13,-68.4393,8.0646], Tmin=(1063.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.95424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CCC(=O)C=CCO(5189)',
    structure = SMILES('[CH]=CCC(=O)C=CCO'),
    E0 = (16.2086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,375,552.5,462.5,1710,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32745,0.0705032,-4.61331e-05,1.40239e-08,-1.76548e-12,2033.73,30.527], Tmin=(100,'K'), Tmax=(1661.56,'K')), NASAPolynomial(coeffs=[12.4885,0.0436342,-2.18767e-05,4.29143e-09,-3.01126e-13,-1675.2,-28.9856], Tmin=(1661.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.2086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C4H5O(93)(92)',
    structure = SMILES('C=CC[C]=O'),
    E0 = (66.8219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,458.926],'cm^-1')),
        HinderedRotor(inertia=(0.0997865,'amu*angstrom^2'), symmetry=1, barrier=(14.9157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.099798,'amu*angstrom^2'), symmetry=1, barrier=(14.9167,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3285.42,'J/mol'), sigma=(5.46087,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=513.18 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51804,0.0238835,1.19491e-05,-2.85418e-08,1.09388e-11,8097.53,17.8098], Tmin=(100,'K'), Tmax=(1083.61,'K')), NASAPolynomial(coeffs=[9.78041,0.0178579,-8.47799e-06,1.72441e-09,-1.27255e-13,5303.46,-23.4402], Tmin=(1083.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.8219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CCC(=O)C1[CH]C1O(5179)',
    structure = SMILES('C=CCC(=O)C1[CH]C1O'),
    E0 = (-12.9392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.288617,0.0661122,-1.74578e-05,-2.31731e-08,1.24383e-11,-1409.35,33.7374], Tmin=(100,'K'), Tmax=(1070.18,'K')), NASAPolynomial(coeffs=[18.29,0.0311364,-1.37185e-05,2.70706e-09,-1.97309e-13,-7112.4,-62.9739], Tmin=(1070.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.9392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCJCO)"""),
)

species(
    label = 'O=C1C=CC(O)C[CH]C1(5190)',
    structure = SMILES('O=C1C=CC(O)C[CH]C1'),
    E0 = (-103.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.974935,0.0153851,0.000189167,-2.90011e-07,1.22209e-10,-12269.9,32.6446], Tmin=(100,'K'), Tmax=(914.403,'K')), NASAPolynomial(coeffs=[37.2285,-0.0102466,1.31083e-05,-2.6369e-09,1.64811e-13,-24458.5,-169.408], Tmin=(914.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + ring(Cycloheptane) + radical(CCJCC=O)"""),
)

species(
    label = 'C=COC(=C)[CH]C=CO(5182)',
    structure = SMILES('C=COC(=C)[CH]C=CO'),
    E0 = (-81.2276,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.941878,0.0879063,-3.7893e-05,-3.39289e-08,2.49775e-11,-9572.4,33.4977], Tmin=(100,'K'), Tmax=(937.367,'K')), NASAPolynomial(coeffs=[28.725,0.013557,-2.52476e-06,3.79629e-10,-3.17615e-14,-17429.5,-119.953], Tmin=(937.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.2276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJCO)"""),
)

species(
    label = 'C=CCC([CH]O)C=C=O(5191)',
    structure = SMILES('C=CCC([CH]O)C=C=O'),
    E0 = (-34.1851,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.176024,0.0958839,-0.000107811,6.74203e-08,-1.72955e-11,-3964.55,34.8691], Tmin=(100,'K'), Tmax=(939.244,'K')), NASAPolynomial(coeffs=[13.1039,0.0393279,-1.74887e-05,3.3097e-09,-2.30922e-13,-6459.14,-28.3662], Tmin=(939.244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.1851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=CC(=O)CC=C(2806)',
    structure = SMILES('[CH]=CC(=O)CC=C'),
    E0 = (205.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,625.128],'cm^-1')),
        HinderedRotor(inertia=(0.0115181,'amu*angstrom^2'), symmetry=1, barrier=(3.3326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.799554,'amu*angstrom^2'), symmetry=1, barrier=(18.3833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800305,'amu*angstrom^2'), symmetry=1, barrier=(18.4006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51193,0.0511228,-2.78333e-05,5.11544e-09,-6.08057e-14,24792.8,24.3932], Tmin=(100,'K'), Tmax=(1699.06,'K')), NASAPolynomial(coeffs=[17.9658,0.0225202,-1.15282e-05,2.22816e-09,-1.52492e-13,17738.9,-68.0134], Tmin=(1699.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1CC(=O)C1C=CO(5173)',
    structure = SMILES('[CH2]C1CC(=O)C1C=CO'),
    E0 = (-79.1351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.171713,0.0643672,6.3723e-06,-6.53693e-08,3.31935e-11,-9361.67,31.5749], Tmin=(100,'K'), Tmax=(932.943,'K')), NASAPolynomial(coeffs=[21.9906,0.0212089,-5.25509e-06,8.33702e-10,-6.08101e-14,-15625.8,-83.9268], Tmin=(932.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.1351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutanone) + radical(Isobutyl)"""),
)

species(
    label = 'C=CCC(=O)C[C]=CO(5174)',
    structure = SMILES('C=CCC(=O)C[C]=CO'),
    E0 = (-34.1929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.392355,0.0788841,-3.70054e-05,-1.20539e-08,1.00132e-11,-3939.22,32.9695], Tmin=(100,'K'), Tmax=(1084.72,'K')), NASAPolynomial(coeffs=[23.2881,0.0264205,-1.26623e-05,2.61184e-09,-1.95076e-13,-11127.4,-92.6542], Tmin=(1084.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.1929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CCC(=O)CC=[C]O(5176)',
    structure = SMILES('C=CCC(=O)CC=[C]O'),
    E0 = (-32.2905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.162463,0.0760693,-3.85621e-05,-4.22513e-09,5.93738e-12,-3720.89,34.8071], Tmin=(100,'K'), Tmax=(1149.95,'K')), NASAPolynomial(coeffs=[21.2217,0.0296383,-1.44579e-05,2.93831e-09,-2.15252e-13,-10487.2,-79.3829], Tmin=(1149.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.2905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]CC(=O)CC=CO(5177)',
    structure = SMILES('C=[C]CC(=O)CC=CO'),
    E0 = (-34.1929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.392355,0.0788841,-3.70054e-05,-1.20539e-08,1.00132e-11,-3939.22,32.9695], Tmin=(100,'K'), Tmax=(1084.72,'K')), NASAPolynomial(coeffs=[23.2881,0.0264205,-1.26623e-05,2.61184e-09,-1.95076e-13,-11127.4,-92.6542], Tmin=(1084.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.1929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CCC(=O)C[CH]C=O(4555)',
    structure = SMILES('C=CCC(=O)CC=C[O]'),
    E0 = (-130.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0836956,0.0668991,-7.1405e-06,-3.76371e-08,1.74907e-11,-15546.8,32.1998], Tmin=(100,'K'), Tmax=(1081.39,'K')), NASAPolynomial(coeffs=[21.2642,0.0305527,-1.49818e-05,3.11231e-09,-2.32963e-13,-22583.4,-82.9957], Tmin=(1081.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CCC(=O)CC=CO(5178)',
    structure = SMILES('[CH]=CCC(=O)CC=CO'),
    E0 = (-24.9385,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,375,552.5,462.5,1710,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.4232,0.0777366,-2.88638e-05,-2.35858e-08,1.47385e-11,-2823.3,33.0067], Tmin=(100,'K'), Tmax=(1054.12,'K')), NASAPolynomial(coeffs=[24.3753,0.0246396,-1.16567e-05,2.43429e-09,-1.84577e-13,-10329.5,-98.7448], Tmin=(1054.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.9385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'O=C1C[CH]CC1C=CO(5180)',
    structure = SMILES('O=C1C[CH]CC1C=CO'),
    E0 = (-150.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.39362,0.0552286,3.61846e-05,-9.67587e-08,4.39955e-11,-17923.1,30.4868], Tmin=(100,'K'), Tmax=(941.892,'K')), NASAPolynomial(coeffs=[22.8413,0.0205522,-5.18524e-06,8.90946e-10,-6.94633e-14,-24842.2,-90.7487], Tmin=(941.892,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-150.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopentanone) + radical(CCJCC=O)"""),
)

species(
    label = 'C=CCC([C]=O)C=CO(5181)',
    structure = SMILES('C=CCC([C]=O)C=CO'),
    E0 = (-99.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1855,455,950,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.52985,0.0875903,-6.67028e-05,1.47131e-08,3.25789e-12,-11838.8,35.7855], Tmin=(100,'K'), Tmax=(981.365,'K')), NASAPolynomial(coeffs=[21.3463,0.023891,-8.26536e-06,1.45837e-09,-1.01868e-13,-17358.9,-75.5909], Tmin=(981.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C=CC[C](O)C1C=CO1(5197)',
    structure = SMILES('C=CC[C](O)C1C=CO1'),
    E0 = (23.1021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.411951,0.0742886,-7.72754e-06,-5.94525e-08,3.24677e-11,2958.17,32.6021], Tmin=(100,'K'), Tmax=(946.494,'K')), NASAPolynomial(coeffs=[26.948,0.0151421,-3.50188e-06,6.17373e-10,-5.12946e-14,-4750.94,-111.254], Tmin=(946.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.1021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1CC(O)=CC=CO1(5198)',
    structure = SMILES('[CH2]C1CC(O)=CC=CO1'),
    E0 = (-105.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.568677,0.0288387,0.000225457,-3.86595e-07,1.72978e-10,-12414.5,28.2292], Tmin=(100,'K'), Tmax=(895.292,'K')), NASAPolynomial(coeffs=[59.6989,-0.0493289,3.62519e-05,-7.29669e-09,4.90015e-13,-30864.6,-298.633], Tmin=(895.292,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=CCC(O)=CC=C=O(3131)',
    structure = SMILES('C=CCC(O)=CC=C=O'),
    E0 = (-139.544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,187.18,188.182,191.451],'cm^-1')),
        HinderedRotor(inertia=(0.558369,'amu*angstrom^2'), symmetry=1, barrier=(14.5703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.567699,'amu*angstrom^2'), symmetry=1, barrier=(14.5901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.591459,'amu*angstrom^2'), symmetry=1, barrier=(14.5634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.565407,'amu*angstrom^2'), symmetry=1, barrier=(14.5777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.30359,0.0864036,-8.19264e-05,3.97243e-08,-7.61744e-12,-16621.4,33.2806], Tmin=(100,'K'), Tmax=(1266.07,'K')), NASAPolynomial(coeffs=[18.9639,0.0255305,-9.80623e-06,1.74862e-09,-1.18749e-13,-21500.2,-64.2196], Tmin=(1266.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-139.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC[C]=CC=C[O](2415)',
    structure = SMILES('C=CC[C]=CC=C[O]'),
    E0 = (311.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.934519,'amu*angstrom^2'), symmetry=1, barrier=(21.4864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.93233,'amu*angstrom^2'), symmetry=1, barrier=(21.4361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.932627,'amu*angstrom^2'), symmetry=1, barrier=(21.4429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.282775,0.0659412,-1.99483e-05,-3.13302e-08,1.95473e-11,37556.9,30.5438], Tmin=(100,'K'), Tmax=(952.309,'K')), NASAPolynomial(coeffs=[20.8062,0.01823,-5.42938e-06,9.51292e-10,-7.0081e-14,31902.5,-76.6313], Tmin=(952.309,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=CC=[C]O(3148)',
    structure = SMILES('[O]C=CC=[C]O'),
    E0 = (57.5428,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17839,'amu*angstrom^2'), symmetry=1, barrier=(27.0934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18425,'amu*angstrom^2'), symmetry=1, barrier=(27.2281,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22215,0.043282,4.97501e-06,-6.05358e-08,3.33983e-11,7037.63,22.2255], Tmin=(100,'K'), Tmax=(897.589,'K')), NASAPolynomial(coeffs=[23.332,-0.00653379,6.81593e-06,-1.4385e-09,9.7429e-14,1106.15,-92.9843], Tmin=(897.589,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.5428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(O)=CC=C[O](5193)',
    structure = SMILES('C=C(O)[CH]C=C[O]'),
    E0 = (-76.2221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.35339,'amu*angstrom^2'), symmetry=1, barrier=(31.117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35269,'amu*angstrom^2'), symmetry=1, barrier=(31.1009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35312,'amu*angstrom^2'), symmetry=1, barrier=(31.1109,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.400843,0.0594199,-7.08619e-06,-5.37493e-08,3.0632e-11,-9019.51,24.2686], Tmin=(100,'K'), Tmax=(924.396,'K')), NASAPolynomial(coeffs=[25.3551,0.00260988,2.06445e-06,-4.65225e-10,2.62503e-14,-15819.3,-105.985], Tmin=(924.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.2221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C=[C]C[C](O)C=CC=O(2815)',
    structure = SMILES('C=[C]CC(O)=CC=C[O]'),
    E0 = (96.4749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.915269,'amu*angstrom^2'), symmetry=1, barrier=(21.0438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.913694,'amu*angstrom^2'), symmetry=1, barrier=(21.0076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.914834,'amu*angstrom^2'), symmetry=1, barrier=(21.0338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.914668,'amu*angstrom^2'), symmetry=1, barrier=(21.03,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.705948,0.0845129,-4.18917e-05,-2.74746e-08,2.27989e-11,11790.3,34.1508], Tmin=(100,'K'), Tmax=(926.87,'K')), NASAPolynomial(coeffs=[27.886,0.0101276,-8.19636e-07,2.82049e-11,-5.58605e-15,4385.01,-112.973], Tmin=(926.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'CHCHO(45)(45)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (245.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06864,0.0187233,-1.21319e-05,-3.33727e-10,2.32882e-12,29739.4,14.7866], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.96288,0.00799899,-4.30606e-06,1.11076e-09,-1.11415e-13,28725.6,-5.17392], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(245.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH]=C(O)CC=C(5194)',
    structure = SMILES('[CH]=C(O)CC=C'),
    E0 = (121.287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.729771,'amu*angstrom^2'), symmetry=1, barrier=(16.7789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.73034,'amu*angstrom^2'), symmetry=1, barrier=(16.792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.728136,'amu*angstrom^2'), symmetry=1, barrier=(16.7413,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16799,0.0502452,-1.51313e-05,-2.48355e-08,1.56226e-11,14700.5,23.4489], Tmin=(100,'K'), Tmax=(942.762,'K')), NASAPolynomial(coeffs=[16.9971,0.0130011,-3.47269e-06,5.79646e-10,-4.26801e-14,10386.4,-59.0356], Tmin=(942.762,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC(O)=[C]C=C[O](3162)',
    structure = SMILES('C=CCC(O)=[C]C=C[O]'),
    E0 = (57.6286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08593,'amu*angstrom^2'), symmetry=1, barrier=(24.9676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08649,'amu*angstrom^2'), symmetry=1, barrier=(24.9805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08582,'amu*angstrom^2'), symmetry=1, barrier=(24.9652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08595,'amu*angstrom^2'), symmetry=1, barrier=(24.9681,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25555,0.103776,-0.000105018,5.03148e-08,-8.88412e-12,7185.47,38.8643], Tmin=(100,'K'), Tmax=(1639.51,'K')), NASAPolynomial(coeffs=[27.8768,0.00934983,5.04373e-07,-3.72783e-10,3.11987e-14,115.418,-112.833], Tmin=(1639.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.6286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC[C](O)C=[C]C=O(2817)',
    structure = SMILES('C=CCC(O)=C[C]=C[O]'),
    E0 = (57.6286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08593,'amu*angstrom^2'), symmetry=1, barrier=(24.9676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08649,'amu*angstrom^2'), symmetry=1, barrier=(24.9805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08582,'amu*angstrom^2'), symmetry=1, barrier=(24.9652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08595,'amu*angstrom^2'), symmetry=1, barrier=(24.9681,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25555,0.103776,-0.000105018,5.03148e-08,-8.88412e-12,7185.47,38.8643], Tmin=(100,'K'), Tmax=(1639.51,'K')), NASAPolynomial(coeffs=[27.8768,0.00934983,5.04373e-07,-3.72783e-10,3.11987e-14,115.418,-112.833], Tmin=(1639.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.6286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CC[C](O)C=CC=O(2823)',
    structure = SMILES('[CH]=CCC(O)=CC=C[O]'),
    E0 = (105.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.957669,'amu*angstrom^2'), symmetry=1, barrier=(22.0187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.956663,'amu*angstrom^2'), symmetry=1, barrier=(21.9956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.956609,'amu*angstrom^2'), symmetry=1, barrier=(21.9943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.955314,'amu*angstrom^2'), symmetry=1, barrier=(21.9646,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.749247,0.083459,-3.38008e-05,-3.93914e-08,2.79072e-11,12906.8,34.236], Tmin=(100,'K'), Tmax=(923.917,'K')), NASAPolynomial(coeffs=[29.3428,0.00774963,5.17611e-07,-2.25536e-10,1.10997e-14,5017.07,-121.164], Tmin=(923.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC[C](O)C=C[C]=O(2824)',
    structure = SMILES('C=CCC(O)=C[CH][C]=O'),
    E0 = (60.9925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1855,455,950,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.136119,0.0712849,-1.8086e-05,-3.70725e-08,2.10986e-11,7502.12,35.2372], Tmin=(100,'K'), Tmax=(993.259,'K')), NASAPolynomial(coeffs=[23.9975,0.0191824,-7.49175e-06,1.51803e-09,-1.17313e-13,483.858,-92.2262], Tmin=(993.259,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.9925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(CCJC=O)"""),
)

species(
    label = 'C=CCC(O)=CC1[CH]O1(5199)',
    structure = SMILES('C=CCC(O)=CC1[CH]O1'),
    E0 = (15.9995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.374364,0.0761622,-1.18792e-05,-6.21385e-08,3.73738e-11,2100.78,33.8156], Tmin=(100,'K'), Tmax=(889.701,'K')), NASAPolynomial(coeffs=[26.4271,0.0120805,1.04696e-06,-5.54572e-10,4.28637e-14,-4901.08,-104.902], Tmin=(889.701,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.9995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCsJO)"""),
)

species(
    label = 'C=CCC1(O)[CH]C=CO1(5200)',
    structure = SMILES('C=CCC1(O)[CH]C=CO1'),
    E0 = (-179.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.395066,0.056537,8.33272e-05,-1.85423e-07,8.65015e-11,-21450.2,25.3143], Tmin=(100,'K'), Tmax=(905.7,'K')), NASAPolynomial(coeffs=[37.8676,-0.00604973,1.07659e-05,-2.30288e-09,1.51421e-13,-32745,-179.582], Tmin=(905.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(C=CCJCO)"""),
)

species(
    label = 'OC1=CC=COC[CH]C1(5201)',
    structure = SMILES('OC1=CC=COC[CH]C1'),
    E0 = (-87.4854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.616717,0.0323709,0.0001338,-2.20368e-07,9.34427e-11,-10361.3,25.2877], Tmin=(100,'K'), Tmax=(925.626,'K')), NASAPolynomial(coeffs=[32.2813,0.0038109,4.61966e-06,-9.54176e-10,4.95597e-14,-20861.6,-150.084], Tmin=(925.626,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.4854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cyclooctane) + radical(CCJCO)"""),
)

species(
    label = 'C=CC[C]=CC=COO(5202)',
    structure = SMILES('C=CC[C]=CC=COO'),
    E0 = (294.066,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1685,370,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.744828,0.091913,-8.45957e-05,3.93214e-08,-7.18314e-12,35549.2,36.8738], Tmin=(100,'K'), Tmax=(1332.58,'K')), NASAPolynomial(coeffs=[21.5653,0.0249442,-9.21257e-06,1.60819e-09,-1.07857e-13,29603.2,-77.1652], Tmin=(1332.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC=C(O)CC=C(5203)',
    structure = SMILES('[CH]=CC=C(O)CC=C'),
    E0 = (173.059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.820322,'amu*angstrom^2'), symmetry=1, barrier=(18.8608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.82151,'amu*angstrom^2'), symmetry=1, barrier=(18.8881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.822981,'amu*angstrom^2'), symmetry=1, barrier=(18.9219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.823453,'amu*angstrom^2'), symmetry=1, barrier=(18.9328,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.130014,0.0747568,-3.30777e-05,-2.38424e-08,1.8352e-11,20977.6,30.7431], Tmin=(100,'K'), Tmax=(942.616,'K')), NASAPolynomial(coeffs=[22.8775,0.0169484,-4.45899e-06,7.37459e-10,-5.42551e-14,14870.9,-88.2799], Tmin=(942.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1CC(O)=CC1C=O(5204)',
    structure = SMILES('[CH2]C1CC(O)=CC1C=O'),
    E0 = (-133.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.259998,0.0653637,-5.96464e-06,-4.54676e-08,2.41103e-11,-15939.5,29.9765], Tmin=(100,'K'), Tmax=(952.959,'K')), NASAPolynomial(coeffs=[19.5959,0.0254449,-8.04869e-06,1.40544e-09,-1.00588e-13,-21497.5,-72.2021], Tmin=(952.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-133.776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopentene) + radical(Isobutyl)"""),
)

species(
    label = 'C=CCC(O)=C=CC=O(2838)',
    structure = SMILES('C=CCC(O)=C=CC=O'),
    E0 = (-108.718,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.788706,'amu*angstrom^2'), symmetry=1, barrier=(18.1339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788832,'amu*angstrom^2'), symmetry=1, barrier=(18.1368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.789678,'amu*angstrom^2'), symmetry=1, barrier=(18.1563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788374,'amu*angstrom^2'), symmetry=1, barrier=(18.1263,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.261835,0.0877121,-8.11087e-05,3.73437e-08,-6.83097e-12,-12917.4,32.1187], Tmin=(100,'K'), Tmax=(1313.29,'K')), NASAPolynomial(coeffs=[19.5727,0.0273001,-1.21077e-05,2.31655e-09,-1.63128e-13,-18127.1,-68.977], Tmin=(1313.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CCC(O)=[C]CC=O(5205)',
    structure = SMILES('C=CCC(O)=[C]CC=O'),
    E0 = (-23.2551,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,1685,370,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.297842,0.0811194,-5.56853e-05,1.25106e-08,9.07085e-13,-2630.97,35.0979], Tmin=(100,'K'), Tmax=(1146.22,'K')), NASAPolynomial(coeffs=[20.4939,0.028328,-1.24666e-05,2.41834e-09,-1.73077e-13,-8695.81,-73.7117], Tmin=(1146.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.2551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CCC(O)=CC[C]=O(5206)',
    structure = SMILES('C=CCC(O)=CC[C]=O'),
    E0 = (-101.136,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1855,455,950,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.273993,0.0788568,-4.43626e-05,-3.58336e-09,7.66997e-12,-11996.9,35.3496], Tmin=(100,'K'), Tmax=(1046.75,'K')), NASAPolynomial(coeffs=[21.3776,0.0257768,-1.07988e-05,2.10803e-09,-1.54012e-13,-18154.5,-77.8573], Tmin=(1046.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]CC(O)=CCC=O(5208)',
    structure = SMILES('C=[C]CC(O)=CCC=O'),
    E0 = (-23.2551,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,1685,370,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.297842,0.0811194,-5.56853e-05,1.25106e-08,9.07085e-13,-2630.97,35.0979], Tmin=(100,'K'), Tmax=(1146.22,'K')), NASAPolynomial(coeffs=[20.4939,0.028328,-1.24666e-05,2.41834e-09,-1.73077e-13,-8695.81,-73.7117], Tmin=(1146.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.2551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CCC(O)=CCC=O(5209)',
    structure = SMILES('[CH]=CCC(O)=CCC=O'),
    E0 = (-14.0007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3120,650,792.5,1650,366.891,367.063],'cm^-1')),
        HinderedRotor(inertia=(0.162994,'amu*angstrom^2'), symmetry=1, barrier=(15.61,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163237,'amu*angstrom^2'), symmetry=1, barrier=(15.6105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163412,'amu*angstrom^2'), symmetry=1, barrier=(15.6116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163306,'amu*angstrom^2'), symmetry=1, barrier=(15.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163462,'amu*angstrom^2'), symmetry=1, barrier=(15.6101,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.301118,0.0796818,-4.66989e-05,1.32734e-10,5.88504e-12,-1516.29,35.0339], Tmin=(100,'K'), Tmax=(1074.71,'K')), NASAPolynomial(coeffs=[21.2187,0.0271243,-1.17783e-05,2.31312e-09,-1.68416e-13,-7732.11,-77.7361], Tmin=(1074.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.0007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC1(O)[CH]C1C=O(5210)',
    structure = SMILES('C=CCC1(O)[CH]C1C=O'),
    E0 = (-4.37804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0541028,0.0730378,-3.3508e-05,-1.30088e-08,1.13316e-11,-372.306,34.3506], Tmin=(100,'K'), Tmax=(991.945,'K')), NASAPolynomial(coeffs=[18.8693,0.0273198,-9.97249e-06,1.81921e-09,-1.29199e-13,-5588.55,-63.7476], Tmin=(991.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.37804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCJCO)"""),
)

species(
    label = 'O=CC1C=C(O)C[CH]C1(5211)',
    structure = SMILES('O=CC1C=C(O)C[CH]C1'),
    E0 = (-154.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.574906,0.0621057,-1.24271e-05,-2.49653e-08,1.2981e-11,-18488.7,28.9042], Tmin=(100,'K'), Tmax=(1030.03,'K')), NASAPolynomial(coeffs=[14.9658,0.0345726,-1.36202e-05,2.53044e-09,-1.78746e-13,-22957.3,-48.2503], Tmin=(1030.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-154.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclohexene) + radical(RCCJCC)"""),
)

species(
    label = 'HCO(14)(15)',
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
    label = '[CH]=C[C](O)CC=C(5195)',
    structure = SMILES('[CH]C=C(O)CC=C'),
    E0 = (208.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.570623,0.0625119,-1.96549e-05,-2.30905e-08,1.452e-11,25253.9,28.1045], Tmin=(100,'K'), Tmax=(968.947,'K')), NASAPolynomial(coeffs=[16.8039,0.0255625,-8.99666e-06,1.5988e-09,-1.12297e-13,20696.8,-56.9822], Tmin=(968.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1CC1(O)C=CC=O(5212)',
    structure = SMILES('[CH2]C1CC1(O)C=CC=O'),
    E0 = (-26.7797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00304922,0.0743264,-3.6117e-05,-1.0122e-08,1.0125e-11,-3064.64,33.9253], Tmin=(100,'K'), Tmax=(1002.72,'K')), NASAPolynomial(coeffs=[18.9761,0.0278895,-1.04416e-05,1.92213e-09,-1.3669e-13,-8342.46,-65.0279], Tmin=(1002.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.7797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclopropane) + radical(Isobutyl)"""),
)

species(
    label = 'C=CCC1(O)C=CC1[O](5213)',
    structure = SMILES('C=CCC1(O)C=CC1[O]'),
    E0 = (60.6165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.163968,0.0637581,6.51644e-06,-6.302e-08,3.10517e-11,7447.28,33.5153], Tmin=(100,'K'), Tmax=(957.38,'K')), NASAPolynomial(coeffs=[22.2774,0.0217944,-6.74467e-06,1.23179e-09,-9.30096e-14,902.058,-84.2755], Tmin=(957.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.6165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C(O)C=CC=O(5192)',
    structure = SMILES('C=C(O)C=CC=O'),
    E0 = (-243.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.939432,'amu*angstrom^2'), symmetry=1, barrier=(21.5994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.941171,'amu*angstrom^2'), symmetry=1, barrier=(21.6394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937363,'amu*angstrom^2'), symmetry=1, barrier=(21.5518,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.402004,0.0694042,-6.75524e-05,3.19153e-08,-5.84626e-12,-29143.7,23.2376], Tmin=(100,'K'), Tmax=(1337.05,'K')), NASAPolynomial(coeffs=[18.8569,0.0141932,-5.61246e-06,1.03134e-09,-7.16018e-14,-34078.7,-71.1572], Tmin=(1337.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=CCC([O])C=CC=O(4553)',
    structure = SMILES('C=CCC([O])C=CC=O'),
    E0 = (-2.23399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,291.374,291.401,291.435,291.442],'cm^-1')),
        HinderedRotor(inertia=(0.00198564,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254875,'amu*angstrom^2'), symmetry=1, barrier=(15.3622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254847,'amu*angstrom^2'), symmetry=1, barrier=(15.3615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254968,'amu*angstrom^2'), symmetry=1, barrier=(15.3621,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.198065,0.0787883,-5.8241e-05,2.08701e-08,-3.01167e-12,-128.439,34.8472], Tmin=(100,'K'), Tmax=(1600.54,'K')), NASAPolynomial(coeffs=[18.6621,0.0326434,-1.49945e-05,2.85664e-09,-1.97992e-13,-6038.88,-62.9153], Tmin=(1600.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.23399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C[CH]C(O)C=CC=O(5214)',
    structure = SMILES('C=C[CH]C(O)C=CC=O'),
    E0 = (-162.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00451375,0.0721248,-2.68885e-05,-1.69865e-08,1.09136e-11,-19380.6,36.5647], Tmin=(100,'K'), Tmax=(1072.28,'K')), NASAPolynomial(coeffs=[19.8188,0.0303591,-1.34828e-05,2.669e-09,-1.94897e-13,-25482,-69.0819], Tmin=(1072.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C=CCC(O)[C]=CC=O(5215)',
    structure = SMILES('C=CCC(O)[C]=CC=O'),
    E0 = (5.24702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.5095,0.0815244,-6.91985e-05,3.08651e-08,-5.8136e-12,752.555,34.3901], Tmin=(100,'K'), Tmax=(1215.72,'K')), NASAPolynomial(coeffs=[12.591,0.0417741,-2.01538e-05,3.97064e-09,-2.83126e-13,-2185.02,-26.2561], Tmin=(1215.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.24702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC(O)C=CC=O(5216)',
    structure = SMILES('C=[C]CC(O)C=CC=O'),
    E0 = (5.24702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.5095,0.0815244,-6.91985e-05,3.08651e-08,-5.8136e-12,752.555,34.3901], Tmin=(100,'K'), Tmax=(1215.72,'K')), NASAPolynomial(coeffs=[12.591,0.0417741,-2.01538e-05,3.97064e-09,-2.83126e-13,-2185.02,-26.2561], Tmin=(1215.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.24702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S)"""),
)

species(
    label = 'C=CCC(O)C=[C]C=O(5217)',
    structure = SMILES('C=CCC(O)C=C=C[O]'),
    E0 = (-35.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,212.005,212.006,212.006,212.006],'cm^-1')),
        HinderedRotor(inertia=(0.574888,'amu*angstrom^2'), symmetry=1, barrier=(18.336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574888,'amu*angstrom^2'), symmetry=1, barrier=(18.336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574892,'amu*angstrom^2'), symmetry=1, barrier=(18.336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574885,'amu*angstrom^2'), symmetry=1, barrier=(18.336,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.443105,0.0827813,-4.85544e-05,-5.12575e-09,1.01687e-11,-4139.45,36.3844], Tmin=(100,'K'), Tmax=(985.177,'K')), NASAPolynomial(coeffs=[22.2338,0.0235566,-8.39336e-06,1.54079e-09,-1.11134e-13,-10201.6,-80.7698], Tmin=(985.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CCC(O)C=CC=O(5218)',
    structure = SMILES('[CH]=CCC(O)C=CC=O'),
    E0 = (14.5014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,329.021,330.113],'cm^-1')),
        HinderedRotor(inertia=(0.163198,'amu*angstrom^2'), symmetry=1, barrier=(12.5449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16394,'amu*angstrom^2'), symmetry=1, barrier=(12.553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162139,'amu*angstrom^2'), symmetry=1, barrier=(12.5488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162595,'amu*angstrom^2'), symmetry=1, barrier=(12.5502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163788,'amu*angstrom^2'), symmetry=1, barrier=(12.5541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281888,0.0826936,-6.90299e-05,2.9306e-08,-5.12468e-12,1876.81,35.132], Tmin=(100,'K'), Tmax=(1324.83,'K')), NASAPolynomial(coeffs=[15.3189,0.0372935,-1.76274e-05,3.43998e-09,-2.43715e-13,-2107.51,-41.6425], Tmin=(1324.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.5014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC(O)C=C[C]=O(5219)',
    structure = SMILES('C=CCC(O)[CH]C=C=O'),
    E0 = (-100.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.185174,0.0898665,-8.22679e-05,3.94415e-08,-7.69107e-12,-11970.5,32.6022], Tmin=(100,'K'), Tmax=(1219.87,'K')), NASAPolynomial(coeffs=[16.3976,0.0354905,-1.54045e-05,2.89999e-09,-2.02214e-13,-16016.2,-50.6959], Tmin=(1219.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(C=CCJCO)"""),
)

species(
    label = 'O=CC=CC1(O)C[CH]C1(5220)',
    structure = SMILES('O=CC=CC1(O)C[CH]C1'),
    E0 = (-39.4015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.352581,0.0673587,-2.71291e-05,-8.95584e-09,6.71318e-12,-4596.81,34.2114], Tmin=(100,'K'), Tmax=(1118.52,'K')), NASAPolynomial(coeffs=[16.3127,0.0342292,-1.4814e-05,2.84446e-09,-2.02364e-13,-9665.1,-51.2706], Tmin=(1118.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.4015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = 'CO(10)(11)',
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
    label = 'C=CC[C](O)C=C(5221)',
    structure = SMILES('[CH2]C=C(O)CC=C'),
    E0 = (-10.3352,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,343.737,343.857],'cm^-1')),
        HinderedRotor(inertia=(0.210039,'amu*angstrom^2'), symmetry=1, barrier=(17.5175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20913,'amu*angstrom^2'), symmetry=1, barrier=(17.5088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207737,'amu*angstrom^2'), symmetry=1, barrier=(17.5089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208345,'amu*angstrom^2'), symmetry=1, barrier=(17.5131,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.574876,0.0600409,-1.13075e-05,-3.61043e-08,2.03683e-11,-1105.84,27.261], Tmin=(100,'K'), Tmax=(954.154,'K')), NASAPolynomial(coeffs=[19.1322,0.0194544,-5.9984e-06,1.05733e-09,-7.71983e-14,-6340.94,-70.2723], Tmin=(954.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.3352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(O)(C=C)C=CC=O(5222)',
    structure = SMILES('[CH2]C(O)(C=C)C=CC=O'),
    E0 = (-50.3973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.419882,0.0830995,-7.01121e-05,3.05854e-08,-5.59447e-12,-5936.45,32.6289], Tmin=(100,'K'), Tmax=(1252.88,'K')), NASAPolynomial(coeffs=[13.5072,0.0413163,-2.00876e-05,3.96701e-09,-2.83033e-13,-9215.83,-33.4606], Tmin=(1252.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.3973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'CHCHCHO(2768)',
    structure = SMILES('[CH]=CC=O'),
    E0 = (171.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,688.575,691.322,692.403,693.509,695.331],'cm^-1')),
        HinderedRotor(inertia=(0.0100227,'amu*angstrom^2'), symmetry=1, barrier=(3.38495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.95468,0.0145151,1.89382e-05,-3.59857e-08,1.47962e-11,20706.8,12.171], Tmin=(100,'K'), Tmax=(981.049,'K')), NASAPolynomial(coeffs=[9.63601,0.00810119,-3.09994e-06,6.30252e-10,-4.9152e-14,18393.6,-25.043], Tmin=(981.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CHCHCHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=CC[C]O(5223)',
    structure = SMILES('C=CC[C]O'),
    E0 = (284.639,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,341.29,341.29,341.29],'cm^-1')),
        HinderedRotor(inertia=(0.142942,'amu*angstrom^2'), symmetry=1, barrier=(11.815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142942,'amu*angstrom^2'), symmetry=1, barrier=(11.815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142942,'amu*angstrom^2'), symmetry=1, barrier=(11.815,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68226,0.0462647,-4.04016e-05,1.81871e-08,-3.27083e-12,34321.5,19.8353], Tmin=(100,'K'), Tmax=(1336.14,'K')), NASAPolynomial(coeffs=[11.8674,0.0157738,-6.17171e-06,1.10824e-09,-7.53014e-14,31599.8,-32.2535], Tmin=(1336.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CH2_triplet)"""),
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
    label = 'Ar(8)',
    structure = SMILES('[Ar]'),
    E0 = (-6.19426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,-745,4.3663], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,-745,4.3663], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.19426,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: BurkeH2O2"""),
)

transitionState(
    label = 'TS1',
    E0 = (87.1773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (105.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (92.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (41.0244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (250.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (154.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (289.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (158.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (292.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (216.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (330.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (220.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (220.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (244.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (256.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (269.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (291.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (328.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (64.9478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (105.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (64.9478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (49.7381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-72.5492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-72.5492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-53.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (97.6228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (79.0555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-74.7972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (28.5471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (88.6472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-30.5546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (30.9738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-46.1502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-167.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-55.2887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-51.0955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (351.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (1.05776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (67.5619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (24.8455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (63.4429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (109.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (-39.5255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (-10.0451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (-40.9509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (36.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (-57.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (-33.5458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (67.5619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (63.4429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (106.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (-39.5255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (-27.9805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (-2.69316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (-12.9582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (-63.1662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (34.4512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (119.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (106.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (-19.1403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (18.1722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (84.4818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (384.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (-22.9817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (50.9953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (22.3748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (109.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (-0.679169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (19.7709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (-33.5797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (-10.0451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (-49.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (321.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (208.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (164.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (207.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (154.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (304.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS79',
    E0 = (380.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS80',
    E0 = (266.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS81',
    E0 = (269.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS82',
    E0 = (306.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS83',
    E0 = (313.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS84',
    E0 = (86.1915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS85',
    E0 = (-88.9882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS86',
    E0 = (-43.2042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS87',
    E0 = (362.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS88',
    E0 = (-1.15796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS89',
    E0 = (412.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS90',
    E0 = (59.3584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS91',
    E0 = (-13.2287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS92',
    E0 = (67.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS93',
    E0 = (118.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS94',
    E0 = (106.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS95',
    E0 = (146.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS96',
    E0 = (110.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS97',
    E0 = (-92.6741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS98',
    E0 = (105.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS99',
    E0 = (110.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS100',
    E0 = (270.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS101',
    E0 = (100.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS102',
    E0 = (-81.0095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS103',
    E0 = (232.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS104',
    E0 = (112.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS105',
    E0 = (411.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS106',
    E0 = (-38.7927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS107',
    E0 = (100.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS108',
    E0 = (-0.676706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS109',
    E0 = (117.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS110',
    E0 = (10.1157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS111',
    E0 = (-15.0936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS112',
    E0 = (12.5867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS113',
    E0 = (-25.0368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS114',
    E0 = (43.1023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS115',
    E0 = (54.6531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS116',
    E0 = (-32.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS117',
    E0 = (104.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS118',
    E0 = (119.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS119',
    E0 = (109.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS120',
    E0 = (19.7709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS121',
    E0 = (106.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS122',
    E0 = (115.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS123',
    E0 = (339.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS124',
    E0 = (208.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS125',
    E0 = (215.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS126',
    E0 = (211.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS127',
    E0 = (158.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS128',
    E0 = (308.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS129',
    E0 = (367.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS130',
    E0 = (273.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS131',
    E0 = (269.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS132',
    E0 = (317.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS133',
    E0 = (272.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS134',
    E0 = (84.5691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS135',
    E0 = (-25.2196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS136',
    E0 = (-87.4854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS137',
    E0 = (14.2969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS138',
    E0 = (359.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS139',
    E0 = (1.74671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS140',
    E0 = (416.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS141',
    E0 = (-60.6241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS142',
    E0 = (119.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS143',
    E0 = (111.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS144',
    E0 = (10.7893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS145',
    E0 = (-7.81365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS146',
    E0 = (-30.7001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS147',
    E0 = (24.5038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS148',
    E0 = (30.3079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS149',
    E0 = (87.2092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS150',
    E0 = (-56.6576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS151',
    E0 = (241.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS152',
    E0 = (-26.7797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS153',
    E0 = (60.6165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS154',
    E0 = (72.1004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS155',
    E0 = (20.6503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS156',
    E0 = (42.9008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS157',
    E0 = (160.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS158',
    E0 = (-18.1681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS159',
    E0 = (187.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS160',
    E0 = (137.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS161',
    E0 = (102.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS162',
    E0 = (58.8099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS163',
    E0 = (3.23185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS164',
    E0 = (25.0731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS165',
    E0 = (175.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS166',
    E0 = (106.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS167',
    E0 = (456.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction109',
    reactants = ['H(3)(3)', 'C=C=CC(O)=CC=CO(3156)'],
    products = ['S(781)(780)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction139',
    reactants = ['OH(5)(5)', 'C=CC=C=CC=CO(5318)'],
    products = ['S(781)(780)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(41610,'cm^3/(mol*s)'), n=2.487, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds;OJ_pri] for rate rule [Ca_Cds-CdH;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from -7.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction149',
    reactants = ['H(3)(3)', 'C=CC=C(O)C=C=CO(3134)'],
    products = ['S(781)(780)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction159',
    reactants = ['H(3)(3)', 'S(778)(777)'],
    products = ['S(781)(780)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(20.92,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2834 used for Od_CO-CdH;HJ
Exact match found for rate rule [Od_CO-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction117',
    reactants = ['OH(5)(5)', '[CH2]C=C[C]=CC=CO(5308)'],
    products = ['S(781)(780)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad/Cd;Y_rad] for rate rule [Cd_rad/Cd;O_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction118',
    reactants = ['H(3)(3)', 'C=C[CH]C(=O)C=C[CH]O(2830)'],
    products = ['S(781)(780)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [O_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction119',
    reactants = ['OH(5)(5)', '[CH]=CC=C(O)C=C[CH2](2092)'],
    products = ['S(781)(780)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.44289e+07,'m^3/(mol*s)'), n=0.0225, Ea=(0.0523,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;Cd_pri_rad] + [O_pri_rad;Y_rad] for rate rule [O_pri_rad;Cd_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction120',
    reactants = ['H(3)(3)', 'S(779)(778)'],
    products = ['S(781)(780)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction121',
    reactants = ['[CH]=C[CH2](891)', 'O[C]=CC=CO(5309)'],
    products = ['S(781)(780)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad/NonDe;Y_rad] for rate rule [Cd_rad/NonDe;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction122',
    reactants = ['H(3)(3)', 'C=C[C]=C(O)C=C[CH]O(3123)'],
    products = ['S(781)(780)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction123',
    reactants = ['CHCHOH(49)(49)', '[CH]=C(O)C=C[CH2](3149)'],
    products = ['S(781)(780)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 89 used for Cd_pri_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction124',
    reactants = ['H(3)(3)', 'C=CC=C(O)[C]=C[CH]O(3113)'],
    products = ['S(781)(780)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction125',
    reactants = ['H(3)(3)', '[CH2]C=CC(O)=C[C]=CO(3147)'],
    products = ['S(781)(780)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction126',
    reactants = ['H(3)(3)', 'C=[C]C=C(O)C=C[CH]O(3126)'],
    products = ['S(781)(780)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction127',
    reactants = ['H(3)(3)', '[CH2]C=CC(O)=CC=[C]O(3145)'],
    products = ['S(781)(780)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)(3)', '[CH]=C[CH]C(O)=CC=CO(3163)'],
    products = ['S(781)(780)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C2H3(28)(29)', '[CH]C(O)=CC=CO(5321)'],
    products = ['S(781)(780)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using an average for rate rule [Cd_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=CC=[C]O(3098)', '[CH]=C[CH]O(3569)'],
    products = ['S(781)(780)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/NonDe] for rate rule [Cd_rad;Cd_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction131',
    reactants = ['S(781)(780)'],
    products = ['[CH2]CC=C(O)C=C=CO(5313)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.92857e+07,'s^-1'), n=1.56, Ea=(259.91,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for 1_3_pentadiene;CH=C_1;CdCJ_2
Exact match found for rate rule [1_3_pentadiene;CH=C_1;CdCJ_2]
Euclidian distance = 0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction132',
    reactants = ['[CH2]C=C=C(O)C=CCO(5314)'],
    products = ['S(781)(780)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.6907e+09,'s^-1'), n=1.02346, Ea=(163.761,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH_end;CddC_2] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH_end;CddC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction166',
    reactants = ['S(781)(780)'],
    products = ['C=C=CC(O)=CC[CH]O(5335)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.92857e+07,'s^-1'), n=1.56, Ea=(259.91,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for 1_3_pentadiene;CH=C_1;CdCJ_2
Exact match found for rate rule [1_3_pentadiene;CH=C_1;CdCJ_2]
Euclidian distance = 0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC=CC(O)=[C]C=CO(5304)'],
    products = ['S(781)(780)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.53605e+09,'s^-1'), n=1.02346, Ea=(163.761,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH_end;CddC_2] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH3_1;CddC_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction135',
    reactants = ['S(781)(780)'],
    products = ['[CH2]C=CC1(O)C=CC1O(5315)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1,3-butadiene_backbone;C=C_1;C=C_2]
Euclidian distance = 0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction24',
    reactants = ['S(781)(780)'],
    products = ['OC=C[CH]C1(O)C=CC1(5299)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;C=C_1;CdH2_2]
Euclidian distance = 1.0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction108',
    reactants = ['S(781)(780)'],
    products = ['O[CH]C1C=C(O)C=CC1(5300)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.21347e+10,'s^-1'), n=0.43, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSM_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R7_SMSM_D;doublebond_intra_HNd_pri;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction110',
    reactants = ['C[C]=CC(O)=CC=CO(5301)'],
    products = ['S(781)(780)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.32e+09,'s^-1'), n=1.12, Ea=(172.799,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 152 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction111',
    reactants = ['CC=[C]C(O)=CC=CO(5302)'],
    products = ['S(781)(780)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 288 used for R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction112',
    reactants = ['CC=CC([O])=CC=CO(5303)'],
    products = ['S(781)(780)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.33e+06,'s^-1'), n=1.64, Ea=(100.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SS(D)MS;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_SS(D)MS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction114',
    reactants = ['CC=CC(O)=C[C]=CO(5305)'],
    products = ['S(781)(780)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6H_SMSMS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction115',
    reactants = ['CC=CC(O)=CC=[C]O(5306)'],
    products = ['S(781)(780)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleNd;Cs_H_out_2H] for rate rule [R7H;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction116',
    reactants = ['CC=CC(O)=CC=C[O](5307)'],
    products = ['S(781)(780)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.81e+06,'s^-1'), n=1.5, Ea=(141.001,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8H_SDSDSD;Y_rad_out;Cs_H_out] for rate rule [R8H_SDSDSD;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction128',
    reactants = ['S(781)(780)'],
    products = ['OC=CC=C(O)C1[CH]C1(5310)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 142 used for R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction129',
    reactants = ['S(781)(780)'],
    products = ['OC=CC1CC=C[C]1O(5311)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.41935e+11,'s^-1'), n=0.104408, Ea=(148.812,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R5_SD_D;doublebond_intra_secNd_HCd;radadd_intra_cs2H]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction130',
    reactants = ['S(781)(780)'],
    products = ['OC1C=CCC(O)[CH]C=1(5312)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_HNd;radadd_intra_cs2H] for rate rule [R7_SDSD_D;doublebond_intra_pri_HNd_O;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction133',
    reactants = ['S(781)(780)'],
    products = ['C=C[CH]C(=O)CC=CO(5175)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction134',
    reactants = ['S(781)(780)'],
    products = ['C=C[CH]C(O)=CCC=O(5207)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction136',
    reactants = ['CH2(17)(18)', '[CH]=CC(O)=CC=CO(5316)'],
    products = ['S(781)(780)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction138',
    reactants = ['S(781)(780)'],
    products = ['C=CC1C(O)=CC1[CH]O(5317)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra_HNd_pri;radadd_intra_cs] for rate rule [R5_SD_D;doublebond_intra_HNd_pri;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction140',
    reactants = ['C=[C]CC(O)=CC=CO(5164)'],
    products = ['S(781)(780)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.3e+10,'s^-1'), n=0.98, Ea=(112.55,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['S(542)(541)'],
    products = ['S(781)(780)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(587605,'s^-1'), n=2.09, Ea=(169.87,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_2Cd;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction141',
    reactants = ['C=CCC(O)=[C]C=CO(5163)'],
    products = ['S(781)(780)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction142',
    reactants = ['[CH]=CCC(O)=CC=CO(5166)'],
    products = ['S(781)(780)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction143',
    reactants = ['C=CCC(O)=C[C]=CO(5165)'],
    products = ['S(781)(780)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction144',
    reactants = ['C=CCC(O)=CC=[C]O(5167)'],
    products = ['S(781)(780)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_single;Cs_H_out] for rate rule [R5H_DSMS;Cd_rad_out_singleNd;Cs_H_out_H/Cd]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['S(568)(567)'],
    products = ['S(781)(780)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(500000,'s^-1'), n=1.95, Ea=(100.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SMSMS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R6H_SMSMS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction145',
    reactants = ['S(781)(780)'],
    products = ['C=CC1[C](O)C1C=CO(5319)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra_secNd_HCd;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction146',
    reactants = ['S(781)(780)'],
    products = ['C=CC1C(O)=C[CH]C1O(5320)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(8.8675e+11,'s^-1'), n=-0.37996, Ea=(137.453,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra_pri;radadd_intra_csHCd] for rate rule [R5_SD_D;doublebond_intra_pri_HNd_O;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction148',
    reactants = ['S(781)(780)'],
    products = ['[CH2]C1C=C(O)C1C=CO(5322)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.17094e+10,'s^-1'), n=0.2405, Ea=(161.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SM;doublebond_intra_2H_pri;radadd_intra_cs] + [R5_SD_D;doublebond_intra;radadd_intra_cs] for rate rule [R5_SD_D;doublebond_intra_2H_pri;radadd_intra_csHCd]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction150',
    reactants = ['C=CC=C(O)C[C]=CO(5323)'],
    products = ['S(781)(780)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.3e+10,'s^-1'), n=0.98, Ea=(112.55,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction151',
    reactants = ['C=C[C]=C(O)CC=CO(5324)'],
    products = ['S(781)(780)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction152',
    reactants = ['C=CC=C(O)CC=[C]O(5325)'],
    products = ['S(781)(780)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleNd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction153',
    reactants = ['C=[C]C=C(O)CC=CO(5326)'],
    products = ['S(781)(780)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction154',
    reactants = ['C=CC=C(O)CC=C[O](5327)'],
    products = ['S(781)(780)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(2.08e+07,'s^-1'), n=1.61, Ea=(113.386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction155',
    reactants = ['[CH]=CC=C(O)CC=CO(5328)'],
    products = ['S(781)(780)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R5H_DSMS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction156',
    reactants = ['S(781)(780)'],
    products = ['C=CC=C(O)C1[CH]C1O(5329)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(2.22857e+12,'s^-1'), n=0, Ea=(182.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;doublebond_intra_pri;radadd_intra_csH(CdCdCd)] for rate rule [R3_D;doublebond_intra_pri_HNd_O;radadd_intra_csH(CdCdCd)]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction158',
    reactants = ['S(781)(780)'],
    products = ['[CH2]C1C=C(O)C=CC1O(5330)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(3.95e+10,'s^-1'), n=0.53, Ea=(131.796,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSM_D;doublebond_intra_2H_pri;radadd_intra_cs] for rate rule [R7_SMSM_D;doublebond_intra_2H_pri;radadd_intra_csHNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction160',
    reactants = ['C=CC=C(O)C=CC[O](5196)'],
    products = ['S(781)(780)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(2.36e+07,'s^-1'), n=1.8, Ea=(66.1072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(Cd-Cd-Cd-Cd-Cd)] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/(Cd-Cd-Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction161',
    reactants = ['C=CC=C(O)C=[C]CO(5331)'],
    products = ['S(781)(780)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.89098e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction162',
    reactants = ['C=CC=C(O)[C]=CCO(5332)'],
    products = ['S(781)(780)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(2.56742e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction163',
    reactants = ['C=C[CH]C(=O)C=CCO(5187)'],
    products = ['S(781)(780)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1'), n=1.64, Ea=(100.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SS(D)MS;Y_rad_out;Cs_H_out] for rate rule [R5H_SS(D)MS;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction164',
    reactants = ['C=[C]C=C(O)C=CCO(5333)'],
    products = ['S(781)(780)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(2.06213e+06,'s^-1'), n=1.71203, Ea=(99.1471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;Cs_H_out_1H] + [R6H;Y_rad_out;Cs_H_out_H/NonDeO] + [R6H_SMSMS;Y_rad_out;Cs_H_out_1H] for rate rule [R6H_SMSMS;Cd_rad_out_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction165',
    reactants = ['[CH]=CC=C(O)C=CCO(5334)'],
    products = ['S(781)(780)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(3.70099e+07,'s^-1'), n=1.485, Ea=(94.7468,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R7H;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction167',
    reactants = ['HCOH(T)(1710)', '[CH]=CC(O)=CC=C(3102)'],
    products = ['S(781)(780)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['S(542)(541)'],
    products = ['[CH2]C1CC(=CC=CO)O1(5161)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 33 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['S(542)(541)'],
    products = ['C=CCC1=CC([CH]O)O1(5162)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R5_SD_D;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction1',
    reactants = ['C3H5(89)(88)', 'O=C=CC=CO(5155)'],
    products = ['S(542)(541)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(0.000357827,'m^3/(mol*s)'), n=2.74787, Ea=(45.8259,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CsJ-CdHH] for rate rule [Ck_O;CsJ-CdHH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=CCC(O)=[C]C=CO(5163)'],
    products = ['S(542)(541)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]CC(O)=CC=CO(5164)'],
    products = ['S(542)(541)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=CCC(O)=C[C]=CO(5165)'],
    products = ['S(542)(541)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(1.11e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['S(542)(541)'],
    products = ['[CH]=CCC(O)=CC=CO(5166)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=CCC(O)=CC=[C]O(5167)'],
    products = ['S(542)(541)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_single;XH_out] for rate rule [R5H_DSMS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['S(568)(567)'],
    products = ['S(542)(541)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(14790.5,'s^-1'), n=1.91656, Ea=(92.3659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6H;O_rad_out;XH_out] + [R6H_SMSMS;Y_rad_out;XH_out] for rate rule [R6H_SMSMS;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['OH(5)(5)', '[CH]=CC=C([O])CC=C(2414)'],
    products = ['S(542)(541)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(3.44289e+07,'m^3/(mol*s)'), n=0.0225, Ea=(0.0523,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;Cd_pri_rad] + [O_pri_rad;Y_rad] for rate rule [O_pri_rad;Cd_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)(3)', 'C=CC[C]([O])C=CC=O(2769)'],
    products = ['S(542)(541)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['C3H5(89)(88)', '[O][C]=CC=CO(5156)'],
    products = ['S(542)(541)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(3.13324e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H3(28)(29)', '[CH2]C([O])=CC=CO(5157)'],
    products = ['S(542)(541)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;Cd_rad] for rate rule [C_rad/H2/Cd;Cd_pri_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)(3)', 'C=C[CH]C(=O)C=C[CH]O(2830)'],
    products = ['S(542)(541)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(4.14e-13,'cm^3/(molecule*s)'), n=0.611, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 53 used for C_rad/H/CdCd;H_rad
Exact match found for rate rule [C_rad/H/CdCd;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.8 to -1.8 kJ/mol.
Ea raised from -1.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)(3)', 'C=[C]CC(=O)C=C[CH]O(2833)'],
    products = ['S(542)(541)'],
    transitionState = 'TS78',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CHCHOH(49)(49)', '[CH]=C([O])CC=C(4101)'],
    products = ['S(542)(541)'],
    transitionState = 'TS79',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using an average for rate rule [Cd_pri_rad;Cd_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)(3)', 'C=CCC(=O)[C]=C[CH]O(2818)'],
    products = ['S(542)(541)'],
    transitionState = 'TS80',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)(3)', 'C=CCC([O])=C[C]=CO(5158)'],
    products = ['S(542)(541)'],
    transitionState = 'TS81',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)(3)', 'C=CCC([O])=CC=[C]O(5159)'],
    products = ['S(542)(541)'],
    transitionState = 'TS82',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)(3)', '[CH]=CCC([O])=CC=CO(5160)'],
    products = ['S(542)(541)'],
    transitionState = 'TS83',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['S(542)(541)'],
    products = ['C=CC[C]1OC1C=CO(5168)'],
    transitionState = 'TS84',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_HCd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['S(542)(541)'],
    products = ['OC=CC=C1C[CH]CO1(5169)'],
    transitionState = 'TS85',
    kinetics = Arrhenius(A=(2.07178e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['S(542)(541)'],
    products = ['C=CCC1=C[CH]C(O)O1(5170)'],
    transitionState = 'TS86',
    kinetics = Arrhenius(A=(2.96016e+10,'s^-1'), n=0.2756, Ea=(101.82,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra_pri_HNd;radadd_intra] for rate rule [R5_SD_D;doublebond_intra_pri_HNd_O;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=CC=C(CC=C)OO(5171)'],
    products = ['S(542)(541)'],
    transitionState = 'TS87',
    kinetics = Arrhenius(A=(9.91772e+09,'s^-1'), n=0, Ea=(65.2704,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OOH;Y_rad_out] for rate rule [R4OOH_DSD;Cd_rad_out_H]
Euclidian distance = 2.2360679775
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['S(542)(541)'],
    products = ['C=CCC(=O)[CH]CC=O(4556)'],
    transitionState = 'TS88',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O(4)(4)', 'C=CC[C]=CC=CO(5172)'],
    products = ['S(542)(541)'],
    transitionState = 'TS89',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction42',
    reactants = ['S(542)(541)'],
    products = ['C=CCC1([O])C=CC1O(5183)'],
    transitionState = 'TS90',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(204.383,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SD_CO;carbonylbond_intra_Nd;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 201.0 to 204.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction43',
    reactants = ['S(542)(541)'],
    products = ['[CH2]C1CC(=O)C=CC1O(5184)'],
    transitionState = 'TS91',
    kinetics = Arrhenius(A=(3.95e+10,'s^-1'), n=0.53, Ea=(131.796,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSR;doublebond_intra_2H_pri;radadd_intra_cs] for rate rule [R7_SMSS_D;doublebond_intra_2H_pri;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)(3)', 'S(459)(458)'],
    products = ['S(542)(541)'],
    transitionState = 'TS92',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(20.92,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2834 used for Od_CO-CdH;HJ
Exact match found for rate rule [Od_CO-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)(3)', 'C=CCC(=O)C=C=CO(2839)'],
    products = ['S(542)(541)'],
    transitionState = 'TS93',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=CCC(=O)C=CC[O](4558)'],
    products = ['S(542)(541)'],
    transitionState = 'TS94',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=CCC(=O)C=[C]CO(5185)'],
    products = ['S(542)(541)'],
    transitionState = 'TS95',
    kinetics = Arrhenius(A=(1.89098e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=CCC(=O)[C]=CCO(5186)'],
    products = ['S(542)(541)'],
    transitionState = 'TS96',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe;Cs_H_out_H/NonDeO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C[CH]C(=O)C=CCO(5187)'],
    products = ['S(542)(541)'],
    transitionState = 'TS97',
    kinetics = Arrhenius(A=(0.0756983,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;C_rad_out_1H;Cs_H_out_1H] for rate rule [R5H_SSMS;C_rad_out_H/Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=[C]CC(=O)C=CCO(5188)'],
    products = ['S(542)(541)'],
    transitionState = 'TS98',
    kinetics = Arrhenius(A=(4.18783e+06,'s^-1'), n=1.59304, Ea=(98.5126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;Cs_H_out_1H] + [R6H;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R6H_RSSMS;Cd_rad_out_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=CCC(=O)C=CCO(5189)'],
    products = ['S(542)(541)'],
    transitionState = 'TS99',
    kinetics = Arrhenius(A=(3.70099e+07,'s^-1'), n=1.485, Ea=(94.7468,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R7H;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C4H5O(93)(92)', '[CH]=C[CH]O(3569)'],
    products = ['S(542)(541)'],
    transitionState = 'TS100',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;CO_rad/NonDe] for rate rule [Cd_rad;CO_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['S(542)(541)'],
    products = ['C=CCC(=O)C1[CH]C1O(5179)'],
    transitionState = 'TS101',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;doublebond_intra_pri_HNd;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra_pri_HNd_O;radadd_intra_csHCO]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction50',
    reactants = ['S(542)(541)'],
    products = ['O=C1C=CC(O)C[CH]C1(5190)'],
    transitionState = 'TS102',
    kinetics = Arrhenius(A=(4.39632e+07,'s^-1'), n=0.77, Ea=(64.0152,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_2H;radadd_intra_csHNd] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_csHO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=COC(=C)[CH]C=CO(5182)'],
    products = ['S(542)(541)'],
    transitionState = 'TS103',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=CCC([CH]O)C=C=O(5191)'],
    products = ['S(542)(541)'],
    transitionState = 'TS104',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=CC(=O)CC=C(2806)', 'HCOH(T)(1710)'],
    products = ['S(542)(541)'],
    transitionState = 'TS105',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['S(542)(541)'],
    products = ['[CH2]C1CC(=O)C1C=CO(5173)'],
    transitionState = 'TS106',
    kinetics = Arrhenius(A=(2.15e+11,'s^-1'), n=0.21, Ea=(106.232,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 159 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_csHCd
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=CCC(=O)C[C]=CO(5174)'],
    products = ['S(542)(541)'],
    transitionState = 'TS107',
    kinetics = Arrhenius(A=(1.82494e+10,'s^-1'), n=0.9, Ea=(134.725,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/OneDe] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C[CH]C(=O)CC=CO(5175)'],
    products = ['S(542)(541)'],
    transitionState = 'TS108',
    kinetics = Arrhenius(A=(3.8e+10,'s^-1'), n=0.87, Ea=(144.348,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_H/Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=CCC(=O)CC=[C]O(5176)'],
    products = ['S(542)(541)'],
    transitionState = 'TS109',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleNd;Cs_H_out_H/CO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=[C]CC(=O)CC=CO(5177)'],
    products = ['S(542)(541)'],
    transitionState = 'TS110',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=CCC(=O)C[CH]C=O(4555)'],
    products = ['S(542)(541)'],
    transitionState = 'TS111',
    kinetics = Arrhenius(A=(2.30755e+06,'s^-1'), n=1.805, Ea=(115.478,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=CCC(=O)CC=CO(5178)'],
    products = ['S(542)(541)'],
    transitionState = 'TS112',
    kinetics = Arrhenius(A=(14044.1,'s^-1'), n=1.88125, Ea=(37.5252,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['S(542)(541)'],
    products = ['O=C1C[CH]CC1C=CO(5180)'],
    transitionState = 'TS113',
    kinetics = Arrhenius(A=(2.16456e+11,'s^-1'), n=-0.0435858, Ea=(119.988,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_linear;doublebond_intra_pri_2H;radadd_intra_csHCd] + [R5_SS_D;doublebond_intra;radadd_intra_csHCd] + [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=CCC([C]=O)C=CO(5181)'],
    products = ['S(542)(541)'],
    transitionState = 'TS114',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction19',
    reactants = ['S(568)(567)'],
    products = ['C=CC[C](O)C1C=CO1(5197)'],
    transitionState = 'TS115',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra] for rate rule [R5_SD_D;doublebond_intra_NdNd_pri;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['S(568)(567)'],
    products = ['[CH2]C1CC(O)=CC=CO1(5198)'],
    transitionState = 'TS116',
    kinetics = Arrhenius(A=(3.19557e+09,'s^-1'), n=0.501475, Ea=(108.906,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra_2H_pri;radadd_intra] for rate rule [R8;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction1',
    reactants = ['H(3)(3)', 'C=CCC(O)=CC=C=O(3131)'],
    products = ['S(568)(567)'],
    transitionState = 'TS117',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=CCC(O)=CC=[C]O(5167)'],
    products = ['S(568)(567)'],
    transitionState = 'TS118',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=CCC(O)=C[C]=CO(5165)'],
    products = ['S(568)(567)'],
    transitionState = 'TS119',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=CCC(O)=[C]C=CO(5163)'],
    products = ['S(568)(567)'],
    transitionState = 'TS120',
    kinetics = Arrhenius(A=(1.11e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C]CC(O)=CC=CO(5164)'],
    products = ['S(568)(567)'],
    transitionState = 'TS121',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R7H;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=CCC(O)=CC=CO(5166)'],
    products = ['S(568)(567)'],
    transitionState = 'TS122',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R8H;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OH(5)(5)', 'C=CC[C]=CC=C[O](2415)'],
    products = ['S(568)(567)'],
    transitionState = 'TS123',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/NonDe] for rate rule [O_pri_rad;Cd_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)(3)', 'C=CC[C]([O])C=CC=O(2769)'],
    products = ['S(568)(567)'],
    transitionState = 'TS124',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['C3H5(89)(88)', '[O]C=CC=[C]O(3148)'],
    products = ['S(568)(567)'],
    transitionState = 'TS125',
    kinetics = Arrhenius(A=(1.28206e+10,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;C_pri_rad] for rate rule [Cd_rad/NonDe;C_rad/H2/Cd]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C2H3(28)(29)', '[CH2]C(O)=CC=C[O](5193)'],
    products = ['S(568)(567)'],
    transitionState = 'TS126',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;Cd_rad] for rate rule [C_rad/H2/Cd;Cd_pri_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)(3)', 'S(779)(778)'],
    products = ['S(568)(567)'],
    transitionState = 'TS127',
    kinetics = Arrhenius(A=(4.14e-13,'cm^3/(molecule*s)'), n=0.611, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 53 used for C_rad/H/CdCd;H_rad
Exact match found for rate rule [C_rad/H/CdCd;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.8 to -1.8 kJ/mol.
Ea raised from -1.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)(3)', 'C=[C]C[C](O)C=CC=O(2815)'],
    products = ['S(568)(567)'],
    transitionState = 'TS128',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CHCHO(45)(45)', '[CH]=C(O)CC=C(5194)'],
    products = ['S(568)(567)'],
    transitionState = 'TS129',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using an average for rate rule [Cd_pri_rad;Cd_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)(3)', 'C=CCC(O)=[C]C=C[O](3162)'],
    products = ['S(568)(567)'],
    transitionState = 'TS130',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)(3)', 'C=CC[C](O)C=[C]C=O(2817)'],
    products = ['S(568)(567)'],
    transitionState = 'TS131',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)(3)', '[CH]=CC[C](O)C=CC=O(2823)'],
    products = ['S(568)(567)'],
    transitionState = 'TS132',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)(3)', 'C=CC[C](O)C=C[C]=O(2824)'],
    products = ['S(568)(567)'],
    transitionState = 'TS133',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['S(568)(567)'],
    products = ['C=CCC(O)=CC1[CH]O1(5199)'],
    transitionState = 'TS134',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri_HCd;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['S(568)(567)'],
    products = ['C=CCC1(O)[CH]C=CO1(5200)'],
    transitionState = 'TS135',
    kinetics = Arrhenius(A=(1.56048e+11,'s^-1'), n=0.0862073, Ea=(116.147,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SD_D;doublebond_intra_pri_NdNd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['S(568)(567)'],
    products = ['OC1=CC=COC[CH]C1(5201)'],
    transitionState = 'TS136',
    kinetics = Arrhenius(A=(1.06679e+09,'s^-1'), n=0.473387, Ea=(53.8815,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_2H;radadd_intra] for rate rule [R8_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 47.9 to 53.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=CC=C(O)C=CC[O](5196)'],
    products = ['S(568)(567)'],
    transitionState = 'TS137',
    kinetics = Arrhenius(A=(4.25e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1_3_pentadiene;CH_end;CdHC_2]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=CC[C]=CC=COO(5202)'],
    products = ['S(568)(567)'],
    transitionState = 'TS138',
    kinetics = Arrhenius(A=(9.91772e+09,'s^-1'), n=0, Ea=(65.2704,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OOH;Y_rad_out] for rate rule [R4OOH_DSD;Cd_rad_out_ND]
Euclidian distance = 2.2360679775
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['S(568)(567)'],
    products = ['C=CCC(=O)C[CH]C=O(4555)'],
    transitionState = 'TS139',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CsC;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O(4)(4)', '[CH]=CC=C(O)CC=C(5203)'],
    products = ['S(568)(567)'],
    transitionState = 'TS140',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['S(568)(567)'],
    products = ['[CH2]C1CC(O)=CC1C=O(5204)'],
    transitionState = 'TS141',
    kinetics = Arrhenius(A=(5.249e+08,'s^-1'), n=0.846, Ea=(80.7428,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS_D;doublebond_intra_2H_pri;radadd_intra_cs] for rate rule [R6_SMS_D;doublebond_intra_2H_pri;radadd_intra_csHDe]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)(3)', 'C=CCC(O)=C=CC=O(2838)'],
    products = ['S(568)(567)'],
    transitionState = 'TS142',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=CCC(O)=[C]CC=O(5205)'],
    products = ['S(568)(567)'],
    transitionState = 'TS143',
    kinetics = Arrhenius(A=(1.82494e+10,'s^-1'), n=0.9, Ea=(134.725,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/OneDe] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=CCC(O)=CC[C]=O(5206)'],
    products = ['S(568)(567)'],
    transitionState = 'TS144',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=CCC(=O)[CH]CC=O(4556)'],
    products = ['S(568)(567)'],
    transitionState = 'TS145',
    kinetics = Arrhenius(A=(2.30755e+06,'s^-1'), n=1.805, Ea=(115.478,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['S(568)(567)'],
    products = ['C=C[CH]C(O)=CCC=O(5207)'],
    transitionState = 'TS146',
    kinetics = Arrhenius(A=(2.86616e+06,'s^-1'), n=1.7775, Ea=(110.667,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H;C_rad_out_H/OneDe;Cs_H_out_H/Cd] + [R4H_SDS;C_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_SDS;C_rad_out_H/OneDe;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C]CC(O)=CCC=O(5208)'],
    products = ['S(568)(567)'],
    transitionState = 'TS147',
    kinetics = Arrhenius(A=(127565,'s^-1'), n=1.89242, Ea=(47.7589,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;Y_rad_out;Cs_H_out_H/OneDe] + [R5H_RSMS;Cd_rad_out;Cs_H_out] for rate rule [R5H_SSMS;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 3.74165738677
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=CCC(O)=CCC=O(5209)'],
    products = ['S(568)(567)'],
    transitionState = 'TS148',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R6H_RSSMS;Cd_rad_out_singleH;Cs_H_out_H/CO]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['S(568)(567)'],
    products = ['C=CCC1(O)[CH]C1C=O(5210)'],
    transitionState = 'TS149',
    kinetics = Arrhenius(A=(1.90962e+10,'s^-1'), n=0.7195, Ea=(228.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3_D;doublebond_intra;radadd_intra_csHDe] + [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri_NdNd;radadd_intra_csHCO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['S(568)(567)'],
    products = ['O=CC1C=C(O)C[CH]C1(5211)'],
    transitionState = 'TS150',
    kinetics = Arrhenius(A=(5.041e+08,'s^-1'), n=0.7, Ea=(84.7093,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SDS_D;doublebond_intra_pri_2H;radadd_intra_cs] for rate rule [R6_SDS_D;doublebond_intra_pri_2H;radadd_intra_csHCO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction90',
    reactants = ['HCO(14)(15)', '[CH]=C[C](O)CC=C(5195)'],
    products = ['S(568)(567)'],
    transitionState = 'TS151',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction43',
    reactants = ['S(568)(567)'],
    products = ['[CH2]C1CC1(O)C=CC=O(5212)'],
    transitionState = 'TS152',
    kinetics = Arrhenius(A=(7.99e+10,'s^-1'), n=0.21, Ea=(114.587,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 104 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_csNdCd
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_csNdCd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 113.2 to 114.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction44',
    reactants = ['S(568)(567)'],
    products = ['C=CCC1(O)C=CC1[O](5213)'],
    transitionState = 'TS153',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(201.983,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SD_CO;carbonylbond_intra_H;radadd_intra_csNdNd]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic
Ea raised from 199.4 to 202.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)(3)', 'S(459)(458)'],
    products = ['S(568)(567)'],
    transitionState = 'TS154',
    kinetics = Arrhenius(A=(6.39e+06,'cm^3/(mol*s)'), n=2.09, Ea=(25.5224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2835 used for Od_CO-CdCs;HJ
Exact match found for rate rule [Od_CO-CdCs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)(3)', 'S(778)(777)'],
    products = ['S(568)(567)'],
    transitionState = 'TS155',
    kinetics = Arrhenius(A=(57.8853,'m^3/(mol*s)'), n=1.77227, Ea=(0.545964,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-OneDe;HJ] for rate rule [Cds-CdH_Cds-CdOs;HJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H3(28)(29)', 'C=C(O)C=CC=O(5192)'],
    products = ['S(568)(567)'],
    transitionState = 'TS156',
    kinetics = Arrhenius(A=(0.0145407,'m^3/(mol*s)'), n=2.409, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-OneDe;CdsJ-H] for rate rule [Cds-HH_Cds-CdOs;CdsJ-H]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from -2.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=CCC([O])C=CC=O(4553)'],
    products = ['S(568)(567)'],
    transitionState = 'TS157',
    kinetics = Arrhenius(A=(1.52488e+09,'s^-1'), n=1.21745, Ea=(162.572,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cd] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['S(568)(567)'],
    products = ['C=C[CH]C(O)C=CC=O(5214)'],
    transitionState = 'TS158',
    kinetics = Arrhenius(A=(11288.6,'s^-1'), n=2.7035, Ea=(123.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_OneDe/O;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=CCC(O)[C]=CC=O(5215)'],
    products = ['S(568)(567)'],
    transitionState = 'TS159',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_NDMustO]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=[C]CC(O)C=CC=O(5216)'],
    products = ['S(568)(567)'],
    transitionState = 'TS160',
    kinetics = Arrhenius(A=(1.1e+10,'s^-1'), n=0.78, Ea=(132.633,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_Cd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=CCC(O)C=[C]C=O(5217)'],
    products = ['S(568)(567)'],
    transitionState = 'TS161',
    kinetics = Arrhenius(A=(2.20504e+10,'s^-1'), n=0.725, Ea=(137.863,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_NonDe] for rate rule [R3H_DS;Cd_rad_out_singleDe;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=CCC(O)C=CC=O(5218)'],
    products = ['S(568)(567)'],
    transitionState = 'TS162',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_Cd]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=CCC(O)C=C[C]=O(5219)'],
    products = ['S(568)(567)'],
    transitionState = 'TS163',
    kinetics = Arrhenius(A=(5.829e+09,'s^-1'), n=0.7425, Ea=(104.027,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Y_rad_out;Cs_H_out_NDMustO] for rate rule [R4H_SDS;CO_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['S(568)(567)'],
    products = ['O=CC=CC1(O)C[CH]C1(5220)'],
    transitionState = 'TS164',
    kinetics = Arrhenius(A=(4.38e+10,'s^-1'), n=0.19, Ea=(166.44,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 212 used for R4_Cs_HH_D;doublebond_intra_pri_2H;radadd_intra_csNdCd
Exact match found for rate rule [R4_Cs_HH_D;doublebond_intra_pri_2H;radadd_intra_csNdCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction53',
    reactants = ['CO(10)(11)', 'C=CC[C](O)C=C(5221)'],
    products = ['S(568)(567)'],
    transitionState = 'TS165',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]C(O)(C=C)C=CC=O(5222)'],
    products = ['S(568)(567)'],
    transitionState = 'TS166',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction55',
    reactants = ['CHCHCHO(2768)', 'C=CC[C]O(5223)'],
    products = ['S(568)(567)'],
    transitionState = 'TS167',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '358',
    isomers = [
        'S(542)(541)',
        'S(568)(567)',
        'S(781)(780)',
    ],
    reactants = [
        ('H(3)(3)', 'S(778)(777)'),
        ('H(3)(3)', 'S(779)(778)'),
        ('H(3)(3)', 'S(459)(458)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '358',
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

