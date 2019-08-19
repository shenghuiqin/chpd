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

network(
    label = '358',
    isomers = [
        'S(781)(780)',
    ],
    reactants = [
        ('H(3)(3)', 'S(778)(777)'),
        ('H(3)(3)', 'S(779)(778)'),
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

