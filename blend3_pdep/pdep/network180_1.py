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
    label = 'C=[C]C=C(O)C=CC=O(3104)',
    structure = SMILES('C=C=CC(O)=CC=C[O]'),
    E0 = (5.04892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.22942,'amu*angstrom^2'), symmetry=1, barrier=(28.2667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23149,'amu*angstrom^2'), symmetry=1, barrier=(28.3144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23253,'amu*angstrom^2'), symmetry=1, barrier=(28.3383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52894,0.10706,-0.00011368,5.55339e-08,-9.892e-12,873.82,35.4661], Tmin=(100,'K'), Tmax=(1638.55,'K')), NASAPolynomial(coeffs=[30.0639,0.00361713,2.87535e-06,-7.81909e-10,5.7262e-14,-6601.69,-128.089], Tmin=(1638.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.04892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'C=CC=C(O)C=C[C]=O(3105)',
    structure = SMILES('[CH2]C=CC(O)=CC=C=O'),
    E0 = (-51.6768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01923,'amu*angstrom^2'), symmetry=1, barrier=(23.4341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01929,'amu*angstrom^2'), symmetry=1, barrier=(23.4354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01898,'amu*angstrom^2'), symmetry=1, barrier=(23.4283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01869,'amu*angstrom^2'), symmetry=1, barrier=(23.4217,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.445392,0.0880923,-8.90038e-05,4.58584e-08,-9.23478e-12,-6047.02,31.5086], Tmin=(100,'K'), Tmax=(1218.81,'K')), NASAPolynomial(coeffs=[19.757,0.0217905,-7.40559e-06,1.22575e-09,-7.98258e-14,-10971.6,-69.9537], Tmin=(1218.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.6768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=C(O)[C]=CC=O(3101)',
    structure = SMILES('[CH2]C=CC(O)=C=CC=O'),
    E0 = (-20.8509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.12168,'amu*angstrom^2'), symmetry=1, barrier=(25.7896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12216,'amu*angstrom^2'), symmetry=1, barrier=(25.8008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12256,'amu*angstrom^2'), symmetry=1, barrier=(25.8099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12184,'amu*angstrom^2'), symmetry=1, barrier=(25.7934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.40369,0.0893876,-8.80957e-05,4.33316e-08,-8.38369e-12,-2342.97,30.3481], Tmin=(100,'K'), Tmax=(1256.58,'K')), NASAPolynomial(coeffs=[20.1814,0.0238604,-9.87518e-06,1.83255e-09,-1.27379e-13,-7516.35,-73.6646], Tmin=(1256.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.8509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=CC=CCJ)"""),
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
    label = 'C=CC=[C]C=CC=O(3097)',
    structure = SMILES('C=CC=C=CC=C[O]'),
    E0 = (218.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29257,'amu*angstrom^2'), symmetry=1, barrier=(29.7187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28492,'amu*angstrom^2'), symmetry=1, barrier=(29.5428,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.143668,0.0655376,-1.15628e-05,-4.94224e-08,2.86896e-11,26489.3,26.6013], Tmin=(100,'K'), Tmax=(929.405,'K')), NASAPolynomial(coeffs=[24.6684,0.00941377,-7.53716e-07,4.43332e-11,-8.03476e-15,19795.9,-101.405], Tmin=(929.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'C=C[CH]C(=O)C=CC=O(2764)',
    structure = SMILES('C=CC=C([O])C=CC=O'),
    E0 = (-53.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.885181,'amu*angstrom^2'), symmetry=1, barrier=(20.3521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.886434,'amu*angstrom^2'), symmetry=1, barrier=(20.3809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.884396,'amu*angstrom^2'), symmetry=1, barrier=(20.334,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.374117,0.0881774,-8.71591e-05,4.28208e-08,-8.24952e-12,-6316.41,30.155], Tmin=(100,'K'), Tmax=(1264.29,'K')), NASAPolynomial(coeffs=[20.4015,0.0224468,-9.17397e-06,1.69887e-09,-1.18107e-13,-11569.7,-74.9477], Tmin=(1264.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C[C]=C(O)C=CC=O(3100)',
    structure = SMILES('C=C[C]=C(O)C=CC=O'),
    E0 = (7.30767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07338,'amu*angstrom^2'), symmetry=1, barrier=(24.6791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07318,'amu*angstrom^2'), symmetry=1, barrier=(24.6745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07534,'amu*angstrom^2'), symmetry=1, barrier=(24.7243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07387,'amu*angstrom^2'), symmetry=1, barrier=(24.6903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.820717,0.0987697,-0.000109251,5.96901e-08,-1.26488e-11,1058.71,29.8509], Tmin=(100,'K'), Tmax=(1161.33,'K')), NASAPolynomial(coeffs=[22.0157,0.0201129,-7.65493e-06,1.36817e-09,-9.36797e-14,-4245.37,-83.7371], Tmin=(1161.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.30767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC=C(O)C=[C]C=O(3103)',
    structure = SMILES('C=CC=C(O)C=C=C[O]'),
    E0 = (5.04892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.22942,'amu*angstrom^2'), symmetry=1, barrier=(28.2667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23149,'amu*angstrom^2'), symmetry=1, barrier=(28.3144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23253,'amu*angstrom^2'), symmetry=1, barrier=(28.3383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52894,0.10706,-0.00011368,5.55339e-08,-9.892e-12,873.82,35.4661], Tmin=(100,'K'), Tmax=(1638.55,'K')), NASAPolynomial(coeffs=[30.0639,0.00361713,2.87535e-06,-7.81909e-10,5.7262e-14,-6601.69,-128.089], Tmin=(1638.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.04892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
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
    label = 'C=C=CC(O)=CC=CO(3156)',
    structure = SMILES('C=C=CC(O)=CC=CO'),
    E0 = (-136.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.19335,0.117059,-0.000127922,6.37793e-08,-1.15109e-11,-16112.3,36.4375], Tmin=(100,'K'), Tmax=(1633.34,'K')), NASAPolynomial(coeffs=[32.3324,0.00127156,4.84942e-06,-1.20321e-09,8.6965e-14,-23877.6,-140.629], Tmin=(1633.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CC=CC(O)=CC=C=O(3157)',
    structure = SMILES('CC=CC(O)=CC=C=O'),
    E0 = (-169.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.389279,0.0930393,-9.84963e-05,5.41118e-08,-1.17817e-11,-20253,29.7498], Tmin=(100,'K'), Tmax=(1118.94,'K')), NASAPolynomial(coeffs=[17.8647,0.0277836,-1.10157e-05,1.98964e-09,-1.36037e-13,-24338,-60.366], Tmin=(1118.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = 'CC=CC(O)=C=CC=O(3130)',
    structure = SMILES('CC=CC(O)=C=CC=O'),
    E0 = (-138.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.887532,'amu*angstrom^2'), symmetry=1, barrier=(20.4061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888904,'amu*angstrom^2'), symmetry=1, barrier=(20.4376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.889438,'amu*angstrom^2'), symmetry=1, barrier=(20.4499,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.890202,'amu*angstrom^2'), symmetry=1, barrier=(20.4675,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.335935,0.0941929,-9.70911e-05,5.09657e-08,-1.06863e-11,-16549.4,28.5481], Tmin=(100,'K'), Tmax=(1150.66,'K')), NASAPolynomial(coeffs=[18.0687,0.0302134,-1.3687e-05,2.6431e-09,-1.874e-13,-20784.9,-62.8263], Tmin=(1150.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C=CC(O)=CCC=O(3129)',
    structure = SMILES('C=C=CC(O)=CCC=O'),
    E0 = (-114.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.955948,'amu*angstrom^2'), symmetry=1, barrier=(21.9791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958828,'amu*angstrom^2'), symmetry=1, barrier=(22.0453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.955361,'amu*angstrom^2'), symmetry=1, barrier=(21.9656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.957617,'amu*angstrom^2'), symmetry=1, barrier=(22.0175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.483644,0.0848869,-6.40415e-05,1.6095e-08,1.00096e-12,-13619.8,30.5018], Tmin=(100,'K'), Tmax=(1083.23,'K')), NASAPolynomial(coeffs=[22.2019,0.0235126,-1.00661e-05,1.96267e-09,-1.42393e-13,-19848.5,-86.8215], Tmin=(1083.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = 'CC=C=C(O)C=CC=O(3181)',
    structure = SMILES('CC=C=C(O)C=CC=O'),
    E0 = (-138.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.335935,0.0941929,-9.70911e-05,5.09657e-08,-1.06863e-11,-16549.4,28.5481], Tmin=(100,'K'), Tmax=(1150.66,'K')), NASAPolynomial(coeffs=[18.0687,0.0302134,-1.3687e-05,2.6431e-09,-1.874e-13,-20784.9,-62.8263], Tmin=(1150.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C=CC(O)C=CC=O(2784)',
    structure = SMILES('C=C=CC(O)C=CC=O'),
    E0 = (-71.9928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676263,0.0785045,-6.70016e-05,2.97846e-08,-5.59314e-12,-8543.8,31.7707], Tmin=(100,'K'), Tmax=(1216.32,'K')), NASAPolynomial(coeffs=[12.2995,0.0402806,-1.98631e-05,3.948e-09,-2.82766e-13,-11371.3,-26.5807], Tmin=(1216.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.9928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CC=CC(=O)C=CC=O(2800)',
    structure = SMILES('CC=CC(=O)C=CC=O'),
    E0 = (-182.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79567,0.0591273,-3.23021e-05,7.05941e-09,-5.7161e-13,-21962.6,22.2864], Tmin=(100,'K'), Tmax=(2822.34,'K')), NASAPolynomial(coeffs=[29.3192,0.0215368,-1.2324e-05,2.3404e-09,-1.53609e-13,-36934.4,-133.194], Tmin=(2822.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-182.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H)"""),
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
    label = 'C=CC=C(O)CC=C=O(3172)',
    structure = SMILES('C=CC=C(O)CC=C=O'),
    E0 = (-160.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.979769,0.0964575,-0.000100184,5.17061e-08,-1.02955e-11,-19115.8,32.6551], Tmin=(100,'K'), Tmax=(1263.63,'K')), NASAPolynomial(coeffs=[23.2447,0.0182864,-5.6232e-06,8.84921e-10,-5.63842e-14,-25119,-89.4116], Tmin=(1263.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-160.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]CC=C(O)C=[C]C=O(3118)',
    structure = SMILES('[CH2]CC=C(O)C=C=C[O]'),
    E0 = (99.8523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.991579,'amu*angstrom^2'), symmetry=1, barrier=(22.7984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.991448,'amu*angstrom^2'), symmetry=1, barrier=(22.7954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990962,'amu*angstrom^2'), symmetry=1, barrier=(22.7842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.99112,'amu*angstrom^2'), symmetry=1, barrier=(22.7878,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93975,0.104448,-0.000107983,5.30248e-08,-9.71357e-12,12245.7,37.207], Tmin=(100,'K'), Tmax=(1529.4,'K')), NASAPolynomial(coeffs=[28.4025,0.0101463,-8.36515e-07,-6.96575e-11,9.77883e-15,4712.53,-116.354], Tmin=(1529.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.8523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(RCCJ) + radical(C=COJ)"""),
)

species(
    label = 'C=C[C]=C(O)C=CC[O](3120)',
    structure = SMILES('C=C[C]=C(O)C=CC[O]'),
    E0 = (167.339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.840089,'amu*angstrom^2'), symmetry=1, barrier=(19.3153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843916,'amu*angstrom^2'), symmetry=1, barrier=(19.4033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842975,'amu*angstrom^2'), symmetry=1, barrier=(19.3816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842685,'amu*angstrom^2'), symmetry=1, barrier=(19.375,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.560433,0.0979597,-0.000110769,6.58183e-08,-1.54552e-11,20292.6,32.2622], Tmin=(100,'K'), Tmax=(1042.82,'K')), NASAPolynomial(coeffs=[17.5672,0.0284248,-1.0747e-05,1.873e-09,-1.24935e-13,16511.9,-55.9525], Tmin=(1042.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CCOJ)"""),
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
    label = 'C=[C]C=C(O)[CH]CC=O(3117)',
    structure = SMILES('[CH2][C]=CC(O)=CCC=O'),
    E0 = (64.6122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,349.258,349.701],'cm^-1')),
        HinderedRotor(inertia=(0.243776,'amu*angstrom^2'), symmetry=1, barrier=(21.0985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243594,'amu*angstrom^2'), symmetry=1, barrier=(21.1002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24368,'amu*angstrom^2'), symmetry=1, barrier=(21.0984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243067,'amu*angstrom^2'), symmetry=1, barrier=(21.1017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24432,'amu*angstrom^2'), symmetry=1, barrier=(21.1,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.275918,0.0809595,-5.66876e-05,1.13341e-08,2.16006e-12,7936.19,32.7331], Tmin=(100,'K'), Tmax=(1074.77,'K')), NASAPolynomial(coeffs=[20.4689,0.0259204,-1.08103e-05,2.06715e-09,-1.48137e-13,2196.69,-74.8007], Tmin=(1074.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.6122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = 'C[CH]C=C(O)[C]=CC=O(3116)',
    structure = SMILES('CC=CC(O)=[C]C=C[O]'),
    E0 = (27.4401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15025,'amu*angstrom^2'), symmetry=1, barrier=(26.4466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14842,'amu*angstrom^2'), symmetry=1, barrier=(26.4045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15074,'amu*angstrom^2'), symmetry=1, barrier=(26.4578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1503,'amu*angstrom^2'), symmetry=1, barrier=(26.4476,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16331,0.108593,-0.000116349,5.92026e-08,-1.11576e-11,3545.32,34.6739], Tmin=(100,'K'), Tmax=(1509.47,'K')), NASAPolynomial(coeffs=[28.5301,0.00909565,5.72302e-07,-4.07448e-10,3.52839e-14,-3651.75,-119.189], Tmin=(1509.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.4401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = 'OC1C=CCOC=CC=1(3158)',
    structure = SMILES('OC1C=CCOC=CC=1'),
    E0 = (-181.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61164,0.0335412,0.000125233,-2.12132e-07,9.10827e-11,-21709.3,20.5265], Tmin=(100,'K'), Tmax=(921.8,'K')), NASAPolynomial(coeffs=[32.5707,0.000816322,6.06772e-06,-1.25302e-09,7.18076e-14,-32102.9,-155.473], Tmin=(921.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-181.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cyclooctane)"""),
)

species(
    label = 'O=CC1C=C(O)C=CC1(3096)',
    structure = SMILES('O=CC1C=C(O)C=CC1'),
    E0 = (-238.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.35857,0.065263,-1.73863e-05,-2.73916e-08,1.59043e-11,-28600.2,25.0542], Tmin=(100,'K'), Tmax=(990.042,'K')), NASAPolynomial(coeffs=[18.2237,0.0266899,-9.86095e-06,1.82668e-09,-1.31384e-13,-33784.7,-69.2738], Tmin=(990.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cyclohexadiene)"""),
)

species(
    label = 'C=CC1OC=CC=C1O(3167)',
    structure = SMILES('C=CC1OC=CC=C1O'),
    E0 = (-207.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00582802,0.054162,6.1312e-05,-1.42793e-07,6.52529e-11,-24829.7,25.4315], Tmin=(100,'K'), Tmax=(930.342,'K')), NASAPolynomial(coeffs=[32.1294,0.00305702,3.34136e-06,-6.67746e-10,3.20244e-14,-34576.7,-147.531], Tmin=(930.342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-207.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(1,3-Cyclohexadiene)"""),
)

species(
    label = 'S(826)(825)',
    structure = SMILES('O=CC=CC1(O)C=CC1'),
    E0 = (-102.851,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4592.46,'J/mol'), sigma=(7.24917,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=717.33 K, Pc=27.35 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.478527,0.0628281,-1.68028e-05,-2.07149e-08,1.09647e-11,-12231,30.0815], Tmin=(100,'K'), Tmax=(1091.02,'K')), NASAPolynomial(coeffs=[17.4585,0.0306381,-1.38792e-05,2.75524e-09,-2.0066e-13,-17725.4,-61.5161], Tmin=(1091.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene)"""),
)

species(
    label = 'S(825)(824)',
    structure = SMILES('C=CC1C(O)=CC1C=O'),
    E0 = (-111.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4470.27,'J/mol'), sigma=(6.92391,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=698.25 K, Pc=30.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0242753,0.0767999,-5.10858e-05,6.77427e-09,4.08457e-12,-13299,28.2186], Tmin=(100,'K'), Tmax=(1016.68,'K')), NASAPolynomial(coeffs=[18.7186,0.0259175,-9.74072e-06,1.77857e-09,-1.25191e-13,-18291.5,-68.325], Tmin=(1016.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = 'C=CC=C(O)C1C=CO1(3173)',
    structure = SMILES('C=CC=C(O)C1C=CO1'),
    E0 = (-100.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.390841,0.0621693,4.66014e-05,-1.35019e-07,6.45966e-11,-11841.8,28.7863], Tmin=(100,'K'), Tmax=(923.426,'K')), NASAPolynomial(coeffs=[35.0079,-0.0020251,6.07641e-06,-1.22227e-09,7.17702e-14,-22180.1,-159.751], Tmin=(923.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = '[O]C=C[CH]C1(O)C=CC1(3140)',
    structure = SMILES('[O]C=C[CH]C1(O)C=CC1'),
    E0 = (5.87588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.332021,0.0477376,7.261e-05,-1.48487e-07,6.58482e-11,869.163,32.5986], Tmin=(100,'K'), Tmax=(931.544,'K')), NASAPolynomial(coeffs=[29.279,0.00773888,1.27791e-06,-2.9496e-10,7.27366e-15,-8181.5,-124.633], Tmin=(931.544,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.87588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=COJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[CH2]C=C[C](O)C1C=CO1(3141)',
    structure = SMILES('[CH2][CH]C=C(O)C1C=CO1'),
    E0 = (135.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.177909,0.0609293,3.88353e-05,-1.17576e-07,5.57805e-11,16521.3,31.4774], Tmin=(100,'K'), Tmax=(932.768,'K')), NASAPolynomial(coeffs=[31.3954,0.00463737,2.15069e-06,-4.38329e-10,1.74479e-14,7189.91,-137.094], Tmin=(932.768,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Allyl_S) + radical(RCCJ)"""),
)

species(
    label = '[O][CH]C1C=C(O)C=CC1(3142)',
    structure = SMILES('[O][CH]C1C=C(O)C=CC1'),
    E0 = (88.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0645229,0.0747738,-4.57712e-05,3.24671e-09,4.75426e-12,10844.5,27.7526], Tmin=(100,'K'), Tmax=(1028.63,'K')), NASAPolynomial(coeffs=[17.8776,0.0281477,-1.07978e-05,1.98036e-09,-1.39149e-13,5982,-64.5106], Tmin=(1028.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]C1OC=CC=C1O(3143)',
    structure = SMILES('[CH2][CH]C1OC=CC=C1O'),
    E0 = (73.2177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0810752,0.0551745,6.09057e-05,-1.45561e-07,6.72665e-11,8985.05,29.6735], Tmin=(100,'K'), Tmax=(925.466,'K')), NASAPolynomial(coeffs=[33.2434,0.00042406,4.93533e-06,-9.98979e-10,5.58409e-14,-1006.57,-149.174], Tmin=(925.466,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.2177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(1,3-Cyclohexadiene) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = 'C[C]=CC(O)=CC=C[O](3144)',
    structure = SMILES('C[C]=CC(O)=CC=C[O]'),
    E0 = (66.2864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01137,'amu*angstrom^2'), symmetry=1, barrier=(23.2534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0113,'amu*angstrom^2'), symmetry=1, barrier=(23.2518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01229,'amu*angstrom^2'), symmetry=1, barrier=(23.2745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01103,'amu*angstrom^2'), symmetry=1, barrier=(23.2455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16313,0.107068,-0.000113018,5.62026e-08,-1.03465e-11,8218.75,35.5579], Tmin=(100,'K'), Tmax=(1543.64,'K')), NASAPolynomial(coeffs=[29.1232,0.00809779,5.46915e-07,-3.5511e-10,2.97847e-14,692.222,-122.056], Tmin=(1543.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.2864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=COJ)"""),
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
    label = 'CC=[C]C(O)=CC=C[O](3146)',
    structure = SMILES('CC=[C]C(O)=CC=C[O]'),
    E0 = (27.4401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15319,'amu*angstrom^2'), symmetry=1, barrier=(26.514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14226,'amu*angstrom^2'), symmetry=1, barrier=(26.2629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14862,'amu*angstrom^2'), symmetry=1, barrier=(26.409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15239,'amu*angstrom^2'), symmetry=1, barrier=(26.4957,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16331,0.108593,-0.000116349,5.92026e-08,-1.11576e-11,3545.32,34.6739], Tmin=(100,'K'), Tmax=(1509.47,'K')), NASAPolynomial(coeffs=[28.5301,0.00909565,5.72302e-07,-4.07448e-10,3.52839e-14,-3651.75,-119.189], Tmin=(1509.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.4401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJC=C)"""),
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
    label = 'CC=C[C]([O])C=CC=O(2794)',
    structure = SMILES('CC=CC([O])=CC=C[O]'),
    E0 = (-33.7506,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.962461,'amu*angstrom^2'), symmetry=1, barrier=(22.1289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.964606,'amu*angstrom^2'), symmetry=1, barrier=(22.1782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.962663,'amu*angstrom^2'), symmetry=1, barrier=(22.1335,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.601334,0.0850063,-4.95027e-05,-1.48209e-08,1.72783e-11,-3878.61,30.9656], Tmin=(100,'K'), Tmax=(925.701,'K')), NASAPolynomial(coeffs=[25.488,0.0141823,-2.64931e-06,3.43118e-10,-2.49728e-14,-10504.4,-102.585], Tmin=(925.701,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.7506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[CH]C=C(O)C=[C]C=O(3121)',
    structure = SMILES('CC=CC(O)=C[C]=C[O]'),
    E0 = (27.4401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15319,'amu*angstrom^2'), symmetry=1, barrier=(26.514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14226,'amu*angstrom^2'), symmetry=1, barrier=(26.2629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14862,'amu*angstrom^2'), symmetry=1, barrier=(26.409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15239,'amu*angstrom^2'), symmetry=1, barrier=(26.4957,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16331,0.108593,-0.000116349,5.92026e-08,-1.11576e-11,3545.32,34.6739], Tmin=(100,'K'), Tmax=(1509.47,'K')), NASAPolynomial(coeffs=[28.5301,0.00909565,5.72302e-07,-4.07448e-10,3.52839e-14,-3651.75,-119.189], Tmin=(1509.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.4401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJC=C)"""),
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
    label = 'C[CH]C=C(O)C=C[C]=O(3127)',
    structure = SMILES('C[CH]C=C(O)[CH]C=C=O'),
    E0 = (-12.9343,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3000,3050,390,425,1340,1360,335,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.792502,0.0957108,-9.55689e-05,4.78715e-08,-9.38681e-12,-1374.9,30.1587], Tmin=(100,'K'), Tmax=(1246.58,'K')), NASAPolynomial(coeffs=[21.6103,0.0238252,-9.06973e-06,1.61216e-09,-1.09576e-13,-6960.3,-82.8598], Tmin=(1246.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.9343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-Cd(CCO)HH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(Allyl_S) + radical(C=CCJCO)"""),
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
    label = '[O]C=CC=C(O)C1[CH]C1(3150)',
    structure = SMILES('[O]C=CC=C(O)C1[CH]C1'),
    E0 = (110.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.139192,0.0639498,2.21101e-05,-9.69616e-08,4.82276e-11,13418.3,32.7976], Tmin=(100,'K'), Tmax=(924.97,'K')), NASAPolynomial(coeffs=[29.33,0.00651973,1.71188e-06,-4.32795e-10,2.17465e-14,4971.86,-123.265], Tmin=(924.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(cyclopropane) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=CC(O)=CC1[CH]O1(3151)',
    structure = SMILES('[CH2]C=CC(O)=CC1[CH]O1'),
    E0 = (103.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.57971,0.101714,-0.000100701,4.84474e-08,-8.469e-12,12766.5,39.4842], Tmin=(100,'K'), Tmax=(1726.49,'K')), NASAPolynomial(coeffs=[23.2273,0.0129689,1.55739e-06,-7.52564e-10,6.18542e-14,8170.68,-86.6153], Tmin=(1726.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(Ethylene_oxide) + radical(CCsJO) + radical(C=CC=CCJ)"""),
)

species(
    label = '[O]C=CC1CC=C[C]1O(3152)',
    structure = SMILES('[O]C=CC1C[CH]C=C1O'),
    E0 = (-60.8842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.196403,0.0574943,3.29355e-05,-9.96625e-08,4.66206e-11,-7161.71,26.176], Tmin=(100,'K'), Tmax=(938.6,'K')), NASAPolynomial(coeffs=[26.1407,0.0131463,-2.0157e-06,3.27449e-10,-3.2648e-14,-14948.8,-112.884], Tmin=(938.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.8842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopentene) + radical(C=COJ) + radical(cyclopentene-allyl)"""),
)

species(
    label = '[CH2]C=CC1(O)[CH]C=CO1(3153)',
    structure = SMILES('[CH2]C=CC1(O)[CH]C=CO1'),
    E0 = (-88.6297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0699264,0.0397398,0.000131053,-2.38452e-07,1.06422e-10,-10474.4,26.4189], Tmin=(100,'K'), Tmax=(906.518,'K')), NASAPolynomial(coeffs=[39.8633,-0.0128168,1.44391e-05,-2.97868e-09,1.94494e-13,-22744.3,-189.538], Tmin=(906.518,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-88.6297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(C=CCJC(O)C=C) + radical(Allyl_P)"""),
)

species(
    label = '[O]C1[CH]C=C(O)C=CC1(3154)',
    structure = SMILES('[O]C1[CH]C=C(O)C=CC1'),
    E0 = (33.7646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03493,0.0947871,-6.86273e-05,5.30964e-09,8.98096e-12,4256.62,12.8668], Tmin=(100,'K'), Tmax=(964.894,'K')), NASAPolynomial(coeffs=[26.3845,0.0176699,-5.56429e-06,9.97212e-10,-7.36039e-14,-2736.25,-127.253], Tmin=(964.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.7646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C1[CH]C(O)=CC=CO1(3155)',
    structure = SMILES('[CH2]C1[CH]C(O)=CC=CO1'),
    E0 = (11.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.786733,0.0327885,0.00021774,-3.82823e-07,1.72631e-10,1655.96,26.2111], Tmin=(100,'K'), Tmax=(895.539,'K')), NASAPolynomial(coeffs=[61.829,-0.0533491,3.78434e-05,-7.57657e-09,5.08393e-13,-17319.9,-312.296], Tmin=(895.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(CJC(C)OC) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C=C[C]=CC=COO(3159)',
    structure = SMILES('[CH2]C=C[C]=CC=COO'),
    E0 = (342.421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1685,370,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.35133,'amu*angstrom^2'), symmetry=1, barrier=(31.0697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35536,'amu*angstrom^2'), symmetry=1, barrier=(31.1623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34717,'amu*angstrom^2'), symmetry=1, barrier=(30.974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35758,'amu*angstrom^2'), symmetry=1, barrier=(31.2134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34581,'amu*angstrom^2'), symmetry=1, barrier=(30.9429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.470958,0.08672,-6.40959e-05,9.50437e-09,6.26842e-12,41355.3,33.5404], Tmin=(100,'K'), Tmax=(939.689,'K')), NASAPolynomial(coeffs=[21.3272,0.0223926,-6.84392e-06,1.1186e-09,-7.57522e-14,36002,-76.9531], Tmin=(939.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C[CH]C(=O)C[CH]C=O(2813)',
    structure = SMILES('C=CC=C([O])CC=C[O]'),
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
    label = '[CH]=CC(O)=CC=C[O](3160)',
    structure = SMILES('[CH]C=C(O)C=CC=O'),
    E0 = (91.1991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.93322,'amu*angstrom^2'), symmetry=1, barrier=(44.4485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93307,'amu*angstrom^2'), symmetry=1, barrier=(44.4451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93296,'amu*angstrom^2'), symmetry=1, barrier=(44.4425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93303,'amu*angstrom^2'), symmetry=1, barrier=(44.4442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0145755,0.0819281,-7.76483e-05,3.70647e-08,-7.03797e-12,11117,27.5375], Tmin=(100,'K'), Tmax=(1269.35,'K')), NASAPolynomial(coeffs=[17.9157,0.0255181,-1.09888e-05,2.05525e-09,-1.4287e-13,6572.35,-63.0949], Tmin=(1269.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.1991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
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
    label = 'C=CC1C(O)=CC1[CH][O](3161)',
    structure = SMILES('C=CC1C(O)=CC1[CH][O]'),
    E0 = (216.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.352446,0.0866627,-8.0427e-05,3.81947e-08,-7.18095e-12,26147.2,31.0424], Tmin=(100,'K'), Tmax=(1289.67,'K')), NASAPolynomial(coeffs=[19.1825,0.0260739,-9.9572e-06,1.76713e-09,-1.19568e-13,21108.4,-68.1721], Tmin=(1289.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C=CC(O)=[C]CC=O(3174)',
    structure = SMILES('[CH2]C=CC(O)=[C]CC=O'),
    E0 = (64.6122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,349.258,349.701],'cm^-1')),
        HinderedRotor(inertia=(0.243776,'amu*angstrom^2'), symmetry=1, barrier=(21.0985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243594,'amu*angstrom^2'), symmetry=1, barrier=(21.1002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24368,'amu*angstrom^2'), symmetry=1, barrier=(21.0984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243067,'amu*angstrom^2'), symmetry=1, barrier=(21.1017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24432,'amu*angstrom^2'), symmetry=1, barrier=(21.1,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.275918,0.0809595,-5.66876e-05,1.13341e-08,2.16006e-12,7936.19,32.7331], Tmin=(100,'K'), Tmax=(1074.77,'K')), NASAPolynomial(coeffs=[20.4689,0.0259204,-1.08103e-05,2.06715e-09,-1.48137e-13,2196.69,-74.8007], Tmin=(1074.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.6122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(Cds_S) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=C(O)[CH]C[C]=O(3109)',
    structure = SMILES('[CH2]C=CC(O)=CC[C]=O'),
    E0 = (-13.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,1855,455,950,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,303.664,303.678],'cm^-1')),
        HinderedRotor(inertia=(0.311206,'amu*angstrom^2'), symmetry=1, barrier=(20.3635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311185,'amu*angstrom^2'), symmetry=1, barrier=(20.364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311215,'amu*angstrom^2'), symmetry=1, barrier=(20.3637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311176,'amu*angstrom^2'), symmetry=1, barrier=(20.3643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311205,'amu*angstrom^2'), symmetry=1, barrier=(20.3639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.262064,0.0787523,-4.52582e-05,-5.34096e-09,9.36399e-12,-1429.23,33.0247], Tmin=(100,'K'), Tmax=(1006.65,'K')), NASAPolynomial(coeffs=[21.8103,0.0226488,-8.75022e-06,1.66808e-09,-1.2195e-13,-7474.29,-81.5614], Tmin=(1006.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C[CH]C(=O)[CH]CC=O(2816)',
    structure = SMILES('[CH2]C=CC([O])=CCC=O'),
    E0 = (-35.4248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,350,440,435,1725,351.009,351.133,351.257],'cm^-1')),
        HinderedRotor(inertia=(0.00136737,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276587,'amu*angstrom^2'), symmetry=1, barrier=(24.1895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276301,'amu*angstrom^2'), symmetry=1, barrier=(24.1884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276243,'amu*angstrom^2'), symmetry=1, barrier=(24.1854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.120753,0.0723982,-3.93744e-05,-1.02366e-09,5.22804e-12,-4109.96,32.3379], Tmin=(100,'K'), Tmax=(1082.72,'K')), NASAPolynomial(coeffs=[18.0518,0.0295664,-1.2471e-05,2.3829e-09,-1.70053e-13,-9365.11,-61.9312], Tmin=(1082.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.4248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C[C]=C(O)[CH]CC=O(3112)',
    structure = SMILES('[CH2]C=[C]C(O)=CCC=O'),
    E0 = (25.7659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,276.924,276.93],'cm^-1')),
        HinderedRotor(inertia=(0.471207,'amu*angstrom^2'), symmetry=1, barrier=(25.6418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470839,'amu*angstrom^2'), symmetry=1, barrier=(25.6415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.471175,'amu*angstrom^2'), symmetry=1, barrier=(25.6421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.471015,'amu*angstrom^2'), symmetry=1, barrier=(25.6416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.471161,'amu*angstrom^2'), symmetry=1, barrier=(25.6414,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.329784,0.0830371,-6.16076e-05,1.59644e-08,8.1615e-13,3265.29,32.0476], Tmin=(100,'K'), Tmax=(1060.88,'K')), NASAPolynomial(coeffs=[19.9963,0.0267002,-1.06564e-05,1.98428e-09,-1.40108e-13,-2189.86,-72.5996], Tmin=(1060.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.7659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=CC1(O)[CH]C1C=O(3175)',
    structure = SMILES('C=C[CH]C1(O)[CH]C1C=O'),
    E0 = (112.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.157168,0.0769094,-4.09632e-05,-9.54989e-09,1.11042e-11,13697.9,32.3081], Tmin=(100,'K'), Tmax=(985.463,'K')), NASAPolynomial(coeffs=[20.9536,0.0233766,-8.425e-06,1.54965e-09,-1.11673e-13,7975.7,-77.1524], Tmin=(985.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C=CCJCO) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1[CH]C(O)=CC1C=O(3176)',
    structure = SMILES('[CH2]C1[CH]C(O)=CC1C=O'),
    E0 = (66.1262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.455628,0.0601365,2.70979e-06,-5.50226e-08,2.79854e-11,8097.04,31.771], Tmin=(100,'K'), Tmax=(943.673,'K')), NASAPolynomial(coeffs=[20.1106,0.021233,-6.04127e-06,1.02831e-09,-7.50588e-14,2410.11,-72.3903], Tmin=(943.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.1262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopentene) + radical(Isobutyl) + radical(CCJCO)"""),
)

species(
    label = 'OC1[CH][CH]COC=CC=1(3166)',
    structure = SMILES('OC1[CH]C=COC[CH]C=1'),
    E0 = (-34.9453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.278101,0.0407488,0.000113923,-2.00536e-07,8.59629e-11,-4031.02,20.9888], Tmin=(100,'K'), Tmax=(933.014,'K')), NASAPolynomial(coeffs=[33.2674,0.00447628,3.17553e-06,-6.02902e-10,2.2911e-14,-14764,-160.407], Tmin=(933.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.9453,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclooctane) + radical(C=CCJCO) + radical(C=CCJCO)"""),
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
    label = 'C=CC1[C](O)C1C=C[O](3164)',
    structure = SMILES('C=CC1[C](O)C1C=C[O]'),
    E0 = (123.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0183951,0.0668857,-2.60654e-06,-5.8628e-08,3.12363e-11,15003.9,32.9981], Tmin=(100,'K'), Tmax=(941.127,'K')), NASAPolynomial(coeffs=[24.2799,0.0156255,-3.557e-06,5.9259e-10,-4.73937e-14,8140.72,-94.7784], Tmin=(941.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C=COJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=CC1C(O)=C[CH]C1[O](3165)',
    structure = SMILES('C=CC1C(O)=C[CH]C1[O]'),
    E0 = (53.9474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.174971,0.0716412,-1.18798e-05,-4.77996e-08,2.64905e-11,6656.97,26.0181], Tmin=(100,'K'), Tmax=(959.274,'K')), NASAPolynomial(coeffs=[24.1855,0.0189139,-5.81918e-06,1.07599e-09,-8.2533e-14,-264.37,-102.21], Tmin=(959.274,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.9474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopentene) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C[CH]C1(O)C=CC1[O](3178)',
    structure = SMILES('C=C[CH]C1(O)C=CC1[O]'),
    E0 = (130.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.427201,0.0494881,5.68221e-05,-1.2377e-07,5.47082e-11,15882.1,33.3691], Tmin=(100,'K'), Tmax=(943.572,'K')), NASAPolynomial(coeffs=[26.2928,0.0132464,-2.26106e-06,4.24492e-10,-4.26772e-14,7733.01,-107.231], Tmin=(943.572,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ)"""),
)

species(
    label = 'S(458)(457)',
    structure = SMILES('C=C[CH]C([O])C=CC=O'),
    E0 = (67.9095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,466.122,466.138,466.658,466.797],'cm^-1')),
        HinderedRotor(inertia=(0.150317,'amu*angstrom^2'), symmetry=1, barrier=(23.1705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150455,'amu*angstrom^2'), symmetry=1, barrier=(23.1786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150011,'amu*angstrom^2'), symmetry=1, barrier=(23.1625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150418,'amu*angstrom^2'), symmetry=1, barrier=(23.174,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4360.96,'J/mol'), sigma=(6.98043,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=681.17 K, Pc=29.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.165628,0.0680478,-2.07106e-05,-2.23407e-08,1.26385e-11,8319.62,35.7619], Tmin=(100,'K'), Tmax=(1066.29,'K')), NASAPolynomial(coeffs=[19.6886,0.0286735,-1.29576e-05,2.59562e-09,-1.91068e-13,2231.15,-68.705], Tmin=(1066.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.9095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C=[C]C(O)C=CC=O(2792)',
    structure = SMILES('[CH2]C=[C]C(O)C=CC=O'),
    E0 = (140.744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,270.321,270.322],'cm^-1')),
        HinderedRotor(inertia=(0.295258,'amu*angstrom^2'), symmetry=1, barrier=(15.3105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00230694,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295249,'amu*angstrom^2'), symmetry=1, barrier=(15.3105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295258,'amu*angstrom^2'), symmetry=1, barrier=(15.3105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295264,'amu*angstrom^2'), symmetry=1, barrier=(15.3105,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.619542,0.0806769,-7.27918e-05,3.5071e-08,-7.18762e-12,17044,33.683], Tmin=(100,'K'), Tmax=(1122.68,'K')), NASAPolynomial(coeffs=[11.5667,0.0416728,-2.06781e-05,4.12463e-09,-2.96362e-13,14586,-20.3975], Tmin=(1122.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[CH]C(O)[C]=CC=O(2771)',
    structure = SMILES('C=C[CH]C(O)[C]=CC=O'),
    E0 = (75.3905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0982422,0.077298,-5.29875e-05,1.29475e-08,1.22714e-13,9225.58,37.3839], Tmin=(100,'K'), Tmax=(1200.62,'K')), NASAPolynomial(coeffs=[19.9041,0.0278345,-1.265e-05,2.46531e-09,-1.75854e-13,3184.53,-67.9288], Tmin=(1200.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.3905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C=[C][CH]C(O)C=CC=O(2774)',
    structure = SMILES('C=[C][CH]C(O)C=CC=O'),
    E0 = (75.3905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0982422,0.077298,-5.29875e-05,1.29475e-08,1.22714e-13,9225.58,37.3839], Tmin=(100,'K'), Tmax=(1200.62,'K')), NASAPolynomial(coeffs=[19.9041,0.0278345,-1.265e-05,2.46531e-09,-1.75854e-13,3184.53,-67.9288], Tmin=(1200.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.3905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C=C[CH]C(O)C=[C]C=O(2775)',
    structure = SMILES('C=C[CH]C(O)C=C=C[O]'),
    E0 = (34.2854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,261.971,261.971,261.971,261.972],'cm^-1')),
        HinderedRotor(inertia=(0.57055,'amu*angstrom^2'), symmetry=1, barrier=(27.7864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570549,'amu*angstrom^2'), symmetry=1, barrier=(27.7864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570554,'amu*angstrom^2'), symmetry=1, barrier=(27.7864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570554,'amu*angstrom^2'), symmetry=1, barrier=(27.7864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.161623,0.0683024,2.44036e-06,-6.66851e-08,3.41219e-11,4294.55,36.1724], Tmin=(100,'K'), Tmax=(953.408,'K')), NASAPolynomial(coeffs=[26.1233,0.0152201,-4.0307e-06,7.61892e-10,-6.31454e-14,-3316.98,-103.016], Tmin=(953.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.2854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C=C[CH]C(O)C=C[C]=O(2778)',
    structure = SMILES('[CH2]C=CC(O)[CH]C=C=O'),
    E0 = (-8.65128,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0558359,0.0701448,-2.20694e-05,-2.71934e-08,1.65279e-11,-883.802,36.0864], Tmin=(100,'K'), Tmax=(997.068,'K')), NASAPolynomial(coeffs=[20.7109,0.0245117,-9.42786e-06,1.80357e-09,-1.32574e-13,-6853.31,-72.7819], Tmin=(997.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.65128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(Allyl_P) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[CH2]C1[CH]C1(O)C=CC=O(3177)',
    structure = SMILES('[CH2]C1[CH]C1(O)C=CC=O'),
    E0 = (173.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.203145,0.0689871,-2.71203e-05,-1.99694e-08,1.40587e-11,20971.5,35.6811], Tmin=(100,'K'), Tmax=(982.358,'K')), NASAPolynomial(coeffs=[19.3778,0.0238648,-8.54025e-06,1.56971e-09,-1.13191e-13,15614.1,-64.5776], Tmin=(982.358,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclopropane) + radical(CCJCO) + radical(Isobutyl)"""),
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
    label = 'C=C[CH][C](O)C=C(3179)',
    structure = SMILES('[CH2]C=CC(O)=C[CH2]'),
    E0 = (44.0886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,350,440,435,1725,307.053],'cm^-1')),
        HinderedRotor(inertia=(0.427031,'amu*angstrom^2'), symmetry=1, barrier=(28.6859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.974192,'amu*angstrom^2'), symmetry=1, barrier=(65.1928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973389,'amu*angstrom^2'), symmetry=1, barrier=(65.1811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.427842,'amu*angstrom^2'), symmetry=1, barrier=(28.6753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.861131,0.0535319,2.84349e-06,-5.15065e-08,2.68292e-11,5430,25.233], Tmin=(100,'K'), Tmax=(922.898,'K')), NASAPolynomial(coeffs=[18.3828,0.0176982,-4.10404e-06,6.01864e-10,-4.21986e-14,487.76,-67.1471], Tmin=(922.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.0886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC=C(O)[CH]CC=O(3122)',
    structure = SMILES('[CH]C=CC(O)=CCC=O'),
    E0 = (79.3991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.526522,0.0852002,-5.46693e-05,7.07036e-09,3.63956e-12,9724.72,32.7113], Tmin=(100,'K'), Tmax=(1082.02,'K')), NASAPolynomial(coeffs=[20.9861,0.0303893,-1.29498e-05,2.47697e-09,-1.76869e-13,3622.4,-79.4567], Tmin=(1082.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.3991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1C=C(O)C1C=C[O](3168)',
    structure = SMILES('[CH2]C1C=C(O)C1C=C[O]'),
    E0 = (107.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.211685,0.0701019,-9.80567e-07,-6.97507e-08,3.81825e-11,13115.2,30.4721], Tmin=(100,'K'), Tmax=(916.732,'K')), NASAPolynomial(coeffs=[27.112,0.0100912,3.27089e-07,-2.45083e-10,1.36846e-14,5617.46,-112.544], Tmin=(916.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Isobutyl) + radical(C=COJ)"""),
)

species(
    label = 'C=CC=C(O)C[C]=C[O](3169)',
    structure = SMILES('C=CC=C(O)C[C]=C[O]'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.705948,0.0845129,-4.18917e-05,-2.74746e-08,2.27989e-11,11790.3,34.1508], Tmin=(100,'K'), Tmax=(926.87,'K')), NASAPolynomial(coeffs=[27.886,0.0101276,-8.19636e-07,2.82049e-11,-5.58605e-15,4385.01,-112.973], Tmin=(926.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C=C[C]=C(O)C[CH]C=O(3110)',
    structure = SMILES('C=C[C]=C(O)CC=C[O]'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25555,0.103776,-0.000105018,5.03148e-08,-8.88412e-12,7185.47,38.8643], Tmin=(100,'K'), Tmax=(1639.51,'K')), NASAPolynomial(coeffs=[27.8768,0.00934983,5.04373e-07,-3.72783e-10,3.11987e-14,115.418,-112.833], Tmin=(1639.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.6286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = 'C=CC=C(O)CC=[C][O](3170)',
    structure = SMILES('C=CC=C(O)C[CH][C]=O'),
    E0 = (41.602,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24979,0.0956637,-9.65664e-05,4.77818e-08,-8.95747e-12,5209.48,36.4779], Tmin=(100,'K'), Tmax=(1457.3,'K')), NASAPolynomial(coeffs=[24.6092,0.0147602,-3.07582e-06,3.39188e-10,-1.69074e-14,-1273.41,-94.3988], Tmin=(1457.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]C=C(O)C[CH]C=O(3115)',
    structure = SMILES('C=[C]C=C(O)CC=C[O]'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25555,0.103776,-0.000105018,5.03148e-08,-8.88412e-12,7185.47,38.8643], Tmin=(100,'K'), Tmax=(1639.51,'K')), NASAPolynomial(coeffs=[27.8768,0.00934983,5.04373e-07,-3.72783e-10,3.11987e-14,115.418,-112.833], Tmin=(1639.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.6286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CC=C(O)C[CH]C=O(3119)',
    structure = SMILES('[CH]=CC=C(O)CC=C[O]'),
    E0 = (105.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.95682,'amu*angstrom^2'), symmetry=1, barrier=(21.9992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.957002,'amu*angstrom^2'), symmetry=1, barrier=(22.0034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.957185,'amu*angstrom^2'), symmetry=1, barrier=(22.0076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.957431,'amu*angstrom^2'), symmetry=1, barrier=(22.0132,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.749182,0.0834582,-3.37977e-05,-3.93956e-08,2.79091e-11,12906.8,34.2358], Tmin=(100,'K'), Tmax=(923.91,'K')), NASAPolynomial(coeffs=[29.3427,0.00774997,5.17413e-07,-2.25489e-10,1.10957e-14,5017.15,-121.163], Tmin=(923.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = 'C=CC=C(O)C1[CH]C1[O](3171)',
    structure = SMILES('C=CC=C(O)C1[CH]C1[O]'),
    E0 = (209.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0646505,0.0631227,2.08569e-05,-9.32472e-08,4.62741e-11,25347,32.8826], Tmin=(100,'K'), Tmax=(927.397,'K')), NASAPolynomial(coeffs=[28.4422,0.0078077,9.23174e-07,-2.73316e-10,1.06576e-14,17150.9,-118.18], Tmin=(927.397,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1C=C(O)C=CC1[O](3180)',
    structure = SMILES('[CH2]C1C=C(O)C=CC1[O]'),
    E0 = (105.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.234005,0.0607047,1.55353e-05,-7.73423e-08,3.80476e-11,12828.3,28.1935], Tmin=(100,'K'), Tmax=(934.01,'K')), NASAPolynomial(coeffs=[23.7736,0.0161304,-3.19414e-06,4.90025e-10,-3.97763e-14,5978.08,-96.8959], Tmin=(934.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(CC(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C[C]=C(O)C=CC=O(3107)',
    structure = SMILES('[CH2]C[C]=C(O)C=CC=O'),
    E0 = (140.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.655137,'amu*angstrom^2'), symmetry=1, barrier=(15.0629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.655316,'amu*angstrom^2'), symmetry=1, barrier=(15.067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.655391,'amu*angstrom^2'), symmetry=1, barrier=(15.0687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.655147,'amu*angstrom^2'), symmetry=1, barrier=(15.0631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.655367,'amu*angstrom^2'), symmetry=1, barrier=(15.0682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.313288,0.0954399,-0.00010246,5.64981e-08,-1.24703e-11,17108.2,32.7825], Tmin=(100,'K'), Tmax=(1094.6,'K')), NASAPolynomial(coeffs=[17.3076,0.0310478,-1.42194e-05,2.75475e-09,-1.95571e-13,13250.6,-53.8205], Tmin=(1094.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C[CH]C(=O)C=CC=O(2810)',
    structure = SMILES('[CH2]CC=C([O])C=CC=O'),
    E0 = (40.9204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,350,440,435,1725,299.002,299.002,299.003],'cm^-1')),
        HinderedRotor(inertia=(0.193211,'amu*angstrom^2'), symmetry=1, barrier=(12.2578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193212,'amu*angstrom^2'), symmetry=1, barrier=(12.2578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193213,'amu*angstrom^2'), symmetry=1, barrier=(12.2578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193215,'amu*angstrom^2'), symmetry=1, barrier=(12.2578,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.117418,0.0865188,-8.40596e-05,4.29285e-08,-8.95424e-12,5060.45,32.2622], Tmin=(100,'K'), Tmax=(1140.68,'K')), NASAPolynomial(coeffs=[14.9103,0.0346451,-1.58457e-05,3.06122e-09,-2.16646e-13,1685.66,-41.052], Tmin=(1140.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.9204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC=C(O)[C]=CC=O(3114)',
    structure = SMILES('[CH2]CC=C(O)[C]=CC=O'),
    E0 = (102.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.778377,'amu*angstrom^2'), symmetry=1, barrier=(17.8964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778076,'amu*angstrom^2'), symmetry=1, barrier=(17.8895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778003,'amu*angstrom^2'), symmetry=1, barrier=(17.8878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777801,'amu*angstrom^2'), symmetry=1, barrier=(17.8832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778755,'amu*angstrom^2'), symmetry=1, barrier=(17.9051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.424362,0.0981498,-0.000109406,6.35202e-08,-1.47444e-11,12439.8,32.3052], Tmin=(100,'K'), Tmax=(1045.25,'K')), NASAPolynomial(coeffs=[16.954,0.031646,-1.39698e-05,2.65083e-09,-1.85896e-13,8806.86,-52.3048], Tmin=(1045.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=C) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC=C(O)C=C[C]=O(3125)',
    structure = SMILES('[CH2]CC=C(O)[CH]C=C=O'),
    E0 = (51.1994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.737969,0.0991636,-0.000106407,5.78228e-08,-1.2346e-11,6332.81,32.8545], Tmin=(100,'K'), Tmax=(1144.74,'K')), NASAPolynomial(coeffs=[20.239,0.0258641,-1.03585e-05,1.88613e-09,-1.29817e-13,1530.22,-71.1828], Tmin=(1144.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.1994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-Cd(CCO)HH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(C=CCJCO) + radical(RCCJ)"""),
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
    label = '[CH]=C[CH]C(O)C=CC=O(2777)',
    structure = SMILES('[CH]=C[CH]C(O)C=CC=O'),
    E0 = (84.6449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,446.521,446.744],'cm^-1')),
        HinderedRotor(inertia=(0.142884,'amu*angstrom^2'), symmetry=1, barrier=(20.2575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143377,'amu*angstrom^2'), symmetry=1, barrier=(20.2638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1428,'amu*angstrom^2'), symmetry=1, barrier=(20.2585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142704,'amu*angstrom^2'), symmetry=1, barrier=(20.2608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14316,'amu*angstrom^2'), symmetry=1, barrier=(20.2633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0966499,0.075838,-4.40812e-05,8.85551e-10,4.8828e-12,10340,37.3001], Tmin=(100,'K'), Tmax=(1103.1,'K')), NASAPolynomial(coeffs=[20.3178,0.0271089,-1.22173e-05,2.41708e-09,-1.75715e-13,4297.08,-70.1669], Tmin=(1103.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.6449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJC(O)C=C) + radical(Cds_P)"""),
)

species(
    label = 'C=CC=C(O)C=[C]C[O](3108)',
    structure = SMILES('C=CC=C(O)C=[C]C[O]'),
    E0 = (206.186,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.681753,'amu*angstrom^2'), symmetry=1, barrier=(15.6748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681967,'amu*angstrom^2'), symmetry=1, barrier=(15.6798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.686184,'amu*angstrom^2'), symmetry=1, barrier=(15.7767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.684375,'amu*angstrom^2'), symmetry=1, barrier=(15.7351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.458974,0.0953433,-0.00010407,5.90263e-08,-1.32509e-11,24961.4,32.7755], Tmin=(100,'K'), Tmax=(1087.19,'K')), NASAPolynomial(coeffs=[17.8638,0.0279299,-1.10592e-05,1.9922e-09,-1.35911e-13,20977.3,-57.1532], Tmin=(1087.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=CC=C(O)[C]=CC[O](3111)',
    structure = SMILES('[CH2]C=CC(O)=C=CC[O]'),
    E0 = (139.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,265.965,266.967,267.002],'cm^-1')),
        HinderedRotor(inertia=(0.272096,'amu*angstrom^2'), symmetry=1, barrier=(13.7052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271475,'amu*angstrom^2'), symmetry=1, barrier=(13.7011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.35669,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.37057,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0762511,0.0878158,-8.71114e-05,4.64709e-08,-1.00332e-11,16888,32.5168], Tmin=(100,'K'), Tmax=(1115.59,'K')), NASAPolynomial(coeffs=[15.2049,0.0330246,-1.34405e-05,2.44606e-09,-1.67456e-13,13478.4,-42.8777], Tmin=(1115.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C[CH]C(=O)C=CC[O](2827)',
    structure = SMILES('C=CC=C([O])C=CC[O]'),
    E0 = (106.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,265.943,266.083,266.125,266.431,266.652],'cm^-1')),
        HinderedRotor(inertia=(0.241144,'amu*angstrom^2'), symmetry=1, barrier=(12.0906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239113,'amu*angstrom^2'), symmetry=1, barrier=(12.0875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240596,'amu*angstrom^2'), symmetry=1, barrier=(12.0923,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0363341,0.086497,-8.58514e-05,4.56095e-08,-9.7745e-12,12914.1,32.2856], Tmin=(100,'K'), Tmax=(1125.54,'K')), NASAPolynomial(coeffs=[15.3964,0.0316516,-1.27594e-05,2.31657e-09,-1.58496e-13,9440.02,-43.9934], Tmin=(1125.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCOJ)"""),
)

species(
    label = 'C=[C]C=C(O)C=CC[O](3124)',
    structure = SMILES('C=C=CC(O)=C[CH]C[O]'),
    E0 = (144.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15343,'amu*angstrom^2'), symmetry=1, barrier=(26.5197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15365,'amu*angstrom^2'), symmetry=1, barrier=(26.5247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1537,'amu*angstrom^2'), symmetry=1, barrier=(26.5259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15362,'amu*angstrom^2'), symmetry=1, barrier=(26.5241,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.669546,0.101146,-0.000110991,6.27169e-08,-1.40624e-11,17576.6,30.3445], Tmin=(100,'K'), Tmax=(1085.54,'K')), NASAPolynomial(coeffs=[18.628,0.0300379,-1.27322e-05,2.37253e-09,-1.64928e-13,13387,-64.3387], Tmin=(1085.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(CCOJ)"""),
)

species(
    label = '[CH]=CC=C(O)C=CC[O](3128)',
    structure = SMILES('[CH]=CC=C(O)C=CC[O]'),
    E0 = (215.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.751648,'amu*angstrom^2'), symmetry=1, barrier=(17.2819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.7521,'amu*angstrom^2'), symmetry=1, barrier=(17.2923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.753503,'amu*angstrom^2'), symmetry=1, barrier=(17.3245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.752948,'amu*angstrom^2'), symmetry=1, barrier=(17.3118,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.656927,0.0961396,-0.000102556,5.57663e-08,-1.18849e-11,26084.5,33.4135], Tmin=(100,'K'), Tmax=(1149.81,'K')), NASAPolynomial(coeffs=[19.9278,0.0245289,-9.13594e-06,1.60073e-09,-1.07844e-13,21350.8,-68.7696], Tmin=(1149.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
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
    E0 = (228.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (192.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (207.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (247.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (188.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (226.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (233.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (434.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (455.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-28.5264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-28.5264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-35.7176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-14.2746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-14.2746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (9.90052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (9.90052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-35.7176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (39.5944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (24.7474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (35.4689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (206.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (331.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (177.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (206.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (191.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-45.4106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-45.9684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-45.9684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-45.2153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-45.9684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-45.9684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (68.9134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (142.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (88.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (73.2177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (239.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (207.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (220.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (197.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (107.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (66.6654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (170.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (38.8663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (135.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (67.9214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (115.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (172.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (172.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (95.3124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (87.8757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (53.1289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (50.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (407.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (86.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (472.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (503.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (216.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (199.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (69.7068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (80.0536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (99.7226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (172.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (83.9536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (-18.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (209.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (166.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (204.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (250.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (101.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (209.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (177.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (83.9536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (142.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (181.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (323.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (238.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (257.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (210.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS79',
    E0 = (130.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS80',
    E0 = (173.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS81',
    E0 = (230.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS82',
    E0 = (123.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS83',
    E0 = (107.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS84',
    E0 = (209.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS85',
    E0 = (204.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS86',
    E0 = (184.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS87',
    E0 = (101.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS88',
    E0 = (138.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS89',
    E0 = (209.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS90',
    E0 = (105.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS91',
    E0 = (303.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS92',
    E0 = (154.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS93',
    E0 = (205.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS94',
    E0 = (163.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS95',
    E0 = (-53.4996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS96',
    E0 = (236.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS97',
    E0 = (368.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS98',
    E0 = (332.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS99',
    E0 = (206.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS100',
    E0 = (276.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS101',
    E0 = (366.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction589',
    reactants = ['H(3)(3)', 'C=[C]C=C(O)C=CC=O(3104)'],
    products = ['S(779)(778)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction590',
    reactants = ['H(3)(3)', 'C=CC=C(O)C=C[C]=O(3105)'],
    products = ['S(779)(778)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction645',
    reactants = ['H(3)(3)', 'C=CC=C(O)[C]=CC=O(3101)'],
    products = ['S(779)(778)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction621',
    reactants = ['OH(5)(5)', 'C=CC=[C]C=CC=O(3097)'],
    products = ['S(779)(778)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(41610,'cm^3/(mol*s)'), n=2.487, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds;OJ_pri] for rate rule [Ca_Cds-CdH;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from -7.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction658',
    reactants = ['H(3)(3)', 'C=C[CH]C(=O)C=CC=O(2764)'],
    products = ['S(779)(778)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.185e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""From training reaction 184 used for Od_R;HJ
Exact match found for rate rule [Od_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction670',
    reactants = ['H(3)(3)', 'C=C[C]=C(O)C=CC=O(3100)'],
    products = ['S(779)(778)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction635',
    reactants = ['H(3)(3)', 'C=CC=C(O)C=[C]C=O(3103)'],
    products = ['S(779)(778)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction603',
    reactants = ['[CH]=C[CH2](891)', '[O]C=CC=[C]O(3148)'],
    products = ['S(779)(778)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad/NonDe;Y_rad] for rate rule [Cd_rad/NonDe;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction604',
    reactants = ['CHCHO(45)(45)', '[CH]=C(O)C=C[CH2](3149)'],
    products = ['S(779)(778)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using an average for rate rule [Cd_pri_rad;Cd_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction611',
    reactants = ['S(779)(778)'],
    products = ['C=C=CC(O)=CC=CO(3156)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction612',
    reactants = ['S(779)(778)'],
    products = ['CC=CC(O)=CC=C=O(3157)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction653',
    reactants = ['S(779)(778)'],
    products = ['CC=CC(O)=C=CC=O(3130)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction654',
    reactants = ['S(779)(778)'],
    products = ['C=C=CC(O)=CCC=O(3129)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction632',
    reactants = ['S(779)(778)'],
    products = ['C=CCC(O)=CC=C=O(3131)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction676',
    reactants = ['S(779)(778)'],
    products = ['CC=C=C(O)C=CC=O(3181)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction677',
    reactants = ['S(779)(778)'],
    products = ['C=C=CC(O)C=CC=O(2784)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction678',
    reactants = ['S(779)(778)'],
    products = ['CC=CC(=O)C=CC=O(2800)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radEndo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction667',
    reactants = ['S(779)(778)'],
    products = ['C=CCC(O)=C=CC=O(2838)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction642',
    reactants = ['S(779)(778)'],
    products = ['C=CC=C(O)C=C=CO(3134)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction643',
    reactants = ['S(779)(778)'],
    products = ['C=CC=C(O)CC=C=O(3172)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction613',
    reactants = ['S(779)(778)'],
    products = ['[CH2]CC=C(O)C=[C]C=O(3118)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.92857e+07,'s^-1'), n=1.56, Ea=(259.91,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for 1_3_pentadiene;CH=C_1;CdCJ_2
Exact match found for rate rule [1_3_pentadiene;CH=C_1;CdCJ_2]
Euclidian distance = 0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction614',
    reactants = ['C=C[C]=C(O)C=CC[O](3120)'],
    products = ['S(779)(778)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.6907e+09,'s^-1'), n=1.02346, Ea=(163.761,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH_end;CddC_2] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH_end;CddC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=CC[C](O)C=C[C]=O(2824)'],
    products = ['S(779)(778)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.94565e+09,'s^-1'), n=0.909333, Ea=(116.555,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH(CJ)_1;unsaturated_end] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH(CJ)_1;unsaturated_end]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction24',
    reactants = ['S(779)(778)'],
    products = ['C=[C]C=C(O)[CH]CC=O(3117)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.92857e+07,'s^-1'), n=1.56, Ea=(259.91,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for 1_3_pentadiene;CH=C_1;CdCJ_2
Exact match found for rate rule [1_3_pentadiene;CH=C_1;CdCJ_2]
Euclidian distance = 0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C[CH]C=C(O)[C]=CC=O(3116)'],
    products = ['S(779)(778)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.53605e+09,'s^-1'), n=1.02346, Ea=(163.761,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH_end;CddC_2] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH3_1;CddC_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction615',
    reactants = ['S(779)(778)'],
    products = ['OC1C=CCOC=CC=1(3158)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Ypri_rad_out] for rate rule [R8;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction655',
    reactants = ['S(779)(778)'],
    products = ['O=CC1C=C(O)C=CC1(3096)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [R6;C_rad_out_H/OneDe;Cpri_rad_out_2H] for rate rule [R6_SDSDS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction633',
    reactants = ['S(779)(778)'],
    products = ['C=CC1OC=CC=C1O(3167)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R6_SDSDS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction679',
    reactants = ['S(779)(778)'],
    products = ['S(826)(825)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SDS;C_rad_out_OneDe/O;Cpri_rad_out_2H]
Euclidian distance = 3.16227766017
family: Birad_recombination"""),
)

reaction(
    label = 'reaction668',
    reactants = ['S(779)(778)'],
    products = ['S(825)(824)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] for rate rule [R4_SDS;C_rad_out_H/OneDe;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction644',
    reactants = ['S(779)(778)'],
    products = ['C=CC=C(O)C1C=CO1(3173)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['S(779)(778)'],
    products = ['[O]C=C[CH]C1(O)C=CC1(3140)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;C=C_1;CdH2_2]
Euclidian distance = 1.0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction586',
    reactants = ['S(779)(778)'],
    products = ['[CH2]C=C[C](O)C1C=CO1(3141)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra] for rate rule [R5_SD_D;doublebond_intra_NdCd_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction587',
    reactants = ['S(779)(778)'],
    products = ['[O][CH]C1C=C(O)C=CC1(3142)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(8.21347e+10,'s^-1'), n=0.43, Ea=(142.403,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7_SMSM_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 141.4 to 142.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction588',
    reactants = ['S(779)(778)'],
    products = ['[CH2][CH]C1OC=CC=C1O(3143)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(9.291e+11,'s^-1'), n=0.234, Ea=(126.717,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 325 used for R7_SMSM_D;doublebond_intra;radadd_intra_O
Exact match found for rate rule [R7_SMSM_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 123.5 to 126.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction591',
    reactants = ['C[C]=CC(O)=CC=C[O](3144)'],
    products = ['S(779)(778)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(7.32e+09,'s^-1'), n=1.12, Ea=(172.799,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 152 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction592',
    reactants = ['[CH2]C=CC(O)=CC=[C]O(3145)'],
    products = ['S(779)(778)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction593',
    reactants = ['CC=[C]C(O)=CC=C[O](3146)'],
    products = ['S(779)(778)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 288 used for R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction594',
    reactants = ['[CH2]C=CC(O)=C[C]=CO(3147)'],
    products = ['S(779)(778)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction595',
    reactants = ['C=CC=C(O)[C]=C[CH]O(3113)'],
    products = ['S(779)(778)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.11e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction596',
    reactants = ['CC=C[C]([O])C=CC=O(2794)'],
    products = ['S(779)(778)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.33e+06,'s^-1'), n=1.64, Ea=(100.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SS(D)MS;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_SS(D)MS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction598',
    reactants = ['C[CH]C=C(O)C=[C]C=O(3121)'],
    products = ['S(779)(778)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6H_SMSMS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction599',
    reactants = ['S(779)(778)'],
    products = ['C=C[CH]C(=O)C=C[CH]O(2830)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(14790.5,'s^-1'), n=1.91656, Ea=(92.3659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6H;O_rad_out;XH_out] + [R6H_SMSMS;Y_rad_out;XH_out] for rate rule [R6H_SMSMS;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction600',
    reactants = ['C=C[C]=C(O)C=C[CH]O(3123)'],
    products = ['S(779)(778)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(280534,'s^-1'), n=2.03283, Ea=(131.895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R6H_SMSMS;Y_rad_out;XH_out] for rate rule [R6H_SMSMS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction601',
    reactants = ['C[CH]C=C(O)C=C[C]=O(3127)'],
    products = ['S(779)(778)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(244756,'s^-1'), n=1.235, Ea=(80.8557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;Y_rad_out;Cs_H_out_2H] for rate rule [R7H;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction602',
    reactants = ['C=[C]C=C(O)C=C[CH]O(3126)'],
    products = ['S(779)(778)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.21847e+06,'s^-1'), n=1.22418, Ea=(82.4275,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;Y_rad_out;XH_out] for rate rule [R7H;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction605',
    reactants = ['S(779)(778)'],
    products = ['[O]C=CC=C(O)C1[CH]C1(3150)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 142 used for R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction606',
    reactants = ['S(779)(778)'],
    products = ['[CH2]C=CC(O)=CC1[CH]O1(3151)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri_HCd;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction607',
    reactants = ['S(779)(778)'],
    products = ['[O]C=CC1CC=C[C]1O(3152)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(5.41935e+11,'s^-1'), n=0.104408, Ea=(148.812,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R5_SD_D;doublebond_intra_secNd_HCd;radadd_intra_cs2H]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction608',
    reactants = ['S(779)(778)'],
    products = ['[CH2]C=CC1(O)[CH]C=CO1(3153)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.25345e+11,'s^-1'), n=0.197304, Ea=(141.375,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra_pri_NdCd;radadd_intra] for rate rule [R5_SD_D;doublebond_intra_pri_NdCd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction609',
    reactants = ['S(779)(778)'],
    products = ['[O]C1[CH]C=C(O)C=CC1(3154)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(4.72071e+06,'s^-1'), n=0.958205, Ea=(106.629,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R7_SDSD_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction610',
    reactants = ['S(779)(778)'],
    products = ['[CH2]C1[CH]C(O)=CC=CO1(3155)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(6.57585e+08,'s^-1'), n=0.579103, Ea=(103.941,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SDSD_D;doublebond_intra_pri;radadd_intra] for rate rule [R7_SDSD_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction616',
    reactants = ['[CH2]C=C[C]=CC=COO(3159)'],
    products = ['S(779)(778)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(9.91772e+09,'s^-1'), n=0, Ea=(65.2704,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OOH;Y_rad_out] for rate rule [R4OOH_DSD;Cd_rad_out_De]
Euclidian distance = 2.2360679775
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction617',
    reactants = ['S(779)(778)'],
    products = ['C=C[CH]C(=O)C[CH]C=O(2813)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction618',
    reactants = ['CH2(17)(18)', '[CH]=CC(O)=CC=C[O](3160)'],
    products = ['S(779)(778)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction619',
    reactants = ['O(4)(4)', '[CH]=CC=C(O)C=C[CH2](2092)'],
    products = ['S(779)(778)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction620',
    reactants = ['S(779)(778)'],
    products = ['C=CC1C(O)=CC1[CH][O](3161)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(269.534,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra_cs] for rate rule [R5_SD_D;doublebond_intra;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 269.5 to 269.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction646',
    reactants = ['[CH2]C=CC(O)=[C]CC=O(3174)'],
    products = ['S(779)(778)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.82494e+10,'s^-1'), n=0.9, Ea=(134.725,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/OneDe] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction647',
    reactants = ['C=CC=C(O)[CH]C[C]=O(3109)'],
    products = ['S(779)(778)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(2.9636e+08,'s^-1'), n=1.49631, Ea=(82.9758,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction648',
    reactants = ['C=C[CH]C(=O)[CH]CC=O(2816)'],
    products = ['S(779)(778)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(2.30755e+06,'s^-1'), n=1.805, Ea=(115.478,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction649',
    reactants = ['C=C[C]=C(O)[CH]CC=O(3112)'],
    products = ['S(779)(778)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(4.05862e+06,'s^-1'), n=1.71075, Ea=(73.9568,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SDS;Cd_rad_out_Cd;Cs_H_out] + [R4H;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction651',
    reactants = ['S(779)(778)'],
    products = ['[CH2]C=CC1(O)[CH]C1C=O(3175)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri_NdCd;radadd_intra_csHCO]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction652',
    reactants = ['S(779)(778)'],
    products = ['[CH2]C1[CH]C(O)=CC1C=O(3176)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(8.8675e+11,'s^-1'), n=-0.37996, Ea=(137.453,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra_pri;radadd_intra_csHDe] for rate rule [R5_SD_D;doublebond_intra_pri;radadd_intra_csHCO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction631',
    reactants = ['S(779)(778)'],
    products = ['OC1[CH][CH]COC=CC=1(3166)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1.06679e+09,'s^-1'), n=0.473387, Ea=(35.0186,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_2H;radadd_intra] for rate rule [R8_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction622',
    reactants = ['C=[C]C[C](O)C=CC=O(2815)'],
    products = ['S(779)(778)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(2.3e+10,'s^-1'), n=0.98, Ea=(112.55,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction623',
    reactants = ['C=CC[C]([O])C=CC=O(2769)'],
    products = ['S(779)(778)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(587605,'s^-1'), n=2.09, Ea=(169.87,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_2Cd;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction624',
    reactants = ['C=CCC(O)=[C]C=C[O](3162)'],
    products = ['S(779)(778)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction625',
    reactants = ['[CH]=CC[C](O)C=CC=O(2823)'],
    products = ['S(779)(778)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction626',
    reactants = ['C=CC[C](O)C=[C]C=O(2817)'],
    products = ['S(779)(778)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction628',
    reactants = ['[CH]=C[CH]C(O)=CC=CO(3163)'],
    products = ['S(779)(778)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R8Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction629',
    reactants = ['S(779)(778)'],
    products = ['C=CC1[C](O)C1C=C[O](3164)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra_secNd_HCd;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction630',
    reactants = ['S(779)(778)'],
    products = ['C=CC1C(O)=C[CH]C1[O](3165)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(8.8675e+11,'s^-1'), n=-0.37996, Ea=(137.453,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SD_D;doublebond_intra_pri;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction657',
    reactants = ['S(779)(778)'],
    products = ['C=C[CH]C1(O)C=CC1[O](3178)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SD_CO;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction153',
    reactants = ['S(458)(457)'],
    products = ['S(779)(778)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction675',
    reactants = ['[CH2]C=[C]C(O)C=CC=O(2792)'],
    products = ['S(779)(778)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_Cd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction659',
    reactants = ['C=C[CH]C(O)[C]=CC=O(2771)'],
    products = ['S(779)(778)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction660',
    reactants = ['C=[C][CH]C(O)C=CC=O(2774)'],
    products = ['S(779)(778)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_Cd] for rate rule [R3HJ;Cd_rad_out_Cd;Cs_H_out_Cd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction661',
    reactants = ['C=C[CH]C(O)C=[C]C=O(2775)'],
    products = ['S(779)(778)'],
    transitionState = 'TS78',
    kinetics = Arrhenius(A=(8.2826e+06,'s^-1'), n=1.67955, Ea=(176.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction662',
    reactants = ['C=C[CH]C(O)C=C[C]=O(2778)'],
    products = ['S(779)(778)'],
    transitionState = 'TS79',
    kinetics = Arrhenius(A=(8.08094e+07,'s^-1'), n=1.1965, Ea=(138.911,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;XH_out] for rate rule [R4H_SDS;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction80',
    reactants = ['S(779)(778)'],
    products = ['[CH2]C1[CH]C1(O)C=CC=O(3177)'],
    transitionState = 'TS80',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(226.622,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_csNdCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 225.9 to 226.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction664',
    reactants = ['CO(10)(11)', 'C=C[CH][C](O)C=C(3179)'],
    products = ['S(779)(778)'],
    transitionState = 'TS81',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction666',
    reactants = ['[CH]=CC=C(O)[CH]CC=O(3122)'],
    products = ['S(779)(778)'],
    transitionState = 'TS82',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_H/CO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction634',
    reactants = ['S(779)(778)'],
    products = ['[CH2]C1C=C(O)C1C=C[O](3168)'],
    transitionState = 'TS83',
    kinetics = Arrhenius(A=(2.17094e+10,'s^-1'), n=0.2405, Ea=(161.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SM;doublebond_intra_2H_pri;radadd_intra_cs] + [R5_SD_D;doublebond_intra;radadd_intra_cs] for rate rule [R5_SD_D;doublebond_intra_2H_pri;radadd_intra_csHCd]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction636',
    reactants = ['C=CC=C(O)C[C]=C[O](3169)'],
    products = ['S(779)(778)'],
    transitionState = 'TS84',
    kinetics = Arrhenius(A=(2.3e+10,'s^-1'), n=0.98, Ea=(112.55,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction637',
    reactants = ['C=C[C]=C(O)C[CH]C=O(3110)'],
    products = ['S(779)(778)'],
    transitionState = 'TS85',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction638',
    reactants = ['C=CC=C(O)CC=[C][O](3170)'],
    products = ['S(779)(778)'],
    transitionState = 'TS86',
    kinetics = Arrhenius(A=(746634,'s^-1'), n=2.07188, Ea=(142.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction639',
    reactants = ['C=[C]C=C(O)C[CH]C=O(3115)'],
    products = ['S(779)(778)'],
    transitionState = 'TS87',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction640',
    reactants = ['[CH]=CC=C(O)C[CH]C=O(3119)'],
    products = ['S(779)(778)'],
    transitionState = 'TS88',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R5H_DSMS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction641',
    reactants = ['S(779)(778)'],
    products = ['C=CC=C(O)C1[CH]C1[O](3171)'],
    transitionState = 'TS89',
    kinetics = Arrhenius(A=(2.22857e+12,'s^-1'), n=0, Ea=(262.828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;doublebond_intra_pri;radadd_intra_csH(CdCdCd)] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_csH(CdCdCd)]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 260.3 to 262.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction669',
    reactants = ['S(779)(778)'],
    products = ['[CH2]C1C=C(O)C=CC1[O](3180)'],
    transitionState = 'TS90',
    kinetics = Arrhenius(A=(8.21347e+10,'s^-1'), n=0.43, Ea=(158.864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSM;multiplebond_intra;radadd_intra_cs] for rate rule [R7_SMSM_CO;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 155.9 to 158.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction671',
    reactants = ['[CH2]C[C]=C(O)C=CC=O(3107)'],
    products = ['S(779)(778)'],
    transitionState = 'TS91',
    kinetics = Arrhenius(A=(9.93038e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction672',
    reactants = ['[CH2]C[CH]C(=O)C=CC=O(2810)'],
    products = ['S(779)(778)'],
    transitionState = 'TS92',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction673',
    reactants = ['[CH2]CC=C(O)[C]=CC=O(3114)'],
    products = ['S(779)(778)'],
    transitionState = 'TS93',
    kinetics = Arrhenius(A=(2.22e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction674',
    reactants = ['[CH2]CC=C(O)C=C[C]=O(3125)'],
    products = ['S(779)(778)'],
    transitionState = 'TS94',
    kinetics = Arrhenius(A=(327.109,'s^-1'), n=2.81, Ea=(112.131,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SMSMS;Y_rad_out;XH_out] for rate rule [R6H_SMSMS;CO_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction665',
    reactants = ['S(779)(778)'],
    products = ['S(778)(777)'],
    transitionState = 'TS95',
    kinetics = Arrhenius(A=(1.58e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Y_12_12] for rate rule [Y_12_12b]
Euclidian distance = 1.0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction663',
    reactants = ['[CH]=C[CH]C(O)C=CC=O(2777)'],
    products = ['S(779)(778)'],
    transitionState = 'TS96',
    kinetics = Arrhenius(A=(3.06452e+09,'s^-1'), n=1.0645, Ea=(151.645,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out;Cs_H_out_Cd] + [RnH;Cd_rad_out_singleH;Cs_H_out_noH] for rate rule [R4HJ_2;Cd_rad_out_singleH;Cs_H_out_Cd]
Euclidian distance = 2.82842712475
family: intra_H_migration"""),
)

reaction(
    label = 'reaction680',
    reactants = ['C=CC=C(O)C=[C]C[O](3108)'],
    products = ['S(779)(778)'],
    transitionState = 'TS97',
    kinetics = Arrhenius(A=(9.93038e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction681',
    reactants = ['C=CC=C(O)[C]=CC[O](3111)'],
    products = ['S(779)(778)'],
    transitionState = 'TS98',
    kinetics = Arrhenius(A=(2.56742e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction682',
    reactants = ['C=C[CH]C(=O)C=CC[O](2827)'],
    products = ['S(779)(778)'],
    transitionState = 'TS99',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1'), n=1.64, Ea=(100.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SS(D)MS;Y_rad_out;XH_out] for rate rule [R5H_SS(D)MS;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction683',
    reactants = ['C=[C]C=C(O)C=CC[O](3124)'],
    products = ['S(779)(778)'],
    transitionState = 'TS100',
    kinetics = Arrhenius(A=(561068,'s^-1'), n=2.03283, Ea=(131.895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R6H_SMSMS;Y_rad_out;XH_out] for rate rule [R6H_SMSMS;Cd_rad_out_Cd;XH_out]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction684',
    reactants = ['[CH]=CC=C(O)C=CC[O](3128)'],
    products = ['S(779)(778)'],
    transitionState = 'TS101',
    kinetics = Arrhenius(A=(3.73886e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7H;Cd_rad_out_singleH;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = '180',
    isomers = [
        'S(779)(778)',
    ],
    reactants = [
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '180',
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

