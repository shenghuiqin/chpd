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
    label = '[CH]=C(O)C=CC=O(3099)',
    structure = SMILES('[CH]=C(O)C=CC=O'),
    E0 = (3.63643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(0.953217,'amu*angstrom^2'), symmetry=1, barrier=(21.9163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.956515,'amu*angstrom^2'), symmetry=1, barrier=(21.9922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.957454,'amu*angstrom^2'), symmetry=1, barrier=(22.0138,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.530708,0.0706358,-7.66542e-05,4.01536e-08,-8.14049e-12,567.094,23.1728], Tmin=(100,'K'), Tmax=(1211.47,'K')), NASAPolynomial(coeffs=[17.9004,0.0132847,-5.64376e-06,1.07669e-09,-7.65193e-14,-3641.47,-63.9582], Tmin=(1211.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.63643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_P)"""),
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
    label = '[CH]=CC=C(O)C=CC=O(3106)',
    structure = SMILES('[CH]=CC=C(O)C=CC=O'),
    E0 = (55.4084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3120,650,792.5,1650,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.979195,'amu*angstrom^2'), symmetry=1, barrier=(22.5136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978818,'amu*angstrom^2'), symmetry=1, barrier=(22.5049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.980639,'amu*angstrom^2'), symmetry=1, barrier=(22.5468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.977629,'amu*angstrom^2'), symmetry=1, barrier=(22.4776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01714,0.0980416,-0.000104471,5.35806e-08,-1.05578e-11,6855.17,31.3666], Tmin=(100,'K'), Tmax=(1253.16,'K')), NASAPolynomial(coeffs=[24.8495,0.0154774,-5.64378e-06,1.00592e-09,-6.94083e-14,372.154,-99.2622], Tmin=(1253.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.4084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
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
    label = 'OC12C=CCOC1C=C2(3132)',
    structure = SMILES('OC12C=CCOC1C=C2'),
    E0 = (-74.711,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33195,0.0359044,6.3145e-05,-1.03948e-07,4.04007e-11,-8869.48,24.2536], Tmin=(100,'K'), Tmax=(1001.07,'K')), NASAPolynomial(coeffs=[16.8079,0.0302494,-1.25651e-05,2.53318e-09,-1.91632e-13,-14783.1,-64.4863], Tmin=(1001.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_6_diene_1_6)"""),
)

species(
    label = 'C=C[C]C(O)C=CC=O(3135)',
    structure = SMILES('C=C[C]C(O)C=CC=O'),
    E0 = (194.136,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.190374,0.101429,-0.000145433,1.27634e-07,-4.64688e-11,23491.1,43.4779], Tmin=(100,'K'), Tmax=(757.936,'K')), NASAPolynomial(coeffs=[7.53275,0.0507824,-2.56303e-05,5.0452e-09,-3.56165e-13,22604.4,10.2325], Tmin=(757.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=CC=C(O)[C]CC=O(3136)',
    structure = SMILES('C=CC=C(O)[C]CC=O'),
    E0 = (71.2682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.886834,0.094438,-9.4403e-05,4.67579e-08,-8.97788e-12,8758.57,36.3653], Tmin=(100,'K'), Tmax=(1280.98,'K')), NASAPolynomial(coeffs=[23.0809,0.0195956,-6.76343e-06,1.14697e-09,-7.62372e-14,2618.17,-85.2001], Tmin=(1280.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.2682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(CsJ2_singlet-(Cds-Cds-Cds-Cds)Cs_6_ring)"""),
)

species(
    label = 'C=CC=C(O)C[C]C=O(3137)',
    structure = SMILES('C=CC=C(O)C[C]C=O'),
    E0 = (143.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.586432,0.0946598,-9.51929e-05,4.90188e-08,-9.98564e-12,17421.2,43.0193], Tmin=(100,'K'), Tmax=(1193.89,'K')), NASAPolynomial(coeffs=[19.4267,0.0276082,-1.09497e-05,1.97763e-09,-1.35256e-13,12642.5,-57.0791], Tmin=(1193.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C[C]C=C(O)C=CC=O(3138)',
    structure = SMILES('C[C]C=C(O)C=CC=O'),
    E0 = (53.2351,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.407743,0.0976345,-0.00010565,5.83108e-08,-1.2858e-11,6560.95,31.755], Tmin=(100,'K'), Tmax=(1096.28,'K')), NASAPolynomial(coeffs=[17.875,0.0309256,-1.43738e-05,2.80367e-09,-1.99796e-13,2552.39,-58.129], Tmin=(1096.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.2351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(CsJ2_singlet-(Cds-Cds-Cds-Cds)Cs_6_ring)"""),
)

species(
    label = '[CH]CC=C(O)C=CC=O(3139)',
    structure = SMILES('[CH]CC=C(O)C=CC=O'),
    E0 = (147.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.285253,0.0988285,-0.000113079,6.78107e-08,-1.64551e-11,17930.9,29.9049], Tmin=(100,'K'), Tmax=(993.727,'K')), NASAPolynomial(coeffs=[15.5277,0.0351786,-1.70026e-05,3.35679e-09,-2.40222e-13,14788.1,-46.2841], Tmin=(993.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(CsJ2_singlet-CsH)"""),
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
    label = 'C=CC=C(O)C=C(3133)',
    structure = SMILES('C=CC=C(O)C=C'),
    E0 = (-68.1997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03669,'amu*angstrom^2'), symmetry=1, barrier=(23.8356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03716,'amu*angstrom^2'), symmetry=1, barrier=(23.8463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03781,'amu*angstrom^2'), symmetry=1, barrier=(23.8614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.329266,0.0643075,-1.91421e-05,-3.71501e-08,2.36962e-11,-8055.14,23.0494], Tmin=(100,'K'), Tmax=(922.561,'K')), NASAPolynomial(coeffs=[22.6109,0.0104644,-1.12933e-06,7.84094e-11,-7.75714e-15,-13986.3,-92.5137], Tmin=(922.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.1997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    E0 = (-27.4031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (251.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (157.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (296.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (289.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (223.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (195.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (211.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (221.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (221.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (160.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (267.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (95.7624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (163.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (103.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (229.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (119.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (14.5839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (74.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (150.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (127.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (228.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (169.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-10.4515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (64.9909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (118.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (43.2582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (109.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (65.8936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (141.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (75.4106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (9.13072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-8.77732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (66.6651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (103.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (95.5105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (139.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (123.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (131.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (206.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (66.6651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (104.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (-32.1842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (43.2582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (183.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (76.1726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (72.3076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (12.0389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (240.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (70.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (24.8541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (45.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (199.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (220.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (93.0109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (160.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (79.4928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (165.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (-69.2748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (-69.2748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (117.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (-52.0142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (7.45306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction524',
    reactants = ['S(778)(777)'],
    products = ['O=CC1C=C(O)C=CC1(3096)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.67327e+12,'s^-1'), n=-0.101958, Ea=(164.285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_5_unsaturated_hexane] for rate rule [linear_1_3_5_hexatriene]
Euclidian distance = 1.0
family: Intra_Diels_alder_monocyclic"""),
)

reaction(
    label = 'reaction525',
    reactants = ['OH(5)(5)', 'C=CC=[C]C=CC=O(3097)'],
    products = ['S(778)(777)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/Cd] for rate rule [O_pri_rad;Cd_rad/Cd]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction526',
    reactants = ['H(3)(3)', 'C=C[CH]C(=O)C=CC=O(2764)'],
    products = ['S(778)(777)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=CC=[C]O(3098)', 'CHCHCHO(2768)'],
    products = ['S(778)(777)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.68887e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_rad;Cd_pri_rad] + [Cd_rad/NonDe;Y_rad] for rate rule [Cd_rad/NonDe;Cd_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction528',
    reactants = ['C2H3(28)(29)', '[CH]=C(O)C=CC=O(3099)'],
    products = ['S(778)(777)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 89 used for Cd_pri_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction529',
    reactants = ['H(3)(3)', 'C=C[C]=C(O)C=CC=O(3100)'],
    products = ['S(778)(777)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction530',
    reactants = ['H(3)(3)', 'C=CC=C(O)[C]=CC=O(3101)'],
    products = ['S(778)(777)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction531',
    reactants = ['HCO(14)(15)', '[CH]=CC(O)=CC=C(3102)'],
    products = ['S(778)(777)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.81e+13,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 91 used for Cd_pri_rad;CO_pri_rad
Exact match found for rate rule [Cd_pri_rad;CO_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction532',
    reactants = ['H(3)(3)', 'C=CC=C(O)C=[C]C=O(3103)'],
    products = ['S(778)(777)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction533',
    reactants = ['H(3)(3)', 'C=[C]C=C(O)C=CC=O(3104)'],
    products = ['S(778)(777)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction534',
    reactants = ['H(3)(3)', 'C=CC=C(O)C=C[C]=O(3105)'],
    products = ['S(778)(777)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_sec_rad;Y_rad] for rate rule [CO_rad/OneDe;H_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction535',
    reactants = ['H(3)(3)', '[CH]=CC=C(O)C=CC=O(3106)'],
    products = ['S(778)(777)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction170',
    reactants = ['S(458)(457)'],
    products = ['S(778)(777)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction536',
    reactants = ['[CH2]C[C]=C(O)C=CC=O(3107)'],
    products = ['S(778)(777)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction537',
    reactants = ['C=C[CH]C(O)[C]=CC=O(2771)'],
    products = ['S(778)(777)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction538',
    reactants = ['C=CC=C(O)C=[C]C[O](3108)'],
    products = ['S(778)(777)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction539',
    reactants = ['C=[C]C[C](O)C=CC=O(2815)'],
    products = ['S(778)(777)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction540',
    reactants = ['C=CC=C(O)[CH]C[C]=O(3109)'],
    products = ['S(778)(777)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction541',
    reactants = ['C=C[CH]C(=O)C[CH]C=O(2813)'],
    products = ['S(778)(777)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction542',
    reactants = ['C=C[C]=C(O)C[CH]C=O(3110)'],
    products = ['S(778)(777)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction543',
    reactants = ['C=C[CH]C(O)C=[C]C=O(2775)'],
    products = ['S(778)(777)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction544',
    reactants = ['C=CC=C(O)[C]=CC[O](3111)'],
    products = ['S(778)(777)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction545',
    reactants = ['[CH]=CC[C](O)C=CC=O(2823)'],
    products = ['S(778)(777)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction546',
    reactants = ['C=C[CH]C(=O)[CH]CC=O(2816)'],
    products = ['S(778)(777)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction547',
    reactants = ['C=C[C]=C(O)[CH]CC=O(3112)'],
    products = ['S(778)(777)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction548',
    reactants = ['C=CC[C](O)C=[C]C=O(2817)'],
    products = ['S(778)(777)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(60.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction549',
    reactants = ['C=CC=C(O)[C]=C[CH]O(3113)'],
    products = ['S(778)(777)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction550',
    reactants = ['[CH]=C[CH]C(O)C=CC=O(2777)'],
    products = ['S(778)(777)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction551',
    reactants = ['[CH2]C[CH]C(=O)C=CC=O(2810)'],
    products = ['S(778)(777)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction552',
    reactants = ['[CH2]CC=C(O)[C]=CC=O(3114)'],
    products = ['S(778)(777)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction553',
    reactants = ['C=[C]C=C(O)C[CH]C=O(3115)'],
    products = ['S(778)(777)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_De] for rate rule [R4radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction554',
    reactants = ['C=C[CH]C(O)C=C[C]=O(2778)'],
    products = ['S(778)(777)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction555',
    reactants = ['CC=C[C]([O])C=CC=O(2794)'],
    products = ['S(778)(777)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction556',
    reactants = ['C[CH]C=C(O)[C]=CC=O(3116)'],
    products = ['S(778)(777)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction557',
    reactants = ['C=[C]C=C(O)[CH]CC=O(3117)'],
    products = ['S(778)(777)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction558',
    reactants = ['C=CC[C](O)C=C[C]=O(2824)'],
    products = ['S(778)(777)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction559',
    reactants = ['[CH2]CC=C(O)C=[C]C=O(3118)'],
    products = ['S(778)(777)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction560',
    reactants = ['[CH]=CC=C(O)C[CH]C=O(3119)'],
    products = ['S(778)(777)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction561',
    reactants = ['C=C[CH]C(=O)C=CC[O](2827)'],
    products = ['S(778)(777)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction562',
    reactants = ['C=C[C]=C(O)C=CC[O](3120)'],
    products = ['S(778)(777)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction563',
    reactants = ['C[CH]C=C(O)C=[C]C=O(3121)'],
    products = ['S(778)(777)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction564',
    reactants = ['[CH]=CC=C(O)[CH]CC=O(3122)'],
    products = ['S(778)(777)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction565',
    reactants = ['C=C[CH]C(=O)C=C[CH]O(2830)'],
    products = ['S(778)(777)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction566',
    reactants = ['C=C[C]=C(O)C=C[CH]O(3123)'],
    products = ['S(778)(777)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction567',
    reactants = ['C=[C]C=C(O)C=CC[O](3124)'],
    products = ['S(778)(777)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction568',
    reactants = ['[CH2]CC=C(O)C=C[C]=O(3125)'],
    products = ['S(778)(777)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction569',
    reactants = ['C=[C]C=C(O)C=C[CH]O(3126)'],
    products = ['S(778)(777)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad_De;XH_Rrad] for rate rule [R7radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction570',
    reactants = ['C[CH]C=C(O)C=C[C]=O(3127)'],
    products = ['S(778)(777)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction571',
    reactants = ['[CH]=CC=C(O)C=CC[O](3128)'],
    products = ['S(778)(777)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction572',
    reactants = ['C=C=CC(O)=CCC=O(3129)'],
    products = ['S(778)(777)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(9.66e+08,'s^-1'), n=1.162, Ea=(185.28,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH2(CO)_1;unsaturated_end] for rate rule [1_3_pentadiene;CH2(CO)_1;CddC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction573',
    reactants = ['CC=CC(O)=C=CC=O(3130)'],
    products = ['S(778)(777)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(2.53605e+09,'s^-1'), n=1.02346, Ea=(163.761,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH_end;CddC_2] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH3_1;CddC_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction574',
    reactants = ['C=CCC(O)=CC=C=O(3131)'],
    products = ['S(778)(777)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(9.66e+08,'s^-1'), n=1.162, Ea=(185.28,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1_3_pentadiene;CH2(C)_1;unsaturated_end]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction575',
    reactants = ['OC12C=CCOC1C=C2(3132)'],
    products = ['S(778)(777)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.24e+10,'s^-1'), n=1.27, Ea=(274.47,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [cyclohexene] for rate rule [cyclohexene_2inring]
Euclidian distance = 1.0
family: Intra_Retro_Diels_alder_bicyclic"""),
)

reaction(
    label = 'reaction578',
    reactants = ['C=C[C]C(O)C=CC=O(3135)'],
    products = ['S(778)(777)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(2.33447e+13,'s^-1'), n=-1.27142, Ea=(26.2576,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [CsJ2-C;CsJ2(C=C);CH]
Euclidian distance = 0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction579',
    reactants = ['C=CC=C(O)[C]CC=O(3136)'],
    products = ['S(778)(777)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(5.65514e+12,'s^-1'), n=-0.428961, Ea=(21.7426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2(C=C);CH2(C)] + [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction580',
    reactants = ['C=CC=C(O)C[C]C=O(3137)'],
    products = ['S(778)(777)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;CsJ2C;CH2(C=C)] for rate rule [CsJ2-C;CsJ2C;CH2(C=C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction581',
    reactants = ['C[C]C=C(O)C=CC=O(3138)'],
    products = ['S(778)(777)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(7.00341e+13,'s^-1'), n=-1.27142, Ea=(26.2576,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction582',
    reactants = ['[CH]CC=C(O)C=CC=O(3139)'],
    products = ['S(778)(777)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;singletcarbene;CH2(C=C)] for rate rule [CsJ2-C;CsJ2H;CH2(C=C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction583',
    reactants = ['S(778)(777)'],
    products = ['S(825)(824)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction584',
    reactants = ['S(778)(777)'],
    products = ['S(826)(825)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;C=C_1;CdH2_2]
Euclidian distance = 1.0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction576',
    reactants = ['CO(10)(11)', 'C=CC=C(O)C=C(3133)'],
    products = ['S(778)(777)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction237',
    reactants = ['S(778)(777)'],
    products = ['S(459)(458)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction577',
    reactants = ['C=CC=C(O)C=C=CO(3134)'],
    products = ['S(778)(777)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

network(
    label = '179',
    isomers = [
        'S(778)(777)',
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
    label = '179',
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

