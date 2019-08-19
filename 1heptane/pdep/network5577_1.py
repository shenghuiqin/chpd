species(
    label = '[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)',
    structure = SMILES('[CH2]C(C=C)C([CH2])(O)C(O)=CC'),
    E0 = (-21.7374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,1558.78,1600,2932.02,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.51634,0.131825,-0.00012362,5.37432e-08,-7.02178e-12,-2369.51,48.7432], Tmin=(100,'K'), Tmax=(967.443,'K')), NASAPolynomial(coeffs=[26.1957,0.0365554,-1.22553e-05,2.0498e-09,-1.36254e-13,-9022.06,-94.4956], Tmin=(967.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.7374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(Isobutyl)"""),
)

species(
    label = 'butadiene13(2459)',
    structure = SMILES('C=CC=C'),
    E0 = (96.4553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.30712,'amu*angstrom^2'), symmetry=1, barrier=(30.0532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80599,0.0102584,6.1726e-05,-9.01643e-08,3.59117e-11,11658.5,12.0621], Tmin=(100,'K'), Tmax=(946.047,'K')), NASAPolynomial(coeffs=[12.4694,0.0100554,-2.41207e-06,4.57077e-10,-3.93161e-14,8010.78,-43.6375], Tmin=(946.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""butadiene13""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2][CH]C=C(2458)',
    structure = SMILES('[CH2]C=C[CH2]'),
    E0 = (274.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.210311,'amu*angstrom^2'), symmetry=1, barrier=(25.2351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779031,'amu*angstrom^2'), symmetry=1, barrier=(93.4717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56318,0.0223429,1.87067e-05,-3.93099e-08,1.63982e-11,33100.5,13.4097], Tmin=(100,'K'), Tmax=(974.264,'K')), NASAPolynomial(coeffs=[9.82995,0.0151966,-5.22272e-06,9.67656e-10,-7.0786e-14,30607.7,-26.985], Tmin=(974.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)C1(O)CC1(O)[CH]C(32559)',
    structure = SMILES('[CH2]C(C=C)C1(O)CC1(O)[CH]C'),
    E0 = (43.2778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.01872,0.115244,-7.33488e-05,4.73457e-11,1.23862e-11,5437.38,46.0467], Tmin=(100,'K'), Tmax=(948.154,'K')), NASAPolynomial(coeffs=[26.7944,0.0347428,-1.09419e-05,1.83377e-09,-1.25604e-13,-1871.82,-101.157], Tmin=(948.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.2778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1CC1C([CH2])(O)C(O)=CC(32560)',
    structure = SMILES('[CH2]C1CC1C([CH2])(O)C(O)=CC'),
    E0 = (2.96562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09909,0.11514,-6.65243e-05,-1.20148e-08,1.80419e-11,593.704,45.3187], Tmin=(100,'K'), Tmax=(935.468,'K')), NASAPolynomial(coeffs=[28.5445,0.0317191,-9.10192e-06,1.46729e-09,-1.0051e-13,-7222.67,-111.609], Tmin=(935.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.96562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(Isobutyl) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C1(O)C(C=C)CC1(O)[CH]C(32561)',
    structure = SMILES('[CH2]C1(O)C(C=C)CC1(O)[CH]C'),
    E0 = (45.7922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07535,0.115301,-7.19191e-05,1.03067e-09,1.06777e-11,5742.59,43.6278], Tmin=(100,'K'), Tmax=(982.242,'K')), NASAPolynomial(coeffs=[26.9872,0.0368095,-1.29243e-05,2.30436e-09,-1.61878e-13,-1889.56,-105.85], Tmin=(982.242,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.7922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCJCO) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1CC(O)(C(O)=CC)C1[CH2](32562)',
    structure = SMILES('[CH2]C1CC(O)(C(O)=CC)C1[CH2]'),
    E0 = (-9.81485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.65759,0.101663,-2.64468e-05,-5.47082e-08,3.38658e-11,-955.686,43.9627], Tmin=(100,'K'), Tmax=(922.659,'K')), NASAPolynomial(coeffs=[28.2705,0.0310145,-7.66944e-06,1.14566e-09,-7.82946e-14,-8993.88,-111.646], Tmin=(922.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.81485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2]C(O)(C(=C)C=C)C(O)=CC(32563)',
    structure = SMILES('[CH2]C(O)(C(=C)C=C)C(O)=CC'),
    E0 = (-127.287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,180,180,180,180,1600,1600.94,2889.49,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150901,'amu*angstrom^2'), symmetry=1, barrier=(3.4695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150901,'amu*angstrom^2'), symmetry=1, barrier=(3.4695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150901,'amu*angstrom^2'), symmetry=1, barrier=(3.4695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150901,'amu*angstrom^2'), symmetry=1, barrier=(3.4695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150901,'amu*angstrom^2'), symmetry=1, barrier=(3.4695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150901,'amu*angstrom^2'), symmetry=1, barrier=(3.4695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150901,'amu*angstrom^2'), symmetry=1, barrier=(3.4695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.53736,0.128883,-0.000116893,4.75933e-08,-5.75984e-12,-15060.7,43.459], Tmin=(100,'K'), Tmax=(1031.66,'K')), NASAPolynomial(coeffs=[28.2812,0.0331057,-1.21133e-05,2.16388e-09,-1.49986e-13,-22681.5,-112.299], Tmin=(1031.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(C=C)C(=C)O(19913)',
    structure = SMILES('[CH2]C(C=C)C(=C)O'),
    E0 = (45.7038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,308.967,562.42],'cm^-1')),
        HinderedRotor(inertia=(0.627096,'amu*angstrom^2'), symmetry=1, barrier=(14.4182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627087,'amu*angstrom^2'), symmetry=1, barrier=(14.418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627088,'amu*angstrom^2'), symmetry=1, barrier=(14.418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0642351,'amu*angstrom^2'), symmetry=1, barrier=(14.4179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0827877,0.0775295,-7.31474e-05,3.58698e-08,-6.86001e-12,5654.3,28.9069], Tmin=(100,'K'), Tmax=(1359.35,'K')), NASAPolynomial(coeffs=[18.099,0.0203041,-5.89137e-06,8.69965e-10,-5.24693e-14,1055.32,-63.1258], Tmin=(1359.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.7038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
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
    label = 'C2H3(30)',
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
    label = '[CH2]C(O)(C=C)C(O)=CC(32267)',
    structure = SMILES('[CH2]C(O)(C=C)C(O)=CC'),
    E0 = (-177.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.21264,0.10183,-9.49263e-05,4.47376e-08,-8.28732e-12,-21149.1,37.2819], Tmin=(100,'K'), Tmax=(1314.41,'K')), NASAPolynomial(coeffs=[23.1597,0.0276596,-1.02824e-05,1.80588e-09,-1.21659e-13,-27556.1,-86.9633], Tmin=(1314.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ)"""),
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
    label = '[CH2]C(C=C)C(=C)C(O)=CC(32564)',
    structure = SMILES('[CH2]C(C=C)C(=C)C(O)=CC'),
    E0 = (60.5507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.199,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.92893,0.122141,-0.000123886,6.62913e-08,-1.4056e-11,7503.1,37.9504], Tmin=(100,'K'), Tmax=(1150.88,'K')), NASAPolynomial(coeffs=[22.4679,0.0373474,-1.33715e-05,2.27412e-09,-1.49956e-13,1887.51,-83.1786], Tmin=(1150.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.5507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(O)([C](C)C=C)C(O)=CC(32565)',
    structure = SMILES('[CH2]C=C(C)C([CH2])(O)C(O)=CC'),
    E0 = (-101.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39742,0.131385,-0.000125465,6.16457e-08,-1.20503e-11,-11918.4,44.3726], Tmin=(100,'K'), Tmax=(1237.34,'K')), NASAPolynomial(coeffs=[25.3593,0.0416547,-1.66864e-05,3.03627e-09,-2.08381e-13,-18787.2,-95.4484], Tmin=(1237.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(C=C)C(C)([O])C(O)=CC(32566)',
    structure = SMILES('[CH2]C(C=C)C(C)([O])C(O)=CC'),
    E0 = (-6.07821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,374.04,762.629,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154507,'amu*angstrom^2'), symmetry=1, barrier=(3.55242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154507,'amu*angstrom^2'), symmetry=1, barrier=(3.55242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154507,'amu*angstrom^2'), symmetry=1, barrier=(3.55242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154507,'amu*angstrom^2'), symmetry=1, barrier=(3.55242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154507,'amu*angstrom^2'), symmetry=1, barrier=(3.55242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154507,'amu*angstrom^2'), symmetry=1, barrier=(3.55242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154507,'amu*angstrom^2'), symmetry=1, barrier=(3.55242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.08254,0.120382,-9.76284e-05,3.45134e-08,-2.51643e-12,-500.35,47.5253], Tmin=(100,'K'), Tmax=(1016.06,'K')), NASAPolynomial(coeffs=[23.9285,0.041126,-1.47911e-05,2.57962e-09,-1.75113e-13,-6980.79,-84.2565], Tmin=(1016.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.07821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2][C](C=C)C(C)(O)C(O)=CC(32567)',
    structure = SMILES('[CH2]C=C([CH2])C(C)(O)C(O)=CC'),
    E0 = (-163.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43772,0.124876,-9.8384e-05,2.96345e-08,-1.41658e-13,-19360.4,43.9208], Tmin=(100,'K'), Tmax=(1033.43,'K')), NASAPolynomial(coeffs=[27.5342,0.0379131,-1.43199e-05,2.60293e-09,-1.82038e-13,-27106.2,-109.167], Tmin=(1033.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-163.019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(O)(C(O)=CC)C(C)[C]=C(25695)',
    structure = SMILES('[CH2]C(O)(C(O)=CC)C(C)[C]=C'),
    E0 = (11.022,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180,1600,1661.6,2837.62,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.61621,0.135288,-0.000137649,7.24423e-08,-1.50272e-11,1572.57,47.7949], Tmin=(100,'K'), Tmax=(1177.66,'K')), NASAPolynomial(coeffs=[25.979,0.0381634,-1.39416e-05,2.41297e-09,-1.61192e-13,-5162.59,-94.8369], Tmin=(1177.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C([O])(C(O)=CC)C(C)C=C(32568)',
    structure = SMILES('[CH2]C([O])(C(O)=CC)C(C)C=C'),
    E0 = (2.28291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.32828,0.125777,-0.000113929,5.33348e-08,-9.92699e-12,513.618,47.508], Tmin=(100,'K'), Tmax=(1299.44,'K')), NASAPolynomial(coeffs=[25.4638,0.0402262,-1.51741e-05,2.66962e-09,-1.79513e-13,-6709.24,-93.8527], Tmin=(1299.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.28291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C([C]=C)C(C)(O)C(O)=CC(32569)',
    structure = SMILES('[CH2]C([C]=C)C(C)(O)C(O)=CC'),
    E0 = (2.66092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180,1600,1617.35,2879.05,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43708,0.130591,-0.000123354,5.556e-08,-8.14888e-12,561.573,48.0567], Tmin=(100,'K'), Tmax=(975.408,'K')), NASAPolynomial(coeffs=[25.3992,0.0375454,-1.27261e-05,2.13348e-09,-1.41517e-13,-5872.85,-90.6911], Tmin=(975.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.66092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C=C)C(C)(O)C(=O)[CH]C(19921)',
    structure = SMILES('[CH2]C(C=C)C(C)(O)C([O])=CC'),
    E0 = (-97.3761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07503,0.122458,-0.000107632,4.53975e-08,-6.07364e-12,-11483.1,47.7843], Tmin=(100,'K'), Tmax=(998.088,'K')), NASAPolynomial(coeffs=[23.0445,0.0410772,-1.43176e-05,2.43235e-09,-1.61999e-13,-17458.2,-78.168], Tmin=(998.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.3761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(C=C)C(C)(O)C(O)=[C]C(32570)',
    structure = SMILES('[CH2]C(C=C)C(C)(O)C(O)=[C]C'),
    E0 = (2.66092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180,1600,1617.35,2879.05,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152095,'amu*angstrom^2'), symmetry=1, barrier=(3.49696,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43708,0.130591,-0.000123354,5.556e-08,-8.14888e-12,561.573,48.0567], Tmin=(100,'K'), Tmax=(975.408,'K')), NASAPolynomial(coeffs=[25.3992,0.0375454,-1.27261e-05,2.13348e-09,-1.41517e-13,-5872.85,-90.6911], Tmin=(975.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.66092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=CC(C)C([CH2])(O)C(O)=CC(32571)',
    structure = SMILES('[CH]=CC(C)C([CH2])(O)C(O)=CC'),
    E0 = (20.2764,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.85987,0.136573,-0.000137632,7.08537e-08,-1.42727e-11,2697.81,48.6004], Tmin=(100,'K'), Tmax=(1217.28,'K')), NASAPolynomial(coeffs=[28.1907,0.0345415,-1.19034e-05,1.99647e-09,-1.31183e-13,-4861.7,-107.306], Tmin=(1217.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.2764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(O)(C(=O)[CH]C)C(C)C=C(19923)',
    structure = SMILES('[CH2]C(O)(C([O])=CC)C(C)C=C'),
    E0 = (-89.0149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.22063,0.126743,-0.000120399,6.01559e-08,-1.19798e-11,-10473.6,47.4032], Tmin=(100,'K'), Tmax=(1217.39,'K')), NASAPolynomial(coeffs=[23.694,0.0415939,-1.54817e-05,2.70081e-09,-1.80829e-13,-16783.2,-82.7174], Tmin=(1217.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.0149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(O)(C(O)=[C]C)C(C)C=C(32572)',
    structure = SMILES('[CH2]C(O)(C(O)=[C]C)C(C)C=C'),
    E0 = (11.022,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180,1600,1661.6,2837.62,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.61621,0.135288,-0.000137649,7.24423e-08,-1.50272e-11,1572.57,47.7949], Tmin=(100,'K'), Tmax=(1177.66,'K')), NASAPolynomial(coeffs=[25.979,0.0381634,-1.39416e-05,2.41297e-09,-1.61192e-13,-5162.59,-94.8369], Tmin=(1177.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH]=CC([CH2])C(C)(O)C(O)=CC(32573)',
    structure = SMILES('[CH]=CC([CH2])C(C)(O)C(O)=CC'),
    E0 = (11.9153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.99074,0.135488,-0.000135732,6.96497e-08,-1.38911e-11,1700.3,49.9768], Tmin=(100,'K'), Tmax=(1279.18,'K')), NASAPolynomial(coeffs=[28.9344,0.031723,-9.4396e-06,1.42572e-09,-8.75752e-14,-6145.36,-110.646], Tmin=(1279.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.9153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=C(O)C(C)(O)C([CH2])C=C(32574)',
    structure = SMILES('[CH2]C=C(O)C(C)(O)C([CH2])C=C'),
    E0 = (-83.6817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.41218,0.135104,-0.00012934,6.21718e-08,-1.14382e-11,-9773.9,51.3803], Tmin=(100,'K'), Tmax=(1472.64,'K')), NASAPolynomial(coeffs=[31.6246,0.0273592,-6.78109e-06,8.89505e-10,-5.00899e-14,-18729.3,-126.582], Tmin=(1472.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.6817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=C(O)C([CH2])(O)C(C)C=C(32575)',
    structure = SMILES('[CH2]C=C(O)C([CH2])(O)C(C)C=C'),
    E0 = (-75.3206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.53296,0.12762,-0.00010244,2.76656e-08,2.7148e-12,-8809.36,47.3019], Tmin=(100,'K'), Tmax=(971.36,'K')), NASAPolynomial(coeffs=[28.6277,0.0339943,-1.14336e-05,1.97441e-09,-1.36134e-13,-16499.7,-110.55], Tmin=(971.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.3206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)C1(O)CC(C)[C]1O(32576)',
    structure = SMILES('[CH2]C(C=C)C1(O)CC(C)[C]1O'),
    E0 = (24.5263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.82835,0.113062,-7.62797e-05,1.24942e-08,5.25028e-12,3173.05,44.2579], Tmin=(100,'K'), Tmax=(991.904,'K')), NASAPolynomial(coeffs=[23.5989,0.0414554,-1.47696e-05,2.59198e-09,-1.77898e-13,-3392.89,-85.8775], Tmin=(991.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.5263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(O)(C(O)=CC)C1[CH]CC1(24042)',
    structure = SMILES('[CH2]C(O)(C(O)=CC)C1[CH]CC1'),
    E0 = (-9.65614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,350,440,435,1725,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.76178,0.108554,-5.97159e-05,-6.68221e-09,1.22139e-11,-937.814,45.6602], Tmin=(100,'K'), Tmax=(990.637,'K')), NASAPolynomial(coeffs=[25.0684,0.0393508,-1.41831e-05,2.55088e-09,-1.79178e-13,-8173.73,-93.2188], Tmin=(990.637,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.65614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(cyclobutane) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C1(O)[C](O)C(C)CC1C=C(32577)',
    structure = SMILES('[CH2]C1(O)[C](O)C(C)CC1C=C'),
    E0 = (-50.3065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77812,0.110366,-6.66608e-05,1.63175e-09,9.16058e-12,-5827.63,41.6339], Tmin=(100,'K'), Tmax=(990.308,'K')), NASAPolynomial(coeffs=[24.252,0.0405824,-1.45125e-05,2.57687e-09,-1.78934e-13,-12716.9,-92.4461], Tmin=(990.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.3065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CJC(C)2O) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1[CH]CCC1(O)C(O)=CC(32578)',
    structure = SMILES('[CH2]C1[CH]CCC1(O)C(O)=CC'),
    E0 = (-93.0802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.192,0.0905941,-4.95833e-06,-6.87078e-08,3.63338e-11,-10986.4,44.2298], Tmin=(100,'K'), Tmax=(940.362,'K')), NASAPolynomial(coeffs=[25.7772,0.0348683,-1.01692e-05,1.69897e-09,-1.20087e-13,-18666.9,-98.0911], Tmin=(940.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-93.0802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopentane) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC(=C)C(C)(O)C(O)=CC(32579)',
    structure = SMILES('C=CC(=C)C(C)(O)C(O)=CC'),
    E0 = (-340.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.37145,0.122424,-8.96983e-05,1.74624e-08,5.16718e-12,-40735.3,41.9876], Tmin=(100,'K'), Tmax=(995.841,'K')), NASAPolynomial(coeffs=[28.3478,0.0351419,-1.26179e-05,2.27154e-09,-1.59938e-13,-48644,-115.077], Tmin=(995.841,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-340.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(C=C)C([CH2])(O)C(=C)O(32580)',
    structure = SMILES('[CH2]C(C=C)C([CH2])(O)C(=C)O'),
    E0 = (14.2882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,217.745,1514.64,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.148237,'amu*angstrom^2'), symmetry=1, barrier=(3.40827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148237,'amu*angstrom^2'), symmetry=1, barrier=(3.40827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148237,'amu*angstrom^2'), symmetry=1, barrier=(3.40827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148237,'amu*angstrom^2'), symmetry=1, barrier=(3.40827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148237,'amu*angstrom^2'), symmetry=1, barrier=(3.40827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148237,'amu*angstrom^2'), symmetry=1, barrier=(3.40827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148237,'amu*angstrom^2'), symmetry=1, barrier=(3.40827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45944,0.123148,-0.000126909,6.5825e-08,-1.30971e-11,1967.18,46.735], Tmin=(100,'K'), Tmax=(1334.03,'K')), NASAPolynomial(coeffs=[27.9172,0.0238481,-6.01456e-06,7.9153e-10,-4.42813e-14,-5406.28,-105.829], Tmin=(1334.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.2882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C=C)CC(O)=C(O)[CH]C(19871)',
    structure = SMILES('[CH2]C(C=C)CC(O)=C(O)[CH]C'),
    E0 = (-50.856,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.2674,0.117985,-7.07607e-05,-1.09281e-08,1.83917e-11,-5872.68,48.6469], Tmin=(100,'K'), Tmax=(935.376,'K')), NASAPolynomial(coeffs=[30.1841,0.0290283,-7.99512e-06,1.27602e-09,-8.82853e-14,-14122.9,-117.394], Tmin=(935.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C(C=C)[C](O)CC(O)=CC(32581)',
    structure = SMILES('[CH2]C(C=C)[C](O)CC(O)=CC'),
    E0 = (-33.1768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.28801,0.133062,-0.000137584,7.557e-08,-1.65533e-11,-3759.32,47.2827], Tmin=(100,'K'), Tmax=(1110.36,'K')), NASAPolynomial(coeffs=[22.5747,0.0434959,-1.65878e-05,2.92319e-09,-1.96699e-13,-9280.64,-75.2684], Tmin=(1110.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.1768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C=CCC([CH2])(O)C(O)=CC(19849)',
    structure = SMILES('[CH2]C=CCC([CH2])(O)C(O)=CC'),
    E0 = (-79.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,1557.28,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.82935,0.134454,-0.000131501,6.56377e-08,-1.28369e-11,-9313.74,48.1778], Tmin=(100,'K'), Tmax=(1251.75,'K')), NASAPolynomial(coeffs=[28.3564,0.0347983,-1.20804e-05,2.03508e-09,-1.33991e-13,-17121,-109.277], Tmin=(1251.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(O)([CH]CC=C)C(O)=CC(32582)',
    structure = SMILES('[CH2]C(O)([CH]CC=C)C(O)=CC'),
    E0 = (-18.9458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,527.679,611.057,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157557,'amu*angstrom^2'), symmetry=1, barrier=(3.62254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157557,'amu*angstrom^2'), symmetry=1, barrier=(3.62254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157557,'amu*angstrom^2'), symmetry=1, barrier=(3.62254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157557,'amu*angstrom^2'), symmetry=1, barrier=(3.62254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157557,'amu*angstrom^2'), symmetry=1, barrier=(3.62254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157557,'amu*angstrom^2'), symmetry=1, barrier=(3.62254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157557,'amu*angstrom^2'), symmetry=1, barrier=(3.62254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157557,'amu*angstrom^2'), symmetry=1, barrier=(3.62254,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.70699,0.131814,-0.000127873,6.34592e-08,-1.23498e-11,-2023.91,50.3931], Tmin=(100,'K'), Tmax=(1257.28,'K')), NASAPolynomial(coeffs=[27.8221,0.0346879,-1.19968e-05,2.01705e-09,-1.32618e-13,-9700.7,-103.882], Tmin=(1257.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.9458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'C=CC1CCC1(O)C(O)=CC(32195)',
    structure = SMILES('C=CC1CCC1(O)C(O)=CC'),
    E0 = (-282.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46236,0.0935432,-3.47843e-06,-7.4671e-08,3.88303e-11,-33706.3,41.1783], Tmin=(100,'K'), Tmax=(952.196,'K')), NASAPolynomial(coeffs=[28.6592,0.0322747,-9.77672e-06,1.72309e-09,-1.26677e-13,-42401.4,-118.202], Tmin=(952.196,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-282.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)',
    structure = SMILES('[CH2]C(C=C)C([CH2])(O)C(=O)CC'),
    E0 = (-22.2369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,428.544,718.295,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4659.56,'J/mol'), sigma=(7.82797,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=727.81 K, Pc=22.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45095,0.14396,-0.000169003,1.07491e-07,-2.73862e-11,-2443.92,45.8821], Tmin=(100,'K'), Tmax=(957.061,'K')), NASAPolynomial(coeffs=[20.3257,0.0487695,-1.98159e-05,3.57396e-09,-2.42303e-13,-6803.79,-63.003], Tmin=(957.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.2369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'C=C[CH]CC(C)[C](O)C(=C)O(32186)',
    structure = SMILES('[CH2]C=CCC(C)C(O)=C([CH2])O'),
    E0 = (-158.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.53355,0.123505,-8.0549e-05,-2.90135e-09,1.56744e-11,-18846.1,45.1873], Tmin=(100,'K'), Tmax=(944.084,'K')), NASAPolynomial(coeffs=[31.454,0.0287737,-8.31977e-06,1.37422e-09,-9.64482e-14,-27459.3,-128.457], Tmin=(944.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
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
    label = '[CH2]C(C=C)[C](O)C(O)=CC(32583)',
    structure = SMILES('[CH2]C(C=C)C(O)=C(O)[CH]C'),
    E0 = (-27.6238,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.55167,0.122816,-0.000125198,6.3671e-08,-1.2366e-11,-3068.34,45.0407], Tmin=(100,'K'), Tmax=(1382.67,'K')), NASAPolynomial(coeffs=[28.8193,0.0221563,-5.25106e-06,6.568e-10,-3.56122e-14,-10796.7,-113.047], Tmin=(1382.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.6238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(O)([CH]C=C)C(O)=CC(32584)',
    structure = SMILES('[CH2]C(O)([CH]C=C)C(O)=CC'),
    E0 = (-124.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58024,0.103204,-5.6109e-05,-1.61781e-08,1.79052e-11,-14806.3,42.4014], Tmin=(100,'K'), Tmax=(954.513,'K')), NASAPolynomial(coeffs=[27.986,0.0259365,-7.96676e-06,1.38094e-09,-9.93704e-14,-22574.9,-109.99], Tmin=(954.513,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(C=CCJC(O)C=C)"""),
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
    E0 = (-21.7374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (43.2778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (15.5002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (45.7922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (37.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (84.5056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-21.7374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (126.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-4.44454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (122.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (95.5265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (102.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (132.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (96.2017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (162.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (61.9426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (46.9695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (61.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (46.9695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (64.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (20.9394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (44.0622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (44.9555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (33.9098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (33.9098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (151.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (103.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (102.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (64.0346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (92.9042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (41.6627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (434.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (135.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (135.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (138.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (226.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-13.4531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (121.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (125.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (353.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (256.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['butadiene13(2459)', 'C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C(C=C)C1(O)CC1(O)[CH]C(32559)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(65.0152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra_cs2H] for rate rule [R4_S_D;doublebond_intra_HNd_secNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 62.9 to 65.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C1CC1C([CH2])(O)C(O)=CC(32560)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C1(O)C(C=C)CC1(O)[CH]C(32561)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.65e+06,'s^-1'), n=1.3, Ea=(67.5296,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_HNd;radadd_intra_cs2H] for rate rule [R5_SS_D;doublebond_intra_HNd_secNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 65.4 to 67.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C1CC(O)(C(O)=CC)C1[CH2](32562)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.07e+06,'s^-1'), n=1.46, Ea=(58.9944,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 343 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH2]C(O)(C(=C)C=C)C(O)=CC(32563)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.88727,'m^3/(mol*s)'), n=1.9687, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 18 used for Cds-CdCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -8.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH]C=C(2458)', 'C=C(O)C(O)=CC(5562)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00784329,'m^3/(mol*s)'), n=2.41519, Ea=(73.439,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-CdOs_Cds;CJ]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 68.1 to 73.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(C=C)C(=C)O(19913)', 'CC=[C]O(31159)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00712612,'m^3/(mol*s)'), n=2.40979, Ea=(7.81798,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CdsJ] for rate rule [Cds-OsCs_Cds;CdsJ-O2s]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['butadiene13(2459)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0136943,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C2H3(30)', '[CH2]C(O)(C=C)C(O)=CC(32267)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6870,'cm^3/(mol*s)'), n=2.41, Ea=(13.7235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CdsJ-H] for rate rule [Cds-Cs\O2s/H_Cds-HH;CdsJ-H]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OH(5)', '[CH2]C(C=C)C(=C)C(O)=CC(32564)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.118912,'m^3/(mol*s)'), n=2.18935, Ea=(6.60377,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdCs_Cds-HH;YJ] for rate rule [Cds-CdCs_Cds-HH;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)([C](C)C=C)C(O)=CC(32565)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C=C)C(C)([O])C(O)=CC(32566)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2][C](C=C)C(C)(O)C(O)=CC(32567)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)(C(O)=CC)C(C)[C]=C(25695)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C([O])(C(O)=CC)C(C)C=C(32568)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([C]=C)C(C)(O)C(O)=CC(32569)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C(C=C)C(C)(O)C(=O)[CH]C(19921)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00963743,'s^-1'), n=3.795, Ea=(83.0524,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_2H;O_H_out] + [R4H_SS(Cd)S;C_rad_out_2H;XH_out] for rate rule [R4H_SS(Cd)S;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(C=C)C(C)(O)C(O)=[C]C(32570)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=CC(C)C([CH2])(O)C(O)=CC(32571)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C(O)(C(=O)[CH]C)C(C)C=C(19923)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);C_rad_out_2H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(O)(C(O)=[C]C)C(C)C=C(32572)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=CC([CH2])C(C)(O)C(O)=CC(32573)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C=C(O)C(C)(O)C([CH2])C=C(32574)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(121000,'s^-1'), n=1.9, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 92 used for R5H_SSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C=C(O)C([CH2])(O)C(C)C=C(32575)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8010,'s^-1'), n=1.94, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 93 used for R6H_RSSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R6H_RSSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]C=C(2458)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C(C=C)C1(O)CC(C)[C]1O(32576)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.82e+08,'s^-1'), n=0.91, Ea=(125.52,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra_secNd;radadd_intra_cs2H] for rate rule [R4_Cs_RR_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C(O)(C(O)=CC)C1[CH]CC1(24042)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C1(O)[C](O)C(C)CC1C=C(32577)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.95277e+09,'s^-1'), n=0.6, Ea=(85.772,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_secNd;radadd_intra_cs2H] for rate rule [R5_SS_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C1[CH]CCC1(O)C(O)=CC(32578)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.12e+09,'s^-1'), n=0.63, Ea=(114.642,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['C=CC(=C)C(C)(O)C(O)=CC(32579)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(S)(23)', '[CH2]C(C=C)C([CH2])(O)C(=C)O(32580)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.97e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C(C=C)CC(O)=C(O)[CH]C(19871)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C(C=C)[C](O)CC(O)=CC(32581)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C=CCC([CH2])(O)C(O)=CC(19849)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(O)([CH]CC=C)C(O)=CC(32582)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['C=CC1CCC1(O)C(O)=CC(32195)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['C=C[CH]CC(C)[C](O)C(=C)O(32186)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CH2(19)', '[CH2]C(C=C)[C](O)C(O)=CC(32583)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction41',
    reactants = ['CH2(19)', '[CH2]C(O)([CH]C=C)C(O)=CC(32584)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '5577',
    isomers = [
        '[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)',
    ],
    reactants = [
        ('butadiene13(2459)', 'C=C(O)C(O)=CC(5562)'),
        ('[CH2][CH]C=C(2458)', 'C=C(O)C(O)=CC(5562)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '5577',
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

