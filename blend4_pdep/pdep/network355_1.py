species(
    label = 'S(244)(243)',
    structure = SMILES('C=C([O])C=CCO'),
    E0 = (-171.328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,269.284,269.284,269.285],'cm^-1')),
        HinderedRotor(inertia=(0.204289,'amu*angstrom^2'), symmetry=1, barrier=(10.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20429,'amu*angstrom^2'), symmetry=1, barrier=(10.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20429,'amu*angstrom^2'), symmetry=1, barrier=(10.5123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4284.46,'J/mol'), sigma=(6.81655,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=669.22 K, Pc=30.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.633106,0.0679166,-6.69863e-05,3.47092e-08,-7.10702e-12,-20479.6,26.3839], Tmin=(100,'K'), Tmax=(1192.86,'K')), NASAPolynomial(coeffs=[14.8693,0.0201783,-6.95604e-06,1.15933e-09,-7.55888e-14,-23875.9,-44.8081], Tmin=(1192.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = 'C5H7O(237)(236)',
    structure = SMILES('C=C=C[CH]CO'),
    E0 = (81.1461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,3025,407.5,1350,352.5,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1976,'amu*angstrom^2'), symmetry=1, barrier=(27.5353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19778,'amu*angstrom^2'), symmetry=1, barrier=(27.5393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19446,'amu*angstrom^2'), symmetry=1, barrier=(27.463,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3721.65,'J/mol'), sigma=(6.20995,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=581.31 K, Pc=35.26 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14012,0.0586746,-4.97917e-05,2.20013e-08,-3.93766e-12,9866.02,21.1068], Tmin=(100,'K'), Tmax=(1327.36,'K')), NASAPolynomial(coeffs=[12.9976,0.0229418,-9.41102e-06,1.71991e-09,-1.17746e-13,6718.22,-39.4564], Tmin=(1327.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.1461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO)"""),
)

species(
    label = 'CH2CO(26)(27)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=CCO(4313)',
    structure = SMILES('[CH]=CCO'),
    E0 = (103.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350],'cm^-1')),
        HinderedRotor(inertia=(0.558245,'amu*angstrom^2'), symmetry=1, barrier=(12.8352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.558882,'amu*angstrom^2'), symmetry=1, barrier=(12.8498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25457,0.0421352,-6.08019e-05,5.22479e-08,-1.79007e-11,12500.9,14.1261], Tmin=(100,'K'), Tmax=(845.742,'K')), NASAPolynomial(coeffs=[5.58243,0.0195677,-8.66594e-06,1.60488e-09,-1.08868e-13,12182.2,0.0724052], Tmin=(845.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_P)"""),
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
    label = '[CH2]C=CC(=C)[O](5186)',
    structure = SMILES('[CH2]C=CC(=C)[O]'),
    E0 = (99.8633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,301.506,301.531],'cm^-1')),
        HinderedRotor(inertia=(0.349175,'amu*angstrom^2'), symmetry=1, barrier=(22.5316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15478,'amu*angstrom^2'), symmetry=1, barrier=(74.524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60392,0.0427776,-7.06827e-06,-2.71088e-08,1.55232e-11,12106.2,20.7036], Tmin=(100,'K'), Tmax=(925.173,'K')), NASAPolynomial(coeffs=[13.8391,0.015431,-4.15959e-06,6.4829e-10,-4.42362e-14,8748.68,-43.2829], Tmin=(925.173,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.8633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
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
    label = 'C=C([O])C=CC[O](5187)',
    structure = SMILES('C=C([O])C=CC[O]'),
    E0 = (54.3769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,228.049,228.049,228.049,2964.37],'cm^-1')),
        HinderedRotor(inertia=(0.387588,'amu*angstrom^2'), symmetry=1, barrier=(14.3039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387588,'amu*angstrom^2'), symmetry=1, barrier=(14.3039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34185,0.0609452,-6.38806e-05,3.89658e-08,-9.95247e-12,6633.63,24.7107], Tmin=(100,'K'), Tmax=(935.551,'K')), NASAPolynomial(coeffs=[8.84339,0.0288723,-1.24575e-05,2.32235e-09,-1.60615e-13,5230,-10.9801], Tmin=(935.551,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.3769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCOJ)"""),
)

species(
    label = 'CH2OH(22)(23)',
    structure = SMILES('[CH2]O'),
    E0 = (-28.7184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3003.59,4000],'cm^-1')),
        HinderedRotor(inertia=(0.057913,'amu*angstrom^2'), symmetry=1, barrier=(25.9304,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (31.0339,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3467.15,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.47834,-0.0013507,2.78485e-05,-3.64869e-08,1.47907e-11,-3500.73,3.30913], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.09314,0.00594761,-2.06497e-06,3.23008e-10,-1.88126e-14,-4034.1,-1.84691], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-28.7184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""CH2OH""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=CC(=C)[O](5188)',
    structure = SMILES('[CH]=CC(=C)[O]'),
    E0 = (264.929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.979636,'amu*angstrom^2'), symmetry=1, barrier=(22.5238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41648,0.0452197,-4.57103e-05,2.25525e-08,-4.15485e-12,31966.6,18.4199], Tmin=(100,'K'), Tmax=(1540.92,'K')), NASAPolynomial(coeffs=[13.6856,0.00610864,-5.68369e-07,-3.63413e-11,6.19928e-15,29047.6,-43.2788], Tmin=(1540.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C([O])C=C[CH]O(5189)',
    structure = SMILES('C=C([O])[CH]C=CO'),
    E0 = (-79.88,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23235,'amu*angstrom^2'), symmetry=1, barrier=(28.3343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24194,'amu*angstrom^2'), symmetry=1, barrier=(28.5545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23521,'amu*angstrom^2'), symmetry=1, barrier=(28.3999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.347664,0.0625043,-1.82904e-05,-4.25449e-08,2.72266e-11,-9459.18,24.5509], Tmin=(100,'K'), Tmax=(911.865,'K')), NASAPolynomial(coeffs=[24.9673,0.00211818,2.72483e-06,-6.50021e-10,4.22106e-14,-15928.6,-102.807], Tmin=(911.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C([O])C=[C]CO(5190)',
    structure = SMILES('C=C([O])C=[C]CO'),
    E0 = (66.5135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1685,370,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,248.404,250.539,257.26],'cm^-1')),
        HinderedRotor(inertia=(0.206828,'amu*angstrom^2'), symmetry=1, barrier=(9.35358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203732,'amu*angstrom^2'), symmetry=1, barrier=(9.35502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213007,'amu*angstrom^2'), symmetry=1, barrier=(9.3875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.812588,0.0699808,-8.28111e-05,5.23006e-08,-1.30966e-11,8114.63,26.2163], Tmin=(100,'K'), Tmax=(978.18,'K')), NASAPolynomial(coeffs=[12.7183,0.021298,-8.16156e-06,1.42646e-09,-9.4944e-14,5785.34,-30.9597], Tmin=(978.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.5135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2][C]=O(1196)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97567e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02477e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'C=C([O])[C]=CCO(5191)',
    structure = SMILES('C=C([O])[C]=CCO'),
    E0 = (27.6672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1685,370,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,252.868,253.076,253.176],'cm^-1')),
        HinderedRotor(inertia=(0.283018,'amu*angstrom^2'), symmetry=1, barrier=(12.8459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282396,'amu*angstrom^2'), symmetry=1, barrier=(12.8483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283096,'amu*angstrom^2'), symmetry=1, barrier=(12.8487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.798664,0.0714463,-8.48341e-05,5.18742e-08,-1.15924e-11,3442.07,25.396], Tmin=(100,'K'), Tmax=(795.488,'K')), NASAPolynomial(coeffs=[12.4266,0.0218011,-7.86112e-06,1.31125e-09,-8.4375e-14,1312.89,-29.7967], Tmin=(795.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.6672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C([O])C=CCO(5192)',
    structure = SMILES('[CH]=C([O])C=CCO'),
    E0 = (75.7679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3615,1277.5,1000,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.484126,'amu*angstrom^2'), symmetry=1, barrier=(11.131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.484489,'amu*angstrom^2'), symmetry=1, barrier=(11.1393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.483731,'amu*angstrom^2'), symmetry=1, barrier=(11.1219,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.67198,0.0701203,-7.90984e-05,4.63345e-08,-1.0641e-11,9235.28,26.6474], Tmin=(100,'K'), Tmax=(1070.2,'K')), NASAPolynomial(coeffs=[14.5378,0.0182956,-6.46091e-06,1.08631e-09,-7.10554e-14,6267.42,-41.188], Tmin=(1070.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.7679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C1OC1[CH]CO(5193)',
    structure = SMILES('C=C1OC1[CH]CO'),
    E0 = (-1.73712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.933465,0.0481514,1.4179e-05,-6.68955e-08,3.30308e-11,-80.7424,24.9615], Tmin=(100,'K'), Tmax=(933.493,'K')), NASAPolynomial(coeffs=[21.7518,0.00812364,-5.24655e-07,4.06714e-11,-9.46511e-15,-6110.21,-85.5189], Tmin=(933.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1.73712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CCJCO)"""),
)

species(
    label = 'C=C(O)[C]=CCO(5194)',
    structure = SMILES('C=C(O)[C]=CCO'),
    E0 = (-110.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3580,3650,1210,1345,900,1100,1685,370,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.775517,'amu*angstrom^2'), symmetry=1, barrier=(17.8307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.775619,'amu*angstrom^2'), symmetry=1, barrier=(17.833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.775299,'amu*angstrom^2'), symmetry=1, barrier=(17.8257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.775464,'amu*angstrom^2'), symmetry=1, barrier=(17.8294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.339754,0.0766289,-8.21415e-05,4.20261e-08,-7.16981e-12,-13111.1,25.5346], Tmin=(100,'K'), Tmax=(899.628,'K')), NASAPolynomial(coeffs=[16.3195,0.0181612,-5.63463e-06,8.77823e-10,-5.53956e-14,-16495.4,-52.698], Tmin=(899.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C(O)C=CCO(5195)',
    structure = SMILES('[CH]=C(O)C=CCO'),
    E0 = (-62.0369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3580,3650,1210,1345,900,1100,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.702054,'amu*angstrom^2'), symmetry=1, barrier=(16.1416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701864,'amu*angstrom^2'), symmetry=1, barrier=(16.1372,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701825,'amu*angstrom^2'), symmetry=1, barrier=(16.1363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.702109,'amu*angstrom^2'), symmetry=1, barrier=(16.1429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00548164,0.077622,-8.38364e-05,4.49828e-08,-9.24828e-12,-7308.76,27.5385], Tmin=(100,'K'), Tmax=(1266.38,'K')), NASAPolynomial(coeffs=[19.2492,0.0133083,-3.47661e-06,4.77266e-10,-2.77257e-14,-11899.6,-68.7282], Tmin=(1266.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.0369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C=[C]CO(5196)',
    structure = SMILES('C=C(O)C=[C]CO'),
    E0 = (-71.2913,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3580,3650,1210,1345,900,1100,1685,370,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.629808,'amu*angstrom^2'), symmetry=1, barrier=(14.4805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.631366,'amu*angstrom^2'), symmetry=1, barrier=(14.5164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.631652,'amu*angstrom^2'), symmetry=1, barrier=(14.5229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.632358,'amu*angstrom^2'), symmetry=1, barrier=(14.5392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.328338,0.0753916,-8.05094e-05,4.22025e-08,-8.13531e-12,-8437.4,26.4498], Tmin=(100,'K'), Tmax=(979.982,'K')), NASAPolynomial(coeffs=[16.8171,0.0173169,-5.74231e-06,9.48226e-10,-6.2295e-14,-12112.3,-55.0265], Tmin=(979.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.2913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=C(O)C=C[CH]O(5197)',
    structure = SMILES('C=C(O)[CH]C=CO'),
    E0 = (-217.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0155296,0.0664323,-1.05458e-05,-6.00699e-08,3.55216e-11,-26016.4,24.3539], Tmin=(100,'K'), Tmax=(911.498,'K')), NASAPolynomial(coeffs=[28.7011,-0.00123313,4.77711e-06,-1.04093e-09,6.75772e-14,-33675.5,-124.824], Tmin=(911.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-217.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C(O)C=CC[O](5198)',
    structure = SMILES('C=C(O)C=CC[O]'),
    E0 = (-83.4279,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,270.824,270.866,271.095],'cm^-1')),
        HinderedRotor(inertia=(0.255396,'amu*angstrom^2'), symmetry=1, barrier=(13.2927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255088,'amu*angstrom^2'), symmetry=1, barrier=(13.2929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255191,'amu*angstrom^2'), symmetry=1, barrier=(13.2927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.799624,0.0671004,-6.44098e-05,3.27439e-08,-6.68785e-12,-9916.05,25.1476], Tmin=(100,'K'), Tmax=(1181.15,'K')), NASAPolynomial(coeffs=[13.5923,0.0237777,-9.39212e-06,1.69072e-09,-1.1519e-13,-12938.1,-38.6994], Tmin=(1181.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.4279,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = 'OCC=C[C]1CO1(5199)',
    structure = SMILES('OC[CH]C=C1CO1'),
    E0 = (-77.2969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.703998,0.0576671,-1.52154e-05,-3.20245e-08,1.93273e-11,-9164.48,21.923], Tmin=(100,'K'), Tmax=(947.198,'K')), NASAPolynomial(coeffs=[19.7295,0.0138753,-3.75088e-06,6.4769e-10,-4.91449e-14,-14408.4,-77.4874], Tmin=(947.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.2969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(methyleneoxirane) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C1[CH]C(CO)O1(5200)',
    structure = SMILES('C=C1[CH]C(CO)O1'),
    E0 = (-134.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.807106,0.0506616,1.69163e-05,-7.5527e-08,3.82563e-11,-16051.9,18.8726], Tmin=(100,'K'), Tmax=(906.539,'K')), NASAPolynomial(coeffs=[22.313,0.00797217,1.17514e-06,-4.29412e-10,2.86792e-14,-22096.1,-94.6016], Tmin=(906.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C=CC(=C)OO(5201)',
    structure = SMILES('[CH2]C=CC(=C)OO'),
    E0 = (86.5533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,710.734],'cm^-1')),
        HinderedRotor(inertia=(0.544229,'amu*angstrom^2'), symmetry=1, barrier=(12.5129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13224,'amu*angstrom^2'), symmetry=1, barrier=(26.0323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.6581,'amu*angstrom^2'), symmetry=1, barrier=(107.099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.65676,'amu*angstrom^2'), symmetry=1, barrier=(107.068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.820249,0.0634873,-5.3216e-05,2.33015e-08,-4.11446e-12,10529.8,26.062], Tmin=(100,'K'), Tmax=(1351.85,'K')), NASAPolynomial(coeffs=[14.1885,0.0239315,-9.3246e-06,1.6561e-09,-1.11484e-13,6915.44,-42.4619], Tmin=(1351.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.5533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'O=C1CC1[CH]CO(5202)',
    structure = SMILES('O=C1CC1[CH]CO'),
    E0 = (-11.4487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08908,0.0568456,-3.89168e-05,9.59843e-09,6.59812e-13,-1265.98,26.5372], Tmin=(100,'K'), Tmax=(1041,'K')), NASAPolynomial(coeffs=[12.8519,0.0238117,-8.84567e-06,1.56559e-09,-1.06791e-13,-4374.13,-33.8503], Tmin=(1041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-11.4487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(CCJCO)"""),
)

species(
    label = 'CC(=O)[C]=CCO(5203)',
    structure = SMILES('CC([O])=C=CCO'),
    E0 = (-119.213,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.325631,'amu*angstrom^2'), symmetry=1, barrier=(7.4869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.325014,'amu*angstrom^2'), symmetry=1, barrier=(7.47272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.323296,'amu*angstrom^2'), symmetry=1, barrier=(7.43322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04684,0.0662527,-6.95082e-05,4.10439e-08,-9.95833e-12,-14232.7,25.7503], Tmin=(100,'K'), Tmax=(991.268,'K')), NASAPolynomial(coeffs=[10.52,0.0280262,-1.16635e-05,2.14112e-09,-1.46973e-13,-16110.8,-19.8694], Tmin=(991.268,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC(=O)C=[C]CO(5204)',
    structure = SMILES('CC(=O)C=[C]CO'),
    E0 = (-77.8133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3615,1277.5,1000,1685,370,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,180,3333.91],'cm^-1')),
        HinderedRotor(inertia=(0.196743,'amu*angstrom^2'), symmetry=1, barrier=(4.52351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194783,'amu*angstrom^2'), symmetry=1, barrier=(4.47845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193934,'amu*angstrom^2'), symmetry=1, barrier=(4.45894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.600699,'amu*angstrom^2'), symmetry=1, barrier=(13.8112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4356,0.0628178,-8.23745e-05,7.51756e-08,-2.90003e-11,-9272.65,26.6464], Tmin=(100,'K'), Tmax=(772.173,'K')), NASAPolynomial(coeffs=[3.61456,0.0411263,-2.00264e-05,3.89723e-09,-2.73728e-13,-9298.98,18.706], Tmin=(772.173,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.8133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + radical(Cds_S)"""),
)

species(
    label = 'S(245)(244)',
    structure = SMILES('CC([O])=CC=CO'),
    E0 = (-227.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4358.71,'J/mol'), sigma=(6.881,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=680.82 K, Pc=30.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.342212,0.0632525,-1.84618e-05,-4.28291e-08,2.77283e-11,-27232.2,24.3496], Tmin=(100,'K'), Tmax=(902.368,'K')), NASAPolynomial(coeffs=[24.4528,0.00343297,2.75282e-06,-7.11476e-10,4.87923e-14,-33499.4,-100.109], Tmin=(902.368,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-227.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC(=O)C=CC[O](5205)',
    structure = SMILES('CC(=O)C=CC[O]'),
    E0 = (-89.9499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,2013.21,2013.56],'cm^-1')),
        HinderedRotor(inertia=(0.247127,'amu*angstrom^2'), symmetry=1, barrier=(5.68193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247068,'amu*angstrom^2'), symmetry=1, barrier=(5.68058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247079,'amu*angstrom^2'), symmetry=1, barrier=(5.68083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80724,0.0558069,-7.13927e-05,7.36168e-08,-3.1659e-11,-10746.9,25.696], Tmin=(100,'K'), Tmax=(774.119,'K')), NASAPolynomial(coeffs=[0.0191784,0.0482043,-2.40274e-05,4.72192e-09,-3.33399e-13,-9965.49,37.124], Tmin=(774.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.9499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + radical(CCOJ)"""),
)

species(
    label = 'O=C1[CH]C(CO)C1(5206)',
    structure = SMILES('O=C1[CH]C(CO)C1'),
    E0 = (-146.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74782,0.0343391,3.16627e-05,-6.42669e-08,2.64523e-11,-17519.2,23.9017], Tmin=(100,'K'), Tmax=(981.687,'K')), NASAPolynomial(coeffs=[13.3984,0.0238778,-8.90418e-06,1.68644e-09,-1.2379e-13,-21590,-41.1736], Tmin=(981.687,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJC=O)"""),
)

species(
    label = 'O=[C]CC=CCO(5207)',
    structure = SMILES('O=[C]CC=CCO'),
    E0 = (-122.34,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,381.448,381.449],'cm^-1')),
        HinderedRotor(inertia=(0.00115862,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0970991,'amu*angstrom^2'), symmetry=1, barrier=(10.0256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0970982,'amu*angstrom^2'), symmetry=1, barrier=(10.0257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0970984,'amu*angstrom^2'), symmetry=1, barrier=(10.0256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38872,0.053032,-3.49004e-05,1.06921e-08,-1.29664e-12,-14616.7,27.4364], Tmin=(100,'K'), Tmax=(1891.16,'K')), NASAPolynomial(coeffs=[16.4144,0.0212515,-9.69346e-06,1.80632e-09,-1.22006e-13,-20300,-54.6281], Tmin=(1891.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'CH2(17)(18)',
    structure = SMILES('[CH2]'),
    E0 = (381.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([971.045,2816.03,3444.23],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71758,0.00127391,2.17347e-06,-3.48858e-09,1.65209e-12,45872.4,1.75298], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.14632,0.00303671,-9.96474e-07,1.50484e-10,-8.57336e-15,46041.3,4.72342], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(381.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=[C]C=CCO(5208)',
    structure = SMILES('O=C=C[CH]CO'),
    E0 = (-144.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3615,1277.5,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,395.212,395.898],'cm^-1')),
        HinderedRotor(inertia=(0.128478,'amu*angstrom^2'), symmetry=1, barrier=(14.0187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09918,'amu*angstrom^2'), symmetry=1, barrier=(25.2723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71738,'amu*angstrom^2'), symmetry=1, barrier=(71.0975,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62616,0.0533103,-5.33187e-05,2.88065e-08,-6.39437e-12,-17306.1,19.2416], Tmin=(100,'K'), Tmax=(1074.61,'K')), NASAPolynomial(coeffs=[9.97609,0.0222293,-9.93385e-06,1.89123e-09,-1.32682e-13,-19100.6,-21.6427], Tmin=(1074.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-144.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH) + radical(C=CCJCO)"""),
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
    E0 = (55.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (128.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (266.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (237.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (134.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (278.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (263.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (243.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (288.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-1.73712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (31.8394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (130.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-26.9827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-70.9123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-10.8273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (15.2076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-48.3187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (154.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (324.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-11.4487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (32.6662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-33.5047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-93.2967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (3.56248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-48.3187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (123.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (236.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CH2CO(26)(27)', '[CH]=CCO(4313)'],
    products = ['S(244)(243)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.0131003,'m^3/(mol*s)'), n=2.40999, Ea=(12.7705,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CdsJ-H] for rate rule [Ck_O;CdsJ-H]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction2',
    reactants = ['OH(5)(5)', '[CH2]C=CC(=C)[O](5186)'],
    products = ['S(244)(243)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.7e+13,'cm^3/(mol*s)','+|-',1e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(400,'K'), comment="""Estimated using template [O_pri_rad;C_pri_rad] for rate rule [O_pri_rad;C_rad/H2/Cd]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)(3)', 'C=C([O])C=CC[O](5187)'],
    products = ['S(244)(243)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2OH(22)(23)', '[CH]=CC(=C)[O](5188)'],
    products = ['S(244)(243)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;Cd_rad] for rate rule [C_rad/H2/O;Cd_pri_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)(3)', 'C=C([O])C=C[CH]O(5189)'],
    products = ['S(244)(243)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;H_rad] for rate rule [C_rad/H/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)(3)', 'C=C([O])C=[C]CO(5190)'],
    products = ['S(244)(243)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][C]=O(1196)', '[CH]=CCO(4313)'],
    products = ['S(244)(243)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.88428e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)(3)', 'C=C([O])[C]=CCO(5191)'],
    products = ['S(244)(243)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)(3)', '[CH]=C([O])C=CCO(5192)'],
    products = ['S(244)(243)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['S(244)(243)'],
    products = ['C=C1OC1[CH]CO(5193)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(169.591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R4_S_(Cd)_D;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 167.5 to 169.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C(O)[C]=CCO(5194)'],
    products = ['S(244)(243)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_2Cd;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(O)C=CCO(5195)'],
    products = ['S(244)(243)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C(O)C=[C]CO(5196)'],
    products = ['S(244)(243)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['S(244)(243)'],
    products = ['C=C(O)C=C[CH]O(5197)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1'), n=1.64, Ea=(100.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SS(D)MS;Y_rad_out;Cs_H_out] for rate rule [R5H_SS(D)MS;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C(O)C=CC[O](5198)'],
    products = ['S(244)(243)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.33753e+06,'s^-1'), n=1.02312, Ea=(72.6006,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;O_rad_out;XH_out] for rate rule [R6H_RSMSR;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['S(244)(243)'],
    products = ['OCC=C[C]1CO1(5199)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(186.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['S(244)(243)'],
    products = ['C=C1[CH]C(CO)O1(5200)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.61e+08,'s^-1'), n=0.96, Ea=(123.01,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [R4_S_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C=CC(=C)OO(5201)'],
    products = ['S(244)(243)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.12e+10,'s^-1'), n=0, Ea=(68.1992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4OOH;C_rad_out_2H] for rate rule [R4OOH_SDS;C_rad_out_2H]
Euclidian distance = 1.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O(4)(4)', 'C5H7O(237)(236)'],
    products = ['S(244)(243)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['S(244)(243)'],
    products = ['O=C1CC1[CH]CO(5202)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(159.88,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd_pri;radadd_intra_cs2H] for rate rule [R4_S_(CO)_D;doublebond_intra_HNd_pri;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 158.5 to 159.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC(=O)[C]=CCO(5203)'],
    products = ['S(244)(243)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC(=O)C=[C]CO(5204)'],
    products = ['S(244)(243)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['S(244)(243)'],
    products = ['S(245)(244)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(423178,'s^-1'), n=1.77, Ea=(78.0316,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;C_rad_out_2H;Cs_H_out] for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CC(=O)C=CC[O](5205)'],
    products = ['S(244)(243)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6H;O_rad_out;Cs_H_out_2H] for rate rule [R6H_RSMSR;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['S(244)(243)'],
    products = ['O=C1[CH]C(CO)C1(5206)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.61e+08,'s^-1'), n=0.96, Ea=(123.01,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=[C]CC=CCO(5207)'],
    products = ['S(244)(243)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(17)(18)', 'O=[C]C=CCO(5208)'],
    products = ['S(244)(243)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

network(
    label = '355',
    isomers = [
        'S(244)(243)',
    ],
    reactants = [
        ('O(4)(4)', 'C5H7O(237)(236)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '355',
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

