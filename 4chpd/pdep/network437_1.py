species(
    label = 'C=C[CH]C(C)C[O](2184)',
    structure = SMILES('[CH2]C=CC(C)C[O]'),
    E0 = (133.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,311.507,2229.83,2229.95],'cm^-1')),
        HinderedRotor(inertia=(0.46564,'amu*angstrom^2'), symmetry=1, barrier=(32.0355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149495,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149551,'amu*angstrom^2'), symmetry=1, barrier=(10.2892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.46554,'amu*angstrom^2'), symmetry=1, barrier=(32.0369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08931,0.0612644,-3.65899e-05,1.07808e-08,-1.30255e-12,16111.9,28.178], Tmin=(100,'K'), Tmax=(1824.38,'K')), NASAPolynomial(coeffs=[13.7084,0.0335965,-1.38414e-05,2.46799e-09,-1.63417e-13,11507.5,-40.2889], Tmin=(1824.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCOJ)"""),
)

species(
    label = 'CH2O(35)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=CC=CC(381)',
    structure = SMILES('C=CC=CC'),
    E0 = (57.8956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.831076,'amu*angstrom^2'), symmetry=1, barrier=(19.1081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833175,'amu*angstrom^2'), symmetry=1, barrier=(19.1563,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3140.68,'J/mol'), sigma=(5.4037,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=490.57 K, Pc=45.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00727,0.0328459,1.55855e-05,-4.25745e-08,1.84259e-11,7044.82,16.9534], Tmin=(100,'K'), Tmax=(972.32,'K')), NASAPolynomial(coeffs=[11.2869,0.0212416,-7.50361e-06,1.3618e-09,-9.72233e-14,3984.25,-34.0139], Tmin=(972.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.8956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][CH]C1OCC1C(3270)',
    structure = SMILES('[CH2][CH]C1OCC1C'),
    E0 = (215.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44859,0.0402828,4.11982e-05,-8.636e-08,3.8665e-11,26061.9,26.3309], Tmin=(100,'K'), Tmax=(899.402,'K')), NASAPolynomial(coeffs=[14.0187,0.0265235,-6.14276e-06,8.3075e-10,-5.26455e-14,22096.2,-42.4563], Tmin=(899.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = 'H(19)',
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
    label = '[CH2]C=CC(C)C=O(3271)',
    structure = SMILES('[CH2]C=CC(C)C=O'),
    E0 = (-14.5231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,1897],'cm^-1')),
        HinderedRotor(inertia=(0.433009,'amu*angstrom^2'), symmetry=1, barrier=(9.95573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777415,'amu*angstrom^2'), symmetry=1, barrier=(17.8743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.432934,'amu*angstrom^2'), symmetry=1, barrier=(9.954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36084,'amu*angstrom^2'), symmetry=1, barrier=(31.2885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10062,0.060669,-4.01461e-05,1.35808e-08,-1.90562e-12,-1640.28,26.0429], Tmin=(100,'K'), Tmax=(1601.76,'K')), NASAPolynomial(coeffs=[12.574,0.0320175,-1.33152e-05,2.41374e-09,-1.62695e-13,-5315.84,-34.7146], Tmin=(1601.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.5231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P)"""),
)

species(
    label = 'C=C=CC(C)C[O](3272)',
    structure = SMILES('C=C=CC(C)C[O]'),
    E0 = (158.182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,256.725,256.942,257.085],'cm^-1')),
        HinderedRotor(inertia=(0.00288235,'amu*angstrom^2'), symmetry=1, barrier=(10.667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513999,'amu*angstrom^2'), symmetry=1, barrier=(24.0226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227677,'amu*angstrom^2'), symmetry=1, barrier=(10.6659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53978,0.0587236,-3.87048e-05,1.38223e-08,-2.20672e-12,19109.2,25.3255], Tmin=(100,'K'), Tmax=(1320.7,'K')), NASAPolynomial(coeffs=[7.53081,0.0405784,-1.8096e-05,3.4192e-09,-2.37456e-13,17526.8,-5.24415], Tmin=(1320.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C=C[CH]C(377)',
    structure = SMILES('[CH2]C=C[CH]C'),
    E0 = (240.064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.0180055,'amu*angstrom^2'), symmetry=1, barrier=(19.7234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34503,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0180001,'amu*angstrom^2'), symmetry=1, barrier=(19.7225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18178,0.0283568,2.70949e-05,-5.14684e-08,2.05693e-11,28948.9,17.5848], Tmin=(100,'K'), Tmax=(990.212,'K')), NASAPolynomial(coeffs=[10.2369,0.0240425,-9.12514e-06,1.70243e-09,-1.22294e-13,25969.9,-28.1844], Tmin=(990.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C=CC(C)[CH]O(3273)',
    structure = SMILES('[CH2]C=CC(C)[CH]O'),
    E0 = (87.6701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.314257,0.0744874,-6.20061e-05,2.73678e-08,-4.90163e-12,10682.4,29.6483], Tmin=(100,'K'), Tmax=(1328.83,'K')), NASAPolynomial(coeffs=[15.1684,0.029774,-1.15332e-05,2.0459e-09,-1.37698e-13,6734.63,-46.2375], Tmin=(1328.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.6701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCsJOH)"""),
)

species(
    label = 'C[C]=CC(C)C[O](3263)',
    structure = SMILES('C[C]=CC(C)C[O]'),
    E0 = (219.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,224.614,224.91,2188.5],'cm^-1')),
        HinderedRotor(inertia=(0.00333599,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238846,'amu*angstrom^2'), symmetry=1, barrier=(8.58033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239262,'amu*angstrom^2'), symmetry=1, barrier=(8.57875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239146,'amu*angstrom^2'), symmetry=1, barrier=(8.57808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22547,0.0492503,3.52168e-05,-1.84294e-07,1.72565e-10,26443.1,24.5755], Tmin=(100,'K'), Tmax=(409.264,'K')), NASAPolynomial(coeffs=[3.84348,0.0492402,-2.26694e-05,4.35167e-09,-3.05519e-13,26178.3,16.5981], Tmin=(409.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C[C](C)CO(3274)',
    structure = SMILES('[CH2]C=C[C](C)CO'),
    E0 = (59.9461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.801747,0.05976,-2.197e-05,-1.09881e-08,7.63114e-12,7334.09,28.5434], Tmin=(100,'K'), Tmax=(1045.26,'K')), NASAPolynomial(coeffs=[13.5275,0.0321336,-1.25651e-05,2.30089e-09,-1.60356e-13,3522.61,-38.9207], Tmin=(1045.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.9461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCJ(C)CO)"""),
)

species(
    label = 'CC=[C]C(C)C[O](3260)',
    structure = SMILES('CC=[C]C(C)C[O]'),
    E0 = (219.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,224.614,224.91,2188.5],'cm^-1')),
        HinderedRotor(inertia=(0.00333599,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238846,'amu*angstrom^2'), symmetry=1, barrier=(8.58033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239262,'amu*angstrom^2'), symmetry=1, barrier=(8.57875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239146,'amu*angstrom^2'), symmetry=1, barrier=(8.57808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22547,0.0492503,3.52168e-05,-1.84294e-07,1.72565e-10,26443.1,24.5755], Tmin=(100,'K'), Tmax=(409.264,'K')), NASAPolynomial(coeffs=[3.84348,0.0492402,-2.26694e-05,4.35167e-09,-3.05519e-13,26178.3,16.5981], Tmin=(409.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC([CH2])CO(3265)',
    structure = SMILES('[CH2]C=CC([CH2])CO'),
    E0 = (112.455,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180,3733.43],'cm^-1')),
        HinderedRotor(inertia=(0.0209959,'amu*angstrom^2'), symmetry=1, barrier=(15.6613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290327,'amu*angstrom^2'), symmetry=1, barrier=(85.403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681178,'amu*angstrom^2'), symmetry=1, barrier=(15.6616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229051,'amu*angstrom^2'), symmetry=1, barrier=(5.26634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.71577,'amu*angstrom^2'), symmetry=1, barrier=(85.4328,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.542244,0.0681927,-4.63e-05,1.20506e-08,4.82621e-13,13656.4,31.0593], Tmin=(100,'K'), Tmax=(1029.65,'K')), NASAPolynomial(coeffs=[13.6419,0.0306952,-1.11835e-05,1.94596e-09,-1.31075e-13,10248.9,-35.9691], Tmin=(1029.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=[C]C(C)CO(3275)',
    structure = SMILES('[CH2]C=[C]C(C)CO'),
    E0 = (145.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0416757,'amu*angstrom^2'), symmetry=1, barrier=(9.53769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.414818,'amu*angstrom^2'), symmetry=1, barrier=(9.53748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0172897,'amu*angstrom^2'), symmetry=1, barrier=(3.95682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0416766,'amu*angstrom^2'), symmetry=1, barrier=(9.53759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109806,'amu*angstrom^2'), symmetry=1, barrier=(25.1289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.59431,0.069848,-5.38944e-05,2.219e-08,-3.76537e-12,17591.8,29.5665], Tmin=(100,'K'), Tmax=(1378,'K')), NASAPolynomial(coeffs=[13.6144,0.0320537,-1.27538e-05,2.28635e-09,-1.54401e-13,14003.5,-37.4224], Tmin=(1378,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=C[C](C)C[O](3257)',
    structure = SMILES('C[CH]C=C(C)C[O]'),
    E0 = (122.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,180,1026.18,1026.69],'cm^-1')),
        HinderedRotor(inertia=(0.154165,'amu*angstrom^2'), symmetry=1, barrier=(3.54455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154394,'amu*angstrom^2'), symmetry=1, barrier=(3.54983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154069,'amu*angstrom^2'), symmetry=1, barrier=(3.54236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154133,'amu*angstrom^2'), symmetry=1, barrier=(3.54383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88619,0.0543366,-2.68033e-05,5.76772e-09,-4.71867e-13,14750.4,23.7058], Tmin=(100,'K'), Tmax=(2860.08,'K')), NASAPolynomial(coeffs=[28.9421,0.0165001,-6.96113e-06,1.14299e-09,-6.76495e-14,-727.167,-135.257], Tmin=(2860.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C]=CC(C)CO(3276)',
    structure = SMILES('[CH2][C]=CC(C)CO'),
    E0 = (145.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0416757,'amu*angstrom^2'), symmetry=1, barrier=(9.53769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.414818,'amu*angstrom^2'), symmetry=1, barrier=(9.53748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0172897,'amu*angstrom^2'), symmetry=1, barrier=(3.95682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0416766,'amu*angstrom^2'), symmetry=1, barrier=(9.53759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109806,'amu*angstrom^2'), symmetry=1, barrier=(25.1289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.59431,0.069848,-5.38944e-05,2.219e-08,-3.76537e-12,17591.8,29.5665], Tmin=(100,'K'), Tmax=(1378,'K')), NASAPolynomial(coeffs=[13.6144,0.0320537,-1.27538e-05,2.28635e-09,-1.54401e-13,14003.5,-37.4224], Tmin=(1378,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=CC)C[O](2186)',
    structure = SMILES('[CH2]C(C=CC)C[O]'),
    E0 = (186.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,1198.17,1198.42],'cm^-1')),
        HinderedRotor(inertia=(0.249587,'amu*angstrom^2'), symmetry=1, barrier=(5.73849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250246,'amu*angstrom^2'), symmetry=1, barrier=(5.75365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249709,'amu*angstrom^2'), symmetry=1, barrier=(5.7413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249535,'amu*angstrom^2'), symmetry=1, barrier=(5.7373,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3875.05,'J/mol'), sigma=(6.66255,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.27 K, Pc=29.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32951,0.0627137,-4.78961e-05,2.42306e-08,-5.88253e-12,22542.8,28.8295], Tmin=(100,'K'), Tmax=(911.733,'K')), NASAPolynomial(coeffs=[5.2164,0.0456614,-1.98423e-05,3.718e-09,-2.58069e-13,21834,10.4366], Tmin=(911.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(CCOJ)"""),
)

species(
    label = 'CC=CC(C)[CH][O](3259)',
    structure = SMILES('CC=CC(C)[CH][O]'),
    E0 = (161.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,337.61,337.61,1757.94],'cm^-1')),
        HinderedRotor(inertia=(0.109695,'amu*angstrom^2'), symmetry=1, barrier=(8.87244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109695,'amu*angstrom^2'), symmetry=1, barrier=(8.87244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109695,'amu*angstrom^2'), symmetry=1, barrier=(8.87244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109695,'amu*angstrom^2'), symmetry=1, barrier=(8.87244,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0986,0.0690427,-6.38247e-05,4.02254e-08,-1.18288e-11,19569,27.43], Tmin=(100,'K'), Tmax=(782.393,'K')), NASAPolynomial(coeffs=[5.42286,0.0469349,-2.14399e-05,4.10997e-09,-2.88759e-13,18892.3,7.62909], Tmin=(782.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2][O](354)',
    structure = SMILES('[CH2][O]'),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88409,-0.00363885,3.28543e-05,-4.13611e-08,1.59631e-11,23210.8,7.47983], Tmin=(100,'K'), Tmax=(933.06,'K')), NASAPolynomial(coeffs=[6.69335,0.000289989,8.61416e-07,-1.56351e-10,7.33778e-15,21991.3,-9.6043], Tmin=(933.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(H3COJ) + radical(CsJOH)"""),
)

species(
    label = '[CH]=C[CH2](321)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,229.711,230.18,230.787],'cm^-1')),
        HinderedRotor(inertia=(1.33306,'amu*angstrom^2'), symmetry=1, barrier=(50.5153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,44095,-3.44606], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C[CH]C[O](3030)',
    structure = SMILES('C[CH]C[O]'),
    E0 = (152.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1353.75,2182.69],'cm^-1')),
        HinderedRotor(inertia=(0.00393585,'amu*angstrom^2'), symmetry=1, barrier=(5.11859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222658,'amu*angstrom^2'), symmetry=1, barrier=(5.11934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.10998,0.016427,1.17435e-05,-1.70304e-08,5.13355e-12,18403.8,17.5215], Tmin=(100,'K'), Tmax=(1242.45,'K')), NASAPolynomial(coeffs=[4.47923,0.0232159,-9.97094e-06,1.87455e-09,-1.29964e-13,17199.3,7.14046], Tmin=(1242.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = 'CC(C[O])C1[CH]C1(3277)',
    structure = SMILES('CC(C[O])C1[CH]C1'),
    E0 = (244.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32621,0.0508077,-1.2977e-05,-8.48881e-09,4.07928e-12,29502,28.6119], Tmin=(100,'K'), Tmax=(1241.4,'K')), NASAPolynomial(coeffs=[9.9593,0.0381347,-1.5963e-05,2.94191e-09,-2.0176e-13,26191.6,-19.6045], Tmin=(1241.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CCOJ) + radical(cyclopropane)"""),
)

species(
    label = '[CH2]C1[CH]C(C)CO1(3176)',
    structure = SMILES('[CH2]C1[CH]C(C)CO1'),
    E0 = (142.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40159,0.0405915,4.31251e-05,-9.05795e-08,4.08234e-11,17260.6,22.2657], Tmin=(100,'K'), Tmax=(896.051,'K')), NASAPolynomial(coeffs=[14.7971,0.0251418,-5.25143e-06,6.47565e-10,-3.97318e-14,13079.6,-50.8243], Tmin=(896.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CJC(C)OC) + radical(CCJCO)"""),
)

species(
    label = 'CC=CC(C)C=O(3268)',
    structure = SMILES('CC=CC(C)C=O'),
    E0 = (-166.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48269,0.0589436,-3.5055e-05,1.05093e-08,-1.3413e-12,-19881.3,24.37], Tmin=(100,'K'), Tmax=(1637.02,'K')), NASAPolynomial(coeffs=[9.59879,0.0391121,-1.68833e-05,3.10895e-09,-2.11144e-13,-22538.6,-18.7858], Tmin=(1637.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'C=C=CC(C)CO(3278)',
    structure = SMILES('C=C=CC(C)CO'),
    E0 = (-67.5227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.543132,0.068764,-5.1161e-05,2.00048e-08,-3.19885e-12,-7990.76,28.0549], Tmin=(100,'K'), Tmax=(1463.99,'K')), NASAPolynomial(coeffs=[14.5946,0.0303714,-1.18238e-05,2.09143e-09,-1.39832e-13,-12105,-45.0912], Tmin=(1463.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.5227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CH2(S)(40)',
    structure = SMILES('[CH2]'),
    E0 = (418.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1358.21,2621.43,3089.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19331,-0.00233105,8.15676e-06,-6.62986e-09,1.93233e-12,50366.2,-0.746734], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.13502,0.00289594,-8.16668e-07,1.13573e-10,-6.36263e-15,50504.1,4.06031], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(418.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=C[CH]CC[O](3242)',
    structure = SMILES('[CH2]C=CCC[O]'),
    E0 = (164.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,245.422,245.638,2444.53],'cm^-1')),
        HinderedRotor(inertia=(0.00245214,'amu*angstrom^2'), symmetry=1, barrier=(10.4,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766442,'amu*angstrom^2'), symmetry=1, barrier=(32.7797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766241,'amu*angstrom^2'), symmetry=1, barrier=(32.778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98027,0.0464601,-2.5237e-05,6.38765e-09,-6.46696e-13,19894.2,22.9639], Tmin=(100,'K'), Tmax=(2131.54,'K')), NASAPolynomial(coeffs=[12.6767,0.0263877,-1.11119e-05,1.9699e-09,-1.28561e-13,15334.1,-36.7357], Tmin=(2131.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(164.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P)"""),
)

species(
    label = 'CC1C=CCOC1(2192)',
    structure = SMILES('CC1C=CCOC1'),
    E0 = (-151.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74582,0.0312494,6.34028e-05,-1.0244e-07,4.1495e-11,-18183.3,18.8564], Tmin=(100,'K'), Tmax=(940.121,'K')), NASAPolynomial(coeffs=[12.8037,0.0305353,-9.38671e-06,1.60239e-09,-1.13271e-13,-22310,-44.6988], Tmin=(940.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-151.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro2hpyran)"""),
)

species(
    label = 'O(20)',
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
    label = '[CH2]C=CC([CH2])C(1156)',
    structure = SMILES('[CH2]C=CC([CH2])C'),
    E0 = (272.282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.00514534,'amu*angstrom^2'), symmetry=1, barrier=(4.602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.674437,'amu*angstrom^2'), symmetry=1, barrier=(15.5066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0173543,'amu*angstrom^2'), symmetry=1, barrier=(15.7571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.3039,'amu*angstrom^2'), symmetry=1, barrier=(75.9633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32272,0.0477992,1.34891e-06,-3.42469e-08,1.64572e-11,32854.3,25.5548], Tmin=(100,'K'), Tmax=(967.428,'K')), NASAPolynomial(coeffs=[12.2803,0.0288316,-1.00792e-05,1.77004e-09,-1.22479e-13,29501.6,-33.3163], Tmin=(967.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Isobutyl)"""),
)

species(
    label = 'CH2(T)(33)',
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
    label = '[CH]=CC(C)C[O](3279)',
    structure = SMILES('[CH]=CC(C)C[O]'),
    E0 = (264.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,313.75,313.763],'cm^-1')),
        HinderedRotor(inertia=(0.00896897,'amu*angstrom^2'), symmetry=1, barrier=(6.88984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0986228,'amu*angstrom^2'), symmetry=1, barrier=(6.8899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0986334,'amu*angstrom^2'), symmetry=1, barrier=(6.89007,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99461,0.0480248,-2.95022e-05,9.49129e-09,-1.3458e-12,31904.3,23.2306], Tmin=(100,'K'), Tmax=(1453.97,'K')), NASAPolynomial(coeffs=[7.15667,0.0338237,-1.48516e-05,2.77384e-09,-1.90794e-13,30403.2,-3.60556], Tmin=(1453.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC=C(C)C[O](3280)',
    structure = SMILES('C=CC=C(C)C[O]'),
    E0 = (91.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,238.949,239.006,239.187],'cm^-1')),
        HinderedRotor(inertia=(0.00295398,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231156,'amu*angstrom^2'), symmetry=1, barrier=(9.37491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231599,'amu*angstrom^2'), symmetry=1, barrier=(9.37432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38968,0.0598288,-4.11097e-05,1.54128e-08,-2.52938e-12,11085.7,24.5377], Tmin=(100,'K'), Tmax=(1330.94,'K')), NASAPolynomial(coeffs=[8.65964,0.0379795,-1.64849e-05,3.07821e-09,-2.12458e-13,9150.5,-12.614], Tmin=(1330.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = 'CH3(34)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=CC=CC[O](451)',
    structure = SMILES('C=CC=CC[O]'),
    E0 = (130.465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,324.17,324.798,326.428],'cm^-1')),
        HinderedRotor(inertia=(0.0834457,'amu*angstrom^2'), symmetry=1, barrier=(6.17438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147528,'amu*angstrom^2'), symmetry=1, barrier=(11.0638,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80098,0.0446979,-2.65813e-05,7.81281e-09,-9.36109e-13,15773.2,22.3196], Tmin=(100,'K'), Tmax=(1862.48,'K')), NASAPolynomial(coeffs=[11.7735,0.0232802,-9.33191e-06,1.63842e-09,-1.07321e-13,12058.5,-31.9938], Tmin=(1862.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = 'C=CC[C](C)C[O](3281)',
    structure = SMILES('C=CC[C](C)C[O]'),
    E0 = (145.328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,360,370,350,311.012,311.013,2590.79,2590.79],'cm^-1')),
        HinderedRotor(inertia=(0.00174279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0810648,'amu*angstrom^2'), symmetry=1, barrier=(5.56436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.562706,'amu*angstrom^2'), symmetry=1, barrier=(38.6247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.562708,'amu*angstrom^2'), symmetry=1, barrier=(38.6247,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26679,0.0574013,-3.02554e-05,7.29651e-09,-6.86864e-13,17579.5,27.8859], Tmin=(100,'K'), Tmax=(2418.73,'K')), NASAPolynomial(coeffs=[21.3339,0.0242134,-9.67258e-06,1.62307e-09,-1.00428e-13,7872.57,-86.6497], Tmin=(2418.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJ(C)CO) + radical(CCOJ)"""),
)

species(
    label = 'C=[C]CC(C)C[O](3282)',
    structure = SMILES('C=[C]CC(C)C[O]'),
    E0 = (230.596,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,336.859,337.196,337.266,1808.63],'cm^-1')),
        HinderedRotor(inertia=(0.0014861,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985831,'amu*angstrom^2'), symmetry=1, barrier=(7.94215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0986145,'amu*angstrom^2'), symmetry=1, barrier=(7.93714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0980394,'amu*angstrom^2'), symmetry=1, barrier=(7.9393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28558,0.0649591,-5.41221e-05,3.15095e-08,-9.01792e-12,27827.3,28.0893], Tmin=(100,'K'), Tmax=(781.025,'K')), NASAPolynomial(coeffs=[4.52858,0.0483499,-2.22225e-05,4.28011e-09,-3.01821e-13,27320.8,13.2453], Tmin=(781.025,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(C[O])CC=C(3283)',
    structure = SMILES('[CH2]C(C[O])CC=C'),
    E0 = (197.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,393.462,393.542,393.693,2040.13],'cm^-1')),
        HinderedRotor(inertia=(0.0705387,'amu*angstrom^2'), symmetry=1, barrier=(7.75979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0704948,'amu*angstrom^2'), symmetry=1, barrier=(7.75925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0705111,'amu*angstrom^2'), symmetry=1, barrier=(7.75952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0230337,'amu*angstrom^2'), symmetry=1, barrier=(68.0198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31297,0.0625336,-4.44251e-05,1.90493e-08,-3.77806e-12,23888,29.2841], Tmin=(100,'K'), Tmax=(1099.58,'K')), NASAPolynomial(coeffs=[6.48314,0.0437259,-1.87687e-05,3.49416e-09,-2.4147e-13,22751,3.85022], Tmin=(1099.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCOJ)"""),
)

species(
    label = 'C=CCC(C)[CH][O](3284)',
    structure = SMILES('C=CCC(C)[CH][O]'),
    E0 = (173.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,332.219,332.761,332.972,3623.15],'cm^-1')),
        HinderedRotor(inertia=(0.00150972,'amu*angstrom^2'), symmetry=1, barrier=(0.119628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120073,'amu*angstrom^2'), symmetry=1, barrier=(9.39972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.537532,'amu*angstrom^2'), symmetry=1, barrier=(42.4009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119907,'amu*angstrom^2'), symmetry=1, barrier=(9.4089,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14176,0.0680456,-5.7002e-05,3.00331e-08,-7.29934e-12,20911.9,27.6786], Tmin=(100,'K'), Tmax=(928.162,'K')), NASAPolynomial(coeffs=[6.35139,0.0455938,-2.07171e-05,3.97039e-09,-2.79235e-13,19944.8,2.93356], Tmin=(928.162,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=CCC(C)C[O](3285)',
    structure = SMILES('[CH]=CCC(C)C[O]'),
    E0 = (239.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3961,0.0624044,-4.20573e-05,1.58099e-08,-2.70504e-12,28936.5,27.6], Tmin=(100,'K'), Tmax=(1234.82,'K')), NASAPolynomial(coeffs=[7.07989,0.0439927,-1.96917e-05,3.73497e-09,-2.60358e-13,27532.8,-1.01987], Tmin=(1234.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C[CH]C(C)CO(3286)',
    structure = SMILES('[CH]C=CC(C)CO'),
    E0 = (126.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.370183,0.0704546,-4.59244e-05,1.54712e-08,-2.1445e-12,15359.7,30.7422], Tmin=(100,'K'), Tmax=(1652.45,'K')), NASAPolynomial(coeffs=[15.2518,0.0344311,-1.32241e-05,2.27846e-09,-1.48549e-13,10441.5,-48.5273], Tmin=(1652.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'CC1[CH][CH]COC1(3193)',
    structure = SMILES('CC1[CH][CH]COC1'),
    E0 = (117.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99674,0.0265344,7.33184e-05,-1.1139e-07,4.50102e-11,14215.5,23.0214], Tmin=(100,'K'), Tmax=(919.833,'K')), NASAPolynomial(coeffs=[11.0657,0.0320621,-9.02144e-06,1.43193e-09,-9.73753e-14,10644.9,-30.3133], Tmin=(919.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxane) + radical(Cs_S) + radical(CCJCO)"""),
)

species(
    label = 'C=CC=C(C)CO(3287)',
    structure = SMILES('C=CC=C(C)CO'),
    E0 = (-134.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.383288,0.0699555,-5.37655e-05,2.17568e-08,-3.56359e-12,-16013.8,27.3044], Tmin=(100,'K'), Tmax=(1449.07,'K')), NASAPolynomial(coeffs=[15.5778,0.0280126,-1.03483e-05,1.78199e-09,-1.17424e-13,-20417.4,-51.6362], Tmin=(1449.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CCC(C)C=O(3288)',
    structure = SMILES('C=CCC(C)C=O'),
    E0 = (-152.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.858167,0.0620306,-3.84701e-05,1.18944e-08,-1.49581e-12,-18239.8,27.7582], Tmin=(100,'K'), Tmax=(1803.09,'K')), NASAPolynomial(coeffs=[15.2148,0.0301816,-1.19746e-05,2.098e-09,-1.37522e-13,-23417,-49.9673], Tmin=(1803.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(C)[CH]C[O](3289)',
    structure = SMILES('C=CC(C)[CH]C[O]'),
    E0 = (193.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17856,0.0589487,-3.38279e-05,9.42266e-09,-1.06476e-12,23403.3,30.5171], Tmin=(100,'K'), Tmax=(1953.03,'K')), NASAPolynomial(coeffs=[14.525,0.0316138,-1.28336e-05,2.25624e-09,-1.47412e-13,18190.1,-42.8058], Tmin=(1953.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = 'C=CC([CH]C)C[O](2183)',
    structure = SMILES('C=CC([CH]C)C[O]'),
    E0 = (188.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,253.191,253.193,1395.2,1395.2],'cm^-1')),
        HinderedRotor(inertia=(0.00262966,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00262966,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134836,'amu*angstrom^2'), symmetry=1, barrier=(6.13384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134835,'amu*angstrom^2'), symmetry=1, barrier=(6.13383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93691,0.0547272,-2.78893e-05,6.27851e-09,-5.45786e-13,22718.8,27.7923], Tmin=(100,'K'), Tmax=(2543.57,'K')), NASAPolynomial(coeffs=[19.3203,0.0273905,-1.17683e-05,2.05324e-09,-1.305e-13,13875.5,-72.301], Tmin=(2543.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cs_S)"""),
)

species(
    label = 'C=CC1OCC1C(2194)',
    structure = SMILES('C=CC1OCC1C'),
    E0 = (-61.8987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56083,0.0339257,6.43495e-05,-1.12762e-07,4.83697e-11,-7338.46,23.0274], Tmin=(100,'K'), Tmax=(908.33,'K')), NASAPolynomial(coeffs=[15.4294,0.0244035,-5.0558e-06,6.59326e-10,-4.38815e-14,-11984.6,-54.2534], Tmin=(908.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.8987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane)"""),
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
    E0 = (133.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (255.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (214.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (381.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (158.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (222.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (340.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (215.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (381.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (196.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (246.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (242.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (189.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (242.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (251.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (438.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (529.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (359.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (190.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (158.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (158.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (583.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (140.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (515.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (645.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (308.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (283.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (252.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (271.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (430.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (317.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (297.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (401.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (261.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (193.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (211.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (222.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (315.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (348.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (140.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['CH2O(35)', 'C=CC=CC(381)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['[CH2][CH]C1OCC1C(3270)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', '[CH2]C=CC(C)C=O(3271)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.6e+09,'cm^3/(mol*s)'), n=0.935, Ea=(17.4473,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2782 used for CO-CsH_O;HJ
Exact match found for rate rule [CO-CsH_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'C=C=CC(C)C[O](3272)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2O(35)', '[CH2]C=C[CH]C(377)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.000319696,'m^3/(mol*s)'), n=2.72689, Ea=(37.1073,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_R;CsJ-CdCsH] + [CO-HH_O;CsJ-OneDe] for rate rule [CO-HH_O;CsJ-CdCsH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['[CH2]C=CC(C)[CH]O(3273)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(153000,'s^-1'), n=2.26, Ea=(88.9937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 325 used for R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['C[C]=CC(C)C[O](3263)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['[CH2]C=C[C](C)CO(3274)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7e-08,'s^-1'), n=6.3, Ea=(82.634,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC=[C]C(C)C[O](3260)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C=CC([CH2])CO(3265)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C=[C]C(C)CO(3275)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['CC=C[C](C)C[O](3257)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.86e+10,'s^-1'), n=0.58, Ea=(109.579,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H;C_rad_out_2H;Cs_H_out_Cs2] for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]=CC(C)CO(3276)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(C=CC)C[O](2186)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(121000,'s^-1'), n=1.9, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 92 used for R5H_SSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CC=CC(C)[CH][O](3259)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.45388e+06,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_SSMS;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C=C[CH]C(377)', '[CH2][O](354)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.71464e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [C_rad/H/CdCs;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C[CH2](321)', 'C[CH]C[O](3030)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.79546e+07,'m^3/(mol*s)'), n=-0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/NonDeC] for rate rule [Cd_rad;C_rad/H/NonDeC]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['CC(C[O])C1[CH]C1(3277)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['[CH2]C1[CH]C(C)CO1(3176)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(57.3526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['CC=CC(C)C=O(3268)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['C=C=CC(C)CO(3278)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(S)(40)', 'C=C[CH]CC[O](3242)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(873476,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H2/OneDeC]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['CC1C=CCOC1(2192)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(20)', '[CH2]C=CC([CH2])C(1156)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/Cs;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(T)(33)', '[CH]=CC(C)C[O](3279)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(19)', 'C=CC=C(C)C[O](3280)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(4.93712,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2570 used for Cds-CsCs_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH3(34)', 'C=CC=CC[O](451)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(11300,'cm^3/(mol*s)'), n=2.41, Ea=(16.736,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 471 used for Cds-CsH_Cds-CdH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-CdH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=CC=CC(381)', '[CH2][O](354)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.00258795,'m^3/(mol*s)'), n=2.41374, Ea=(2.11244,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-CdH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=CC[C](C)C[O](3281)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.94e+08,'s^-1'), n=1.27, Ea=(125.938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 162 used for R2H_S;C_rad_out_Cs2;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_Cs2;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]CC(C)C[O](3282)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C[O])CC=C(3283)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(25000,'s^-1'), n=2.28, Ea=(119.244,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 85 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=CCC(C)[CH][O](3284)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['[CH]=CCC(C)C[O](3285)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(8.32e+10,'s^-1'), n=0.77, Ea=(268.194,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 195 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['[CH]=C[CH]C(C)CO(3286)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.427,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cd_H_out_singleH] for rate rule [R6HJ_3;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['CC1[CH][CH]COC1(3193)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(60.8864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['C=CC=C(C)CO(3287)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['C=CCC(C)C=O(3288)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['C=CC(C)[CH]C[O](3289)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(6.59e+09,'s^-1'), n=0.99, Ea=(182.004,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for cCs(-HC)CJ;CsJ-CdH;CH3
Exact match found for rate rule [cCs(-HC)CJ;CsJ-CdH;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=CC([CH]C)C[O](2183)'],
    products = ['C=C[CH]C(C)C[O](2184)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C[CH]C(C)C[O](2184)'],
    products = ['C=CC1OCC1C(2194)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '437',
    isomers = [
        'C=C[CH]C(C)C[O](2184)',
    ],
    reactants = [
        ('CH2O(35)', 'C=CC=CC(381)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '437',
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

