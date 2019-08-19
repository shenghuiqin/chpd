species(
    label = '[O]CO[CH]C(OO)C1CCCC1(2208)',
    structure = SMILES('[O]CO[CH]C(OO)C1CCCC1'),
    E0 = (-172.426,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10841,0.113782,-5.75702e-05,-9.25569e-09,1.19734e-11,-20500.3,46.7377], Tmin=(100,'K'), Tmax=(1036.49,'K')), NASAPolynomial(coeffs=[26.0794,0.0457987,-1.82278e-05,3.416e-09,-2.42966e-13,-28535,-100.834], Tmin=(1036.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(OCOJ) + radical(CCsJOCs)"""),
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
    label = 'O=CC(OO)C1CCCC1(1993)',
    structure = SMILES('O=CC(OO)C1CCCC1'),
    E0 = (-329.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4681.46,'J/mol'), sigma=(7.80337,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=731.23 K, Pc=22.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0321018,0.063013,3.67252e-05,-9.13792e-08,3.86252e-11,-39459.5,37.6681], Tmin=(100,'K'), Tmax=(989.421,'K')), NASAPolynomial(coeffs=[20.4293,0.0387808,-1.48147e-05,2.82869e-09,-2.07746e-13,-46345.9,-74.9226], Tmin=(989.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-329.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cyclopentane)"""),
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
    label = '[O]COC=C(OO)C1CCCC1(3850)',
    structure = SMILES('[O]COC=C(OO)C1CCCC1'),
    E0 = (-264.155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,350,440,435,1725,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (173.186,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59603,0.0907801,2.0325e-05,-1.04637e-07,4.98209e-11,-31539.8,44.3914], Tmin=(100,'K'), Tmax=(960.776,'K')), NASAPolynomial(coeffs=[31.6765,0.0311994,-9.90374e-06,1.85803e-09,-1.42699e-13,-41576.9,-133.76], Tmin=(960.776,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-264.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(OCOJ)"""),
)

species(
    label = 'O=CO[CH]C(OO)C1CCCC1(3851)',
    structure = SMILES('O=CO[CH]C(OO)C1CCCC1'),
    E0 = (-365.437,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (173.186,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67982,0.104623,-4.43208e-05,-1.79661e-08,1.40656e-11,-43729.9,45.289], Tmin=(100,'K'), Tmax=(1039.45,'K')), NASAPolynomial(coeffs=[24.379,0.0454879,-1.83578e-05,3.46247e-09,-2.46955e-13,-51370,-92.1299], Tmin=(1039.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdOsH) + ring(Cyclopentane) + radical(CCsJOC(O)H)"""),
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
    label = '[CH]1CCCC1(1948)',
    structure = SMILES('[CH]1CCCC1'),
    E0 = (94.2072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,262.657,800.537,800.537,800.537,800.537,800.537,800.537,800.537,800.537,1736.13,1736.13,1736.13,1736.13,1736.13,1736.13,1736.13,1736.13,1736.13],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17949,0.00227247,9.52802e-05,-1.14305e-07,4.10425e-11,11374.2,17.1178], Tmin=(100,'K'), Tmax=(968.374,'K')), NASAPolynomial(coeffs=[6.21397,0.0302738,-1.08828e-05,2.00838e-09,-1.45034e-13,8885.9,-7.23767], Tmin=(968.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.2072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclopentane) + radical(cyclopentane)"""),
)

species(
    label = '[O]COC=COO(3852)',
    structure = SMILES('[O]COC=COO'),
    E0 = (-196.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16619,'amu*angstrom^2'), symmetry=1, barrier=(26.8129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16685,'amu*angstrom^2'), symmetry=1, barrier=(26.8282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16655,'amu*angstrom^2'), symmetry=1, barrier=(26.8212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1662,'amu*angstrom^2'), symmetry=1, barrier=(26.8132,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.529553,0.0578288,-1.30862e-05,-4.37994e-08,2.61229e-11,-23529.5,25.5978], Tmin=(100,'K'), Tmax=(933.607,'K')), NASAPolynomial(coeffs=[24.6718,0.00149538,1.74394e-06,-3.48517e-10,1.66834e-14,-30090.1,-100.21], Tmin=(933.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(OCOJ)"""),
)

species(
    label = 'HO2(26)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263985,1.5223e-05,-1.71671e-08,6.26738e-12,322.677,4.84428], Tmin=(100,'K'), Tmax=(923.913,'K')), NASAPolynomial(coeffs=[4.15133,0.00191146,-4.11274e-07,6.34957e-11,-4.86385e-15,83.4208,3.09341], Tmin=(923.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[O]COC=CC1CCCC1(2198)',
    structure = SMILES('[O]COC=CC1CCCC1'),
    E0 = (-174.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.240107,0.0559763,9.4682e-05,-1.75284e-07,7.43995e-11,-20751.6,35.4243], Tmin=(100,'K'), Tmax=(945.323,'K')), NASAPolynomial(coeffs=[29.1972,0.025487,-6.20566e-06,1.13174e-09,-9.42134e-14,-30520.4,-127.17], Tmin=(945.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(544.598,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(OCOJ)"""),
)

species(
    label = '[O]COC[C](OO)C1CCCC1(3853)',
    structure = SMILES('[O]COC[C](OO)C1CCCC1'),
    E0 = (-165.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80084,0.110523,-6.62598e-05,1.30786e-08,1.05794e-12,-19740.2,46.5006], Tmin=(100,'K'), Tmax=(1207.58,'K')), NASAPolynomial(coeffs=[22.6052,0.0522507,-2.19132e-05,4.07426e-09,-2.82303e-13,-27280.3,-82.6619], Tmin=(1207.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(OCOJ) + radical(C2CsJOOH)"""),
)

species(
    label = 'O[CH]O[CH]C(OO)C1CCCC1(3854)',
    structure = SMILES('O[CH]O[CH]C(OO)C1CCCC1'),
    E0 = (-210.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.95485,0.135752,-0.000122952,5.67465e-08,-1.03473e-11,-25096.2,49.0979], Tmin=(100,'K'), Tmax=(1330.45,'K')), NASAPolynomial(coeffs=[28.9167,0.03993,-1.49196e-05,2.61329e-09,-1.7534e-13,-33576.9,-113.764], Tmin=(1330.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(602.799,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(CCsJOCs) + radical(OCJO)"""),
)

species(
    label = '[O][CH]OCC(OO)C1CCCC1(3855)',
    structure = SMILES('[O][CH]OCC(OO)C1CCCC1'),
    E0 = (-163.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.865233,0.111623,-7.92172e-05,2.92133e-08,-4.59414e-12,-19544.9,42.6033], Tmin=(100,'K'), Tmax=(1403.63,'K')), NASAPolynomial(coeffs=[15.4515,0.0651239,-2.95256e-05,5.61179e-09,-3.90468e-13,-24125.4,-41.648], Tmin=(1403.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-163.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(OCJO) + radical(OCOJ)"""),
)

species(
    label = '[O]COCC(OO)[C]1CCCC1(3856)',
    structure = SMILES('[O]COCC(OO)[C]1CCCC1'),
    E0 = (-161.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,3615,1310,387.5,850,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.2221,0.0934267,-1.52056e-05,-4.13462e-08,2.03653e-11,-19274.2,48.9591], Tmin=(100,'K'), Tmax=(1050.86,'K')), NASAPolynomial(coeffs=[21.5936,0.0522674,-2.16681e-05,4.12522e-09,-2.94738e-13,-26592,-74.2478], Tmin=(1050.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(OCOJ) + radical(C2CJCOOH)"""),
)

species(
    label = '[O]COCC(O[O])C1CCCC1(3857)',
    structure = SMILES('[O]COCC(O[O])C1CCCC1'),
    E0 = (-200.877,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,492.5,1135,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.54359,0.103137,-4.27736e-05,-1.44115e-08,1.15686e-11,-23944.3,46.0262], Tmin=(100,'K'), Tmax=(1063.8,'K')), NASAPolynomial(coeffs=[21.9445,0.051153,-2.07072e-05,3.86698e-09,-2.72417e-13,-30997.5,-78.4061], Tmin=(1063.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(611.114,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(OCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]COCC(OO)C1[CH]CCC1(3858)',
    structure = SMILES('[O]COCC(OO)C1[CH]CCC1'),
    E0 = (-158.339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,3615,1310,387.5,850,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.66369,0.105948,-5.04036e-05,-5.23533e-09,7.73541e-12,-18824.3,48.2616], Tmin=(100,'K'), Tmax=(1114.74,'K')), NASAPolynomial(coeffs=[22.5337,0.0519482,-2.1913e-05,4.14257e-09,-2.92094e-13,-26258.7,-80.2535], Tmin=(1114.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(Cs_S) + radical(OCOJ)"""),
)

species(
    label = 'OCO[CH][C](OO)C1CCCC1(3859)',
    structure = SMILES('OCO[CH][C](OO)C1CCCC1'),
    E0 = (-212.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.44816,0.118347,-5.60549e-05,-2.52899e-08,2.17931e-11,-25355.6,47.7712], Tmin=(100,'K'), Tmax=(969.818,'K')), NASAPolynomial(coeffs=[30.7488,0.0356786,-1.21036e-05,2.17851e-09,-1.56761e-13,-34345.9,-124.521], Tmin=(969.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-212.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(602.799,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(C2CsJOOH) + radical(CCsJOCs)"""),
)

species(
    label = '[O]COCC(OO)C1C[CH]CC1(3860)',
    structure = SMILES('[O]COCC(OO)C1C[CH]CC1'),
    E0 = (-167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67414,0.108334,-6.45235e-05,1.34769e-08,5.7551e-13,-19867.7,48.3177], Tmin=(100,'K'), Tmax=(1230.93,'K')), NASAPolynomial(coeffs=[21.7851,0.0527539,-2.19614e-05,4.05601e-09,-2.79468e-13,-27207.7,-76.0893], Tmin=(1230.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(cyclopentane) + radical(OCOJ)"""),
)

species(
    label = 'OCO[CH]C(OO)[C]1CCCC1(3861)',
    structure = SMILES('OCO[CH]C(OO)[C]1CCCC1'),
    E0 = (-208.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99139,0.102451,-8.02563e-06,-7.75746e-08,4.09919e-11,-24884.1,50.6824], Tmin=(100,'K'), Tmax=(956.708,'K')), NASAPolynomial(coeffs=[31.8609,0.0323123,-9.99861e-06,1.80551e-09,-1.34974e-13,-34628.9,-128.214], Tmin=(956.708,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(602.799,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(C2CJCOOH) + radical(CCsJOCs)"""),
)

species(
    label = '[O]OC([CH]OCO)C1CCCC1(3862)',
    structure = SMILES('[O]OC([CH]OCO)C1CCCC1'),
    E0 = (-247.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.27829,0.111743,-3.4084e-05,-5.2669e-08,3.30948e-11,-29555.7,47.6261], Tmin=(100,'K'), Tmax=(948.559,'K')), NASAPolynomial(coeffs=[32.0821,0.0314161,-9.16276e-06,1.57664e-09,-1.1508e-13,-38979,-131.64], Tmin=(948.559,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-247.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(CCsJOCs) + radical(ROOJ)"""),
)

species(
    label = 'OCO[CH]C(OO)C1[CH]CCC1(3863)',
    structure = SMILES('OCO[CH]C(OO)C1[CH]CCC1'),
    E0 = (-205.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.3834,0.114488,-4.20424e-05,-4.21817e-08,2.82847e-11,-24436.5,49.8009], Tmin=(100,'K'), Tmax=(965.441,'K')), NASAPolynomial(coeffs=[31.9633,0.0333439,-1.09929e-05,1.99484e-09,-1.46284e-13,-33918.7,-129.456], Tmin=(965.441,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-205.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(602.799,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(Cs_S) + radical(CCsJOCs)"""),
)

species(
    label = 'OCO[CH]C(OO)C1C[CH]CC1(3864)',
    structure = SMILES('OCO[CH]C(OO)C1C[CH]CC1'),
    E0 = (-213.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.30628,0.116006,-5.39066e-05,-2.52629e-08,2.14048e-11,-25483.9,49.532], Tmin=(100,'K'), Tmax=(968.464,'K')), NASAPolynomial(coeffs=[29.642,0.0366264,-1.23914e-05,2.21405e-09,-1.58217e-13,-34137.6,-116.305], Tmin=(968.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(602.799,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(CCsJOCs) + radical(cyclopentane)"""),
)

species(
    label = '[O][CH]C(OO)C1CCCC1(2032)',
    structure = SMILES('[O][CH]C(OO)C1CCCC1'),
    E0 = (-3.06575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.720628,0.0925877,-6.32804e-05,2.17495e-08,-3.03318e-12,-189.696,39.0194], Tmin=(100,'K'), Tmax=(1656.29,'K')), NASAPolynomial(coeffs=[20.6219,0.0410454,-1.66022e-05,2.96149e-09,-1.97353e-13,-7259.65,-74.715], Tmin=(1656.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.06575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'O=COCC(OO)C1CCCC1(3865)',
    structure = SMILES('O=COCC(OO)C1CCCC1'),
    E0 = (-560.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29743,0.0957251,-2.08478e-05,-3.51414e-08,1.80549e-11,-67157.9,44.3704], Tmin=(100,'K'), Tmax=(1060.48,'K')), NASAPolynomial(coeffs=[21.4084,0.0533833,-2.22056e-05,4.21555e-09,-3.00018e-13,-74408.6,-77.9855], Tmin=(1060.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-560.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdOsH) + ring(Cyclopentane)"""),
)

species(
    label = 'OCOC=C(OO)C1CCCC1(3866)',
    structure = SMILES('OCOC=C(OO)C1CCCC1'),
    E0 = (-491.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.11776,0.0932782,4.57257e-05,-1.54653e-07,7.40229e-11,-58861.8,45.7887], Tmin=(100,'K'), Tmax=(932.102,'K')), NASAPolynomial(coeffs=[40.2756,0.0167927,-8.70144e-07,3.47449e-11,-1.64325e-14,-71345.2,-180.324], Tmin=(932.102,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-491.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane)"""),
)

species(
    label = 'OH(21)',
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
    label = '[O]COC1OC1C1CCCC1(3867)',
    structure = SMILES('[O]COC1OC1C1CCCC1'),
    E0 = (-275.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (157.187,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.358703,0.0673178,5.73866e-05,-1.29212e-07,5.69349e-11,-32919,37.0501], Tmin=(100,'K'), Tmax=(935.658,'K')), NASAPolynomial(coeffs=[23.5712,0.0373306,-1.04703e-05,1.7384e-09,-1.24632e-13,-40562.4,-93.7216], Tmin=(935.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-275.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsOsHH) + ring(Ethylene_oxide) + ring(Cyclopentane) + radical(OCOJ)"""),
)

species(
    label = '[O]COC([CH]OO)C1CCCC1(3868)',
    structure = SMILES('[O]COC([CH]OO)C1CCCC1'),
    E0 = (-162.207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.86155,0.110078,-5.66635e-05,-2.05544e-09,7.40024e-12,-19282,46.5619], Tmin=(100,'K'), Tmax=(1096.11,'K')), NASAPolynomial(coeffs=[23.4612,0.0508844,-2.11141e-05,3.96977e-09,-2.79683e-13,-26828.7,-87.0316], Tmin=(1096.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(CCsJOOH) + radical(OCOJ)"""),
)

species(
    label = 'OOC(C1CCCC1)C1OCO1(2212)',
    structure = SMILES('OOC(C1CCCC1)C1OCO1'),
    E0 = (-431.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24276,0.0738776,8.52509e-05,-1.76485e-07,7.57455e-11,-51700.9,42.4083], Tmin=(100,'K'), Tmax=(957.744,'K')), NASAPolynomial(coeffs=[33.6741,0.0311787,-9.39595e-06,1.8286e-09,-1.47683e-13,-63119.1,-149.231], Tmin=(957.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-431.747,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsOsHH) + ring(Cyclobutane) + ring(Cyclopentane)"""),
)

species(
    label = '[O]COC(O)C([O])C1CCCC1(3869)',
    structure = SMILES('[O]COC(O)C([O])C1CCCC1'),
    E0 = (-394.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.55502,0.0997422,-2.49271e-05,-3.87045e-08,2.13571e-11,-47198.6,47.5109], Tmin=(100,'K'), Tmax=(1023.93,'K')), NASAPolynomial(coeffs=[24.3458,0.0470582,-1.87939e-05,3.55879e-09,-2.5558e-13,-55045.1,-90.4726], Tmin=(1023.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-394.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(611.114,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsOsHH) + ring(Cyclopentane) + radical(OCOJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]C(OO)C1CCCC1(2038)',
    structure = SMILES('[CH]C(OO)C1CCCC1'),
    E0 = (198.974,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0237545,0.0689259,7.51208e-06,-5.80683e-08,2.66659e-11,24093.2,34.3848], Tmin=(100,'K'), Tmax=(994.323,'K')), NASAPolynomial(coeffs=[19.0121,0.0378969,-1.43923e-05,2.68862e-09,-1.93431e-13,18056,-68.6663], Tmin=(994.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O]C[O](3870)',
    structure = SMILES('[O]C[O]'),
    E0 = (43.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1404.84,1404.86],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.27779,-0.00180121,2.07568e-05,-1.64549e-08,3.76041e-12,5255.97,7.17787], Tmin=(100,'K'), Tmax=(1502.36,'K')), NASAPolynomial(coeffs=[2.16884,0.0155415,-8.26787e-06,1.62055e-09,-1.12064e-13,4566.12,13.806], Tmin=(1502.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCOJ)"""),
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
    label = '[CH2]O[CH]C(OO)C1CCCC1(2039)',
    structure = SMILES('[CH2]O[CH]C(OO)C1CCCC1'),
    E0 = (-20.0954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57447,0.10489,-5.31414e-05,-7.71822e-09,1.08919e-11,-2200.76,44.1225], Tmin=(100,'K'), Tmax=(1019.54,'K')), NASAPolynomial(coeffs=[22.9803,0.0444107,-1.69147e-05,3.0888e-09,-2.1655e-13,-9071.26,-83.9544], Tmin=(1019.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.0954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(CCsJOCs) + radical(CsJOCC)"""),
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
    E0 = (-172.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-41.0732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-136.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-100.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-85.7195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-141.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-15.3434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-48.6526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1.59754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1.62318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-90.9209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-112.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-101.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-147.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-101.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-96.0368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-65.5246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-34.0709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (189.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-109.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-147.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-126.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-2.27196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-164.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-70.7546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (242.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (222.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['CH2O(35)', 'O=CC(OO)C1CCCC1(1993)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(19)', '[O]COC=C(OO)C1CCCC1(3850)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', 'O=CO[CH]C(OO)C1CCCC1(3851)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.6e+09,'cm^3/(mol*s)'), n=0.935, Ea=(17.4473,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [CO-NdH_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][O](354)', 'O=CC(OO)C1CCCC1(1993)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]1CCCC1(1948)', '[O]COC=COO(3852)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00116219,'m^3/(mol*s)'), n=2.41, Ea=(16.8889,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CsJ-CsCsH] for rate rule [Cds-OsH_Cds;CsJ-CsCsH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['HO2(26)', '[O]COC=CC1CCCC1(2198)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.49939e-05,'m^3/(mol*s)'), n=3.11, Ea=(30.2643,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds;OJ-O2s] for rate rule [Cds-CsH_Cds-OsH;OJ-O2s]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]COC[C](OO)C1CCCC1(3853)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.81094e+06,'s^-1'), n=1.92799, Ea=(150.639,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_1H] for rate rule [R2H_S;C_rad_out_NDMustO;Cs_H_out_H/NonDeO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['O[CH]O[CH]C(OO)C1CCCC1(3854)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(779620,'s^-1'), n=2.1225, Ea=(123.773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] + [R2H_S;O_rad_out;Cs_H_out_1H] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][CH]OCC(OO)C1CCCC1(3855)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.53164e+09,'s^-1'), n=1.09, Ea=(165.519,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]COCC(OO)[C]1CCCC1(3856)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.5e+08,'s^-1'), n=1.01, Ea=(163.594,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_Cs2_cy5;Cs_H_out_1H] for rate rule [R3H_SS_Cs;C_rad_out_Cs2_cy5;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]COCC(O[O])C1CCCC1(3857)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.84e+09,'s^-1'), n=0.82, Ea=(109.956,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 273 used for R4H_SSS_O(Cs)Cs;O_rad_out;Cs_H_out_H/NonDeO
Exact match found for rate rule [R4H_SSS_O(Cs)Cs;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]COCC(OO)C1[CH]CCC1(3858)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.00310895,'s^-1'), n=4.065, Ea=(46.2541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['OCO[CH][C](OO)C1CCCC1(3859)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(57.9,'s^-1'), n=2.9, Ea=(71.128,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_OOH/Cs] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_OOH/Cs]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['[O]COCC(OO)C1C[CH]CC1(3860)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.1128,'s^-1'), n=3.28, Ea=(24.7274,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_1H;Cs_H_out_H/NonDeC] for rate rule [R5H_CCC;C_rad_out_H/NonDeO;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['OCO[CH]C(OO)[C]1CCCC1(3861)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.56e+09,'s^-1'), n=0, Ea=(71.3372,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6Hall;O_rad_out;Cs_H_out_Cs2] for rate rule [R6HJ_3;O_rad_out;Cs_H_out_Cs2_cy5]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['[O]OC([CH]OCO)C1CCCC1(3862)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(970995,'s^-1'), n=0.905106, Ea=(76.3888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;O_rad_out;XH_out] for rate rule [R7HJ_3;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['OCO[CH]C(OO)C1[CH]CCC1(3863)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.8e+08,'s^-1'), n=0, Ea=(106.901,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R7HJ_3;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['OCO[CH]C(OO)C1C[CH]CC1(3864)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.18651e+08,'s^-1'), n=1.46287, Ea=(138.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R8Hall;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][O](354)', '[O][CH]C(OO)C1CCCC1(2032)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.35017e+06,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['O=COCC(OO)C1CCCC1(3865)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['OCOC=C(OO)C1CCCC1(3866)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['OH(21)', '[O]COC1OC1C1CCCC1(3867)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.38e+12,'s^-1','*|/',1.2), n=0, Ea=(45.6433,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2OO_S;C_sec_rad_intra;OOH] for rate rule [R2OO_S;C_rad/H/NonDeO_intra;OOH]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]COC([CH]OO)C1CCCC1(3868)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['OOC(C1CCCC1)C1OCO1(2212)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    products = ['[O]COC(O)C([O])C1CCCC1(3869)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.69e+11,'s^-1'), n=0, Ea=(101.671,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2OOH_S;C_rad_out_1H] for rate rule [R2OOH_S;C_rad_out_H/NonDeO]
Euclidian distance = 1.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C(OO)C1CCCC1(2038)', '[O]C[O](3870)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(109477,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['O(20)', '[CH2]O[CH]C(OO)C1CCCC1(2039)'],
    products = ['[O]CO[CH]C(OO)C1CCCC1(2208)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/O;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '461',
    isomers = [
        '[O]CO[CH]C(OO)C1CCCC1(2208)',
    ],
    reactants = [
        ('CH2O(35)', 'O=CC(OO)C1CCCC1(1993)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '461',
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

