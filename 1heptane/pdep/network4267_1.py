species(
    label = 'C=C([CH]C)C(=C)[CH]C(24182)',
    structure = SMILES('[CH2]C(=CC)C([CH2])=CC'),
    E0 = (249.687,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.735277,'amu*angstrom^2'), symmetry=1, barrier=(16.9055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0632434,'amu*angstrom^2'), symmetry=1, barrier=(29.514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.737545,'amu*angstrom^2'), symmetry=1, barrier=(16.9576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.732781,'amu*angstrom^2'), symmetry=1, barrier=(16.8481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.739219,'amu*angstrom^2'), symmetry=1, barrier=(16.9961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.384005,0.0840749,-5.09991e-05,5.50851e-09,4.14197e-12,30198.9,28.4131], Tmin=(100,'K'), Tmax=(1039.09,'K')), NASAPolynomial(coeffs=[18.1326,0.0354522,-1.35159e-05,2.44392e-09,-1.69358e-13,25127.7,-67.5143], Tmin=(1039.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'CH3CHCCH2(18175)',
    structure = SMILES('C=C=CC'),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.759584,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74635,0.0218189,8.22353e-06,-2.14768e-08,8.55624e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.6,'K')), NASAPolynomial(coeffs=[6.82078,0.0192338,-7.45622e-06,1.36536e-09,-9.53195e-14,16028,-10.4333], Tmin=(1025.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C1([CH]C)CC1=CC(25275)',
    structure = SMILES('[CH2]C1([CH]C)CC1=CC'),
    E0 = (462.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.263258,0.0692237,-2.26363e-05,-1.35463e-08,8.13734e-12,55737.7,31.4039], Tmin=(100,'K'), Tmax=(1105.46,'K')), NASAPolynomial(coeffs=[15.171,0.0400578,-1.66801e-05,3.13624e-09,-2.2049e-13,50927.8,-48.8594], Tmin=(1105.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Methylene_cyclopropane) + radical(Neopentyl) + radical(Cs_S)"""),
)

species(
    label = 'C=[C][CH]C(18176)',
    structure = SMILES('[CH2][C]=CC'),
    E0 = (361.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.352622,'amu*angstrom^2'), symmetry=1, barrier=(8.10748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828631,'amu*angstrom^2'), symmetry=1, barrier=(19.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12013e-13,43485.7,14.8304], Tmin=(100,'K'), Tmax=(2065.83,'K')), NASAPolynomial(coeffs=[10.7464,0.014324,-5.20136e-06,8.69079e-10,-5.48385e-14,40045.6,-31.3799], Tmin=(2065.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=CC)C(C)=[C]C(25412)',
    structure = SMILES('[CH2]C(=CC)C(C)=[C]C'),
    E0 = (336.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,1685,370,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,222.04],'cm^-1')),
        HinderedRotor(inertia=(0.395973,'amu*angstrom^2'), symmetry=1, barrier=(13.8694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.396086,'amu*angstrom^2'), symmetry=1, barrier=(13.8683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395737,'amu*angstrom^2'), symmetry=1, barrier=(13.8691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395039,'amu*angstrom^2'), symmetry=1, barrier=(13.8689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395901,'amu*angstrom^2'), symmetry=1, barrier=(13.8689,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.116365,0.0876489,-7.20737e-05,3.21805e-08,-5.96317e-12,40565.5,28.3373], Tmin=(100,'K'), Tmax=(1264.63,'K')), NASAPolynomial(coeffs=[14.5979,0.041109,-1.68732e-05,3.08148e-09,-2.10818e-13,36843.8,-46.1055], Tmin=(1264.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=[C]C)C(C)=CC(25413)',
    structure = SMILES('[CH2]C(=[C]C)C(C)=CC'),
    E0 = (336.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,1685,370,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,222.04],'cm^-1')),
        HinderedRotor(inertia=(0.395973,'amu*angstrom^2'), symmetry=1, barrier=(13.8694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.396086,'amu*angstrom^2'), symmetry=1, barrier=(13.8683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395737,'amu*angstrom^2'), symmetry=1, barrier=(13.8691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395039,'amu*angstrom^2'), symmetry=1, barrier=(13.8689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395901,'amu*angstrom^2'), symmetry=1, barrier=(13.8689,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.116365,0.0876489,-7.20737e-05,3.21805e-08,-5.96317e-12,40565.5,28.3373], Tmin=(100,'K'), Tmax=(1264.63,'K')), NASAPolynomial(coeffs=[14.5979,0.041109,-1.68732e-05,3.08148e-09,-2.10818e-13,36843.8,-46.1055], Tmin=(1264.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=CC)[C](C)C=C(24605)',
    structure = SMILES('[CH2]C=C(C)C([CH2])=CC'),
    E0 = (216.244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.712083,'amu*angstrom^2'), symmetry=1, barrier=(16.3722,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.555659,'amu*angstrom^2'), symmetry=1, barrier=(96.3851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0202512,'amu*angstrom^2'), symmetry=1, barrier=(16.3711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712008,'amu*angstrom^2'), symmetry=1, barrier=(16.3705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.19211,'amu*angstrom^2'), symmetry=1, barrier=(96.3849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0883175,0.0775021,-3.58132e-05,-7.55711e-09,8.27771e-12,26166.1,29.3215], Tmin=(100,'K'), Tmax=(1017.17,'K')), NASAPolynomial(coeffs=[16.4341,0.0376674,-1.41425e-05,2.53759e-09,-1.75328e-13,21504.4,-57.0638], Tmin=(1017.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][C](C=C)C(C)=CC(24606)',
    structure = SMILES('[CH2]C=C([CH2])C(C)=CC'),
    E0 = (216.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0883175,0.0775021,-3.58132e-05,-7.55711e-09,8.27771e-12,26166.1,29.3215], Tmin=(100,'K'), Tmax=(1017.17,'K')), NASAPolynomial(coeffs=[16.4341,0.0376674,-1.41425e-05,2.53759e-09,-1.75328e-13,21504.4,-57.0638], Tmin=(1017.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C(=CC)[C]1CC1C(25414)',
    structure = SMILES('[CH2]C(=CC)[C]1CC1C'),
    E0 = (289.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.71289,0.0520158,3.84829e-05,-8.55933e-08,3.61457e-11,35003.5,26.4903], Tmin=(100,'K'), Tmax=(968.714,'K')), NASAPolynomial(coeffs=[16.7686,0.0352996,-1.24057e-05,2.26286e-09,-1.62921e-13,29566.5,-62.466], Tmin=(968.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(Allyl_T) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]1C(=CC)CC1C(25415)',
    structure = SMILES('[CH2]C1=C([CH]C)CC1C'),
    E0 = (304.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.583091,0.0531885,4.0938e-05,-9.08388e-08,3.83549e-11,36774.2,26.4705], Tmin=(100,'K'), Tmax=(972.301,'K')), NASAPolynomial(coeffs=[18.2947,0.0339462,-1.21014e-05,2.24934e-09,-1.64353e-13,30795.4,-71.5147], Tmin=(972.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(Cyclobutene) + radical(Allyl_P) + radical(Allyl_S)"""),
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
    label = '[CH2]C(=C)C([CH2])=CC(25416)',
    structure = SMILES('[CH2]C(=C)C([CH2])=CC'),
    E0 = (285.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,311.383],'cm^-1')),
        HinderedRotor(inertia=(0.327475,'amu*angstrom^2'), symmetry=1, barrier=(22.5291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.327466,'amu*angstrom^2'), symmetry=1, barrier=(22.5294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.327318,'amu*angstrom^2'), symmetry=1, barrier=(22.5272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.327483,'amu*angstrom^2'), symmetry=1, barrier=(22.5297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.335271,0.0676667,-2.76626e-05,-1.62749e-08,1.21982e-11,34506.8,24.024], Tmin=(100,'K'), Tmax=(980.594,'K')), NASAPolynomial(coeffs=[17.5531,0.0266059,-9.47854e-06,1.70194e-09,-1.19937e-13,29727.4,-65.8563], Tmin=(980.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C=C([CH]C)C[C]=CC(24184)',
    structure = SMILES('[CH2]C(=CC)C[C]=CC'),
    E0 = (366.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1685,370,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,579.702],'cm^-1')),
        HinderedRotor(inertia=(0.147406,'amu*angstrom^2'), symmetry=1, barrier=(3.38916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64226,'amu*angstrom^2'), symmetry=1, barrier=(14.7668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64164,'amu*angstrom^2'), symmetry=1, barrier=(14.7526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.643937,'amu*angstrom^2'), symmetry=1, barrier=(14.8054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145327,'amu*angstrom^2'), symmetry=1, barrier=(3.34136,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3683.66,'J/mol'), sigma=(6.4482,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=575.38 K, Pc=31.18 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.29648,0.0786067,-5.42868e-05,1.96375e-08,-2.97459e-12,44273.2,31.2372], Tmin=(100,'K'), Tmax=(1490.43,'K')), NASAPolynomial(coeffs=[13.9025,0.0420909,-1.75363e-05,3.199e-09,-2.17227e-13,40217.5,-39.8334], Tmin=(1490.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'CC=C1CCC1=CC(25269)',
    structure = SMILES('CC=C1CCC1=CC'),
    E0 = (114.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.677799,0.0585738,5.80411e-06,-4.1598e-08,1.78951e-11,13856,25.5085], Tmin=(100,'K'), Tmax=(1034.79,'K')), NASAPolynomial(coeffs=[13.4814,0.0415234,-1.65073e-05,3.07348e-09,-2.16896e-13,9469.28,-45.0922], Tmin=(1034.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(12methylenecyclobutane)"""),
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
    label = '[CH2]C([C]=CC)=CC(25417)',
    structure = SMILES('[CH2]C([C]=CC)=CC'),
    E0 = (334.774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.7606,'amu*angstrom^2'), symmetry=1, barrier=(17.4877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.760854,'amu*angstrom^2'), symmetry=1, barrier=(17.4935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.760586,'amu*angstrom^2'), symmetry=1, barrier=(17.4874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15146,'amu*angstrom^2'), symmetry=1, barrier=(49.4663,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.352604,0.0734369,-5.91187e-05,2.57941e-08,-4.60694e-12,40400.9,25.1788], Tmin=(100,'K'), Tmax=(1327.42,'K')), NASAPolynomial(coeffs=[14.2321,0.0316126,-1.18565e-05,2.05761e-09,-1.36512e-13,36716.1,-45.7131], Tmin=(1327.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=CJC=C) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1([CH]C)C(=C)C1C(25296)',
    structure = SMILES('[CH2]C1([CH]C)C(=C)C1C'),
    E0 = (466.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.29276,0.0655305,-4.50464e-06,-3.74661e-08,1.7759e-11,56253.7,30.0992], Tmin=(100,'K'), Tmax=(1027.4,'K')), NASAPolynomial(coeffs=[16.6435,0.0372633,-1.49065e-05,2.81296e-09,-2.01072e-13,51026,-58.316], Tmin=(1027.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Neopentyl) + radical(Cs_S)"""),
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
    label = '[CH2]C(=CC)C(=C)C=C(24604)',
    structure = SMILES('[CH2]C(=CC)C(=C)C=C'),
    E0 = (242.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,181.962,683.313],'cm^-1')),
        HinderedRotor(inertia=(0.669842,'amu*angstrom^2'), symmetry=1, barrier=(19.1337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0582339,'amu*angstrom^2'), symmetry=1, barrier=(19.1767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.83204,'amu*angstrom^2'), symmetry=1, barrier=(19.1302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.52237,'amu*angstrom^2'), symmetry=1, barrier=(104.569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.293043,0.0682771,-2.00337e-05,-2.05401e-08,1.21516e-11,29332.3,27.0261], Tmin=(100,'K'), Tmax=(1018.57,'K')), NASAPolynomial(coeffs=[15.7386,0.0358123,-1.37404e-05,2.51366e-09,-1.76142e-13,24723.4,-54.9529], Tmin=(1018.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC(=C)C([CH2])=CC(25418)',
    structure = SMILES('[CH2]CC(=C)C([CH2])=CC'),
    E0 = (316.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0368535,'amu*angstrom^2'), symmetry=1, barrier=(17.9864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00736317,'amu*angstrom^2'), symmetry=1, barrier=(3.60618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781153,'amu*angstrom^2'), symmetry=1, barrier=(17.9602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779478,'amu*angstrom^2'), symmetry=1, barrier=(17.9217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781104,'amu*angstrom^2'), symmetry=1, barrier=(17.9591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.348925,0.0836004,-5.1879e-05,7.14877e-09,3.44908e-12,38270.9,31.5928], Tmin=(100,'K'), Tmax=(1044.14,'K')), NASAPolynomial(coeffs=[17.9255,0.0352115,-1.34219e-05,2.42456e-09,-1.67785e-13,33276.3,-63.0036], Tmin=(1044.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(CC)C([CH2])=CC(25419)',
    structure = SMILES('[CH]=C(CC)C([CH2])=CC'),
    E0 = (358.664,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.701639,'amu*angstrom^2'), symmetry=1, barrier=(16.1321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344302,'amu*angstrom^2'), symmetry=1, barrier=(16.1602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0492932,'amu*angstrom^2'), symmetry=1, barrier=(16.1378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.702005,'amu*angstrom^2'), symmetry=1, barrier=(16.1405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.702379,'amu*angstrom^2'), symmetry=1, barrier=(16.1491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.468616,0.0864938,-5.84569e-05,1.27697e-08,1.75707e-12,43308.4,30.6389], Tmin=(100,'K'), Tmax=(1047.28,'K')), NASAPolynomial(coeffs=[18.4195,0.034593,-1.31104e-05,2.35762e-09,-1.62637e-13,38242.2,-66.6572], Tmin=(1047.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(=[C]C)C(=C)CC(25420)',
    structure = SMILES('[CH2]C(=[C]C)C(=C)CC'),
    E0 = (349.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.159905,'amu*angstrom^2'), symmetry=1, barrier=(15.9368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.693159,'amu*angstrom^2'), symmetry=1, barrier=(15.9371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.693127,'amu*angstrom^2'), symmetry=1, barrier=(15.9364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.693165,'amu*angstrom^2'), symmetry=1, barrier=(15.9372,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0150632,'amu*angstrom^2'), symmetry=1, barrier=(15.9371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.583231,0.089245,-7.16619e-05,3.00631e-08,-5.07891e-12,42198.9,31.1306], Tmin=(100,'K'), Tmax=(1412.15,'K')), NASAPolynomial(coeffs=[19.0319,0.0336833,-1.2643e-05,2.20036e-09,-1.46165e-13,36659.1,-70.2702], Tmin=(1412.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([CH]C)C(C)=CC(25421)',
    structure = SMILES('[CH]C(=CC)C(C)=CC'),
    E0 = (317.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.247945,0.0873521,-6.16843e-05,2.31486e-08,-3.62747e-12,38328.8,29.1665], Tmin=(100,'K'), Tmax=(1460.93,'K')), NASAPolynomial(coeffs=[15.297,0.0447902,-1.7984e-05,3.20673e-09,-2.14924e-13,33786.8,-51.7212], Tmin=(1460.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(317.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C](C=C)C(=C)CC(24623)',
    structure = SMILES('[CH2]C(C=C)=C([CH2])CC'),
    E0 = (228.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0497728,0.0733281,-1.6094e-05,-3.35123e-08,1.88363e-11,27601.1,30.4448], Tmin=(100,'K'), Tmax=(975.095,'K')), NASAPolynomial(coeffs=[18.3695,0.0342638,-1.21408e-05,2.16747e-09,-1.52112e-13,22274,-66.8493], Tmin=(975.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Allyl_P)"""),
)

species(
    label = 'C[CH][C]1CCC1=CC(25422)',
    structure = SMILES('C[CH]C1CCC=1[CH]C'),
    E0 = (303.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.788866,0.0500701,4.22235e-05,-8.64809e-08,3.53174e-11,36611.5,25.2586], Tmin=(100,'K'), Tmax=(987.239,'K')), NASAPolynomial(coeffs=[16.2187,0.0373502,-1.4111e-05,2.65357e-09,-1.92503e-13,31138.2,-61.2734], Tmin=(987.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(Cyclobutene) + radical(Allyl_S) + radical(Allyl_S)"""),
)

species(
    label = '[CH2][C]1C(=C)C(C)C1C(25423)',
    structure = SMILES('[CH2]C1=C([CH2])C(C)C1C'),
    E0 = (305.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.377097,0.0563026,3.9705e-05,-9.53284e-08,4.14811e-11,36937,26.2973], Tmin=(100,'K'), Tmax=(959.735,'K')), NASAPolynomial(coeffs=[20.4056,0.0304853,-1.006e-05,1.83774e-09,-1.35603e-13,30437.2,-83.3398], Tmin=(959.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(Cyclobutene) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C=CC(=C)C(C)=CC(24616)',
    structure = SMILES('C=CC(=C)C(C)=CC'),
    E0 = (91.1774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.236638,0.0713806,-3.04205e-05,-5.26762e-09,5.54498e-12,11111.2,26.9518], Tmin=(100,'K'), Tmax=(1093.32,'K')), NASAPolynomial(coeffs=[14.1536,0.040705,-1.6104e-05,2.93544e-09,-2.02595e-13,6858.32,-46.9636], Tmin=(1093.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.1774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=[C]C(C)C(=C)[CH]C(24183)',
    structure = SMILES('[CH2]C(=CC)C(C)[C]=C'),
    E0 = (369.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,345.333,347.343],'cm^-1')),
        HinderedRotor(inertia=(0.119405,'amu*angstrom^2'), symmetry=1, barrier=(9.93037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281457,'amu*angstrom^2'), symmetry=1, barrier=(24.022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116909,'amu*angstrom^2'), symmetry=1, barrier=(9.94809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117447,'amu*angstrom^2'), symmetry=1, barrier=(9.9744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116555,'amu*angstrom^2'), symmetry=1, barrier=(9.93684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.33,'J/mol'), sigma=(6.4092,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=566.27 K, Pc=31.24 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.299693,0.0839308,-6.74533e-05,3.06742e-08,-6.02582e-12,44564.4,29.0122], Tmin=(100,'K'), Tmax=(1163.73,'K')), NASAPolynomial(coeffs=[10.857,0.0476425,-2.06788e-05,3.8782e-09,-2.69295e-13,42107.3,-23.5217], Tmin=(1163.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C1C(=CC)CC1C(25265)',
    structure = SMILES('C=C1C(=CC)CC1C'),
    E0 = (118.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.689924,0.0550304,2.3689e-05,-6.56265e-08,2.77602e-11,14372.8,24.9628], Tmin=(100,'K'), Tmax=(993.204,'K')), NASAPolynomial(coeffs=[15.3775,0.0380508,-1.43595e-05,2.66472e-09,-1.90565e-13,9375.16,-56.2678], Tmin=(993.204,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane)"""),
)

species(
    label = 'CHCH3(T)(95)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438699,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82363,-0.000909515,3.2138e-05,-3.7348e-08,1.3309e-11,41371.4,7.10948], Tmin=(100,'K'), Tmax=(960.812,'K')), NASAPolynomial(coeffs=[4.30487,0.00943069,-3.27559e-06,5.95121e-10,-4.27307e-14,40709.1,1.84202], Tmin=(960.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C([C]=C)=CC(24774)',
    structure = SMILES('[CH2]C([C]=C)=CC'),
    E0 = (370.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.17315,'amu*angstrom^2'), symmetry=1, barrier=(26.9731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17496,'amu*angstrom^2'), symmetry=1, barrier=(27.0146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1727,'amu*angstrom^2'), symmetry=1, barrier=(26.9626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0818,0.0569416,-3.56598e-05,4.1841e-09,3.20998e-12,44708.4,20.7527], Tmin=(100,'K'), Tmax=(982.69,'K')), NASAPolynomial(coeffs=[12.9204,0.0239405,-8.46845e-06,1.46434e-09,-9.91425e-14,41648.3,-39.886], Tmin=(982.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C([CH]C)C(=C)CC(25424)',
    structure = SMILES('[CH]C(=CC)C(=C)CC'),
    E0 = (330.753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.442166,0.0858934,-5.1432e-05,9.5936e-09,1.54315e-12,39950.3,30.9724], Tmin=(100,'K'), Tmax=(1106.5,'K')), NASAPolynomial(coeffs=[16.3579,0.0427111,-1.66841e-05,2.99222e-09,-2.04007e-13,35158.1,-56.633], Tmin=(1106.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC(=C)C(=C)CC(24630)',
    structure = SMILES('C=CC(=C)C(=C)CC'),
    E0 = (104.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.296747,0.0670054,-1.0269e-05,-3.13536e-08,1.59568e-11,12721.3,27.8384], Tmin=(100,'K'), Tmax=(1010.3,'K')), NASAPolynomial(coeffs=[15.6889,0.0379462,-1.44599e-05,2.64736e-09,-1.86033e-13,7984.11,-54.6302], Tmin=(1010.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C1C(=C)C(C)C1C(25274)',
    structure = SMILES('C=C1C(=C)C(C)C1C'),
    E0 = (122.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.691732,0.0515838,4.13669e-05,-8.96066e-08,3.77135e-11,14890,23.0693], Tmin=(100,'K'), Tmax=(969.873,'K')), NASAPolynomial(coeffs=[17.4573,0.0342784,-1.20439e-05,2.21718e-09,-1.61071e-13,9199.74,-69.8715], Tmin=(969.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12methylenecyclobutane)"""),
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
    E0 = (291.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (462.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (538.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (497.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (380.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (399.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (350.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (722.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (343.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (380.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (705.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (537.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (257.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (716.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (466.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (454.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (430.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (503.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (393.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (361.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (350.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (380.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (375.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (274.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (463.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (257.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (714.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (375.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (258.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (257.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['CH3CHCCH2(18175)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(41.5431,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 41.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['[CH2]C1([CH]C)CC1=CC(25275)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.36e+09,'s^-1'), n=0.84, Ea=(212.534,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra_cs2H] for rate rule [R4_S_(Cd)_D;doublebond_intra_HNd;radadd_intra_cs2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 210.2 to 212.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH3CHCCH2(18175)', 'C=[C][CH]C(18176)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(=CC)C(C)=[C]C(25412)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(=[C]C)C(C)=CC(25413)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['[CH2]C(=CC)[C](C)C=C(24605)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.6e+06,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['[CH2][C](C=C)C(C)=CC(24606)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.66e+06,'s^-1'), n=1.64, Ea=(100.416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 96 used for R5H_SS(D)MS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_SS(D)MS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C][CH]C(18176)', 'C=[C][CH]C(18176)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.73038e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['[CH2]C(=CC)[C]1CC1C(25414)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.36786e+12,'s^-1'), n=-0.105173, Ea=(93.5715,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe_HNd;radadd_intra_cs2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['[CH2][C]1C(=CC)CC1C(25415)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.43734e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH2(S)(23)', '[CH2]C(=C)C([CH2])=CC(25416)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([CH]C)C[C]=CC(24184)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['CC=C1CCC1=CC(25269)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CH2(19)', '[CH2]C([C]=CC)=CC(25417)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['[CH2]C1([CH]C)C(=C)C1C(25296)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.72658e+10,'s^-1'), n=0.535608, Ea=(216.807,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S_D;doublebond_intra;radadd_intra_csHNd] + [R4_S_D;doublebond_intra_HNd;radadd_intra_cs] for rate rule [R4_S_(Cd)_D;doublebond_intra_HNd;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 214.2 to 216.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2]C(=CC)C(=C)C=C(24604)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC(=C)C([CH2])=CC(25418)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(CC)C([CH2])=CC(25419)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(=[C]C)C(=C)CC(25420)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_H/NonDeC]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C([CH]C)C(C)=CC(25421)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['[CH2][C](C=C)C(=C)CC(24623)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.66e+06,'s^-1'), n=1.64, Ea=(100.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SS(D)MS;C_rad_out_single;Cs_H_out_2H] for rate rule [R5H_SS(D)MS;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['C[CH][C]1CCC1=CC(25422)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.21867e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['[CH2][C]1C(=C)C(C)C1C(25423)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.16207e+08,'s^-1'), n=0.911389, Ea=(125.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['C=CC(=C)C(C)=CC(24616)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.27566e+10,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]C(C)C(=C)[CH]C(24183)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['C=C1C(=CC)CC1C(25265)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CHCH3(T)(95)', '[CH2]C([C]=C)=CC(24774)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C([CH]C)C(=C)CC(25424)'],
    products = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['C=CC(=C)C(=C)CC(24630)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.926e+10,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C([CH]C)C(=C)[CH]C(24182)'],
    products = ['C=C1C(=C)C(C)C1C(25274)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

network(
    label = '4267',
    isomers = [
        'C=C([CH]C)C(=C)[CH]C(24182)',
    ],
    reactants = [
        ('CH3CHCCH2(18175)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4267',
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

