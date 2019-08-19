species(
    label = '[O]CC1[CH]C(CCCC1)OO(1365)',
    structure = SMILES('[O]CC1[CH]C(CCCC1)OO'),
    E0 = (-2.97398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11172,0.0669519,0.000106148,-2.14136e-07,9.54693e-11,-131.44,43.5085], Tmin=(100,'K'), Tmax=(921.921,'K')), NASAPolynomial(coeffs=[38.4548,0.0129801,2.46224e-06,-6.79309e-10,3.39735e-14,-12428.7,-171.288], Tmin=(921.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.97398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCOJ) + radical(CCJCOOH)"""),
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
    label = 'OOC1C=CCCCC1(259)',
    structure = SMILES('OOC1C=CCCCC1'),
    E0 = (-124.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,3615,1310,387.5,850,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4468.79,'J/mol'), sigma=(7.53736,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=698.02 K, Pc=23.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.372801,0.0206158,0.000217263,-3.35136e-07,1.4146e-10,-14763.2,32.5035], Tmin=(100,'K'), Tmax=(913.852,'K')), NASAPolynomial(coeffs=[42.4349,-0.00964692,1.44126e-05,-2.93401e-09,1.83612e-13,-28874.9,-201.781], Tmin=(913.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene)"""),
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
    label = '[O]CC1=CC(CCCC1)OO(1640)',
    structure = SMILES('[O]CC1=CC(CCCC1)OO'),
    E0 = (-90.7785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.187,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.674513,0.0526417,0.000142826,-2.53564e-07,1.10016e-10,-10703.2,41.6336], Tmin=(100,'K'), Tmax=(921.444,'K')), NASAPolynomial(coeffs=[39.6659,0.00732407,5.29804e-06,-1.1861e-09,6.57522e-14,-23647.9,-179.586], Tmin=(921.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.7785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCOJ)"""),
)

species(
    label = '[O]CC1C=C(CCCC1)OO(1641)',
    structure = SMILES('[O]CC1C=C(CCCC1)OO'),
    E0 = (-78.9294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.187,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.607643,0.0529289,0.000136715,-2.45077e-07,1.06744e-10,-9282.11,41.656], Tmin=(100,'K'), Tmax=(919.515,'K')), NASAPolynomial(coeffs=[38.5505,0.00810886,5.06609e-06,-1.17128e-09,6.66171e-14,-21789.9,-172.828], Tmin=(919.515,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.9294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCOJ)"""),
)

species(
    label = 'O=CC1[CH]C(CCCC1)OO(1642)',
    structure = SMILES('O=CC1[CH]C(CCCC1)OO'),
    E0 = (-148.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3615,1310,387.5,850,1000,2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.187,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02816,0.0622168,0.000119237,-2.32747e-07,1.03748e-10,-17618.7,41.819], Tmin=(100,'K'), Tmax=(917.963,'K')), NASAPolynomial(coeffs=[40.8526,0.00553285,6.28078e-06,-1.41075e-09,8.37324e-14,-30608.5,-185.519], Tmin=(917.963,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cycloheptane) + radical(CCJCOOH)"""),
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
    label = 'OOC1[CH][CH]CCCC1(365)',
    structure = SMILES('OOC1[CH][CH]CCCC1'),
    E0 = (133.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.113571,0.0375133,0.0001521,-2.5566e-07,1.10283e-10,16228.8,36.6232], Tmin=(100,'K'), Tmax=(914.695,'K')), NASAPolynomial(coeffs=[37.2426,4.92732e-05,8.71038e-06,-1.86634e-09,1.15238e-13,4211.41,-167.754], Tmin=(914.695,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(cycloheptane)"""),
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
    label = '[O]CC1C=CCCCC1(1643)',
    structure = SMILES('[O]CC1C=CCCCC1'),
    E0 = (11.1344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728492,0.0183427,0.000210411,-3.15065e-07,1.31155e-10,1506.93,32.7608], Tmin=(100,'K'), Tmax=(914.629,'K')), NASAPolynomial(coeffs=[36.2479,0.00209846,8.93501e-06,-1.93775e-09,1.18425e-13,-10808.4,-167.235], Tmin=(914.629,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.1344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCOJ)"""),
)

species(
    label = '[O]C[C]1CCCCC(C1)OO(1644)',
    structure = SMILES('[O]C[C]1CCCCC(C1)OO'),
    E0 = (-50.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.868308,0.058318,0.000133228,-2.42754e-07,1.05652e-10,-5891.13,41.218], Tmin=(100,'K'), Tmax=(922.858,'K')), NASAPolynomial(coeffs=[38.7449,0.0125549,2.91837e-06,-7.51261e-10,3.68424e-14,-18565.3,-175.768], Tmin=(922.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCJ(C)CO) + radical(CCOJ)"""),
)

species(
    label = '[O]CC1CCCC[C](C1)OO(1645)',
    structure = SMILES('[O]CC1CCCC[C](C1)OO'),
    E0 = (-16.4898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00132,0.0674612,9.6957e-05,-1.98493e-07,8.83292e-11,-1763.8,41.7941], Tmin=(100,'K'), Tmax=(924.524,'K')), NASAPolynomial(coeffs=[35.8937,0.0175392,-4.02846e-08,-1.99892e-10,1.85644e-15,-13274.4,-158.664], Tmin=(924.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.4898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(C2CsJOOH) + radical(CCOJ)"""),
)

species(
    label = 'O[CH]C1[CH]C(CCCC1)OO(1646)',
    structure = SMILES('O[CH]C1[CH]C(CCCC1)OO'),
    E0 = (-48.3812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68017,0.0776775,8.99774e-05,-2.10608e-07,9.80118e-11,-5570.3,44.2386], Tmin=(100,'K'), Tmax=(913.18,'K')), NASAPolynomial(coeffs=[43.4889,0.00384102,7.54918e-06,-1.71124e-09,1.07508e-13,-18990.7,-197.887], Tmin=(913.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.3812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCsJOH) + radical(CCJCOOH)"""),
)

species(
    label = '[O]CC1[CH]CCCC(C1)OO(1647)',
    structure = SMILES('[O]CC1[CH]CCCC(C1)OO'),
    E0 = (-33.9601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.967572,0.0654673,0.000103892,-2.06944e-07,9.16369e-11,-3865.04,43.0331], Tmin=(100,'K'), Tmax=(924.403,'K')), NASAPolynomial(coeffs=[36.5433,0.01611,6.91811e-07,-3.30651e-10,1.00282e-14,-15626.2,-161.091], Tmin=(924.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.9601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCOJ) + radical(cycloheptane)"""),
)

species(
    label = '[O][CH]C1CCCCC(C1)OO(1648)',
    structure = SMILES('[O][CH]C1CCCCC(C1)OO'),
    E0 = (-23.0911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.36867,0.0732132,9.248e-05,-2.03214e-07,9.24771e-11,-2542.22,41.6755], Tmin=(100,'K'), Tmax=(919.819,'K')), NASAPolynomial(coeffs=[39.6416,0.0113893,3.28955e-06,-8.54801e-10,4.72188e-14,-15015.7,-179.541], Tmin=(919.819,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.0911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O]CC1CCC[CH]C(C1)OO(1649)',
    structure = SMILES('[O]CC1CCC[CH]C(C1)OO'),
    E0 = (-2.97398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11164,0.0669509,0.000106151,-2.14141e-07,9.54717e-11,-131.444,43.5082], Tmin=(100,'K'), Tmax=(921.917,'K')), NASAPolynomial(coeffs=[38.4546,0.0129805,2.462e-06,-6.7925e-10,3.39686e-14,-12428.6,-171.286], Tmin=(921.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.97398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(CCOJ)"""),
)

species(
    label = 'OC[C]1[CH]C(CCCC1)OO(1650)',
    structure = SMILES('OC[C]1[CH]C(CCCC1)OO'),
    E0 = (-76.1052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18099,0.0627991,0.00013065,-2.50023e-07,1.11118e-10,-8919.16,43.7852], Tmin=(100,'K'), Tmax=(916.535,'K')), NASAPolynomial(coeffs=[42.5872,0.00501474,7.17346e-06,-1.60665e-09,9.70465e-14,-22538.1,-194.084], Tmin=(916.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.1052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(CCJ(C)CO)"""),
)

species(
    label = '[O]CC1CCCCC(C1)O[O](1651)',
    structure = SMILES('[O]CC1CCCCC(C1)O[O]'),
    E0 = (-51.3842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.853226,0.0610974,0.000118196,-2.25131e-07,9.94332e-11,-5962.9,41.035], Tmin=(100,'K'), Tmax=(918.101,'K')), NASAPolynomial(coeffs=[37.4173,0.0129568,3.08366e-06,-8.44783e-10,4.7091e-14,-17988.5,-167.549], Tmin=(918.101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.3842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = 'OCC1[CH]CCCC([CH]1)OO(1652)',
    structure = SMILES('OCC1[CH]CCCC([CH]1)OO'),
    E0 = (-59.2503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,3615,1310,387.5,850,1000,3615,1277.5,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27919,0.0699359,0.000101357,-2.14265e-07,9.71232e-11,-6893.11,45.5965], Tmin=(100,'K'), Tmax=(917.242,'K')), NASAPolynomial(coeffs=[40.3795,0.00858026,4.94088e-06,-1.18461e-09,7.01136e-14,-19596.5,-179.373], Tmin=(917.242,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.2503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(cycloheptane)"""),
)

species(
    label = '[O]CC1C[CH]CCC(C1)OO(1653)',
    structure = SMILES('[O]CC1C[CH]CCC(C1)OO'),
    E0 = (-33.9601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.967572,0.0654673,0.000103892,-2.06944e-07,9.16369e-11,-3865.04,43.0331], Tmin=(100,'K'), Tmax=(924.403,'K')), NASAPolynomial(coeffs=[36.5433,0.01611,6.91811e-07,-3.30651e-10,1.00282e-14,-15626.2,-161.091], Tmin=(924.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.9601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(cycloheptane) + radical(CCOJ)"""),
)

species(
    label = '[O]CC1CC[CH]CC(C1)OO(1654)',
    structure = SMILES('[O]CC1CC[CH]CC(C1)OO'),
    E0 = (-33.9601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.967572,0.0654673,0.000103892,-2.06944e-07,9.16369e-11,-3865.04,43.0331], Tmin=(100,'K'), Tmax=(924.403,'K')), NASAPolynomial(coeffs=[36.5433,0.01611,6.91811e-07,-3.30651e-10,1.00282e-14,-15626.2,-161.091], Tmin=(924.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.9601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCOJ) + radical(cycloheptane)"""),
)

species(
    label = 'OCC1[CH]C(CC[CH]C1)OO(1655)',
    structure = SMILES('OCC1[CH]C(CC[CH]C1)OO'),
    E0 = (-59.2503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27908,0.0699344,0.000101363,-2.14273e-07,9.71267e-11,-6893.11,45.5961], Tmin=(100,'K'), Tmax=(917.237,'K')), NASAPolynomial(coeffs=[40.3792,0.00858084,4.94053e-06,-1.18453e-09,7.01067e-14,-19596.4,-179.371], Tmin=(917.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.2503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(cycloheptane) + radical(CCJCOOH)"""),
)

species(
    label = 'OCC1[CH][C](CCCC1)OO(1656)',
    structure = SMILES('OCC1[CH][C](CCCC1)OO'),
    E0 = (-41.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31252,0.0719246,9.44417e-05,-2.05842e-07,9.3828e-11,-4791.89,44.356], Tmin=(100,'K'), Tmax=(917.137,'K')), NASAPolynomial(coeffs=[39.7288,0.0100117,4.20751e-06,-1.05355e-09,6.19164e-14,-17244.2,-176.94], Tmin=(917.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(C2CsJOOH)"""),
)

species(
    label = 'OCC1[CH]C(C[CH]CC1)OO(1657)',
    structure = SMILES('OCC1[CH]C(C[CH]CC1)OO'),
    E0 = (-59.2503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27908,0.0699344,0.000101363,-2.14273e-07,9.71267e-11,-6893.11,45.5961], Tmin=(100,'K'), Tmax=(917.237,'K')), NASAPolynomial(coeffs=[40.3792,0.00858084,4.94053e-06,-1.18453e-09,7.01067e-14,-19596.4,-179.371], Tmin=(917.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.2503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(cycloheptane) + radical(CCJCOOH)"""),
)

species(
    label = 'OCC1[CH]C([CH]CCC1)OO(1658)',
    structure = SMILES('OCC1[CH]C([CH]CCC1)OO'),
    E0 = (-28.2642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.42392,0.0714261,0.0001036,-2.21453e-07,1.00962e-10,-3159.48,46.0741], Tmin=(100,'K'), Tmax=(915.229,'K')), NASAPolynomial(coeffs=[42.2987,0.00543746,6.71873e-06,-1.53502e-09,9.42034e-14,-16402.2,-189.613], Tmin=(915.229,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.2642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(CCJCOOH)"""),
)

species(
    label = '[O]OC1[CH]C(CO)CCCC1(1659)',
    structure = SMILES('[O]OC1[CH]C(CO)CCCC1'),
    E0 = (-76.6744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16552,0.0655705,0.000115666,-2.32497e-07,1.0496e-10,-8990.94,43.601], Tmin=(100,'K'), Tmax=(911.723,'K')), NASAPolynomial(coeffs=[41.2709,0.00539778,7.34948e-06,-1.70268e-09,1.07501e-13,-21966.2,-185.93], Tmin=(911.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.6744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(ROOJ)"""),
)

species(
    label = '[CH2]CCCC(C=COO)C[O](1660)',
    structure = SMILES('[CH2]CCCC(C=COO)C[O]'),
    E0 = (67.2115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,200,800,1000,1200,1400,1600],'cm^-1')),
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
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13688,0.119967,-0.000114529,6.31807e-08,-1.49463e-11,8262.83,46.5245], Tmin=(100,'K'), Tmax=(989.549,'K')), NASAPolynomial(coeffs=[12.8172,0.0635601,-2.90247e-05,5.57523e-09,-3.9266e-13,5501.21,-20.649], Tmin=(989.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.2115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]CCCC(C=CC[O])OO(1661)',
    structure = SMILES('[CH2]CCCC(C=CC[O])OO'),
    E0 = (52.6259,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,200,800,1000,1200,1400,1600],'cm^-1')),
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
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.876069,0.116921,-0.000111106,6.39303e-08,-1.63367e-11,6496.8,46.5862], Tmin=(100,'K'), Tmax=(904.853,'K')), NASAPolynomial(coeffs=[9.76777,0.0698689,-3.31067e-05,6.46255e-09,-4.59074e-13,4570.58,-3.69999], Tmin=(904.853,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.6259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(RCCJ)"""),
)

species(
    label = 'OCC1=CC(CCCC1)OO(1662)',
    structure = SMILES('OCC1=CC(CCCC1)OO'),
    E0 = (-316.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01508,0.0551288,0.000156042,-2.79852e-07,1.22633e-10,-37832.2,41.9958], Tmin=(100,'K'), Tmax=(916.182,'K')), NASAPolynomial(coeffs=[44.2714,0.00105368,9.39806e-06,-2.0174e-09,1.2322e-13,-52159,-205.422], Tmin=(916.182,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-316.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cycloheptene)"""),
)

species(
    label = 'O=CC1CCCCC(C1)OO(1663)',
    structure = SMILES('O=CC1CCCCC(C1)OO'),
    E0 = (-348.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05529,0.060214,0.00013506,-2.51802e-07,1.10911e-10,-41719.7,39.6113], Tmin=(100,'K'), Tmax=(918.668,'K')), NASAPolynomial(coeffs=[41.611,0.00681067,6.12104e-06,-1.38555e-09,8.09043e-14,-55144.7,-193.011], Tmin=(918.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-348.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cycloheptane)"""),
)

species(
    label = 'OCC1C=C(CCCC1)OO(1664)',
    structure = SMILES('OCC1C=C(CCCC1)OO'),
    E0 = (-304.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.948545,0.0554191,0.000149925,-2.71369e-07,1.19367e-10,-36411.1,42.0196], Tmin=(100,'K'), Tmax=(914.313,'K')), NASAPolynomial(coeffs=[43.1613,0.00182951,9.17126e-06,-2.00379e-09,1.24185e-13,-50303.2,-198.694], Tmin=(914.313,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-304.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cycloheptene)"""),
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
    label = '[O]CC1CCCCC2OC12(1665)',
    structure = SMILES('[O]CC1CCCCC2OC12'),
    E0 = (-114.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.571335,0.0547928,4.85451e-05,-8.47621e-08,3.07568e-11,-13684.4,29.8271], Tmin=(100,'K'), Tmax=(1073.53,'K')), NASAPolynomial(coeffs=[13.6123,0.0568374,-2.50626e-05,4.88502e-09,-3.51678e-13,-19402.2,-47.6032], Tmin=(1073.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + polycyclic(s2_3_7_ane) + radical(CCOJ)"""),
)

species(
    label = '[O]C[CH]C1CCCCC1OO(1612)',
    structure = SMILES('[O]C[CH]C1CCCCC1OO'),
    E0 = (-33.6555,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01735,0.0954836,-4.0897e-05,-1.96985e-09,3.74851e-12,-3855.19,39.6772], Tmin=(100,'K'), Tmax=(1312.29,'K')), NASAPolynomial(coeffs=[19.9734,0.0567231,-2.54214e-05,4.81396e-09,-3.3394e-13,-11536.1,-75.5705], Tmin=(1312.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.6555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclohexane) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = '[O]CC1[CH]CCCCC1OO(1366)',
    structure = SMILES('[O]CC1[CH]CCCCC1OO'),
    E0 = (-33.9601,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(5082.61,'J/mol'), sigma=(8.64278,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=793.89 K, Pc=17.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.967572,0.0654673,0.000103892,-2.06944e-07,9.16369e-11,-3865.04,43.0331], Tmin=(100,'K'), Tmax=(924.403,'K')), NASAPolynomial(coeffs=[36.5433,0.01611,6.91811e-07,-3.30651e-10,1.00282e-14,-15626.2,-161.091], Tmin=(924.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.9601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(cycloheptane) + radical(CCOJ)"""),
)

species(
    label = '[O]CC1CCCCC1[CH]OO(1666)',
    structure = SMILES('[O]CC1CCCCC1[CH]OO'),
    E0 = (-41.3736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00535,0.100445,-5.62146e-05,1.35281e-08,-1.14389e-12,-4787.6,37.3852], Tmin=(100,'K'), Tmax=(1825.03,'K')), NASAPolynomial(coeffs=[29.0334,0.0455644,-2.01126e-05,3.62982e-09,-2.38577e-13,-17576.5,-130.604], Tmin=(1825.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.3736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(Cyclohexane) + radical(CCsJOOH) + radical(CCOJ)"""),
)

species(
    label = 'OOC1CCCCC2COC21(1370)',
    structure = SMILES('OOC1CCCCC2COC21'),
    E0 = (-401.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.85642,0.128014,-0.000111373,5.05198e-08,-9.42033e-12,-48064.9,4.78736], Tmin=(100,'K'), Tmax=(1258.46,'K')), NASAPolynomial(coeffs=[20.9276,0.0555956,-2.50563e-05,4.79395e-09,-3.36717e-13,-53799.5,-110.371], Tmin=(1258.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-401.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + polycyclic(PolycyclicRing)"""),
)

species(
    label = '[O]CC1CCCCC([O])C1O(1667)',
    structure = SMILES('[O]CC1CCCCC([O])C1O'),
    E0 = (-217.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.990198,0.061139,0.000126191,-2.40015e-07,1.06441e-10,-25874.5,41.4793], Tmin=(100,'K'), Tmax=(915.937,'K')), NASAPolynomial(coeffs=[40.1991,0.00795803,5.7963e-06,-1.3648e-09,8.20592e-14,-38734.5,-182.63], Tmin=(915.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-217.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCOJ) + radical(CC(C)OJ)"""),
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
    label = '[CH2]C1[CH]C(CCCC1)OO(774)',
    structure = SMILES('[CH2]C1[CH]C(CCCC1)OO'),
    E0 = (136.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.706848,0.0513185,0.000152698,-2.72367e-07,1.19945e-10,16603.1,40.2755], Tmin=(100,'K'), Tmax=(909.296,'K')), NASAPolynomial(coeffs=[41.9192,0.000752258,1.02033e-05,-2.26429e-09,1.45335e-13,3189.68,-192.448], Tmin=(909.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(Isobutyl)"""),
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
    E0 = (-2.97398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (130.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (144.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (80.8694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (85.0093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (73.2266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (43.9356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (151.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (129.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (86.0197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (145.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (141.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (141.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (72.3067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (81.1244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (24.8481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (43.7613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (43.7613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (22.7576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (68.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (76.7312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (76.7312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (73.4148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (326.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (94.4075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (79.8219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (60.4262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (60.4262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (21.9993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (42.6693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (156.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (156.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (156.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (5.31034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (98.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (379.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['CH2O(35)', 'OOC1C=CCCCC1(259)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(19)', '[O]CC1=CC(CCCC1)OO(1640)'],
    products = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.72e+07,'cm^3/(mol*s)'), n=1.64, Ea=(9.07928,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2568 used for Cds-CsCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', '[O]CC1C=C(CCCC1)OO(1641)'],
    products = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'O=CC1[CH]C(CCCC1)OO(1642)'],
    products = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.6e+09,'cm^3/(mol*s)'), n=0.935, Ea=(17.4473,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2782 used for CO-CsH_O;HJ
Exact match found for rate rule [CO-CsH_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OOC1C=CCCCC1(259)', '[CH2][O](354)'],
    products = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00429395,'m^3/(mol*s)'), n=2.41179, Ea=(16.3987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-CsH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OOC1[CH][CH]CCCC1(365)', 'CH2O(35)'],
    products = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(225.36,'m^3/(mol*s)'), n=0.996465, Ea=(58.8821,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-HH_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['HO2(26)', '[O]CC1C=CCCCC1(1643)'],
    products = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(46.2,'cm^3/(mol*s)'), n=3.09, Ea=(30.1248,'kJ/mol'), T0=(1,'K'), Tmin=(400,'K'), Tmax=(1100,'K'), comment="""From training reaction 2773 used for Cds-CsH_Cds-CsH;OJ-O2s
Exact match found for rate rule [Cds-CsH_Cds-CsH;OJ-O2s]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['[O]C[C]1CCCCC(C1)OO(1644)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 150 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['[O]CC1CCCC[C](C1)OO(1645)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1557.48,'s^-1'), n=2.79033, Ea=(132.312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['O[CH]C1[CH]C(CCCC1)OO(1646)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(153000,'s^-1'), n=2.26, Ea=(88.9937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 325 used for R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['[O]CC1[CH]CCCC(C1)OO(1647)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.84148e+07,'s^-1'), n=1.625, Ea=(148.323,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)] + [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['[O][CH]C1CCCCC(C1)OO(1648)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.17661e+06,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]CC1CCC[CH]C(C1)OO(1649)'],
    products = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(734000,'s^-1'), n=2.24, Ea=(144.348,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['OC[C]1[CH]C(CCCC1)OO(1650)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['[O]CC1CCCCC(C1)O[O](1651)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 323 used for R4H_SSS;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['OCC1[CH]CCCC([CH]1)OO(1652)'],
    products = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 323 used for R4H_SSS;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['[O]CC1C[CH]CCC(C1)OO(1653)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00228,'s^-1'), n=3.95, Ea=(46.7353,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['[O]CC1CC[CH]CC(C1)OO(1654)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00228,'s^-1'), n=3.95, Ea=(46.7353,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['OCC1[CH]C(CC[CH]C1)OO(1655)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8e+10,'s^-1'), n=0, Ea=(25.7316,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""Estimated using template [R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_CCC;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['OCC1[CH][C](CCCC1)OO(1656)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(57.9,'s^-1'), n=2.9, Ea=(71.128,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_OOH/Cs] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_OOH/Cs]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['OCC1[CH]C(C[CH]CC1)OO(1657)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.12e+09,'s^-1'), n=0, Ea=(79.7052,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6H_SSSSS;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6H_SSSSS;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['OCC1[CH]C([CH]CCC1)OO(1658)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.12e+09,'s^-1'), n=0, Ea=(79.7052,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6HJ_3;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['[O]OC1[CH]C(CO)CCCC1(1659)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(970995,'s^-1'), n=0.905106, Ea=(76.3888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;O_rad_out;XH_out] for rate rule [R7HJ_3;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['OOC1[CH][CH]CCCC1(365)', '[CH2][O](354)'],
    products = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.35017e+06,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]CCCC(C=COO)C[O](1660)'],
    products = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_HNd;radadd_intra_cs2H] for rate rule [R7_linear;doublebond_intra_pri_HNd_O;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]CCCC(C=CC[O])OO(1661)'],
    products = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 49 used for R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['OCC1=CC(CCCC1)OO(1662)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['O=CC1CCCCC(C1)OO(1663)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['OCC1C=C(CCCC1)OO(1664)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['OH(21)', '[O]CC1CCCCC2OC12(1665)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.38e+12,'s^-1','*|/',1.2), n=0, Ea=(45.6433,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 3 used for R2OO_S;C_rad/H/NonDeC_intra;OOH
Exact match found for rate rule [R2OO_S;C_rad/H/NonDeC_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['[O]C[CH]C1CCCCC1OO(1612)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['[O]CC1[CH]CCCCC1OO(1366)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['[O]CC1CCCCC1[CH]OO(1666)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['OOC1CCCCC2COC21(1370)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    products = ['[O]CC1CCCCC([O])C1O(1667)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.69e+11,'s^-1'), n=0, Ea=(101.671,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for R2OOH_S;C_rad_out_H/NonDeC
Exact match found for rate rule [R2OOH_S;C_rad_out_H/NonDeC]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(20)', '[CH2]C1[CH]C(CCCC1)OO(774)'],
    products = ['[O]CC1[CH]C(CCCC1)OO(1365)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/Cs;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '255',
    isomers = [
        '[O]CC1[CH]C(CCCC1)OO(1365)',
    ],
    reactants = [
        ('CH2O(35)', 'OOC1C=CCCCC1(259)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '255',
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

