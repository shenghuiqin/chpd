species(
    label = 'CCC(O)O[O](11243)',
    structure = SMILES('CCC(O)O[O]'),
    E0 = (-257.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,492.5,1135,1000,3615,1277.5,1000,278.968],'cm^-1')),
        HinderedRotor(inertia=(0.13352,'amu*angstrom^2'), symmetry=1, barrier=(7.37361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133519,'amu*angstrom^2'), symmetry=1, barrier=(7.37361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13352,'amu*angstrom^2'), symmetry=1, barrier=(7.37361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133519,'amu*angstrom^2'), symmetry=1, barrier=(7.37361,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0859,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5709,0.0564872,-5.56998e-05,3.20723e-08,-7.87793e-12,-30833.8,23.4723], Tmin=(100,'K'), Tmax=(959.77,'K')), NASAPolynomial(coeffs=[8.107,0.0292471,-1.31273e-05,2.50124e-09,-1.75337e-13,-32088.4,-7.79233], Tmin=(959.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-257.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ)"""),
)

species(
    label = 'C2H5CHO(70)',
    structure = SMILES('CCC=O'),
    E0 = (-204.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.207559,'amu*angstrom^2'), symmetry=1, barrier=(4.77219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208362,'amu*angstrom^2'), symmetry=1, barrier=(4.79065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90578,0.0240644,-7.06356e-06,-9.81837e-10,5.55825e-13,-24535.9,13.5806], Tmin=(100,'K'), Tmax=(1712.49,'K')), NASAPolynomial(coeffs=[7.69109,0.0189242,-7.84934e-06,1.38273e-09,-8.99057e-14,-27060.1,-14.6647], Tmin=(1712.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""propanal""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'HO2(10)',
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
    label = 'CH3OH(22)',
    structure = SMILES('CO'),
    E0 = (-211.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1964.49,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0117392,'amu*angstrom^2'), symmetry=1, barrier=(6.09573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.0419,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4005.93,'J/mol'), sigma=(3.626,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.65851,-0.0162983,6.91938e-05,-7.58373e-08,2.80428e-11,-25612,-0.897331], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.52727,0.0103179,-3.62893e-06,5.77448e-10,-3.42183e-14,-26002.9,5.16759], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-211.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH3OH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CH2CHOO(46)',
    structure = SMILES('C=CO[O]'),
    E0 = (92.9866,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,492.5,1135,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.968919,'amu*angstrom^2'), symmetry=1, barrier=(22.2774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3284.22,'J/mol'), sigma=(4.037,'angstroms'), dipoleMoment=(1.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41461,0.0297546,-2.58394e-05,1.12221e-08,-1.86781e-12,11371.2,19.0226], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[1.41603,0.0299337,-2.63852e-05,1.17707e-08,-2.05116e-12,11361.7,18.9695], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), E0=(92.9866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHOO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CC[C](O)OO(11591)',
    structure = SMILES('CC[C](O)OO'),
    E0 = (-203.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (91.0859,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40022,0.0601665,-5.92307e-05,3.21473e-08,-7.30907e-12,-24424,24.0127], Tmin=(100,'K'), Tmax=(1038.01,'K')), NASAPolynomial(coeffs=[9.67581,0.0282763,-1.3147e-05,2.54985e-09,-1.80648e-13,-26142.1,-16.2209], Tmin=(1038.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(Cs_P)"""),
)

species(
    label = 'CCC([O])OO(11242)',
    structure = SMILES('CCC([O])OO'),
    E0 = (-183.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,180,670.234],'cm^-1')),
        HinderedRotor(inertia=(0.168006,'amu*angstrom^2'), symmetry=1, barrier=(3.86279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167914,'amu*angstrom^2'), symmetry=1, barrier=(3.86067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167055,'amu*angstrom^2'), symmetry=1, barrier=(3.84092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97004,'amu*angstrom^2'), symmetry=1, barrier=(45.2952,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0859,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24776,0.0472998,-5.47295e-06,-8.50833e-08,8.92061e-11,-22000.6,21.3121], Tmin=(100,'K'), Tmax=(440.495,'K')), NASAPolynomial(coeffs=[4.76732,0.0377872,-1.85973e-05,3.66781e-09,-2.61102e-13,-22352.3,9.75042], Tmin=(440.495,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ)"""),
)

species(
    label = 'C[CH]C(O)OO(11592)',
    structure = SMILES('C[CH]C(O)OO'),
    E0 = (-208.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3615,1310,387.5,850,1000,3615,1277.5,1000,357.26],'cm^-1')),
        HinderedRotor(inertia=(0.110601,'amu*angstrom^2'), symmetry=1, barrier=(10.0174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110601,'amu*angstrom^2'), symmetry=1, barrier=(10.0174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310749,'amu*angstrom^2'), symmetry=1, barrier=(28.1451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0028879,'amu*angstrom^2'), symmetry=1, barrier=(10.0174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310748,'amu*angstrom^2'), symmetry=1, barrier=(28.1451,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0859,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39179,0.0613508,-6.40226e-05,3.78556e-08,-9.44386e-12,-25005.7,25.6646], Tmin=(100,'K'), Tmax=(949.447,'K')), NASAPolynomial(coeffs=[8.908,0.0296842,-1.39923e-05,2.7251e-09,-1.9335e-13,-26432.9,-10.2066], Tmin=(949.447,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]CC(O)OO(11593)',
    structure = SMILES('[CH2]CC(O)OO'),
    E0 = (-203.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (91.0859,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40022,0.0601665,-5.92307e-05,3.21473e-08,-7.30907e-12,-24424,25.1113], Tmin=(100,'K'), Tmax=(1038.01,'K')), NASAPolynomial(coeffs=[9.67581,0.0282763,-1.3147e-05,2.54985e-09,-1.80648e-13,-26142.1,-15.1223], Tmin=(1038.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(RCCJ)"""),
)

species(
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1483.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53764,-0.00122828,5.36759e-06,-4.93128e-09,1.45955e-12,-1037.99,4.6718], Tmin=(100,'K'), Tmax=(1087.71,'K')), NASAPolynomial(coeffs=[3.16427,0.00169454,-8.00335e-07,1.5903e-10,-1.14891e-14,-1048.45,6.08303], Tmin=(1087.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'HOCHC2H5(1012)',
    structure = SMILES('CC[CH]O'),
    E0 = (-91.8987,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.00178298,'amu*angstrom^2'), symmetry=1, barrier=(20.2441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250308,'amu*angstrom^2'), symmetry=1, barrier=(5.75507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0237095,'amu*angstrom^2'), symmetry=1, barrier=(20.299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.0871,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50157,0.0263381,4.92316e-06,-2.16134e-08,9.24954e-12,-10993.1,16.8679], Tmin=(100,'K'), Tmax=(1011.24,'K')), NASAPolynomial(coeffs=[8.06358,0.0195007,-7.42704e-06,1.35669e-09,-9.50157e-14,-12893.3,-13.8611], Tmin=(1011.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.8987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""HOCHC2H5""", comment="""Thermo library: CBS_QB3_1dHR"""),
)

species(
    label = 'CC[CH]O[O](238)',
    structure = SMILES('CC[CH]O[O]'),
    E0 = (139.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,208.348],'cm^-1')),
        HinderedRotor(inertia=(0.243143,'amu*angstrom^2'), symmetry=1, barrier=(7.00807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245015,'amu*angstrom^2'), symmetry=1, barrier=(7.01593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00088615,'amu*angstrom^2'), symmetry=1, barrier=(7.02015,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20785,0.0420296,-4.3379e-05,3.01324e-08,-9.24759e-12,16864.4,19.4406], Tmin=(100,'K'), Tmax=(771.405,'K')), NASAPolynomial(coeffs=[5.43298,0.0253065,-1.08618e-05,2.03087e-09,-1.40556e-13,16366.8,4.71818], Tmin=(771.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCsJOOH) + radical(ROOJ)"""),
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
    label = 'CCC([O])O[O](10578)',
    structure = SMILES('CCC([O])O[O]'),
    E0 = (-31.3675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,492.5,1135,1000,1380,1390,370,380,2900,435,180,1976.03],'cm^-1')),
        HinderedRotor(inertia=(0.279223,'amu*angstrom^2'), symmetry=1, barrier=(6.41989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279209,'amu*angstrom^2'), symmetry=1, barrier=(6.41956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279357,'amu*angstrom^2'), symmetry=1, barrier=(6.42298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0779,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71392,0.0565161,-7.87295e-05,7.28942e-08,-2.76894e-11,-3696.37,23.8074], Tmin=(100,'K'), Tmax=(798.825,'K')), NASAPolynomial(coeffs=[3.64684,0.0352595,-1.70745e-05,3.29576e-09,-2.29714e-13,-3635.78,17.2285], Tmin=(798.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.3675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = 'C2H5(29)',
    structure = SMILES('C[CH2]'),
    E0 = (107.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1190.6,1642.82,1642.96,3622.23,3622.39],'cm^-1')),
        HinderedRotor(inertia=(0.866817,'amu*angstrom^2'), symmetry=1, barrier=(19.9298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24186,-0.00356905,4.82667e-05,-5.85401e-08,2.25805e-11,12969,4.44704], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.32196,0.0123931,-4.39681e-06,7.0352e-10,-4.18435e-14,12175.9,0.171104], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(107.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H5""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[O]O[CH]O(4538)',
    structure = SMILES('[O]O[CH]O'),
    E0 = (9.28844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,492.5,1135,1000],'cm^-1')),
        HinderedRotor(inertia=(0.156407,'amu*angstrom^2'), symmetry=1, barrier=(3.59612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156804,'amu*angstrom^2'), symmetry=1, barrier=(3.60523,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58857,0.0452524,-0.00011197,1.22556e-07,-4.65581e-11,1154.51,14.2974], Tmin=(100,'K'), Tmax=(891.504,'K')), NASAPolynomial(coeffs=[-0.0459,0.0212294,-1.12418e-05,2.13324e-09,-1.41938e-13,3048.61,34.6933], Tmin=(891.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.28844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = '[CH2]C(O)O[O](159)',
    structure = SMILES('[CH2]C(O)O[O]'),
    E0 = (-19.3298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,492.5,1135,1000,3000,3100,440,815,1455,1000,3615,1277.5,1000],'cm^-1')),
        HinderedRotor(inertia=(0.143819,'amu*angstrom^2'), symmetry=1, barrier=(6.57764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143926,'amu*angstrom^2'), symmetry=1, barrier=(6.57118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14444,'amu*angstrom^2'), symmetry=1, barrier=(6.57847,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3937.33,'J/mol'), sigma=(6.43989,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=615.00 K, Pc=33.45 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0219,0.0460188,-6.08592e-05,4.45311e-08,-1.31672e-11,-2255.9,21.5986], Tmin=(100,'K'), Tmax=(824.704,'K')), NASAPolynomial(coeffs=[8.06888,0.016691,-7.51954e-06,1.41494e-09,-9.77081e-14,-3253.34,-6.40945], Tmin=(824.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.3298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CJCOOH)"""),
)

species(
    label = 'CH3(17)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C[CH]C(O)O[O](11594)',
    structure = SMILES('C[CH]C(O)O[O]'),
    E0 = (-56.6576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3615,1277.5,1000,492.5,1135,1000,240.519],'cm^-1')),
        HinderedRotor(inertia=(0.124986,'amu*angstrom^2'), symmetry=1, barrier=(5.13216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124982,'amu*angstrom^2'), symmetry=1, barrier=(5.13208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12499,'amu*angstrom^2'), symmetry=1, barrier=(5.13212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125009,'amu*angstrom^2'), symmetry=1, barrier=(5.13221,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0779,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43509,0.060511,-7.91265e-05,6.19411e-08,-2.018e-11,-6725.83,26.2582], Tmin=(100,'K'), Tmax=(755.497,'K')), NASAPolynomial(coeffs=[7.53848,0.0276459,-1.27814e-05,2.43219e-09,-1.68881e-13,-7632.34,-1.37182], Tmin=(755.497,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.6576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCJCOOH) + radical(ROOJ)"""),
)

species(
    label = 'CC[C](O)O[O](11595)',
    structure = SMILES('CC[C](O)O[O]'),
    E0 = (-51.8262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3615,1277.5,1000,492.5,1135,1000,201.998],'cm^-1')),
        HinderedRotor(inertia=(0.180652,'amu*angstrom^2'), symmetry=1, barrier=(5.22495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180645,'amu*angstrom^2'), symmetry=1, barrier=(5.22506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18044,'amu*angstrom^2'), symmetry=1, barrier=(5.22507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180731,'amu*angstrom^2'), symmetry=1, barrier=(5.22458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0779,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52255,0.0582352,-6.97416e-05,4.91052e-08,-1.44428e-11,-6147.41,24.3332], Tmin=(100,'K'), Tmax=(818.557,'K')), NASAPolynomial(coeffs=[7.93242,0.026912,-1.23415e-05,2.35565e-09,-1.64554e-13,-7196.77,-5.30729], Tmin=(818.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.8262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_P)"""),
)

species(
    label = '[CH2]CC(O)O[O](11596)',
    structure = SMILES('[CH2]CC(O)O[O]'),
    E0 = (-51.8262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3615,1277.5,1000,3000,3100,440,815,1455,1000,492.5,1135,1000,199.898],'cm^-1')),
        HinderedRotor(inertia=(0.00421224,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179958,'amu*angstrom^2'), symmetry=1, barrier=(5.10971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180121,'amu*angstrom^2'), symmetry=1, barrier=(5.10985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180005,'amu*angstrom^2'), symmetry=1, barrier=(5.10984,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0779,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52255,0.0582352,-6.97416e-05,4.91052e-08,-1.44428e-11,-6147.41,25.4318], Tmin=(100,'K'), Tmax=(818.557,'K')), NASAPolynomial(coeffs=[7.93242,0.026912,-1.23415e-05,2.35565e-09,-1.64554e-13,-7196.77,-4.20868], Tmin=(818.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.8262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(RCCJ) + radical(ROOJ)"""),
)

species(
    label = 'H2O(7)',
    structure = SMILES('O'),
    E0 = (-251.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1591.57,3573.8,3573.81],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.0153,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(6727.26,'J/mol'), sigma=(2.641,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.99882,-0.000554831,2.76774e-06,-1.55665e-09,3.0233e-13,-30274.6,-0.0308942], Tmin=(100,'K'), Tmax=(1281.44,'K')), NASAPolynomial(coeffs=[3.1956,0.0019524,-1.67118e-07,-2.97938e-11,4.45139e-15,-30068.7,4.04335], Tmin=(1281.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'propen1peroxy(5702)',
    structure = SMILES('CC=CO[O]'),
    E0 = (62.3012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.791668,'amu*angstrom^2'), symmetry=1, barrier=(18.202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788245,'amu*angstrom^2'), symmetry=1, barrier=(18.1233,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05109,0.0384282,-2.89472e-05,1.11693e-08,-1.74495e-12,7566.86,16.6482], Tmin=(100,'K'), Tmax=(1503.71,'K')), NASAPolynomial(coeffs=[10.5824,0.0157341,-6.30895e-06,1.13262e-09,-7.62989e-14,5001.14,-27.9908], Tmin=(1503.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.3012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""propen1peroxy""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'CC(O)O[O](689)',
    structure = SMILES('CC(O)O[O]'),
    E0 = (-233.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3615,1277.5,1000,492.5,1135,1000],'cm^-1')),
        HinderedRotor(inertia=(0.349917,'amu*angstrom^2'), symmetry=1, barrier=(8.04527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.350493,'amu*angstrom^2'), symmetry=1, barrier=(8.05852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.349449,'amu*angstrom^2'), symmetry=1, barrier=(8.03451,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0593,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3937.33,'J/mol'), sigma=(6.43989,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=615.00 K, Pc=33.45 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1881,0.0420711,-4.32163e-05,2.58456e-08,-6.53772e-12,-27995.2,19.0097], Tmin=(100,'K'), Tmax=(938.25,'K')), NASAPolynomial(coeffs=[7.15873,0.0208806,-9.33938e-06,1.77516e-09,-1.2422e-13,-28928,-4.65405], Tmin=(938.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ)"""),
)

species(
    label = 'propen1ol(10572)',
    structure = SMILES('CC=CO'),
    E0 = (-166.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.635471,'amu*angstrom^2'), symmetry=1, barrier=(14.6107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.638617,'amu*angstrom^2'), symmetry=1, barrier=(14.6831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43518,0.0321429,-2.11365e-05,7.43426e-09,-1.09824e-12,-20020.5,13.0756], Tmin=(100,'K'), Tmax=(1529.28,'K')), NASAPolynomial(coeffs=[8.00322,0.0175792,-6.85177e-06,1.2071e-09,-8.02562e-14,-21723.5,-16.1523], Tmin=(1529.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""propen1ol""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O(4)',
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
    label = 'CCC([O])O(11590)',
    structure = SMILES('CCC([O])O'),
    E0 = (-254.877,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,210.647,2506.05],'cm^-1')),
        HinderedRotor(inertia=(0.0950942,'amu*angstrom^2'), symmetry=1, barrier=(2.99275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0950896,'amu*angstrom^2'), symmetry=1, barrier=(2.9916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.35528,'amu*angstrom^2'), symmetry=1, barrier=(11.1764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0865,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88177,0.0309076,2.40949e-05,-1.14402e-07,1.04569e-10,-30621.1,18.2709], Tmin=(100,'K'), Tmax=(404.624,'K')), NASAPolynomial(coeffs=[3.529,0.0330502,-1.55109e-05,3.02197e-09,-2.14668e-13,-30743.4,14.87], Tmin=(404.624,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-254.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ)"""),
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
    E0 = (-133.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-102.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-77.7307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-124.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-169.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-100.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (168.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (180.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (117.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (116.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (155.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (159.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (159.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (72.4982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (210.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (186.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-133.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-11.8724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CCC(O)O[O](11243)'],
    products = ['C2H5CHO(70)', 'HO2(10)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(6.813e+10,'s^-1','*|/',10), n=0.493, Ea=(123.219,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 13 used for R2OO_O_HNd
Exact match found for rate rule [R2OO_O_HNd]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CCC(O)O[O](11243)'],
    products = ['CC[C](O)OO(11591)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3e+08,'s^-1'), n=1.23, Ea=(154.18,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 270 used for R3H_SS_O;O_rad_out;Cs_H_out_NDMustO
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CCC([O])OO(11242)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(94.0113,'s^-1'), n=2.81534, Ea=(105.641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C[CH]C(O)OO(11592)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 323 used for R4H_SSS;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CCC(O)O[O](11243)'],
    products = ['[CH2]CC(O)OO(11593)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.07e+06,'s^-1'), n=1.55, Ea=(87.9477,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 255 used for R5H_SSSS_OCC;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R5H_SSSS_OCC;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O2(2)', 'HOCHC2H5(1012)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(8.36e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for C_rad/H/CsO;O2_birad
Exact match found for rate rule [C_rad/H/CsO;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC[CH]O[O](238)', 'OH(5)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/CsO] for rate rule [O_pri_rad;C_rad/H/CsO]
Euclidian distance = 2.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'CCC([O])O[O](10578)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C2H5(29)', '[O]O[CH]O(4538)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.15e+14,'cm^3/(mol*s)','*|/',2), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [C_rad/H/NonDe;C_rad/H2/Cs] for rate rule [C_rad/H/O2;C_rad/H2/Cs]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)O[O](159)', 'CH3(17)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', 'C[CH]C(O)O[O](11594)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', 'CC[C](O)O[O](11595)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.68468e+06,'m^3/(mol*s)'), n=-0.0135, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/NonDe;Y_rad] for rate rule [C_rad/NonDeCO;H_rad]
Euclidian distance = 2.2360679775
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH2]CC(O)O[O](11596)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['H2O(7)', 'propen1peroxy(5702)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.238e-05,'m^3/(mol*s)'), n=3.25, Ea=(261.918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd/monosub_Cd/monosub;H_OH] for rate rule [Cd/H/Nd_Cd/H/Nd;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CH3OH(22)', 'CH2CHOO(46)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.79e-05,'cm^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cd_Cd;CH3OH] for rate rule [Cd/H2_Cd/H/Nd;CH3OH]
Euclidian distance = 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CH2(S)(23)', 'CC(O)O[O](689)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['CCC(O)O[O](11243)'],
    products = ['HO2(10)', 'propen1ol(10572)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.3e+09,'s^-1'), n=1.01, Ea=(123.846,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 6 used for R2OO_HNd_HNd
Exact match found for rate rule [R2OO_HNd_HNd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O(4)', 'CCC([O])O(11590)'],
    products = ['CCC(O)O[O](11243)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '2060',
    isomers = [
        'CCC(O)O[O](11243)',
    ],
    reactants = [
        ('C2H5CHO(70)', 'HO2(10)'),
        ('CH3OH(22)', 'CH2CHOO(46)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '2060',
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

