species(
    label = 'O[C]=CC[CH]O(1255)',
    structure = SMILES('O[C]=CC[CH]O'),
    E0 = (34.9002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,1685,370,3010,987.5,1337.5,450,1655,311.231,311.254],'cm^-1')),
        HinderedRotor(inertia=(0.181118,'amu*angstrom^2'), symmetry=1, barrier=(12.4497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181105,'amu*angstrom^2'), symmetry=1, barrier=(12.4497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181108,'amu*angstrom^2'), symmetry=1, barrier=(12.4497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181104,'amu*angstrom^2'), symmetry=1, barrier=(12.4496,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.86678,0.0658295,-7.55763e-05,4.25469e-08,-8.80929e-12,4313.16,26.5646], Tmin=(100,'K'), Tmax=(926.301,'K')), NASAPolynomial(coeffs=[14.768,0.0138416,-4.41148e-06,7.00613e-10,-4.46589e-14,1392.83,-41.2987], Tmin=(926.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.9002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CCsJOH)"""),
)

species(
    label = 'HCCOH(50)',
    structure = SMILES('C#CO'),
    E0 = (80.0402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3615,1277.5,1000,287.694,287.7,1588.42],'cm^-1')),
        HinderedRotor(inertia=(0.002037,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05541,0.0252003,-3.80822e-05,3.09891e-08,-9.898e-12,9768.72,12.2272], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.3751,0.00549429,-1.88137e-06,2.93804e-10,-1.71772e-14,8932.78,-8.24498], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(80.0402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CH2CHOH(42)',
    structure = SMILES('C=CO'),
    E0 = (-138.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.72808,'amu*angstrom^2'), symmetry=1, barrier=(39.7321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28758,0.0197013,1.96383e-06,-1.9439e-08,1.02617e-11,-16537.3,14.1333], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[7.49818,0.0103957,-3.66891e-06,5.85206e-10,-3.47374e-14,-18164.3,-13.8388], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-138.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'O=CCC=[C]O(8483)',
    structure = SMILES('O=CCC=[C]O'),
    E0 = (-62.1872,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1685,370,3615,1277.5,1000,3010,987.5,1337.5,450,1655,296.192],'cm^-1')),
        HinderedRotor(inertia=(0.285524,'amu*angstrom^2'), symmetry=1, barrier=(17.7584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284403,'amu*angstrom^2'), symmetry=1, barrier=(17.7551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285329,'amu*angstrom^2'), symmetry=1, barrier=(17.7607,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72707,0.0389403,-6.5296e-06,-2.17987e-08,1.11795e-11,-7387.92,22.8136], Tmin=(100,'K'), Tmax=(1027.39,'K')), NASAPolynomial(coeffs=[14.3933,0.0141827,-6.23627e-06,1.27579e-09,-9.64687e-14,-11286.5,-44.9429], Tmin=(1027.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.1872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=CJO)"""),
)

species(
    label = 'O=C=CC[CH]O(5758)',
    structure = SMILES('O=C=CC[CH]O'),
    E0 = (-81.213,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.584872,'amu*angstrom^2'), symmetry=1, barrier=(13.4474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.584038,'amu*angstrom^2'), symmetry=1, barrier=(13.4282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00431203,'amu*angstrom^2'), symmetry=1, barrier=(13.4545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46062,0.0594856,-8.19004e-05,6.44597e-08,-2.06057e-11,-9679.67,22.1835], Tmin=(100,'K'), Tmax=(799.065,'K')), NASAPolynomial(coeffs=[8.30796,0.0233336,-1.05159e-05,1.96607e-09,-1.34688e-13,-10714.1,-8.9404], Tmin=(799.065,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH) + radical(CCsJOH)"""),
)

species(
    label = 'O[C]=CC=CO(8484)',
    structure = SMILES('O[C]=CC=CO'),
    E0 = (-83.9199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3580,3650,1210,1345,900,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.09553,'amu*angstrom^2'), symmetry=1, barrier=(25.1883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09555,'amu*angstrom^2'), symmetry=1, barrier=(25.1888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09679,'amu*angstrom^2'), symmetry=1, barrier=(25.2175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803519,0.0503097,1.529e-06,-6.70039e-08,3.84255e-11,-9959.14,22.3195], Tmin=(100,'K'), Tmax=(888.056,'K')), NASAPolynomial(coeffs=[26.7325,-0.0104688,9.5814e-06,-2.02664e-09,1.39785e-13,-16773.1,-112.13], Tmin=(888.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.9199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = 'OC#CC[CH]O(8485)',
    structure = SMILES('OC#CC[CH]O'),
    E0 = (28.5863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2100,2250,500,550,3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,317.4,317.444],'cm^-1')),
        HinderedRotor(inertia=(0.137625,'amu*angstrom^2'), symmetry=1, barrier=(9.84344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00167338,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174558,'amu*angstrom^2'), symmetry=1, barrier=(12.4795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.911408,'amu*angstrom^2'), symmetry=1, barrier=(65.2147,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38122,0.0600055,-8.12374e-05,6.09271e-08,-1.82154e-11,3530.15,22.9118], Tmin=(100,'K'), Tmax=(866.862,'K')), NASAPolynomial(coeffs=[9.42301,0.0203721,-8.28598e-06,1.46201e-09,-9.64606e-14,2230.83,-14.1888], Tmin=(866.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.5863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CtH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=[C]O(172)',
    structure = SMILES('[CH]=[C]O'),
    E0 = (348.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,1058.91,1059.8],'cm^-1')),
        HinderedRotor(inertia=(0.315045,'amu*angstrom^2'), symmetry=1, barrier=(7.24351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22339,0.0198737,-3.18466e-05,3.10436e-08,-1.17427e-11,41917.3,13.7606], Tmin=(100,'K'), Tmax=(829.31,'K')), NASAPolynomial(coeffs=[3.78101,0.0111997,-5.33323e-06,1.02849e-09,-7.13893e-14,42030.6,12.4155], Tmin=(829.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH]O(284)',
    structure = SMILES('[CH2][CH]O'),
    E0 = (143.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.217851,'amu*angstrom^2'), symmetry=1, barrier=(5.00882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0197382,'amu*angstrom^2'), symmetry=1, barrier=(14.867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84769,0.0229973,-1.85068e-05,7.77211e-09,-1.30725e-12,17300.5,13.0245], Tmin=(100,'K'), Tmax=(1418.34,'K')), NASAPolynomial(coeffs=[7.97636,0.00853337,-3.21015e-06,5.82141e-10,-3.99268e-14,15845.7,-13.5107], Tmin=(1418.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CJCO) + radical(CCsJOH)"""),
)

species(
    label = '[O]CCC=[C]O(8486)',
    structure = SMILES('[O]CCC=[C]O'),
    E0 = (80.3074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6172,0.0528827,-5.12501e-05,2.78003e-08,-6.23906e-12,9744.18,25.1865], Tmin=(100,'K'), Tmax=(1062.72,'K')), NASAPolynomial(coeffs=[9.41615,0.023528,-9.81653e-06,1.80808e-09,-1.24481e-13,8086.57,-12.9134], Tmin=(1062.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.3074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=CJO)"""),
)

species(
    label = 'O=C[CH]C[CH]O(5761)',
    structure = SMILES('[O]C=CC[CH]O'),
    E0 = (-63.3814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,333.351,333.353,333.355],'cm^-1')),
        HinderedRotor(inertia=(0.200499,'amu*angstrom^2'), symmetry=1, barrier=(15.8107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200498,'amu*angstrom^2'), symmetry=1, barrier=(15.8106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200499,'amu*angstrom^2'), symmetry=1, barrier=(15.8106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.663175,0.0617871,-6.12456e-05,3.00647e-08,-5.64521e-12,-7493.01,25.5823], Tmin=(100,'K'), Tmax=(1420.22,'K')), NASAPolynomial(coeffs=[17.0256,0.0111206,-2.89316e-06,4.01608e-10,-2.37239e-14,-11678.5,-57.4699], Tmin=(1420.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.3814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCsJOH) + radical(C=COJ)"""),
)

species(
    label = 'O[C]=C[CH]CO(8487)',
    structure = SMILES('O[C]=C[CH]CO'),
    E0 = (-28.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.509826,0.0658157,-6.86097e-05,3.51964e-08,-6.88772e-12,-3290.48,25.4976], Tmin=(100,'K'), Tmax=(1358.65,'K')), NASAPolynomial(coeffs=[17.7086,0.0106142,-2.62361e-06,3.44349e-10,-1.95331e-14,-7542.43,-61.1968], Tmin=(1358.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJCO) + radical(C=CJO)"""),
)

species(
    label = 'O[CH]C[C]=CO(8488)',
    structure = SMILES('O[CH]C[C]=CO'),
    E0 = (32.9978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.687318,'amu*angstrom^2'), symmetry=1, barrier=(15.8028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.687614,'amu*angstrom^2'), symmetry=1, barrier=(15.8096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.687447,'amu*angstrom^2'), symmetry=1, barrier=(15.8058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.686798,'amu*angstrom^2'), symmetry=1, barrier=(15.7908,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.64352,0.0684197,-7.24963e-05,3.15807e-08,-2.84543e-12,4094.68,24.7124], Tmin=(100,'K'), Tmax=(894.972,'K')), NASAPolynomial(coeffs=[17.6777,0.00928822,-1.88503e-06,2.08199e-10,-1.11281e-14,364.762,-59.3817], Tmin=(894.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.9978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = 'O[C]=[C]CCO(8489)',
    structure = SMILES('O[C]=[C]CCO'),
    E0 = (92.4441,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3580,3650,1210,1345,900,1100,1670,1700,300,440,322.313,322.455],'cm^-1')),
        HinderedRotor(inertia=(0.148009,'amu*angstrom^2'), symmetry=1, barrier=(10.8936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146726,'amu*angstrom^2'), symmetry=1, barrier=(10.8994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147223,'amu*angstrom^2'), symmetry=1, barrier=(10.8985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147939,'amu*angstrom^2'), symmetry=1, barrier=(10.8932,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05278,0.0622997,-7.1373e-05,4.25191e-08,-9.91969e-12,11226.8,26.8205], Tmin=(100,'K'), Tmax=(1054.61,'K')), NASAPolynomial(coeffs=[13.2492,0.0160405,-5.57748e-06,9.26856e-10,-6.01058e-14,8654.27,-32.6686], Tmin=(1054.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.4441,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'O[CH][CH]C=CO(8490)',
    structure = SMILES('O[CH]C=C[CH]O'),
    E0 = (-101.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80733,0.0328991,2.45073e-05,-6.10662e-08,2.71877e-11,-12124.7,24.8527], Tmin=(100,'K'), Tmax=(952.542,'K')), NASAPolynomial(coeffs=[15.7273,0.0126021,-3.61763e-06,6.71966e-10,-5.31136e-14,-16507.7,-50.7131], Tmin=(952.542,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJO) + radical(C=CCJO)"""),
)

species(
    label = 'O=CC[CH][CH]O(5763)',
    structure = SMILES('O=CC[CH][CH]O'),
    E0 = (-6.73096,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,348.246,1839.05],'cm^-1')),
        HinderedRotor(inertia=(0.00138966,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928456,'amu*angstrom^2'), symmetry=1, barrier=(7.99073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928395,'amu*angstrom^2'), symmetry=1, barrier=(7.99069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.092824,'amu*angstrom^2'), symmetry=1, barrier=(7.99058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42645,0.0612384,-8.47733e-05,7.03001e-08,-2.37783e-11,-721.324,25.6232], Tmin=(100,'K'), Tmax=(813.775,'K')), NASAPolynomial(coeffs=[6.94305,0.0281764,-1.28716e-05,2.41773e-09,-1.6585e-13,-1422.3,1.35524], Tmin=(813.775,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.73096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCJCO) + radical(CCsJOH)"""),
)

species(
    label = 'O=[C][CH]CCO(5762)',
    structure = SMILES('O=[C][CH]CCO'),
    E0 = (-59.4409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61391,0.0540039,-5.83541e-05,3.68306e-08,-9.62936e-12,-7064.5,23.4882], Tmin=(100,'K'), Tmax=(921.149,'K')), NASAPolynomial(coeffs=[8.45852,0.0242815,-9.95344e-06,1.80093e-09,-1.22187e-13,-8325.46,-8.97079], Tmin=(921.149,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.4409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'OC=CC=CO(8491)',
    structure = SMILES('OC=CC=CO'),
    E0 = (-323.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.583761,0.0486673,2.80231e-05,-1.05328e-07,5.44794e-11,-38779.7,19.2865], Tmin=(100,'K'), Tmax=(892.479,'K')), NASAPolynomial(coeffs=[30.9038,-0.0146483,1.24596e-05,-2.58582e-09,1.76132e-13,-47082.1,-139.734], Tmin=(892.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-323.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH)"""),
)

species(
    label = 'O=CCC=CO(5765)',
    structure = SMILES('O=CCC=CO'),
    E0 = (-301.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52612,0.037118,2.036e-05,-6.02228e-08,2.70604e-11,-36209.3,20.4039], Tmin=(100,'K'), Tmax=(972.884,'K')), NASAPolynomial(coeffs=[18.2525,0.010519,-3.64943e-06,7.84402e-10,-6.56802e-14,-41459.6,-70.0882], Tmin=(972.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-301.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'O=C=CCCO(5766)',
    structure = SMILES('O=C=CCCO'),
    E0 = (-261.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86998,0.0490307,-4.43395e-05,2.32664e-08,-5.24646e-12,-31377.6,21.1693], Tmin=(100,'K'), Tmax=(1033.74,'K')), NASAPolynomial(coeffs=[7.68523,0.0265292,-1.16894e-05,2.21036e-09,-1.54347e-13,-32579.9,-7.07903], Tmin=(1033.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-261.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = '[CH2]C(O)C=[C]O(1257)',
    structure = SMILES('[CH2]C(O)C=[C]O'),
    E0 = (53.2121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,340.613],'cm^-1')),
        HinderedRotor(inertia=(0.167235,'amu*angstrom^2'), symmetry=1, barrier=(13.7673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167187,'amu*angstrom^2'), symmetry=1, barrier=(13.7671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167209,'amu*angstrom^2'), symmetry=1, barrier=(13.7671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167247,'amu*angstrom^2'), symmetry=1, barrier=(13.7672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.372397,0.0645561,-6.7258e-05,3.36046e-08,-6.27118e-12,6543.59,29.9932], Tmin=(100,'K'), Tmax=(1500.67,'K')), NASAPolynomial(coeffs=[18.5905,0.0069173,-5.70202e-07,-5.26759e-11,7.46619e-15,2097.98,-61.8877], Tmin=(1500.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.2121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CJCO)"""),
)

species(
    label = 'OC1=CCC1O(8480)',
    structure = SMILES('OC1=CCC1O'),
    E0 = (-242.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70002,0.0298393,4.56077e-05,-9.38361e-08,4.18318e-11,-29111.9,20.1871], Tmin=(100,'K'), Tmax=(928.648,'K')), NASAPolynomial(coeffs=[20.1069,0.00386597,1.45045e-06,-3.18036e-10,1.40431e-14,-34829.3,-79.6294], Tmin=(928.648,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-242.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'HCOH(T)(285)',
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
    label = '[CH2]C=[C]O(6303)',
    structure = SMILES('[CH2]C=[C]O'),
    E0 = (224.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.206965,'amu*angstrom^2'), symmetry=1, barrier=(4.75853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.69007,'amu*angstrom^2'), symmetry=1, barrier=(15.8661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.71481,0.0306825,-3.73315e-05,3.08792e-08,-1.09328e-11,27020.3,15.471], Tmin=(100,'K'), Tmax=(771.374,'K')), NASAPolynomial(coeffs=[4.68419,0.0180544,-8.07735e-06,1.53601e-09,-1.0692e-13,26788.3,6.94699], Tmin=(771.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=CJO) + radical(Allyl_P)"""),
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
    E0 = (34.9002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (185.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (161.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (135.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (255.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (241.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (246.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (164.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (197.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (151.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (267.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (221.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (185.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (67.9404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (61.7824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (491.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (98.3003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (43.2682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (43.2682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (211.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (42.8079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (430.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O[C]=CC[CH]O(1255)'],
    products = ['HCCOH(50)', 'CH2CHOH(42)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'O=CCC=[C]O(8483)'],
    products = ['O[C]=CC[CH]O(1255)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2818 used for Od_CO-CsH;HJ
Exact match found for rate rule [Od_CO-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'O=C=CC[CH]O(5758)'],
    products = ['O[C]=CC[CH]O(1255)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.185e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""Estimated using template [Od_R;HJ] for rate rule [Od_Cdd;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'O[C]=CC=CO(8484)'],
    products = ['O[C]=CC[CH]O(1255)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(105.863,'m^3/(mol*s)'), n=1.66278, Ea=(8.08939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'OC#CC[CH]O(8485)'],
    products = ['O[C]=CC[CH]O(1255)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]O(172)', 'CH2CHOH(42)'],
    products = ['O[C]=CC[CH]O(1255)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00668046,'m^3/(mol*s)'), n=2.5095, Ea=(31.5264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-OsH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['HCCOH(50)', '[CH2][CH]O(284)'],
    products = ['O[C]=CC[CH]O(1255)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O[C]=CC[CH]O(1255)'],
    products = ['[O]CCC=[C]O(8486)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4500,'s^-1'), n=2.62, Ea=(129.286,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 322 used for R2H_S;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O[C]=CC[CH]O(1255)'],
    products = ['O=C[CH]C[CH]O(5761)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O[C]=CC[CH]O(1255)'],
    products = ['O[C]=C[CH]CO(8487)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.58526e-05,'s^-1'), n=4.94167, Ea=(116.966,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_1H;Cs_H_out_H/OneDe] + [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_1H] for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O[CH]C[C]=CO(8488)'],
    products = ['O[C]=CC[CH]O(1255)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O[C]=[C]CCO(8489)'],
    products = ['O[C]=CC[CH]O(1255)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.00524938,'s^-1'), n=4.55, Ea=(129.286,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/NonDeO] + [R3H_SS_Cs;Cd_rad_out;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O[C]=CC[CH]O(1255)'],
    products = ['O[CH][CH]C=CO(8490)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O[C]=CC[CH]O(1255)'],
    products = ['O=CC[CH][CH]O(5763)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O[C]=CC[CH]O(1255)'],
    products = ['O=[C][CH]CCO(5762)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_1H;XH_out] for rate rule [R5HJ_3;C_rad_out_H/NonDeO;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]O(284)', '[CH]=[C]O(172)'],
    products = ['O[C]=CC[CH]O(1255)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['O[C]=CC[CH]O(1255)'],
    products = ['OC=CC=CO(8491)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O[C]=CC[CH]O(1255)'],
    products = ['O=CCC=CO(5765)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O[C]=CC[CH]O(1255)'],
    products = ['O=C=CCCO(5766)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)C=[C]O(1257)'],
    products = ['O[C]=CC[CH]O(1255)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O[C]=CC[CH]O(1255)'],
    products = ['OC1=CCC1O(8480)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeO;CdsinglepriND_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HCOH(T)(285)', '[CH2]C=[C]O(6303)'],
    products = ['O[C]=CC[CH]O(1255)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '840',
    isomers = [
        'O[C]=CC[CH]O(1255)',
    ],
    reactants = [
        ('HCCOH(50)', 'CH2CHOH(42)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '840',
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

