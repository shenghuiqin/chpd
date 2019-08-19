species(
    label = 'S(557)(556)',
    structure = SMILES('[O]C([CH]C=C=O)C=CC=O'),
    E0 = (-17.2521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,458.251,458.251,458.251,458.251],'cm^-1')),
        HinderedRotor(inertia=(0.158464,'amu*angstrom^2'), symmetry=1, barrier=(23.6137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158464,'amu*angstrom^2'), symmetry=1, barrier=(23.6137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158464,'amu*angstrom^2'), symmetry=1, barrier=(23.6137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158464,'amu*angstrom^2'), symmetry=1, barrier=(23.6137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4320.36,'J/mol'), sigma=(7.26806,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=674.83 K, Pc=25.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.350069,0.0821835,-6.70099e-05,2.57063e-08,-3.85463e-12,-1907.11,38.3751], Tmin=(100,'K'), Tmax=(1598.08,'K')), NASAPolynomial(coeffs=[23.9129,0.0214514,-1.00036e-05,1.9245e-09,-1.34168e-13,-9661.75,-90.0539], Tmin=(1598.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.2521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ)"""),
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
    label = 'O=C=C[CH]C(=O)C=CC=O(7144)',
    structure = SMILES('[O]C(C=CC=O)=CC=C=O'),
    E0 = (-119.39,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2120,512.5,787.5,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.710751,'amu*angstrom^2'), symmetry=1, barrier=(16.3416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.71268,'amu*angstrom^2'), symmetry=1, barrier=(16.3859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.713579,'amu*angstrom^2'), symmetry=1, barrier=(16.4066,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0940242,0.0969957,-0.000130077,9.44439e-08,-2.78069e-11,-14218.1,31.4797], Tmin=(100,'K'), Tmax=(825.529,'K')), NASAPolynomial(coeffs=[12.7085,0.0349636,-1.73658e-05,3.42369e-09,-2.43197e-13,-16331.9,-27.8307], Tmin=(825.529,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CCO)H) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C(C=C=C=O)C=CC=O(7145)',
    structure = SMILES('[O]C(C=C=C=O)C=CC=O'),
    E0 = (135.135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2120,512.5,787.5,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000585947,'amu*angstrom^2'), symmetry=1, barrier=(6.65287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288376,'amu*angstrom^2'), symmetry=1, barrier=(6.63032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285983,'amu*angstrom^2'), symmetry=1, barrier=(6.57531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36954,0.0612954,-4.96352e-05,1.95053e-08,-3.14206e-12,16344.1,17.5524], Tmin=(100,'K'), Tmax=(1414.25,'K')), NASAPolynomial(coeffs=[13.0996,0.0281188,-1.44471e-05,2.91797e-09,-2.09896e-13,13026.3,-43.104], Tmin=(1414.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(CC(C)OJ)"""),
)

species(
    label = 'S(525)(524)',
    structure = SMILES('O=CC=CC=O'),
    E0 = (-204.827,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.770505,'amu*angstrom^2'), symmetry=1, barrier=(17.7154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.772944,'amu*angstrom^2'), symmetry=1, barrier=(17.7715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3741.66,'J/mol'), sigma=(5.7541,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=584.44 K, Pc=44.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03135,0.0345588,-1.9969e-05,4.287e-09,-3.16505e-13,-24613.1,14.2636], Tmin=(100,'K'), Tmax=(2498.65,'K')), NASAPolynomial(coeffs=[27.6611,0.000737979,-3.0321e-06,6.66313e-10,-4.41138e-14,-38671.9,-130.618], Tmin=(2498.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.827,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = '[CH]C=C=O(7112)',
    structure = SMILES('[CH]C=C=O'),
    E0 = (292.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3010,987.5,1337.5,450,1655,245.068,247.168,247.199],'cm^-1')),
        HinderedRotor(inertia=(1.15951,'amu*angstrom^2'), symmetry=1, barrier=(50.7417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00496,0.0183622,-5.86278e-07,-8.1995e-09,3.3864e-12,35273.9,13.1075], Tmin=(100,'K'), Tmax=(1106.37,'K')), NASAPolynomial(coeffs=[5.96289,0.0150316,-6.0541e-06,1.11098e-09,-7.67696e-14,34168.7,-3.49842], Tmin=(1106.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
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
    label = 'O=C=CC=CC=CC=O(7146)',
    structure = SMILES('O=C=CC=CC=CC=O'),
    E0 = (-43.3021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.819906,'amu*angstrom^2'), symmetry=1, barrier=(18.8512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.819661,'amu*angstrom^2'), symmetry=1, barrier=(18.8456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.819637,'amu*angstrom^2'), symmetry=1, barrier=(18.8451,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.844097,0.0751875,-7.38281e-05,3.94745e-08,-8.92999e-12,-5099.18,27.366], Tmin=(100,'K'), Tmax=(1033.29,'K')), NASAPolynomial(coeffs=[10.7249,0.0369379,-1.83024e-05,3.65007e-09,-2.62457e-13,-7141.13,-20.6266], Tmin=(1033.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.3021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CCO)H) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'S(780)(779)',
    structure = SMILES('O=[C]C=CC=O'),
    E0 = (-44.2053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.882295,'amu*angstrom^2'), symmetry=1, barrier=(20.2857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.879664,'amu*angstrom^2'), symmetry=1, barrier=(20.2252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3741.66,'J/mol'), sigma=(5.7541,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=584.44 K, Pc=44.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70479,0.0574491,-9.74937e-05,9.58172e-08,-3.8062e-11,-5240.81,19.2682], Tmin=(100,'K'), Tmax=(736.114,'K')), NASAPolynomial(coeffs=[5.25854,0.0278824,-1.6346e-05,3.39827e-09,-2.46552e-13,-5486.14,5.0995], Tmin=(736.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.2053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
)

species(
    label = 'CHCHCHO(1679)',
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
    label = '[O][CH]C=CC=O(5602)',
    structure = SMILES('[O]C=CC=C[O]'),
    E0 = (-40.7388,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.42008,'amu*angstrom^2'), symmetry=1, barrier=(32.6505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42346,0.0345738,3.51007e-05,-9.27285e-08,4.46226e-11,-4786.29,19.0904], Tmin=(100,'K'), Tmax=(909.432,'K')), NASAPolynomial(coeffs=[24.1228,-0.00681144,6.94748e-06,-1.41392e-09,9.1781e-14,-11332.3,-101.556], Tmin=(909.432,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.7388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'O=C=CC=C(O)C=CC=O(7162)',
    structure = SMILES('O=C=CC=C(O)C=CC=O'),
    E0 = (-257.195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2782.5,750,1395,475,1775,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.843107,'amu*angstrom^2'), symmetry=1, barrier=(19.3847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.841343,'amu*angstrom^2'), symmetry=1, barrier=(19.3441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.841714,'amu*angstrom^2'), symmetry=1, barrier=(19.3527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.845833,'amu*angstrom^2'), symmetry=1, barrier=(19.4474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.369783,0.0999539,-0.000119194,7.30845e-08,-1.79172e-11,-30779.2,30.9648], Tmin=(100,'K'), Tmax=(990.016,'K')), NASAPolynomial(coeffs=[16.6382,0.0312358,-1.50769e-05,2.97337e-09,-2.12674e-13,-34146.8,-50.9183], Tmin=(990.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-257.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CCO)H) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'S(558)(557)',
    structure = SMILES('O=C=CCC(=O)C=CC=O'),
    E0 = (-254.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0341412,'amu*angstrom^2'), symmetry=1, barrier=(13.3068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202629,'amu*angstrom^2'), symmetry=1, barrier=(4.65883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.578584,'amu*angstrom^2'), symmetry=1, barrier=(13.3028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205588,'amu*angstrom^2'), symmetry=1, barrier=(4.72687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4173.58,'J/mol'), sigma=(6.77968,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=651.90 K, Pc=30.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.33415,0.0896355,-0.000125339,1.14301e-07,-4.47385e-11,-30518.1,33.3045], Tmin=(100,'K'), Tmax=(703.695,'K')), NASAPolynomial(coeffs=[5.60267,0.0517251,-2.75555e-05,5.58284e-09,-4.01819e-13,-31062.4,11.1392], Tmin=(703.695,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-254.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'O=C=C=CC(O)C=CC=O(7163)',
    structure = SMILES('O=C=C=CC(O)C=CC=O'),
    E0 = (-95.2263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24461,0.0648131,-5.38133e-05,2.23181e-08,-3.83931e-12,-11357.8,18.1957], Tmin=(100,'K'), Tmax=(1326.01,'K')), NASAPolynomial(coeffs=[12.4103,0.0311308,-1.57109e-05,3.16139e-09,-2.2756e-13,-14319,-38.8227], Tmin=(1326.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.2263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
)

species(
    label = 'O=CC=CC(=O)C=CC=O(7177)',
    structure = SMILES('O=CC=CC(=O)C=CC=O'),
    E0 = (-270.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.240641,0.0944848,-0.000147571,1.48204e-07,-6.12111e-11,-32365.5,32.2285], Tmin=(100,'K'), Tmax=(736.087,'K')), NASAPolynomial(coeffs=[3.21232,0.0581519,-3.24001e-05,6.64198e-09,-4.79469e-13,-32256.1,22.5168], Tmin=(736.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-270.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'O=C=CC1OC1C=CC=O(7166)',
    structure = SMILES('O=C=CC1OC1C=CC=O'),
    E0 = (-149.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.297315,0.0928019,-0.000104051,6.25571e-08,-1.49609e-11,-17789.9,31.5139], Tmin=(100,'K'), Tmax=(1022.02,'K')), NASAPolynomial(coeffs=[15.8872,0.0294582,-1.10824e-05,1.9129e-09,-1.26473e-13,-21098.1,-46.9198], Tmin=(1022.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(Ethylene_oxide)"""),
)

species(
    label = 'O=CC=CC1C=CC(=O)O1(7178)',
    structure = SMILES('O=CC=CC1C=CC(=O)O1'),
    E0 = (-244.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15678,0.0514005,-6.51669e-06,-1.95016e-08,8.08176e-12,-29314.6,32.8569], Tmin=(100,'K'), Tmax=(1201.69,'K')), NASAPolynomial(coeffs=[13.63,0.0346955,-1.66389e-05,3.29753e-09,-2.36272e-13,-34104,-37.0657], Tmin=(1201.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclopentane)"""),
)

species(
    label = 'O=C=C[CH]C1OC1[CH]C=O(7139)',
    structure = SMILES('[O]C=CC1OC1[CH]C=C=O'),
    E0 = (8.45882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.20569,0.100597,-9.93819e-05,4.73169e-08,-8.28989e-12,1271.81,37.3713], Tmin=(100,'K'), Tmax=(1676.58,'K')), NASAPolynomial(coeffs=[25.2666,0.0122104,1.33681e-07,-3.80714e-10,3.42399e-14,-4729.51,-99.7881], Tmin=(1676.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.45882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + ring(Ethylene_oxide) + radical(C=CCJCO) + radical(C=COJ)"""),
)

species(
    label = '[O]C1C([CH]C=O)C1C=C=O(7140)',
    structure = SMILES('[O]C=CC1C([O])C1C=C=O'),
    E0 = (90.5293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.157094,0.0709811,-1.41219e-05,-4.80316e-08,2.76175e-11,11056.5,33.8909], Tmin=(100,'K'), Tmax=(947.759,'K')), NASAPolynomial(coeffs=[25.43,0.0135403,-3.21463e-06,5.71366e-10,-4.72887e-14,3936.14,-100.156], Tmin=(947.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.5293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=[C]C1[CH]C(C=CC=O)O1(7141)',
    structure = SMILES('O=[C]C1[CH]C(C=CC=O)O1'),
    E0 = (90.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.584964,0.0573071,6.31517e-06,-5.77653e-08,2.88558e-11,11043.6,37.3423], Tmin=(100,'K'), Tmax=(941.689,'K')), NASAPolynomial(coeffs=[19.7846,0.020297,-5.68535e-06,9.61536e-10,-7.03943e-14,5452.56,-64.6178], Tmin=(941.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Oxetane) + radical(CCCJ=O) + radical(CCJCO)"""),
)

species(
    label = '[O]C1C=CC([CH]C=C=O)O1(7142)',
    structure = SMILES('[O]C1C=CC([CH]C=C=O)O1'),
    E0 = (-37.9156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.578496,0.0566002,8.17932e-06,-5.16218e-08,2.30159e-11,-4420.72,32.9799], Tmin=(100,'K'), Tmax=(1023.37,'K')), NASAPolynomial(coeffs=[19.0073,0.0274351,-1.19039e-05,2.39294e-09,-1.78727e-13,-10437.3,-67.321], Tmin=(1023.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-37.9156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + ring(25dihydrofuran) + radical(CCOJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[O]C1C=CC([O])C1C=C=O(7143)',
    structure = SMILES('[O]C1C=CC([O])C1C=C=O'),
    E0 = (122.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.310038,0.0615318,1.65495e-07,-5.15213e-08,2.54846e-11,14914,33.2709], Tmin=(100,'K'), Tmax=(984.229,'K')), NASAPolynomial(coeffs=[21.8151,0.0203912,-7.63368e-06,1.51386e-09,-1.16065e-13,8440.3,-81.5189], Tmin=(984.229,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.752,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + ring(Cyclopentene) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=C=C[CH][C](O)C=CC=O(7147)',
    structure = SMILES('[O]C=CC=C(O)[CH]C=C=O'),
    E0 = (-110.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.34017,0.115261,-0.000130136,6.79807e-08,-1.31793e-11,-13093.3,36.0085], Tmin=(100,'K'), Tmax=(1432.36,'K')), NASAPolynomial(coeffs=[31.7909,0.00495184,1.08477e-06,-4.0261e-10,3.17294e-14,-21332.7,-135.549], Tmin=(1432.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(C=CCJCO) + radical(C=COJ)"""),
)

species(
    label = '[O][C](C=CC=O)CC=C=O(7148)',
    structure = SMILES('[O]C=CC=C([O])CC=C=O'),
    E0 = (-90.0459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2120,512.5,787.5,350,440,435,1725,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.830871,'amu*angstrom^2'), symmetry=1, barrier=(19.1034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.832428,'amu*angstrom^2'), symmetry=1, barrier=(19.1392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.832653,'amu*angstrom^2'), symmetry=1, barrier=(19.1443,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35738,0.102716,-0.000114105,6.11803e-08,-1.24414e-11,-10624.2,36.778], Tmin=(100,'K'), Tmax=(1300.84,'K')), NASAPolynomial(coeffs=[25.8767,0.0125592,-2.74939e-06,3.21665e-10,-1.69417e-14,-17167,-99.6875], Tmin=(1300.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.0459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C(C=CC=O)C[C]=C=O(7149)',
    structure = SMILES('[O]C(C=CC=O)C[C]=C=O'),
    E0 = (114.885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,1685,370,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,243.922,243.925,243.929,243.931],'cm^-1')),
        HinderedRotor(inertia=(0.450028,'amu*angstrom^2'), symmetry=1, barrier=(19.0031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.45007,'amu*angstrom^2'), symmetry=1, barrier=(19.0031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00283318,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.450045,'amu*angstrom^2'), symmetry=1, barrier=(19.0031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.380899,0.0859057,-8.55372e-05,4.49265e-08,-9.85303e-12,13942.6,34.6444], Tmin=(100,'K'), Tmax=(1069.68,'K')), NASAPolynomial(coeffs=[12.9029,0.0390818,-1.98783e-05,4.00636e-09,-2.89653e-13,11263.6,-26.6108], Tmin=(1069.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + radical(CCCJ=C=O) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=C=C[CH]C(O)[C]=CC=O(7150)',
    structure = SMILES('O=C=C[CH]C(O)[C]=CC=O'),
    E0 = (-9.77113,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2120,512.5,787.5,3025,407.5,1350,352.5,3615,1277.5,1000,1685,370,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.110025,0.0835179,-7.43816e-05,3.24161e-08,-5.67408e-12,-1033.79,37.3566], Tmin=(100,'K'), Tmax=(1349.56,'K')), NASAPolynomial(coeffs=[18.2432,0.0297734,-1.46471e-05,2.90843e-09,-2.08016e-13,-5928.24,-55.5617], Tmin=(1349.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.77113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[O]C([C]=CC=O)CC=C=O(7151)',
    structure = SMILES('[O]C([C]=CC=O)CC=C=O'),
    E0 = (150.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,1685,370,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,274.42,274.42,274.421,274.424],'cm^-1')),
        HinderedRotor(inertia=(0.00223858,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247846,'amu*angstrom^2'), symmetry=1, barrier=(13.2445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247831,'amu*angstrom^2'), symmetry=1, barrier=(13.2445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247837,'amu*angstrom^2'), symmetry=1, barrier=(13.2445,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.291046,0.0900871,-0.000104314,6.78572e-08,-1.86525e-11,18220.8,35.7465], Tmin=(100,'K'), Tmax=(864.511,'K')), NASAPolynomial(coeffs=[10.5291,0.0427167,-2.21225e-05,4.4755e-09,-3.23799e-13,16450.6,-12.1558], Tmin=(864.511,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=C=C[CH]C(O)C=[C]C=O(7152)',
    structure = SMILES('[O]C=C=CC(O)[CH]C=C=O'),
    E0 = (-50.8762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.51515,0.0809203,-4.03821e-05,-2.05757e-08,1.73824e-11,-5939.84,38.1772], Tmin=(100,'K'), Tmax=(967.338,'K')), NASAPolynomial(coeffs=[25.8869,0.0148381,-4.73244e-06,9.0661e-10,-7.10356e-14,-13063.9,-98.742], Tmin=(967.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.8762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[O]C(C=[C]C=O)CC=C=O(7153)',
    structure = SMILES('[O]C=C=CC([O])CC=C=O'),
    E0 = (109.341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,184.665,186.511,186.788,186.985,189.071],'cm^-1')),
        HinderedRotor(inertia=(0.861515,'amu*angstrom^2'), symmetry=1, barrier=(21.3432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.865256,'amu*angstrom^2'), symmetry=1, barrier=(21.3378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.864523,'amu*angstrom^2'), symmetry=1, barrier=(21.3398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.867969,0.0940248,-9.39136e-05,4.57734e-08,-8.63516e-12,13337,38.4618], Tmin=(100,'K'), Tmax=(1302.41,'K')), NASAPolynomial(coeffs=[23.7366,0.0184574,-6.88068e-06,1.22312e-09,-8.35458e-14,6928,-86.7418], Tmin=(1302.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=C=[C][CH]C(O)C=CC=O(7154)',
    structure = SMILES('O=C=[C][CH]C(O)C=CC=O'),
    E0 = (-45.3321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.317521,0.0849979,-7.34959e-05,3.03453e-08,-4.92754e-12,-5288.66,38.1414], Tmin=(100,'K'), Tmax=(1475.49,'K')), NASAPolynomial(coeffs=[22.3126,0.0236484,-1.11271e-05,2.16535e-09,-1.5285e-13,-11966.7,-79.8386], Tmin=(1475.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.3321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJC(O)C=C) + radical(CCCJ=C=O)"""),
)

species(
    label = 'O=[C]C=CC(O)[CH]C=C=O(7155)',
    structure = SMILES('O=[C]C=CC(O)[CH]C=C=O'),
    E0 = (-86.9915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.258553,0.0910723,-8.9233e-05,4.28802e-08,-8.19159e-12,-10307.3,37.0393], Tmin=(100,'K'), Tmax=(1256.33,'K')), NASAPolynomial(coeffs=[19.6674,0.0276306,-1.34865e-05,2.68558e-09,-1.93164e-13,-15314,-63.6387], Tmin=(1256.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.9915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJC(O)C=C) + radical(C=CCJ=O)"""),
)

species(
    label = '[O]C(C=C[C]=O)CC=C=O(7156)',
    structure = SMILES('[O]C([CH]C=C=O)CC=C=O'),
    E0 = (44.4036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2110,2130,495,530,650,925,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.136877,0.0954806,-0.000107904,6.50904e-08,-1.59992e-11,5485.78,32.9878], Tmin=(100,'K'), Tmax=(979.084,'K')), NASAPolynomial(coeffs=[14.3448,0.0363162,-1.72614e-05,3.37045e-09,-2.39523e-13,2650.04,-36.5717], Tmin=(979.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.4036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-(Cdd-O2d)CsH) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O]C(C=CC=O)C1[CH]C1=O(7157)',
    structure = SMILES('[O]C(C=CC=O)C1[CH]C1=O'),
    E0 = (148.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.297775,0.0813594,-6.11345e-05,1.69115e-08,-1.43818e-13,18047.1,35.5809], Tmin=(100,'K'), Tmax=(1131.58,'K')), NASAPolynomial(coeffs=[21.1772,0.0244388,-1.08556e-05,2.12095e-09,-1.52776e-13,11971.1,-76.0505], Tmin=(1131.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(cyclopropanone) + radical(CC(C)OJ) + radical(CCJC=O)"""),
)

species(
    label = 'O=[C]C=CC1[CH]C(C=O)O1(7130)',
    structure = SMILES('O=C=C[CH]C1[CH]C(C=O)O1'),
    E0 = (62.5025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.3736,0.062027,2.06252e-06,-5.86167e-08,3.08716e-11,7664.33,33.4769], Tmin=(100,'K'), Tmax=(919.936,'K')), NASAPolynomial(coeffs=[20.6676,0.019915,-4.48732e-06,6.37691e-10,-4.42298e-14,1978.58,-73.3455], Tmin=(919.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.5025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + ring(Oxetane) + radical(C=CCJCO) + radical(CCJCO)"""),
)

species(
    label = '[O]C1[CH]C(C=O)C1C=C=O(7158)',
    structure = SMILES('[O]C1[CH]C(C=O)C1C=C=O'),
    E0 = (144.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.284188,0.06916,-3.65164e-05,-3.79738e-09,6.51926e-12,17580.7,35.2168], Tmin=(100,'K'), Tmax=(1048.93,'K')), NASAPolynomial(coeffs=[17.71,0.0270171,-1.10131e-05,2.08751e-09,-1.4921e-13,12587.7,-56.0594], Tmin=(1048.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C(C=CC=O)C1C=[C]O1(7159)',
    structure = SMILES('[O]C(C=CC=O)C1C=[C]O1'),
    E0 = (216.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.18565,0.0722885,-2.28394e-05,-3.21465e-08,1.9185e-11,26211.3,37.0466], Tmin=(100,'K'), Tmax=(1005.24,'K')), NASAPolynomial(coeffs=[24.7203,0.0175727,-7.42903e-06,1.5601e-09,-1.22092e-13,18961.3,-94.3959], Tmin=(1005.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(C=CJO) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=CC=CC1[CH][CH]C(=O)O1(7160)',
    structure = SMILES('[O]C1=C[CH]C(C=CC=O)O1'),
    E0 = (-130.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.477796,0.0525942,3.53812e-05,-9.40062e-08,4.21227e-11,-15522.1,28.8874], Tmin=(100,'K'), Tmax=(958.961,'K')), NASAPolynomial(coeffs=[23.7663,0.0169582,-5.08221e-06,1.00518e-09,-8.25555e-14,-22816.6,-97.2353], Tmin=(958.961,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(2,3-Dihydrofuran) + radical(C=COJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'S(556)(555)',
    structure = SMILES('O=C=C[CH]C1C=C[CH]OO1'),
    E0 = (108.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4040.83,'J/mol'), sigma=(7.10062,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=631.17 K, Pc=25.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.761823,0.0362626,9.75141e-05,-1.65729e-07,6.87308e-11,13214.9,30.2471], Tmin=(100,'K'), Tmax=(955.002,'K')), NASAPolynomial(coeffs=[27.5213,0.0121235,-2.70014e-06,6.53558e-10,-6.67665e-14,4093.54,-118.616], Tmin=(955.002,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + ring(36dihydro12dioxin) + radical(C=CCJC(O)C=C) + radical(C=CCJO)"""),
)

species(
    label = '[O]C1C=C[CH]OC1C=C=O(7161)',
    structure = SMILES('[O]C1[CH]C=COC1C=C=O'),
    E0 = (-5.53929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.351042,0.0781959,-3.09208e-05,-3.25926e-08,2.32559e-11,-493.367,25.1076], Tmin=(100,'K'), Tmax=(924.216,'K')), NASAPolynomial(coeffs=[24.5613,0.0157135,-3.09529e-06,4.13883e-10,-2.99368e-14,-7034.57,-103.592], Tmin=(924.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.53929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = 'O=CC=CC1[CH]C=[C]OO1(7123)',
    structure = SMILES('O=CC=CC1[CH]C=[C]OO1'),
    E0 = (215.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.651598,0.0541178,1.6894e-05,-6.30542e-08,2.80225e-11,26100,32.8985], Tmin=(100,'K'), Tmax=(993.33,'K')), NASAPolynomial(coeffs=[19.1774,0.0257792,-1.0172e-05,1.99654e-09,-1.4949e-13,20137.2,-67.8423], Tmin=(993.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + ring(34dihydro12dioxin) + radical(C=CCJC(O)C=C) + radical(C=CJO)"""),
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
    label = 'C=CC([O])[CH]C=C=O(7164)',
    structure = SMILES('C=CC([O])[CH]C=C=O'),
    E0 = (106.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,520.166,520.383,520.527,520.549],'cm^-1')),
        HinderedRotor(inertia=(0.142769,'amu*angstrom^2'), symmetry=1, barrier=(27.4303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142724,'amu*angstrom^2'), symmetry=1, barrier=(27.4302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142777,'amu*angstrom^2'), symmetry=1, barrier=(27.4367,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.887913,0.0528099,-3.28189e-06,-3.83465e-08,1.92745e-11,12903.2,30.8213], Tmin=(100,'K'), Tmax=(990.799,'K')), NASAPolynomial(coeffs=[18.2209,0.0191812,-7.39757e-06,1.44824e-09,-1.0885e-13,7684.41,-61.6434], Tmin=(990.799,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[O][CH]C(C=C=O)C=CC=O(7165)',
    structure = SMILES('[O][CH]C(C=C=O)C=CC=O'),
    E0 = (90.8767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,245.738,245.866,245.927,246.118],'cm^-1')),
        HinderedRotor(inertia=(0.19732,'amu*angstrom^2'), symmetry=1, barrier=(8.44505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343764,'amu*angstrom^2'), symmetry=1, barrier=(14.6599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19741,'amu*angstrom^2'), symmetry=1, barrier=(8.44519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199907,'amu*angstrom^2'), symmetry=1, barrier=(8.44802,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.61936,0.115257,-0.000195173,1.81492e-07,-6.58431e-11,11083.1,36.2849], Tmin=(100,'K'), Tmax=(819.633,'K')), NASAPolynomial(coeffs=[7.93291,0.0478876,-2.49707e-05,4.89971e-09,-3.4187e-13,10542.1,1.97854], Tmin=(819.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.8767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCd(CCO)H) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O]C([CH]C=C=O)C=C=CO(7167)',
    structure = SMILES('[O]C([CH]C=C=O)C=C=CO'),
    E0 = (38.022,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,540,610,2055,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1380,1390,370,380,2900,435,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.27514,'amu*angstrom^2'), symmetry=1, barrier=(29.3179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27337,'amu*angstrom^2'), symmetry=1, barrier=(29.2773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27251,'amu*angstrom^2'), symmetry=1, barrier=(29.2575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27288,'amu*angstrom^2'), symmetry=1, barrier=(29.266,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.76494,0.0839171,-3.79901e-05,-3.16202e-08,2.36128e-11,4763.64,37.4714], Tmin=(100,'K'), Tmax=(946.07,'K')), NASAPolynomial(coeffs=[29.0447,0.00940811,-1.55126e-06,2.70882e-10,-2.69863e-14,-3182.68,-116.877], Tmin=(946.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=C=C[CH][CH]C=CC=O(7168)',
    structure = SMILES('[O]C=CC=C[CH]C=C=O'),
    E0 = (87.8386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3025,407.5,1350,352.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.30356,'amu*angstrom^2'), symmetry=1, barrier=(29.9714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30244,'amu*angstrom^2'), symmetry=1, barrier=(29.9456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3034,'amu*angstrom^2'), symmetry=1, barrier=(29.9677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.111807,0.0743674,-3.73509e-05,-1.98856e-08,1.71558e-11,10727.3,29.2678], Tmin=(100,'K'), Tmax=(943.897,'K')), NASAPolynomial(coeffs=[23.6618,0.0132699,-3.26615e-06,5.43064e-10,-4.17879e-14,4473.07,-93.4096], Tmin=(943.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.8386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CCJC=C=O)"""),
)

species(
    label = 'O=[C][CH]C1OC1C=CC=O(7169)',
    structure = SMILES('O=[C][CH]C1OC1C=CC=O'),
    E0 = (93.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.137603,0.0854106,-8.62634e-05,4.61795e-08,-9.82285e-12,11412.1,37.1111], Tmin=(100,'K'), Tmax=(1145.91,'K')), NASAPolynomial(coeffs=[16.6121,0.0269431,-9.72944e-06,1.65374e-09,-1.08818e-13,7573.34,-45.9777], Tmin=(1145.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(93.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Ethylene_oxide) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1C=CC(=O)C1[CH]C=O(7170)',
    structure = SMILES('[O]C=CC1C(=O)C=CC1[O]'),
    E0 = (-23.0596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.659188,0.052815,2.10865e-05,-7.00534e-08,3.13401e-11,-2634.67,32.5307], Tmin=(100,'K'), Tmax=(980.498,'K')), NASAPolynomial(coeffs=[20.3305,0.0222871,-8.27735e-06,1.63168e-09,-1.24713e-13,-8882.28,-74.1726], Tmin=(980.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.0596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1C=CC(=O)C([O])C=C1(7073)',
    structure = SMILES('[O]C1C=CC(=O)C([O])C=C1'),
    E0 = (112.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04146,0.0119057,0.000194301,-2.97783e-07,1.25625e-10,13719.3,34.4917], Tmin=(100,'K'), Tmax=(915.671,'K')), NASAPolynomial(coeffs=[39.3213,-0.0172093,1.57581e-05,-3.07737e-09,1.92224e-13,919.169,-178.429], Tmin=(915.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + ring(Cycloheptane) + radical(C=OCOJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C(C=[C]C=O)C=CC=O(7171)',
    structure = SMILES('[O]C=C=CC([O])C=CC=O'),
    E0 = (91.038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.973344,'amu*angstrom^2'), symmetry=1, barrier=(22.3791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.974427,'amu*angstrom^2'), symmetry=1, barrier=(22.404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975368,'amu*angstrom^2'), symmetry=1, barrier=(22.4256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.347211,0.089868,-8.49321e-05,3.89589e-08,-7.0589e-12,11110.5,36.5233], Tmin=(100,'K'), Tmax=(1327.03,'K')), NASAPolynomial(coeffs=[21.0101,0.0254916,-1.21643e-05,2.40207e-09,-1.71902e-13,5442.15,-72.5559], Tmin=(1327.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]C=[C]C(O)C=CC=O(7172)',
    structure = SMILES('O=[C]C=[C]C(O)C=CC=O'),
    E0 = (62.4037,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1855,455,950,3615,1277.5,1000,1685,370,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.62777,'amu*angstrom^2'), symmetry=1, barrier=(14.4337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627768,'amu*angstrom^2'), symmetry=1, barrier=(14.4336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627018,'amu*angstrom^2'), symmetry=1, barrier=(14.4164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627565,'amu*angstrom^2'), symmetry=1, barrier=(14.4289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.626525,'amu*angstrom^2'), symmetry=1, barrier=(14.405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.337184,0.0972547,-9.72112e-05,-1.8738e-08,7.67354e-11,7620.2,35.4234], Tmin=(100,'K'), Tmax=(479.059,'K')), NASAPolynomial(coeffs=[9.89252,0.0458938,-2.53903e-05,5.1648e-09,-3.70336e-13,6378.53,-7.04764], Tmin=(479.059,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.4037,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(C=CCJ=O)"""),
)

species(
    label = '[O]C([C]=CC=O)C=CC=O(7173)',
    structure = SMILES('[O]C([C]=CC=O)C=CC=O'),
    E0 = (132.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,1380,1390,370,380,2900,435,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,243.976,243.984,243.985],'cm^-1')),
        HinderedRotor(inertia=(0.30731,'amu*angstrom^2'), symmetry=1, barrier=(12.9826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307336,'amu*angstrom^2'), symmetry=1, barrier=(12.9826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.3074,'amu*angstrom^2'), symmetry=1, barrier=(12.9826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307304,'amu*angstrom^2'), symmetry=1, barrier=(12.9826,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.273099,0.0922641,-0.000117681,9.10692e-08,-3.07514e-11,16017.9,35.744], Tmin=(100,'K'), Tmax=(703.469,'K')), NASAPolynomial(coeffs=[7.70275,0.0500198,-2.76062e-05,5.70999e-09,-4.17301e-13,14972.6,2.51326], Tmin=(703.469,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[O][C](C=CC=O)C=CC=O(7174)',
    structure = SMILES('[O]C=CC=C([O])C=CC=O'),
    E0 = (-121.213,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,2782.5,750,1395,475,1775,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0386,'amu*angstrom^2'), symmetry=1, barrier=(23.8795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03479,'amu*angstrom^2'), symmetry=1, barrier=(23.7919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.034,'amu*angstrom^2'), symmetry=1, barrier=(23.7737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53508,0.103178,-0.000109572,5.50845e-08,-1.04951e-11,-14363.5,35.5983], Tmin=(100,'K'), Tmax=(1378.76,'K')), NASAPolynomial(coeffs=[28.212,0.0110201,-2.93863e-06,4.43196e-10,-2.87746e-14,-22009.6,-115.45], Tmin=(1378.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-121.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C1[CH]C(C=O)C(=O)C=C1(7175)',
    structure = SMILES('[O]C1[CH]C(C=O)C(=O)C=C1'),
    E0 = (8.61165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30897,0.0394118,4.71885e-05,-8.03332e-08,2.96922e-11,1149.52,30.353], Tmin=(100,'K'), Tmax=(1062.23,'K')), NASAPolynomial(coeffs=[15.1761,0.036152,-1.73454e-05,3.56048e-09,-2.64978e-13,-4558.59,-50.3866], Tmin=(1062.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.61165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + ring(Cyclohexane) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1C=C[CH]OC(=O)C=C1(7176)',
    structure = SMILES('[O]C1[CH]C=COC(=O)C=C1'),
    E0 = (-6.60796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25167,0.0307857,9.18187e-05,-1.45302e-07,5.79494e-11,-669.095,28.8312], Tmin=(100,'K'), Tmax=(971.043,'K')), NASAPolynomial(coeffs=[21.6818,0.0208454,-7.47147e-06,1.57514e-09,-1.28934e-13,-8135.84,-87.1489], Tmin=(971.043,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.60796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)O2s) + ring(Cyclooctane) + radical(CC(C)OJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[CH]=CC([O])C=CC=O(5640)',
    structure = SMILES('[CH]=CC([O])C=CC=O'),
    E0 = (264.885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,242.235,242.348],'cm^-1')),
        HinderedRotor(inertia=(0.358141,'amu*angstrom^2'), symmetry=1, barrier=(14.9268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357158,'amu*angstrom^2'), symmetry=1, barrier=(14.9276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357621,'amu*angstrom^2'), symmetry=1, barrier=(14.9288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26288,0.0641803,-5.31865e-05,2.20801e-08,-3.81796e-12,31953.2,29.0106], Tmin=(100,'K'), Tmax=(1316.78,'K')), NASAPolynomial(coeffs=[12.0587,0.0313856,-1.58285e-05,3.16621e-09,-2.27018e-13,29110,-26.0439], Tmin=(1316.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[C]=O(1191)',
    structure = SMILES('[C]=O'),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08916,0.00200416,-1.61661e-05,2.55058e-08,-1.16424e-11,52802.7,4.52505], Tmin=(100,'K'), Tmax=(856.11,'K')), NASAPolynomial(coeffs=[0.961625,0.00569045,-3.48044e-06,7.19202e-10,-5.08041e-14,53738.7,21.4663], Tmin=(856.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CdCdJ2_triplet)"""),
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
    E0 = (113.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (362.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (107.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (207.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (140.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (252.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (10.6008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (5.60925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (7.72112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (7.72112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-9.39179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-10.1393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (29.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (90.9039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (104.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (178.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (178.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (96.7199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (178.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (295.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (132.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (291.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (402.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (153.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (128.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (77.3605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (132.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (178.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (110.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (144.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (216.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (59.4309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (108.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (58.4547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (215.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (292.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (233.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (181.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (330.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (94.1078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (62.0914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (112.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (253.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (204.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (282.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (121.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (67.4571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (33.6337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (165.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (703.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction179',
    reactants = ['H(3)(3)', 'O=C=C[CH]C(=O)C=CC=O(7144)'],
    products = ['S(557)(556)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction180',
    reactants = ['H(3)(3)', '[O]C(C=C=C=O)C=CC=O(7145)'],
    products = ['S(557)(556)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction181',
    reactants = ['S(525)(524)', '[CH]C=C=O(7112)'],
    products = ['S(557)(556)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.5e+07,'cm^3/(mol*s)'), n=2.16, Ea=(19.2464,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdH_O;YJ] for rate rule [CO-CdH_O;Y_1centerbirad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction182',
    reactants = ['O(4)(4)', 'O=C=CC=CC=CC=O(7146)'],
    products = ['S(557)(556)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.338459,'m^3/(mol*s)'), n=2.24862, Ea=(7.77083,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-CdH;YJ] for rate rule [Cds-CdH_Cds-CdH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction183',
    reactants = ['S(780)(779)', 'CHCHCHO(1679)'],
    products = ['S(557)(556)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0131003,'m^3/(mol*s)'), n=2.40999, Ea=(12.7705,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CdsJ-H] for rate rule [CO_O;CdsJ-H]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O][CH]C=CC=O(5602)', '[CH]C=C=O(7112)'],
    products = ['S(557)(556)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.83595e+06,'m^3/(mol*s)'), n=0.223047, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction201',
    reactants = ['S(557)(556)'],
    products = ['O=C=CC=C(O)C=CC=O(7162)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction202',
    reactants = ['S(557)(556)'],
    products = ['S(558)(557)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction203',
    reactants = ['S(557)(556)'],
    products = ['O=C=C=CC(O)C=CC=O(7163)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction219',
    reactants = ['S(557)(556)'],
    products = ['O=CC=CC(=O)C=CC=O(7177)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction206',
    reactants = ['S(557)(556)'],
    products = ['O=C=CC1OC1C=CC=O(7166)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction220',
    reactants = ['S(557)(556)'],
    products = ['O=CC=CC1C=CC(=O)O1(7178)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSDS;O_rad;Ypri_rad_out]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction174',
    reactants = ['S(557)(556)'],
    products = ['O=C=C[CH]C1OC1[CH]C=O(7139)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(46.3011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra_HDe_pri;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction175',
    reactants = ['S(557)(556)'],
    products = ['[O]C1C([CH]C=O)C1C=C=O(7140)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_csHCd] for rate rule [R4_S_D;doublebond_intra_HDe_pri;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction176',
    reactants = ['S(557)(556)'],
    products = ['O=[C]C1[CH]C(C=CC=O)O1(7141)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction177',
    reactants = ['S(557)(556)'],
    products = ['[O]C1C=CC([CH]C=C=O)O1(7142)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSM;multiplebond_intra;radadd_intra] for rate rule [R6_SSM_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction178',
    reactants = ['S(557)(556)'],
    products = ['[O]C1C=CC([O])C1C=C=O(7143)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSM;multiplebond_intra;radadd_intra_cs] for rate rule [R6_SSM_CO;carbonylbond_intra_H;radadd_intra_csHCd]
Euclidian distance = 3.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction184',
    reactants = ['S(557)(556)'],
    products = ['O=C=C[CH][C](O)C=CC=O(7147)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction185',
    reactants = ['S(557)(556)'],
    products = ['[O][C](C=CC=O)CC=C=O(7148)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.00613e+10,'s^-1'), n=0.845153, Ea=(195.403,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction186',
    reactants = ['[O]C(C=CC=O)C[C]=C=O(7149)'],
    products = ['S(557)(556)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.62853e+10,'s^-1'), n=0.944167, Ea=(180.133,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction187',
    reactants = ['O=C=C[CH]C(O)[C]=CC=O(7150)'],
    products = ['S(557)(556)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction188',
    reactants = ['[O]C([C]=CC=O)CC=C=O(7151)'],
    products = ['S(557)(556)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction189',
    reactants = ['S(557)(556)'],
    products = ['O=C=C[CH]C(O)C=[C]C=O(7152)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.71035e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_singleDe]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction190',
    reactants = ['[O]C(C=[C]C=O)CC=C=O(7153)'],
    products = ['S(557)(556)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_DSS;Cd_rad_out_singleDe;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction191',
    reactants = ['S(557)(556)'],
    products = ['O=C=[C][CH]C(O)C=CC=O(7154)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4Hall;O_rad_out;Cd_H_out_doubleC] for rate rule [R4HJ_2;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction192',
    reactants = ['S(557)(556)'],
    products = ['O=[C]C=CC(O)[CH]C=C=O(7155)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(722272,'s^-1'), n=1.6737, Ea=(94.6126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;O_rad_out;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction193',
    reactants = ['[O]C(C=C[C]=O)CC=C=O(7156)'],
    products = ['S(557)(556)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(252000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5H_SMSS;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction194',
    reactants = ['S(557)(556)'],
    products = ['[O]C(C=CC=O)C1[CH]C1=O(7157)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.61353e+11,'s^-1'), n=0.395207, Ea=(195.485,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs] for rate rule [R3_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction195',
    reactants = ['S(557)(556)'],
    products = ['O=[C]C=CC1[CH]C(C=O)O1(7130)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.65487e+07,'s^-1'), n=1.11281, Ea=(128.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_pri;radadd_intra] for rate rule [R4_S_D;doublebond_intra_pri_HCO;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction196',
    reactants = ['S(557)(556)'],
    products = ['[O]C1[CH]C(C=O)C1C=C=O(7158)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(8.81207e+07,'s^-1'), n=1.08774, Ea=(162.224,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R4_S_D;doublebond_intra_pri_HCO;radadd_intra_csHCd]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic
Ea raised from 160.6 to 162.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction197',
    reactants = ['S(557)(556)'],
    products = ['[O]C(C=CC=O)C1C=[C]O1(7159)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.80882e+13,'s^-1'), n=0.063734, Ea=(233.787,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4;multiplebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 232.9 to 233.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction198',
    reactants = ['S(557)(556)'],
    products = ['O=CC=CC1[CH][CH]C(=O)O1(7160)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(5.04724e+10,'s^-1'), n=0.246651, Ea=(76.6831,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction137',
    reactants = ['S(557)(556)'],
    products = ['S(556)(555)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.00901e+09,'s^-1'), n=0.463766, Ea=(125.893,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSM;multiplebond_intra;radadd_intra] for rate rule [R6_SSM_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 120.7 to 125.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction199',
    reactants = ['S(557)(556)'],
    products = ['[O]C1C=C[CH]OC1C=C=O(7161)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.00901e+09,'s^-1'), n=0.463766, Ea=(75.7068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSM;multiplebond_intra;radadd_intra_csHCd] for rate rule [R6_SSM_CO;carbonyl_intra_H;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction200',
    reactants = ['S(557)(556)'],
    products = ['O=CC=CC1[CH]C=[C]OO1(7123)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.41956e+10,'s^-1'), n=0.267163, Ea=(233.113,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;multiplebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 229.8 to 233.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction204',
    reactants = ['CO(10)(11)', 'C=CC([O])[CH]C=C=O(7164)'],
    products = ['S(557)(556)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction205',
    reactants = ['[O][CH]C(C=C=O)C=CC=O(7165)'],
    products = ['S(557)(556)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction207',
    reactants = ['[O]C([CH]C=C=O)C=C=CO(7167)'],
    products = ['S(557)(556)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction208',
    reactants = ['O(4)(4)', 'O=C=C[CH][CH]C=CC=O(7168)'],
    products = ['S(557)(556)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction210',
    reactants = ['S(557)(556)'],
    products = ['O=[C][CH]C1OC1C=CC=O(7169)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(111.36,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction211',
    reactants = ['S(557)(556)'],
    products = ['[O]C1C=CC(=O)C1[CH]C=O(7170)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(79.3435,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS_D;doublebond_intra;radadd_intra] for rate rule [R6_SMS_D;doublebond_intra_HDe_pri;radadd_intra_CO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction212',
    reactants = ['S(557)(556)'],
    products = ['[O]C1C=CC(=O)C([O])C=C1(7073)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.19e+11,'s^-1'), n=0.08, Ea=(130.023,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8;multiplebond_intra;radadd_intra] for rate rule [R8;carbonylbond_intra_H;radadd_intra_CO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 121.8 to 130.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O]C(C=[C]C=O)C=CC=O(7171)'],
    products = ['S(557)(556)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction214',
    reactants = ['O=[C]C=[C]C(O)C=CC=O(7172)'],
    products = ['S(557)(556)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[O]C([C]=CC=O)C=CC=O(7173)'],
    products = ['S(557)(556)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction216',
    reactants = ['S(557)(556)'],
    products = ['[O][C](C=CC=O)C=CC=O(7174)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(8.08094e+07,'s^-1'), n=1.1965, Ea=(138.911,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;XH_out] for rate rule [R4H_SDS;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction217',
    reactants = ['S(557)(556)'],
    products = ['[O]C1[CH]C(C=O)C(=O)C=C1(7175)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(5.041e+08,'s^-1'), n=0.7, Ea=(84.7093,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SDS_D;doublebond_intra_pri;radadd_intra] for rate rule [R6_SDS_D;doublebond_intra_pri_HCO;radadd_intra_CO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction218',
    reactants = ['S(557)(556)'],
    products = ['[O]C1C=C[CH]OC(=O)C=C1(7176)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(50.8858,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra_CO] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_CO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction221',
    reactants = ['CO(10)(11)', '[CH]=CC([O])C=CC=O(5640)'],
    products = ['S(557)(556)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.51e+11,'cm^3/(mol*s)','*|/',5), n=0, Ea=(20.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 7 used for COm;Cd_pri_rad
Exact match found for rate rule [COm;Cd_pri_rad]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction222',
    reactants = ['[C]=O(1191)', '[CH]=CC([O])C=CC=O(5640)'],
    products = ['S(557)(556)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

network(
    label = '424',
    isomers = [
        'S(557)(556)',
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
    label = '424',
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

