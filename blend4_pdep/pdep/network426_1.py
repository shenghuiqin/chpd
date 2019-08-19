species(
    label = 'S(560)(559)',
    structure = SMILES('C=CC(=O)C(C=O)C=C=O'),
    E0 = (-245.483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,375,552.5,462.5,1710,260.436,260.442,260.514],'cm^-1')),
        HinderedRotor(inertia=(0.00248555,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00248383,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140938,'amu*angstrom^2'), symmetry=1, barrier=(6.77909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.488228,'amu*angstrom^2'), symmetry=1, barrier=(23.5031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4117.78,'J/mol'), sigma=(6.74437,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=643.19 K, Pc=30.46 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.794687,0.0782275,-7.89878e-05,4.74386e-08,-1.25977e-11,-29415.9,33.3538], Tmin=(100,'K'), Tmax=(873.798,'K')), NASAPolynomial(coeffs=[7.94813,0.0454814,-2.27749e-05,4.55131e-09,-3.27424e-13,-30666.1,-0.192424], Tmin=(873.798,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-245.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
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
    label = 'CH2CHCO(2811)',
    structure = SMILES('C=C[C]=O'),
    E0 = (83.3963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.27992,'amu*angstrom^2'), symmetry=1, barrier=(29.4278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08334,0.0153204,6.54286e-06,-1.77538e-08,7.39385e-12,10067.5,11.895], Tmin=(100,'K'), Tmax=(1002.99,'K')), NASAPolynomial(coeffs=[6.97837,0.0111885,-4.32951e-06,8.06755e-10,-5.75269e-14,8712.68,-9.76679], Tmin=(1002.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.3963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'HCCO(47)(47)',
    structure = SMILES('C#C[O]'),
    E0 = (166.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,180],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1247.18,'J/mol'), sigma=(2.5,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87608,0.0221205,-3.58869e-05,3.05403e-08,-1.01281e-11,20163.4,13.6968], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91479,0.00371409,-1.30137e-06,2.06473e-10,-1.21477e-14,19359.6,-5.50567], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(166.705,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=CC([O])=CC=O(7240)',
    structure = SMILES('C=CC([O])=CC=O'),
    E0 = (-105.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.804458,'amu*angstrom^2'), symmetry=1, barrier=(18.4961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805941,'amu*angstrom^2'), symmetry=1, barrier=(18.5302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16867,0.0608176,-5.94523e-05,2.9485e-08,-5.85638e-12,-12604.2,21.9805], Tmin=(100,'K'), Tmax=(1208.64,'K')), NASAPolynomial(coeffs=[13.3959,0.0203512,-9.23038e-06,1.78305e-09,-1.26335e-13,-15559.8,-39.3255], Tmin=(1208.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=CC(=O)[CH]C=C=O(7278)',
    structure = SMILES('C=CC([O])=CC=C=O'),
    E0 = (4.09784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2120,512.5,787.5,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.751708,'amu*angstrom^2'), symmetry=1, barrier=(17.2833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.750015,'amu*angstrom^2'), symmetry=1, barrier=(17.2443,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.665792,0.072967,-8.35439e-05,5.02662e-08,-1.20063e-11,613.296,25.6599], Tmin=(100,'K'), Tmax=(1022.36,'K')), NASAPolynomial(coeffs=[13.6615,0.0221208,-8.94219e-06,1.61912e-09,-1.1046e-13,-2043.94,-37.3241], Tmin=(1022.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.09784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cd-Cd(CCO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = 'C=CC(=O)[C](C=O)C=C=O(7279)',
    structure = SMILES('C=CC([O])=C(C=O)C=C=O'),
    E0 = (-125.737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2120,512.5,787.5,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.719354,'amu*angstrom^2'), symmetry=1, barrier=(16.5394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.716991,'amu*angstrom^2'), symmetry=1, barrier=(16.485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.717262,'amu*angstrom^2'), symmetry=1, barrier=(16.4913,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.566655,0.107468,-0.000155118,1.16927e-07,-3.50765e-11,-14964.6,32.0384], Tmin=(100,'K'), Tmax=(816.573,'K')), NASAPolynomial(coeffs=[14.9308,0.0315509,-1.56588e-05,3.06586e-09,-2.16142e-13,-17495.5,-39.587], Tmin=(816.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-125.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-CdCs(CCO)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=[C]C(C=O)C=C=O(7280)',
    structure = SMILES('O=[C]C(C=O)C=C=O'),
    E0 = (-138.207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,1855,455,950,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.543857,'amu*angstrom^2'), symmetry=1, barrier=(12.5043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00578688,'amu*angstrom^2'), symmetry=1, barrier=(3.29876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.968631,'amu*angstrom^2'), symmetry=1, barrier=(22.2707,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68922,0.0520081,-4.80553e-05,2.32768e-08,-4.6526e-12,-16540.3,27.9244], Tmin=(100,'K'), Tmax=(1176.78,'K')), NASAPolynomial(coeffs=[10.33,0.0226368,-1.06165e-05,2.06692e-09,-1.46665e-13,-18573.9,-15.1691], Tmin=(1176.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C=CC(=O)C([C]=C=O)C=O(7281)',
    structure = SMILES('C=CC(=O)C([C]=C=O)C=O'),
    E0 = (-39.6909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2120,512.5,787.5,1685,370,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,188.167,188.168,2121.37],'cm^-1')),
        HinderedRotor(inertia=(0.00476146,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.43165,'amu*angstrom^2'), symmetry=1, barrier=(10.8448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24265,'amu*angstrom^2'), symmetry=1, barrier=(31.2224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24269,'amu*angstrom^2'), symmetry=1, barrier=(31.2224,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907858,0.0787578,-7.68016e-05,1.76939e-08,1.98832e-11,-4672.9,32.8568], Tmin=(100,'K'), Tmax=(524.221,'K')), NASAPolynomial(coeffs=[7.37951,0.0446332,-2.28127e-05,4.55185e-09,-3.25479e-13,-5561.05,3.81503], Tmin=(524.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.6909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)CJ=C=O)"""),
)

species(
    label = 'C=[C]C(=O)C(C=O)C=C=O(7282)',
    structure = SMILES('C=C=C([O])C(C=O)C=C=O'),
    E0 = (-49.1941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,540,610,2055,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,256.794,257.154,258.682],'cm^-1')),
        HinderedRotor(inertia=(0.282792,'amu*angstrom^2'), symmetry=1, barrier=(13.4725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283915,'amu*angstrom^2'), symmetry=1, barrier=(13.4696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284112,'amu*angstrom^2'), symmetry=1, barrier=(13.4564,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0426603,0.09047,-0.000108786,6.95235e-08,-1.7864e-11,-5777.06,34.2459], Tmin=(100,'K'), Tmax=(945.528,'K')), NASAPolynomial(coeffs=[14.1704,0.0307025,-1.39689e-05,2.66911e-09,-1.87218e-13,-8448.65,-33.1208], Tmin=(945.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.1941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cds-CdsCsOs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC(=O)C([C]=O)C=C=O(7283)',
    structure = SMILES('C=CC(=O)C([C]=O)C=C=O'),
    E0 = (-86.7941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,2120,512.5,787.5,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,211.925,268.436,1219.64],'cm^-1')),
        HinderedRotor(inertia=(0.17516,'amu*angstrom^2'), symmetry=1, barrier=(8.71212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383909,'amu*angstrom^2'), symmetry=1, barrier=(17.6547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217906,'amu*angstrom^2'), symmetry=1, barrier=(8.62198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016492,'amu*angstrom^2'), symmetry=1, barrier=(17.703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.809136,0.0774398,-8.29569e-05,5.13264e-08,-1.36926e-11,-10330.2,34.4561], Tmin=(100,'K'), Tmax=(880.974,'K')), NASAPolynomial(coeffs=[8.87775,0.0408043,-2.05781e-05,4.12126e-09,-2.9673e-13,-11751.8,-3.44763], Tmin=(880.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.7941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH]=CC(=O)C(C=O)C=C=O(7284)',
    structure = SMILES('[CH]=CC(=O)C(C=O)C=C=O'),
    E0 = (1.61343,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,1380,1390,370,380,2900,435,3120,650,792.5,1650,268.828,980.812],'cm^-1')),
        HinderedRotor(inertia=(0.145956,'amu*angstrom^2'), symmetry=1, barrier=(7.48545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145947,'amu*angstrom^2'), symmetry=1, barrier=(7.485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145974,'amu*angstrom^2'), symmetry=1, barrier=(7.48528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.502114,'amu*angstrom^2'), symmetry=1, barrier=(25.7505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.51509,0.084294,-0.000105142,7.81748e-08,-2.47752e-11,312.598,34.7525], Tmin=(100,'K'), Tmax=(754.919,'K')), NASAPolynomial(coeffs=[8.46153,0.0421885,-2.14784e-05,4.29113e-09,-3.07505e-13,-887.17,-1.35025], Tmin=(754.919,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.61343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CC([O])[C](C=O)C=C=O(7285)',
    structure = SMILES('C=CC([O])C(C=O)=C[C]=O'),
    E0 = (47.6049,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.681551,'amu*angstrom^2'), symmetry=1, barrier=(15.6702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681478,'amu*angstrom^2'), symmetry=1, barrier=(15.6685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681425,'amu*angstrom^2'), symmetry=1, barrier=(15.6673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.68133,'amu*angstrom^2'), symmetry=1, barrier=(15.6651,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0426419,0.0995826,-0.000129219,9.47905e-08,-2.93472e-11,5861.6,35.4815], Tmin=(100,'K'), Tmax=(772.801,'K')), NASAPolynomial(coeffs=[10.2915,0.0460976,-2.54137e-05,5.24945e-09,-3.83241e-13,4264.21,-11.7123], Tmin=(772.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.6049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC(=O)[C](C=O)C[C]=O(7286)',
    structure = SMILES('C=CC([O])=C(C=O)C[C]=O'),
    E0 = (-88.3002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,325,375,415,465,420,450,1700,1750,1855,455,950,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,295.567,295.57,295.583],'cm^-1')),
        HinderedRotor(inertia=(0.236863,'amu*angstrom^2'), symmetry=1, barrier=(14.6817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236799,'amu*angstrom^2'), symmetry=1, barrier=(14.6816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236867,'amu*angstrom^2'), symmetry=1, barrier=(14.6819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236833,'amu*angstrom^2'), symmetry=1, barrier=(14.6818,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.126258,0.0878599,-8.57005e-05,4.15551e-08,-8.16247e-12,-10482.9,32.397], Tmin=(100,'K'), Tmax=(1205.08,'K')), NASAPolynomial(coeffs=[16.8223,0.0324427,-1.67228e-05,3.39671e-09,-2.46529e-13,-14507,-51.2666], Tmin=(1205.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-88.3002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC(=O)[C](C=C=O)C[O](7287)',
    structure = SMILES('C=CC([O])=C(C=C=O)C[O]'),
    E0 = (34.6121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2120,512.5,787.5,325,375,415,465,420,450,1700,1750,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.303678,'amu*angstrom^2'), symmetry=1, barrier=(6.98216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303712,'amu*angstrom^2'), symmetry=1, barrier=(6.98294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303713,'amu*angstrom^2'), symmetry=1, barrier=(6.98295,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.66692,0.112049,-0.00017674,1.50277e-07,-4.99469e-11,4321.92,35.5933], Tmin=(100,'K'), Tmax=(853.205,'K')), NASAPolynomial(coeffs=[11.5712,0.0381344,-1.7715e-05,3.29961e-09,-2.23137e-13,2835.61,-17.9777], Tmin=(853.205,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.6121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cd-CdCs(CCO)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCOJ)"""),
)

species(
    label = 'C=CC(=O)C([C]=C[O])C=O(7288)',
    structure = SMILES('C=CC(=O)C([C]=C[O])C=O'),
    E0 = (9.42799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1685,370,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,706.868],'cm^-1')),
        HinderedRotor(inertia=(0.0842462,'amu*angstrom^2'), symmetry=1, barrier=(1.93698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.764663,'amu*angstrom^2'), symmetry=1, barrier=(17.5811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0495648,'amu*angstrom^2'), symmetry=1, barrier=(17.5809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0495743,'amu*angstrom^2'), symmetry=1, barrier=(17.581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.450287,0.0804903,-7.36859e-05,3.51488e-08,-6.93107e-12,1259.62,35.0325], Tmin=(100,'K'), Tmax=(1188.41,'K')), NASAPolynomial(coeffs=[13.8025,0.0355489,-1.69613e-05,3.32784e-09,-2.37055e-13,-1913.97,-31.6891], Tmin=(1188.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.42799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C([O])C(C=O)C=C=O(7289)',
    structure = SMILES('C=[C]C([O])C(C=O)C=C=O'),
    E0 = (159.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2120,512.5,787.5,1685,370,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,387.55,387.551,387.551,387.551],'cm^-1')),
        HinderedRotor(inertia=(0.00112239,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13636,'amu*angstrom^2'), symmetry=1, barrier=(14.5335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13636,'amu*angstrom^2'), symmetry=1, barrier=(14.5335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136362,'amu*angstrom^2'), symmetry=1, barrier=(14.5335,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.20343,0.0860055,-8.80851e-05,4.76795e-08,-1.05466e-11,19345.9,37.7084], Tmin=(100,'K'), Tmax=(1079.37,'K')), NASAPolynomial(coeffs=[14.2198,0.0340622,-1.58989e-05,3.09378e-09,-2.19767e-13,16320.2,-30.9831], Tmin=(1079.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=CC([O])C([C]=C=O)C=O(7290)',
    structure = SMILES('C=CC([O])C([C]=C=O)C=O'),
    E0 = (127.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2120,512.5,787.5,1685,370,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,180,180,464.085,464.085],'cm^-1')),
        HinderedRotor(inertia=(0.0155059,'amu*angstrom^2'), symmetry=1, barrier=(2.36985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1304,'amu*angstrom^2'), symmetry=1, barrier=(19.9297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1304,'amu*angstrom^2'), symmetry=1, barrier=(19.9297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.866811,'amu*angstrom^2'), symmetry=1, barrier=(19.9297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.164667,0.0856711,-8.46243e-05,4.37108e-08,-9.20504e-12,15493.5,37.2013], Tmin=(100,'K'), Tmax=(1131.45,'K')), NASAPolynomial(coeffs=[14.8997,0.0335786,-1.55639e-05,3.01956e-09,-2.14138e-13,12159.1,-35.7067], Tmin=(1131.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)CJ=C=O)"""),
)

species(
    label = 'C=CC(=O)C([C]=C=O)C[O](7291)',
    structure = SMILES('C=CC(=O)C([C]=C=O)C[O]'),
    E0 = (107.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,2120,512.5,787.5,1685,370,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.030548,0.101809,-0.000165918,1.589e-07,-5.98174e-11,13025,36.7651], Tmin=(100,'K'), Tmax=(813.81,'K')), NASAPolynomial(coeffs=[4.52033,0.0529937,-2.71944e-05,5.32504e-09,-3.72086e-13,13160.1,21.1281], Tmin=(813.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(CC(C)CJ=C=O) + radical(CCOJ)"""),
)

species(
    label = 'C=CC([O])C([C]=O)C=C=O(7292)',
    structure = SMILES('C=CC([O])C([C]=O)C=C=O'),
    E0 = (80.5782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1855,455,950,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,387.836,387.87,387.879,387.987],'cm^-1')),
        HinderedRotor(inertia=(0.00112061,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150428,'amu*angstrom^2'), symmetry=1, barrier=(16.065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150481,'amu*angstrom^2'), symmetry=1, barrier=(16.0649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150514,'amu*angstrom^2'), symmetry=1, barrier=(16.0657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0324131,0.0839962,-7.9153e-05,3.76486e-08,-7.14158e-12,9840.35,39.231], Tmin=(100,'K'), Tmax=(1267.09,'K')), NASAPolynomial(coeffs=[17.9677,0.0271719,-1.18825e-05,2.25437e-09,-1.58095e-13,5278.88,-51.8697], Tmin=(1267.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.5782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)CJ=O) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC(=O)C([C]=O)C[C]=O(7293)',
    structure = SMILES('C=CC(=O)C([C]=O)C[C]=O'),
    E0 = (-51.8814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1850,1860,440,470,900,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.761098,0.0810516,-8.41032e-05,3.60478e-08,4.54349e-12,-6132.66,36.4515], Tmin=(100,'K'), Tmax=(554.03,'K')), NASAPolynomial(coeffs=[7.50079,0.0453843,-2.27117e-05,4.5017e-09,-3.21049e-13,-7078.85,6.11695], Tmin=(554.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.8814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)CJ=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]CC(=O)[C](C=O)C=C=O(7294)',
    structure = SMILES('[CH2]CC([O])=C(C=O)C=C=O'),
    E0 = (-31.5993,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,2120,512.5,787.5,3010,987.5,1337.5,450,1655,325,375,415,465,420,450,1700,1750,257.199,257.218,257.239],'cm^-1')),
        HinderedRotor(inertia=(0.173525,'amu*angstrom^2'), symmetry=1, barrier=(8.14393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173695,'amu*angstrom^2'), symmetry=1, barrier=(8.14421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173489,'amu*angstrom^2'), symmetry=1, barrier=(8.14391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173631,'amu*angstrom^2'), symmetry=1, barrier=(8.14439,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.532827,0.10779,-0.000157847,1.26405e-07,-4.09146e-11,-3644.86,36.5957], Tmin=(100,'K'), Tmax=(754.364,'K')), NASAPolynomial(coeffs=[12.6343,0.0379754,-1.90323e-05,3.73525e-09,-2.63279e-13,-5631.52,-23.2173], Tmin=(754.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.5993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=CC(=O)[C](C=O)C=C[O](7295)',
    structure = SMILES('[CH2]C=C([O])C(C=O)=CC=O'),
    E0 = (-137.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.902929,'amu*angstrom^2'), symmetry=1, barrier=(20.7601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.901484,'amu*angstrom^2'), symmetry=1, barrier=(20.7269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.90299,'amu*angstrom^2'), symmetry=1, barrier=(20.7615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.901839,'amu*angstrom^2'), symmetry=1, barrier=(20.735,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.165622,0.0929548,-9.56513e-05,4.96975e-08,-1.03892e-11,-16397.9,31.3766], Tmin=(100,'K'), Tmax=(1146.06,'K')), NASAPolynomial(coeffs=[17.3721,0.0317427,-1.55331e-05,3.09131e-09,-2.22365e-13,-20417.7,-55.6234], Tmin=(1146.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-137.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=CC([O])C(C=O)C=C=O(7296)',
    structure = SMILES('[CH]=CC([O])C(C=O)C=C=O'),
    E0 = (168.986,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,334.96,335.418,336.204],'cm^-1')),
        HinderedRotor(inertia=(0.00150206,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203943,'amu*angstrom^2'), symmetry=1, barrier=(16.3279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204823,'amu*angstrom^2'), symmetry=1, barrier=(16.3255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205189,'amu*angstrom^2'), symmetry=1, barrier=(16.3294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0281005,0.0865829,-8.59959e-05,4.38771e-08,-9.01344e-12,20468,38.2616], Tmin=(100,'K'), Tmax=(1167.46,'K')), NASAPolynomial(coeffs=[16.4227,0.0304119,-1.38262e-05,2.66596e-09,-1.88619e-13,16639.9,-43.3716], Tmin=(1167.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C[C](O)C([C]=C=O)C=O(7297)',
    structure = SMILES('[CH2]C=C(O)C([C]=C=O)C=O'),
    E0 = (-6.31224,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,350,440,435,1725,2120,512.5,787.5,3615,1277.5,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,1685,370,321.168,321.193],'cm^-1')),
        HinderedRotor(inertia=(0.232359,'amu*angstrom^2'), symmetry=1, barrier=(16.9984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232362,'amu*angstrom^2'), symmetry=1, barrier=(16.9976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232311,'amu*angstrom^2'), symmetry=1, barrier=(16.9987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.509015,'amu*angstrom^2'), symmetry=1, barrier=(37.2575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232072,'amu*angstrom^2'), symmetry=1, barrier=(16.9966,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.530865,0.0977774,-0.000109532,6.22018e-08,-1.393e-11,-594.334,35.8522], Tmin=(100,'K'), Tmax=(1089.82,'K')), NASAPolynomial(coeffs=[18.9104,0.0264215,-1.132e-05,2.12332e-09,-1.48288e-13,-4831.85,-59.6131], Tmin=(1089.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.31224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(CC(C)CJ=C=O) + radical(Allyl_P)"""),
)

species(
    label = 'C=CC(=O)C([C]=C=O)[CH]O(7298)',
    structure = SMILES('C=CC(=O)C([C]=C=O)[CH]O'),
    E0 = (61.7885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,2120,512.5,787.5,3025,407.5,1350,352.5,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.615157,0.112687,-0.000182365,1.62253e-07,-5.68688e-11,7586.85,37.5551], Tmin=(100,'K'), Tmax=(819.783,'K')), NASAPolynomial(coeffs=[9.78476,0.0434592,-2.18778e-05,4.23857e-09,-2.9401e-13,6502.77,-6.76397], Tmin=(819.783,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.7885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(CC(C)CJ=C=O) + radical(CCsJOH)"""),
)

species(
    label = 'C=C[C](O)C([C]=O)C=C=O(7299)',
    structure = SMILES('[CH2]C=C(O)C([C]=O)C=C=O'),
    E0 = (-53.4154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,2120,512.5,787.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,350,440,435,1725,245.261,415.637],'cm^-1')),
        HinderedRotor(inertia=(0.735744,'amu*angstrom^2'), symmetry=1, barrier=(16.9162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138001,'amu*angstrom^2'), symmetry=1, barrier=(16.9161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137997,'amu*angstrom^2'), symmetry=1, barrier=(16.9161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137991,'amu*angstrom^2'), symmetry=1, barrier=(16.9161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137998,'amu*angstrom^2'), symmetry=1, barrier=(16.9161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.702452,0.0957969,-0.000103001,5.48275e-08,-1.13492e-11,-6248.56,37.7913], Tmin=(100,'K'), Tmax=(1186.34,'K')), NASAPolynomial(coeffs=[21.5807,0.0206638,-8.00213e-06,1.44219e-09,-9.91073e-14,-11535.6,-73.5193], Tmin=(1186.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.4154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C=CC(=O)C([C]=O)[CH]C=O(7300)',
    structure = SMILES('C=CC(=O)C([C]=O)C=C[O]'),
    E0 = (-69.7252,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1855,455,950,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,375,552.5,462.5,1710,469.806,469.808,469.815,469.818],'cm^-1')),
        HinderedRotor(inertia=(0.0144806,'amu*angstrom^2'), symmetry=1, barrier=(2.26817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812032,'amu*angstrom^2'), symmetry=1, barrier=(18.6702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812034,'amu*angstrom^2'), symmetry=1, barrier=(18.6703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119199,'amu*angstrom^2'), symmetry=1, barrier=(18.6702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.101978,0.0796834,-6.84046e-05,2.91313e-08,-4.96038e-12,-8240.86,36.967], Tmin=(100,'K'), Tmax=(1395.95,'K')), NASAPolynomial(coeffs=[18.3125,0.0275024,-1.23343e-05,2.35366e-09,-1.6479e-13,-13325.1,-56.9628], Tmin=(1395.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.7252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C[CH]C(=O)[C](C=O)C=C=O(7301)',
    structure = SMILES('CC=C([O])C(C=C=O)=C[O]'),
    E0 = (-97.7525,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2750,2800,2850,1350,1500,750,1050,1375,1000,2120,512.5,787.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.74452,'amu*angstrom^2'), symmetry=1, barrier=(17.118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743002,'amu*angstrom^2'), symmetry=1, barrier=(17.0831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.74365,'amu*angstrom^2'), symmetry=1, barrier=(17.098,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.909413,0.108999,-0.000141829,9.36212e-08,-2.41351e-11,-11580.9,33.7338], Tmin=(100,'K'), Tmax=(955.429,'K')), NASAPolynomial(coeffs=[19.0989,0.0252342,-1.03248e-05,1.86444e-09,-1.26427e-13,-15404.4,-61.883], Tmin=(955.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.7525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC(=O)[C](C=O)C=[C]O(7302)',
    structure = SMILES('C=CC(=O)[C](C=O)C=[C]O'),
    E0 = (-10.4317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,360,370,350,2950,3100,1380,975,1025,1650,1685,370,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0752526,0.0775456,-6.29425e-05,2.50799e-08,-3.95119e-12,-1097.87,37.6397], Tmin=(100,'K'), Tmax=(1519.54,'K')), NASAPolynomial(coeffs=[20.3297,0.0238318,-9.91914e-06,1.8169e-09,-1.23855e-13,-7299.07,-69.34], Tmin=(1519.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.4317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CCJ(C)C=O)"""),
)

species(
    label = '[CH]=C[C](O)C(C=O)C=C=O(7303)',
    structure = SMILES('[CH]C=C(O)C(C=O)C=C=O'),
    E0 = (7.08122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.491527,0.0963684,-9.8086e-05,5.20471e-08,-1.10797e-11,1015.5,36.7075], Tmin=(100,'K'), Tmax=(1134.49,'K')), NASAPolynomial(coeffs=[17.5666,0.0326984,-1.39024e-05,2.57758e-09,-1.78319e-13,-3081.83,-52.6911], Tmin=(1134.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.08122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]CC(=O)C([C]=C=O)C=O(7304)',
    structure = SMILES('[CH2]CC(=O)C([C]=C=O)C=O'),
    E0 = (40.1979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,2120,512.5,787.5,1685,370,1380,1390,370,380,2900,435,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.316795,0.103736,-0.000155639,1.33505e-07,-4.61533e-11,4981.69,37.3514], Tmin=(100,'K'), Tmax=(798.512,'K')), NASAPolynomial(coeffs=[9.64265,0.0429633,-2.10347e-05,4.05727e-09,-2.82028e-13,3738.09,-6.28378], Tmin=(798.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.1979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(CJCC=O) + radical(CC(C)CJ=C=O)"""),
)

species(
    label = '[CH2]CC(=O)C([C]=O)C=C=O(7305)',
    structure = SMILES('[CH2]CC(=O)C([C]=O)C=C=O'),
    E0 = (-6.90527,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,2120,512.5,787.5,1855,455,950,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0803582,0.0967117,-0.000130311,9.99299e-08,-3.14815e-11,-690.012,37.8416], Tmin=(100,'K'), Tmax=(770.536,'K')), NASAPolynomial(coeffs=[10.9246,0.0395816,-1.90941e-05,3.70278e-09,-2.60024e-13,-2385.92,-12.382], Tmin=(770.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.90527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(CC(C)CJ=O) + radical(CJCC=O)"""),
)

species(
    label = 'C=[C]C(=O)C(C=O)C[C]=O(7306)',
    structure = SMILES('C=C=C([O])C(C=O)C[C]=O'),
    E0 = (-12.7908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,350,440,435,1725,540,610,2055,1855,455,950,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,187.748,188.529,188.562],'cm^-1')),
        HinderedRotor(inertia=(0.292518,'amu*angstrom^2'), symmetry=1, barrier=(7.39112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294783,'amu*angstrom^2'), symmetry=1, barrier=(7.39506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294332,'amu*angstrom^2'), symmetry=1, barrier=(7.39138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.293699,'amu*angstrom^2'), symmetry=1, barrier=(7.39382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.533863,0.107698,-0.000164113,1.36388e-07,-4.4811e-11,-1382.75,37.3364], Tmin=(100,'K'), Tmax=(837.072,'K')), NASAPolynomial(coeffs=[11.8999,0.0369815,-1.71407e-05,3.20742e-09,-2.18263e-13,-3068.43,-18.0733], Tmin=(837.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.7908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]C(=O)C(C=C=O)C[O](7307)',
    structure = SMILES('C=C=C([O])C(C=C=O)C[O]'),
    E0 = (97.8918,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2120,512.5,787.5,540,610,2055,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,270.982,270.985,271,271.013,271.018],'cm^-1')),
        HinderedRotor(inertia=(0.167765,'amu*angstrom^2'), symmetry=1, barrier=(8.74605,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167839,'amu*angstrom^2'), symmetry=1, barrier=(8.74599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167834,'amu*angstrom^2'), symmetry=1, barrier=(8.74657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.637427,0.109843,-0.000166354,1.37289e-07,-4.49145e-11,11933.2,36.1473], Tmin=(100,'K'), Tmax=(830.08,'K')), NASAPolynomial(coeffs=[12.3685,0.0372765,-1.73434e-05,3.2551e-09,-2.22093e-13,10114.8,-22.1232], Tmin=(830.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.8918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCd(CCO)H) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(CCOJ)"""),
)

species(
    label = 'C[CH]C(=O)C([C]=C=O)C=O(7308)',
    structure = SMILES('CC=C([O])C([C]=C=O)C=O'),
    E0 = (-20.0066,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,343.91,343.915,343.919],'cm^-1')),
        HinderedRotor(inertia=(0.114706,'amu*angstrom^2'), symmetry=1, barrier=(9.62729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114708,'amu*angstrom^2'), symmetry=1, barrier=(9.62742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194124,'amu*angstrom^2'), symmetry=1, barrier=(16.2932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114699,'amu*angstrom^2'), symmetry=1, barrier=(9.62725,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0569721,0.0948866,-0.000120096,8.49506e-08,-2.45472e-11,-2265.2,35.3828], Tmin=(100,'K'), Tmax=(839.39,'K')), NASAPolynomial(coeffs=[12.0146,0.0373615,-1.72985e-05,3.30669e-09,-2.30959e-13,-4291.77,-20.7422], Tmin=(839.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.0066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(CC(C)CJ=C=O) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[CH]C(=O)C([C]=O)C=C=O(7309)',
    structure = SMILES('CC=C([O])C([C]=O)C=C=O'),
    E0 = (-67.1098,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2120,512.5,787.5,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,258.703,258.71,258.72],'cm^-1')),
        HinderedRotor(inertia=(0.00251852,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218503,'amu*angstrom^2'), symmetry=1, barrier=(10.3772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218511,'amu*angstrom^2'), symmetry=1, barrier=(10.3773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218492,'amu*angstrom^2'), symmetry=1, barrier=(10.3772,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0670558,0.0894428,-0.000101565,6.21387e-08,-1.54052e-11,-7932.24,36.2602], Tmin=(100,'K'), Tmax=(975.314,'K')), NASAPolynomial(coeffs=[13.828,0.0330051,-1.47648e-05,2.80652e-09,-1.96516e-13,-10616.5,-29.7841], Tmin=(975.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.1098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(CC(C)CJ=O) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC(=O)C([C]=O)C=[C]O(7310)',
    structure = SMILES('C=CC(=O)C([C]=O)C=[C]O'),
    E0 = (28.5564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,1855,455,950,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.223603,0.0846145,-8.56265e-05,4.52951e-08,-9.74357e-12,3569.27,38.2506], Tmin=(100,'K'), Tmax=(1110.97,'K')), NASAPolynomial(coeffs=[14.7502,0.0323113,-1.5007e-05,2.91729e-09,-2.07192e-13,341.605,-33.3603], Tmin=(1110.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.5564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C=[C]C(=O)C([CH]C=O)C=O(7311)',
    structure = SMILES('C=C=C([O])C(C=O)C=C[O]'),
    E0 = (-30.6346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.747094,'amu*angstrom^2'), symmetry=1, barrier=(17.1772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.746779,'amu*angstrom^2'), symmetry=1, barrier=(17.1699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.747225,'amu*angstrom^2'), symmetry=1, barrier=(17.1802,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.578984,0.0977675,-0.000110446,6.32784e-08,-1.42314e-11,-3517,35.7245], Tmin=(100,'K'), Tmax=(1088.88,'K')), NASAPolynomial(coeffs=[19.248,0.0249328,-1.01108e-05,1.84814e-09,-1.27339e-13,-7834.81,-61.6173], Tmin=(1088.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.6346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C]C(=O)C([CH]O)C=C=O(7312)',
    structure = SMILES('C=C=C([O])C([CH]O)C=C=O'),
    E0 = (52.4845,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2120,512.5,787.5,3025,407.5,1350,352.5,540,610,2055,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.550515,'amu*angstrom^2'), symmetry=1, barrier=(12.6574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.550002,'amu*angstrom^2'), symmetry=1, barrier=(12.6456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.550276,'amu*angstrom^2'), symmetry=1, barrier=(12.6519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549995,'amu*angstrom^2'), symmetry=1, barrier=(12.6455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.2522,0.121118,-0.00018442,1.43137e-07,-4.32455e-11,6496.26,37.043], Tmin=(100,'K'), Tmax=(860.841,'K')), NASAPolynomial(coeffs=[17.6729,0.0276689,-1.19825e-05,2.1578e-09,-1.43094e-13,3442.18,-50.2376], Tmin=(860.841,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.4845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCd(CCO)H) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=CC(=O)C(C=O)C[C]=O(7313)',
    structure = SMILES('[CH]=CC(=O)C(C=O)C[C]=O'),
    E0 = (36.5261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1855,455,950,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,1380,1390,370,380,2900,435,375,552.5,462.5,1710,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0930538,'amu*angstrom^2'), symmetry=1, barrier=(2.13949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0930445,'amu*angstrom^2'), symmetry=1, barrier=(2.13928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0930375,'amu*angstrom^2'), symmetry=1, barrier=(2.13911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0930424,'amu*angstrom^2'), symmetry=1, barrier=(2.13923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.468879,'amu*angstrom^2'), symmetry=1, barrier=(10.7804,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11439,0.0771946,-4.63234e-05,-7.41038e-08,1.02473e-10,4482.74,34.583], Tmin=(100,'K'), Tmax=(463.828,'K')), NASAPolynomial(coeffs=[7.54036,0.0458891,-2.30583e-05,4.53233e-09,-3.19726e-13,3627.27,5.72218], Tmin=(463.828,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.5261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(=O)C(C=C=O)C[O](7314)',
    structure = SMILES('[CH]=CC(=O)C(C=C=O)C[O]'),
    E0 = (148.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,375,552.5,462.5,1710,180,180,857.011,857.08],'cm^-1')),
        HinderedRotor(inertia=(0.238272,'amu*angstrom^2'), symmetry=1, barrier=(5.47834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238312,'amu*angstrom^2'), symmetry=1, barrier=(5.47926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238193,'amu*angstrom^2'), symmetry=1, barrier=(5.47652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23817,'amu*angstrom^2'), symmetry=1, barrier=(5.47601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0616825,0.101439,-0.000162634,1.5281e-07,-5.68682e-11,17994.9,37.45], Tmin=(100,'K'), Tmax=(807.888,'K')), NASAPolynomial(coeffs=[5.61778,0.0505493,-2.58727e-05,5.06961e-09,-3.54706e-13,17820.3,15.8603], Tmin=(807.888,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C(=O)C(C=O)C=[C]O(7315)',
    structure = SMILES('C=C=C([O])C(C=O)C=[C]O'),
    E0 = (67.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,540,610,2055,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.604707,'amu*angstrom^2'), symmetry=1, barrier=(13.9034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.604748,'amu*angstrom^2'), symmetry=1, barrier=(13.9043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.605414,'amu*angstrom^2'), symmetry=1, barrier=(13.9197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.605806,'amu*angstrom^2'), symmetry=1, barrier=(13.9287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.743573,0.106002,-0.000138777,9.30676e-08,-2.44505e-11,8305.44,38.038], Tmin=(100,'K'), Tmax=(936.808,'K')), NASAPolynomial(coeffs=[17.9875,0.0260239,-1.07181e-05,1.93713e-09,-1.31189e-13,4795.93,-51.1062], Tmin=(936.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=CC(=O)C([CH]C=O)C=O(7316)',
    structure = SMILES('[CH]=CC(=O)C(C=O)C=C[O]'),
    E0 = (18.6824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3120,650,792.5,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,375,552.5,462.5,1710,522.007,522.532,523.419],'cm^-1')),
        HinderedRotor(inertia=(0.0621163,'amu*angstrom^2'), symmetry=1, barrier=(1.58824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.76098,'amu*angstrom^2'), symmetry=1, barrier=(19.0437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.750697,'amu*angstrom^2'), symmetry=1, barrier=(19.0392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.746328,'amu*angstrom^2'), symmetry=1, barrier=(19.0398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.222273,0.0816459,-7.34123e-05,3.34151e-08,-6.16407e-12,2383.99,35.7775], Tmin=(100,'K'), Tmax=(1281.86,'K')), NASAPolynomial(coeffs=[16.322,0.0314074,-1.46245e-05,2.84094e-09,-2.01213e-13,-1743.52,-45.8919], Tmin=(1281.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.6824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CC(=O)C([CH]O)C=C=O(7317)',
    structure = SMILES('[CH]=CC(=O)C([CH]O)C=C=O'),
    E0 = (103.093,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,2120,512.5,787.5,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,3120,650,792.5,1650,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.64463,0.112293,-0.000178985,1.56006e-07,-5.38339e-11,12556.7,38.2343], Tmin=(100,'K'), Tmax=(812.343,'K')), NASAPolynomial(coeffs=[10.8817,0.0410158,-2.05568e-05,3.98333e-09,-2.76646e-13,11163.2,-12.0289], Tmin=(812.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=CC(=O)C(C=O)C=[C]O(7318)',
    structure = SMILES('[CH]=CC(=O)C(C=O)C=[C]O'),
    E0 = (116.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1685,370,375,552.5,462.5,1710,3120,650,792.5,1650,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,269.234,269.261],'cm^-1')),
        HinderedRotor(inertia=(0.343906,'amu*angstrom^2'), symmetry=1, barrier=(17.6901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486816,'amu*angstrom^2'), symmetry=1, barrier=(25.0419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244721,'amu*angstrom^2'), symmetry=1, barrier=(12.5882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486821,'amu*angstrom^2'), symmetry=1, barrier=(25.042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0011087,'amu*angstrom^2'), symmetry=1, barrier=(12.5882,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.162624,0.088533,-9.66671e-05,5.63626e-08,-1.34428e-11,14202.3,37.7238], Tmin=(100,'K'), Tmax=(1005.05,'K')), NASAPolynomial(coeffs=[13.6428,0.0348832,-1.65969e-05,3.25081e-09,-2.31584e-13,11492.6,-27.3783], Tmin=(1005.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(=O)C([C]C=O)C=O(7326)',
    structure = SMILES('C=CC(=O)C([C]C=O)C=O'),
    E0 = (59.9236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53806,0.0676866,-4.54339e-06,-1.47938e-07,1.58785e-10,7281.54,42.0782], Tmin=(100,'K'), Tmax=(427.462,'K')), NASAPolynomial(coeffs=[5.51323,0.05222,-2.65261e-05,5.27511e-09,-3.76314e-13,6743.15,23.9565], Tmin=(427.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.9236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C[C]C(=O)C(C=O)C=C=O(7327)',
    structure = SMILES('C[C]C(=O)C(C=O)C=C=O'),
    E0 = (53.2712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,2120,512.5,787.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.392603,0.0888746,-9.87191e-05,5.30417e-08,-4.79643e-12,6527.97,44.5841], Tmin=(100,'K'), Tmax=(585.549,'K')), NASAPolynomial(coeffs=[8.57679,0.0460432,-2.24957e-05,4.39723e-09,-3.10821e-13,5345.34,7.5661], Tmin=(585.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.2712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]CC(=O)C(C=O)C=C=O(7328)',
    structure = SMILES('[CH]CC(=O)C(C=O)C=C=O'),
    E0 = (75.7072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,375,552.5,462.5,1710,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.194869,0.091177,-0.000116518,8.79739e-08,-2.79514e-11,9235.64,36.0883], Tmin=(100,'K'), Tmax=(757.243,'K')), NASAPolynomial(coeffs=[9.29,0.0431373,-2.13649e-05,4.20888e-09,-2.98896e-13,7858.09,-5.26206], Tmin=(757.243,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.7072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=C[C]1OC(=O)[CH]C1C=O(7329)',
    structure = SMILES('[CH2]C=C1OC(=O)[CH]C1C=O'),
    E0 = (-114.051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.687845,0.0592606,-1.2895e-05,-2.38759e-08,1.24687e-11,-13586.2,32.9019], Tmin=(100,'K'), Tmax=(1044.08,'K')), NASAPolynomial(coeffs=[15.8536,0.0300711,-1.24969e-05,2.39263e-09,-1.71929e-13,-18328.9,-48.465], Tmin=(1044.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + ring(Cyclopentane) + radical(Allyl_P) + radical(CCJCO)"""),
)

species(
    label = 'C=C[C]1OO[CH]C1C=C=O(7330)',
    structure = SMILES('[CH2]C=C1OO[CH]C1C=C=O'),
    E0 = (266.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2120,512.5,787.5,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.146115,0.0719254,-3.94559e-05,-2.27045e-09,6.20509e-12,32188.8,33.1935], Tmin=(100,'K'), Tmax=(1049.88,'K')), NASAPolynomial(coeffs=[18.2165,0.027647,-1.12965e-05,2.13872e-09,-1.52676e-13,27040.4,-61.3137], Tmin=(1049.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCd(CCO)H) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + ring(Cyclopentane) + radical(CCsJOOC) + radical(Allyl_P)"""),
)

species(
    label = 'C=CC(=O)C1[CH]OC(=O)[CH]1(7331)',
    structure = SMILES('C=CC(=O)C1[CH]OC(=O)[CH]1'),
    E0 = (-49.2025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.589526,0.0648733,-3.22062e-05,-3.41906e-09,5.7658e-12,-5786.12,34.4052], Tmin=(100,'K'), Tmax=(1027.83,'K')), NASAPolynomial(coeffs=[14.3595,0.0311024,-1.18431e-05,2.13232e-09,-1.473e-13,-9663.56,-37.4973], Tmin=(1027.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.2025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-OdCsOs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + ring(butyrolactone) + radical(CCsJOC(O)) + radical(CCJCO)"""),
)

species(
    label = 'S(559)(558)',
    structure = SMILES('[O]C1=CC(C=O)C([O])=CC1'),
    E0 = (-175.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(5198.54,'J/mol'), sigma=(7.85497,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=812.00 K, Pc=24.34 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0726473,0.0761439,-5.6658e-05,1.49843e-08,8.44305e-13,-20965.8,31.9487], Tmin=(100,'K'), Tmax=(1029.31,'K')), NASAPolynomial(coeffs=[17.7739,0.0256043,-9.60138e-06,1.73113e-09,-1.20277e-13,-25576.5,-58.6568], Tmin=(1029.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-175.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cyclohexadiene) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=C=CC1[CH]OC[CH]C1=O(7332)',
    structure = SMILES('[O]C1=CCO[CH]C1C=C=O'),
    E0 = (-31.8328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,2120,512.5,787.5,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.993459,0.0925296,-9.44368e-05,4.84772e-08,-9.44776e-12,-3634.03,31.1332], Tmin=(100,'K'), Tmax=(1415.48,'K')), NASAPolynomial(coeffs=[21.9875,0.0170605,-3.30538e-06,3.01655e-10,-1.10708e-14,-9085.2,-83.9962], Tmin=(1415.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.8328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCd(CCO)H) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + ring(36dihydro2hpyran) + radical(CCsJOCs) + radical(C=C(C)OJ)"""),
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
    label = 'C=CC(=O)CC=C=O(7210)',
    structure = SMILES('C=CC(=O)CC=C=O'),
    E0 = (-131.28,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,807.891],'cm^-1')),
        HinderedRotor(inertia=(0.133365,'amu*angstrom^2'), symmetry=1, barrier=(3.06632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130124,'amu*angstrom^2'), symmetry=1, barrier=(2.99182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627718,'amu*angstrom^2'), symmetry=1, barrier=(14.4325,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63965,0.0584724,-5.00744e-05,2.61513e-08,-6.411e-12,-15709.9,25.5698], Tmin=(100,'K'), Tmax=(904.781,'K')), NASAPolynomial(coeffs=[5.74413,0.0403268,-1.9992e-05,3.98615e-09,-2.86619e-13,-16452.6,6.17864], Tmin=(904.781,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(=CC=O)OC=C=O(7319)',
    structure = SMILES('C=CC(=CC=O)OC=C=O'),
    E0 = (-134.319,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2120,512.5,787.5,2782.5,750,1395,475,1775,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.805197,'amu*angstrom^2'), symmetry=1, barrier=(18.5131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805068,'amu*angstrom^2'), symmetry=1, barrier=(18.5101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.806242,'amu*angstrom^2'), symmetry=1, barrier=(18.5371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805754,'amu*angstrom^2'), symmetry=1, barrier=(18.5259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.148678,0.0948826,-0.000106066,6.14051e-08,-1.43487e-11,-16008.4,31.7042], Tmin=(100,'K'), Tmax=(1031.01,'K')), NASAPolynomial(coeffs=[15.8089,0.0329732,-1.59967e-05,3.16568e-09,-2.26957e-13,-19299,-45.7693], Tmin=(1031.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)OsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(=CC=C=O)OC=O(7320)',
    structure = SMILES('C=CC(=CC=C=O)OC=O'),
    E0 = (-261.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.857046,'amu*angstrom^2'), symmetry=1, barrier=(19.7052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.857093,'amu*angstrom^2'), symmetry=1, barrier=(19.7062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.857102,'amu*angstrom^2'), symmetry=1, barrier=(19.7065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.857104,'amu*angstrom^2'), symmetry=1, barrier=(19.7065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0631758,0.0908337,-9.84299e-05,5.6709e-08,-1.33676e-11,-31284.4,30.2413], Tmin=(100,'K'), Tmax=(1015.79,'K')), NASAPolynomial(coeffs=[14.039,0.0357994,-1.71614e-05,3.3721e-09,-2.40648e-13,-34123.7,-37.4027], Tmin=(1015.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-261.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cd-Cd(CCO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + group(Cds-OdOsH)"""),
)

species(
    label = 'C=CC(O)=C(C=O)C=C=O(7321)',
    structure = SMILES('C=CC(O)=C(C=O)C=C=O'),
    E0 = (-263.542,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.833242,'amu*angstrom^2'), symmetry=1, barrier=(19.1579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834038,'amu*angstrom^2'), symmetry=1, barrier=(19.1762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833755,'amu*angstrom^2'), symmetry=1, barrier=(19.1697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833747,'amu*angstrom^2'), symmetry=1, barrier=(19.1695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.821206,0.110105,-0.000142781,9.32288e-08,-2.39959e-11,-31526.5,31.4523], Tmin=(100,'K'), Tmax=(951.02,'K')), NASAPolynomial(coeffs=[18.6511,0.028203,-1.35995e-05,2.67162e-09,-1.9039e-13,-35230.2,-61.5122], Tmin=(951.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-CdCs(CCO)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C(O)C(C=O)C=C=O(7322)',
    structure = SMILES('C=C=C(O)C(C=O)C=C=O'),
    E0 = (-186.999,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,540,610,2055,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.784198,'amu*angstrom^2'), symmetry=1, barrier=(18.0302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784202,'amu*angstrom^2'), symmetry=1, barrier=(18.0304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784532,'amu*angstrom^2'), symmetry=1, barrier=(18.0379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784479,'amu*angstrom^2'), symmetry=1, barrier=(18.0367,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.449132,0.0960405,-0.000107332,6.0836e-08,-1.35997e-11,-22328.9,34.5013], Tmin=(100,'K'), Tmax=(1091.57,'K')), NASAPolynomial(coeffs=[18.6518,0.0260469,-1.11501e-05,2.09459e-09,-1.4646e-13,-26499,-59.3237], Tmin=(1091.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-186.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cds-CdsCsOs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC(=O)OC=CC=C=O(7323)',
    structure = SMILES('C=CC(=O)OC=CC=C=O'),
    E0 = (-174.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2120,512.5,787.5,180,180,180,437.675,714.075,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.166427,'amu*angstrom^2'), symmetry=1, barrier=(3.82649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166427,'amu*angstrom^2'), symmetry=1, barrier=(3.82649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166427,'amu*angstrom^2'), symmetry=1, barrier=(3.82649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166427,'amu*angstrom^2'), symmetry=1, barrier=(3.82649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.270416,0.101587,-0.00014796,1.21659e-07,-4.07318e-11,-20787.6,34.3227], Tmin=(100,'K'), Tmax=(756.086,'K')), NASAPolynomial(coeffs=[11.0112,0.0393322,-1.93539e-05,3.76619e-09,-2.64013e-13,-22420.1,-16.4643], Tmin=(756.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cd-Cd(CCO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(=O)C=COC=C=O(7324)',
    structure = SMILES('C=CC(=O)C=COC=C=O'),
    E0 = (-122.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2120,512.5,787.5,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.190454,0.0865968,-8.89853e-05,4.83705e-08,-1.07568e-11,-14635.5,33.3163], Tmin=(100,'K'), Tmax=(1073.11,'K')), NASAPolynomial(coeffs=[14.1574,0.0345356,-1.62146e-05,3.16228e-09,-2.24879e-13,-17633.1,-35.052], Tmin=(1073.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(=O)C(C=C=O)=CO(7325)',
    structure = SMILES('C=CC(=O)C(C=C=O)=CO'),
    E0 = (-252.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,375,552.5,462.5,1710,2120,512.5,787.5,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.745947,'amu*angstrom^2'), symmetry=1, barrier=(17.1508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744502,'amu*angstrom^2'), symmetry=1, barrier=(17.1176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.746988,'amu*angstrom^2'), symmetry=1, barrier=(17.1747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744603,'amu*angstrom^2'), symmetry=1, barrier=(17.1199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.470239,0.101656,-0.000125039,7.92404e-08,-1.99679e-11,-30154,33.0237], Tmin=(100,'K'), Tmax=(967.521,'K')), NASAPolynomial(coeffs=[16.8601,0.0300082,-1.39596e-05,2.70224e-09,-1.9116e-13,-33507.5,-50.0128], Tmin=(967.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-252.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-CdCs(CCO)) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
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
    E0 = (39.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (61.0496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (36.5761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (86.2685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (148.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (172.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (166.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (124.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (213.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (70.4663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-65.4388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (57.4734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (31.7289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (182.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (216.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (196.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (143.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (11.5188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (57.3693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-48.6097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (232.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (11.4698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (101.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-35.6334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-44.7519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-58.5275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (28.7933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (41.5992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (79.4229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (18.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (26.4342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (137.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (19.2184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-42.1365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (53.5297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (8.59042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (91.7095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (61.4994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (173.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (106.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (43.6556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (128.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (141.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (147.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (110.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (112.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (36.2611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (474.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (101.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (-175.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (-31.8328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (93.4263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (179.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (52.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (-123.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (-43.8853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (139.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (190.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (-108.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (-98.4333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction338',
    reactants = ['S(780)(779)', 'CH2CHCO(2811)'],
    products = ['S(560)(559)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(6.64e+07,'m^3/(mol*s)'), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_sec_rad;C_sec_rad] for rate rule [CO_rad/OneDe;C_rad/H/TwoDe]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction339',
    reactants = ['HCCO(47)(47)', 'C=CC([O])=CC=O(7240)'],
    products = ['S(560)(559)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.8959e+06,'m^3/(mol*s)'), n=-0.00978389, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cs_rad;Cd_allenic] + [C_rad/H/TwoDe;Y_rad] for rate rule [C_rad/H/TwoDe;Cd_allenic]
Euclidian distance = 2.0
family: R_Recombination
Ea raised from -6.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction340',
    reactants = ['HCO(14)(15)', 'C=CC(=O)[CH]C=C=O(7278)'],
    products = ['S(560)(559)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.64e+07,'m^3/(mol*s)'), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_sec_rad;CO_rad] for rate rule [C_rad/H/TwoDe;CO_pri_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction341',
    reactants = ['H(3)(3)', 'C=CC(=O)[C](C=O)C=C=O(7279)'],
    products = ['S(560)(559)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.25122e+07,'m^3/(mol*s)'), n=0.204, Ea=(0.213384,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_ter_rad;H_rad] for rate rule [C_rad/ThreeDe;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction342',
    reactants = ['C2H3(28)(29)', 'O=[C]C(C=O)C=C=O(7280)'],
    products = ['S(560)(559)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.6647e+07,'m^3/(mol*s)'), n=-0.0666667, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;CO_rad/NonDe] + [Cd_pri_rad;CO_rad] for rate rule [Cd_pri_rad;CO_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction343',
    reactants = ['H(3)(3)', 'C=CC(=O)C([C]=C=O)C=O(7281)'],
    products = ['S(560)(559)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction344',
    reactants = ['H(3)(3)', 'C=[C]C(=O)C(C=O)C=C=O(7282)'],
    products = ['S(560)(559)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction345',
    reactants = ['H(3)(3)', 'C=CC(=O)C([C]=O)C=C=O(7283)'],
    products = ['S(560)(559)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_rad/NonDe;Y_rad] for rate rule [CO_rad/NonDe;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction346',
    reactants = ['H(3)(3)', '[CH]=CC(=O)C(C=O)C=C=O(7284)'],
    products = ['S(560)(559)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction347',
    reactants = ['C=CC([O])[C](C=O)C=C=O(7285)'],
    products = ['S(560)(559)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction348',
    reactants = ['C=CC(=O)[C](C=O)C[C]=O(7286)'],
    products = ['S(560)(559)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction349',
    reactants = ['C=CC(=O)[C](C=C=O)C[O](7287)'],
    products = ['S(560)(559)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction350',
    reactants = ['C=CC(=O)C([C]=C[O])C=O(7288)'],
    products = ['S(560)(559)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction351',
    reactants = ['C=[C]C([O])C(C=O)C=C=O(7289)'],
    products = ['S(560)(559)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction352',
    reactants = ['C=CC([O])C([C]=C=O)C=O(7290)'],
    products = ['S(560)(559)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction353',
    reactants = ['C=CC(=O)C([C]=C=O)C[O](7291)'],
    products = ['S(560)(559)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction354',
    reactants = ['C=CC([O])C([C]=O)C=C=O(7292)'],
    products = ['S(560)(559)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction355',
    reactants = ['C=CC(=O)C([C]=O)C[C]=O(7293)'],
    products = ['S(560)(559)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction356',
    reactants = ['[CH2]CC(=O)[C](C=O)C=C=O(7294)'],
    products = ['S(560)(559)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction357',
    reactants = ['C=CC(=O)[C](C=O)C=C[O](7295)'],
    products = ['S(560)(559)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction358',
    reactants = ['[CH]=CC([O])C(C=O)C=C=O(7296)'],
    products = ['S(560)(559)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction359',
    reactants = ['C=C[C](O)C([C]=C=O)C=O(7297)'],
    products = ['S(560)(559)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_De] for rate rule [R4radEndo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction360',
    reactants = ['C=CC(=O)C([C]=C=O)[CH]O(7298)'],
    products = ['S(560)(559)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction361',
    reactants = ['C=C[C](O)C([C]=O)C=C=O(7299)'],
    products = ['S(560)(559)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radEndo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction362',
    reactants = ['C=CC(=O)C([C]=O)[CH]C=O(7300)'],
    products = ['S(560)(559)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction363',
    reactants = ['C[CH]C(=O)[C](C=O)C=C=O(7301)'],
    products = ['S(560)(559)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.34494e+09,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction364',
    reactants = ['C=CC(=O)[C](C=O)C=[C]O(7302)'],
    products = ['S(560)(559)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction365',
    reactants = ['[CH]=C[C](O)C(C=O)C=C=O(7303)'],
    products = ['S(560)(559)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction366',
    reactants = ['[CH2]CC(=O)C([C]=C=O)C=O(7304)'],
    products = ['S(560)(559)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction367',
    reactants = ['[CH2]CC(=O)C([C]=O)C=C=O(7305)'],
    products = ['S(560)(559)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction368',
    reactants = ['C=[C]C(=O)C(C=O)C[C]=O(7306)'],
    products = ['S(560)(559)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction369',
    reactants = ['C=[C]C(=O)C(C=C=O)C[O](7307)'],
    products = ['S(560)(559)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction370',
    reactants = ['C[CH]C(=O)C([C]=C=O)C=O(7308)'],
    products = ['S(560)(559)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction371',
    reactants = ['C[CH]C(=O)C([C]=O)C=C=O(7309)'],
    products = ['S(560)(559)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction372',
    reactants = ['C=CC(=O)C([C]=O)C=[C]O(7310)'],
    products = ['S(560)(559)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction373',
    reactants = ['C=[C]C(=O)C([CH]C=O)C=O(7311)'],
    products = ['S(560)(559)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction374',
    reactants = ['C=[C]C(=O)C([CH]O)C=C=O(7312)'],
    products = ['S(560)(559)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction375',
    reactants = ['[CH]=CC(=O)C(C=O)C[C]=O(7313)'],
    products = ['S(560)(559)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction376',
    reactants = ['[CH]=CC(=O)C(C=C=O)C[O](7314)'],
    products = ['S(560)(559)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction377',
    reactants = ['C=[C]C(=O)C(C=O)C=[C]O(7315)'],
    products = ['S(560)(559)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction378',
    reactants = ['[CH]=CC(=O)C([CH]C=O)C=O(7316)'],
    products = ['S(560)(559)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction379',
    reactants = ['[CH]=CC(=O)C([CH]O)C=C=O(7317)'],
    products = ['S(560)(559)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction380',
    reactants = ['[CH]=CC(=O)C(C=O)C=[C]O(7318)'],
    products = ['S(560)(559)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction389',
    reactants = ['C=CC(=O)C([C]C=O)C=O(7326)'],
    products = ['S(560)(559)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.61832e+16,'s^-1'), n=-0.885455, Ea=(87.4392,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [CsJ2-C;CsJ2(CsC);CH]
Euclidian distance = 0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction390',
    reactants = ['C[C]C(=O)C(C=O)C=C=O(7327)'],
    products = ['S(560)(559)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.84394e+15,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CsJ2-C;CsJ2C;CH] for rate rule [CsJ2-C;CsJ2C;CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction391',
    reactants = ['[CH]CC(=O)C(C=O)C=C=O(7328)'],
    products = ['S(560)(559)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;singletcarbene;CH2(C)] + [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH2(C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction392',
    reactants = ['C=C[C]1OC(=O)[CH]C1C=O(7329)'],
    products = ['S(560)(559)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction393',
    reactants = ['C=C[C]1OO[CH]C1C=C=O(7330)'],
    products = ['S(560)(559)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(6.8031e+16,'s^-1'), n=-1.26274, Ea=(207.863,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 0 used for R5JJ_Cd
Exact match found for rate rule [R5JJ_Cd]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction394',
    reactants = ['C=CC(=O)C1[CH]OC(=O)[CH]1(7331)'],
    products = ['S(560)(559)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction337',
    reactants = ['S(559)(558)'],
    products = ['S(560)(559)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction395',
    reactants = ['O=C=CC1[CH]OC[CH]C1=O(7332)'],
    products = ['S(560)(559)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction381',
    reactants = ['CO(10)(11)', 'C=CC(=O)CC=C=O(7210)'],
    products = ['S(560)(559)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.532e+06,'cm^3/(mol*s)'), n=2.07, Ea=(343.925,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO;C_sec] for rate rule [CO;C/H2/TwoDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction382',
    reactants = ['C=CC(=CC=O)OC=C=O(7319)'],
    products = ['S(560)(559)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction383',
    reactants = ['C=CC(=CC=C=O)OC=O(7320)'],
    products = ['S(560)(559)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction384',
    reactants = ['C=CC(O)=C(C=O)C=C=O(7321)'],
    products = ['S(560)(559)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction385',
    reactants = ['C=C=C(O)C(C=O)C=C=O(7322)'],
    products = ['S(560)(559)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction386',
    reactants = ['C=CC(=O)OC=CC=C=O(7323)'],
    products = ['S(560)(559)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction387',
    reactants = ['C=CC(=O)C=COC=C=O(7324)'],
    products = ['S(560)(559)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction388',
    reactants = ['C=CC(=O)C(C=C=O)=CO(7325)'],
    products = ['S(560)(559)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction270',
    reactants = ['S(560)(559)'],
    products = ['S(558)(557)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

network(
    label = '426',
    isomers = [
        'S(560)(559)',
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
    label = '426',
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

