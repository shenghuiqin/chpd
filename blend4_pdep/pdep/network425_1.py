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
    label = '[CH2]C(=O)C=CC=O(6095)',
    structure = SMILES('C=C([O])C=CC=O'),
    E0 = (-105.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.804578,'amu*angstrom^2'), symmetry=1, barrier=(18.4988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805879,'amu*angstrom^2'), symmetry=1, barrier=(18.5287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16866,0.0608176,-5.94524e-05,2.9485e-08,-5.85639e-12,-12604.2,21.9805], Tmin=(100,'K'), Tmax=(1208.63,'K')), NASAPolynomial(coeffs=[13.3958,0.0203512,-9.23039e-06,1.78305e-09,-1.26336e-13,-15559.8,-39.3255], Tmin=(1208.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ)"""),
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
    label = 'O=[C]CC=C=O(7179)',
    structure = SMILES('O=[C]CC=C=O'),
    E0 = (-22.7318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2120,512.5,787.5,1855,455,950,328.891],'cm^-1')),
        HinderedRotor(inertia=(0.109925,'amu*angstrom^2'), symmetry=1, barrier=(8.4527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109874,'amu*angstrom^2'), symmetry=1, barrier=(8.43667,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35721,0.0370312,-3.36965e-05,1.68264e-08,-3.51704e-12,-2675.61,20.4424], Tmin=(100,'K'), Tmax=(1123.22,'K')), NASAPolynomial(coeffs=[7.78866,0.017689,-7.8661e-06,1.49527e-09,-1.04757e-13,-3895.76,-6.39236], Tmin=(1123.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.7318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C=CCC(=O)[C]=CC=O(7180)',
    structure = SMILES('[O]C(=C=CC=O)CC=C=O'),
    E0 = (-57.3971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,350,440,435,1725,2120,512.5,787.5,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.609369,'amu*angstrom^2'), symmetry=1, barrier=(14.0106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608809,'amu*angstrom^2'), symmetry=1, barrier=(13.9977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.610204,'amu*angstrom^2'), symmetry=1, barrier=(14.0298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.150723,0.100644,-0.000140997,1.02467e-07,-2.76085e-11,-6762.54,32.5151], Tmin=(100,'K'), Tmax=(638.925,'K')), NASAPolynomial(coeffs=[11.9445,0.036719,-1.86179e-05,3.67323e-09,-2.59815e-13,-8548.92,-22.3034], Tmin=(638.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.3971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsOs) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
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
    label = '[CH]=CC(=O)CC=C=O(7181)',
    structure = SMILES('[CH]=CC(=O)CC=C=O'),
    E0 = (115.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,375,552.5,462.5,1710,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0941364,'amu*angstrom^2'), symmetry=1, barrier=(2.16438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0892541,'amu*angstrom^2'), symmetry=1, barrier=(2.05213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.644702,'amu*angstrom^2'), symmetry=1, barrier=(14.823,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30799,0.0651965,-7.87847e-05,6.07152e-08,-2.05364e-11,14020.9,27.1527], Tmin=(100,'K'), Tmax=(702.892,'K')), NASAPolynomial(coeffs=[6.25483,0.037044,-1.87038e-05,3.7284e-09,-2.66934e-13,13325.5,5.03125], Tmin=(702.892,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'O=C=CCC(=O)C=[C]C=O(7182)',
    structure = SMILES('[O]C=C=CC(=O)CC=C=O'),
    E0 = (-58.0311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,264.227,264.308,264.393,264.874],'cm^-1')),
        HinderedRotor(inertia=(0.374027,'amu*angstrom^2'), symmetry=1, barrier=(18.4828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.372682,'amu*angstrom^2'), symmetry=1, barrier=(18.4837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.372605,'amu*angstrom^2'), symmetry=1, barrier=(18.4767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.325678,0.0836565,-8.58322e-05,4.57563e-08,-9.94612e-12,-6849.55,32.3997], Tmin=(100,'K'), Tmax=(1096.83,'K')), NASAPolynomial(coeffs=[14.4002,0.0323285,-1.56372e-05,3.0909e-09,-2.214e-13,-9937.01,-36.8024], Tmin=(1096.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.0311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'O=C=[C]CC(=O)C=CC=O(7183)',
    structure = SMILES('O=C=[C]CC(=O)C=CC=O'),
    E0 = (-52.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,2120,512.5,787.5,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,180,180,921.992],'cm^-1')),
        HinderedRotor(inertia=(0.274406,'amu*angstrom^2'), symmetry=1, barrier=(6.30913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274706,'amu*angstrom^2'), symmetry=1, barrier=(6.31602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.933178,'amu*angstrom^2'), symmetry=1, barrier=(21.4556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.932925,'amu*angstrom^2'), symmetry=1, barrier=(21.4498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.275082,0.0912236,-0.000134568,1.23696e-07,-4.77445e-11,-6187.52,33.224], Tmin=(100,'K'), Tmax=(720.126,'K')), NASAPolynomial(coeffs=[6.27546,0.0485389,-2.6171e-05,5.30651e-09,-3.81337e-13,-6809.16,7.92997], Tmin=(720.126,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + radical(CCCJ=C=O)"""),
)

species(
    label = 'O=[C]C=CC(=O)CC=C=O(7184)',
    structure = SMILES('O=C=C[CH]C(=O)CC=C=O'),
    E0 = (-95.0472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2110,2130,495,530,650,925,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,323.406,323.413,323.496,1347.4],'cm^-1')),
        HinderedRotor(inertia=(0.110275,'amu*angstrom^2'), symmetry=1, barrier=(8.18572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00161184,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.521841,'amu*angstrom^2'), symmetry=1, barrier=(38.7263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11029,'amu*angstrom^2'), symmetry=1, barrier=(8.18585,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.337623,0.0876572,-0.000114881,8.88788e-08,-2.88284e-11,-11306.2,34.8755], Tmin=(100,'K'), Tmax=(743.688,'K')), NASAPolynomial(coeffs=[9.06665,0.0407028,-2.01658e-05,3.96533e-09,-2.81015e-13,-12604.4,-4.65121], Tmin=(743.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.0472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-(Cdd-O2d)CsH) + radical(CCJC(C)=C=O)"""),
)

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
    label = 'O=[C]C[CH]C(=O)C=CC=O(7185)',
    structure = SMILES('[O]C(C=CC=O)=CC[C]=O'),
    E0 = (-80.9824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.647047,'amu*angstrom^2'), symmetry=1, barrier=(14.8769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647022,'amu*angstrom^2'), symmetry=1, barrier=(14.8763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.646703,'amu*angstrom^2'), symmetry=1, barrier=(14.869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.646864,'amu*angstrom^2'), symmetry=1, barrier=(14.8727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.187874,0.0866802,-8.35544e-05,4.00186e-08,-7.77844e-12,-9605.19,32.6273], Tmin=(100,'K'), Tmax=(1214.82,'K')), NASAPolynomial(coeffs=[16.5871,0.0326834,-1.68822e-05,3.43055e-09,-2.48956e-13,-13589.6,-49.6806], Tmin=(1214.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.9824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + radical(CCCJ=O) + radical(C=C(C)OJ)"""),
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
    label = '[O]C[C]=CC(=O)CC=C=O(7186)',
    structure = SMILES('[O]C[C]=CC(=O)CC=C=O'),
    E0 = (143.106,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,375,552.5,462.5,1710,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277147,0.096797,-0.000163197,1.63436e-07,-6.32674e-11,17331.2,37.324], Tmin=(100,'K'), Tmax=(821.933,'K')), NASAPolynomial(coeffs=[1.78527,0.0563823,-2.908e-05,5.69612e-09,-3.9756e-13,18200.5,37.1403], Tmin=(821.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=[C]CC(=O)C=CC=O(7187)',
    structure = SMILES('[O]C=[C]CC(=O)C=CC=O'),
    E0 = (5.29756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,359.955,359.955,359.955,359.955],'cm^-1')),
        HinderedRotor(inertia=(0.184564,'amu*angstrom^2'), symmetry=1, barrier=(16.9695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184564,'amu*angstrom^2'), symmetry=1, barrier=(16.9695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184564,'amu*angstrom^2'), symmetry=1, barrier=(16.9695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184564,'amu*angstrom^2'), symmetry=1, barrier=(16.9695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.899011,0.0774738,-6.16382e-05,2.29222e-08,-3.51256e-12,739.834,30.9464], Tmin=(100,'K'), Tmax=(1450.87,'K')), NASAPolynomial(coeffs=[15.062,0.0384269,-2.12693e-05,4.37303e-09,-3.16361e-13,-3369.93,-42.6531], Tmin=(1450.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.29756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'O=[C]C[CH]C(=O)CC=C=O(7188)',
    structure = SMILES('[O]C(=CC[C]=O)CC=C=O'),
    E0 = (-49.8152,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,2120,512.5,787.5,350,440,435,1725,278.091,279.579,279.963,280.799],'cm^-1')),
        HinderedRotor(inertia=(0.0021564,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200529,'amu*angstrom^2'), symmetry=1, barrier=(10.9553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19555,'amu*angstrom^2'), symmetry=1, barrier=(10.9674,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198235,'amu*angstrom^2'), symmetry=1, barrier=(10.9795,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.176703,0.0881747,-9.38127e-05,5.22263e-08,-1.18692e-11,-5857.05,34.5045], Tmin=(100,'K'), Tmax=(1051.31,'K')), NASAPolynomial(coeffs=[14.3854,0.0341142,-1.66802e-05,3.31468e-09,-2.38134e-13,-8844.62,-34.7553], Tmin=(1051.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.8152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(CCCJ=O)"""),
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
    label = 'O=C=C[CH]C(=O)C[CH]C=O(7189)',
    structure = SMILES('[O]C=CCC(=O)[CH]C=C=O'),
    E0 = (-72.8236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,375,552.5,462.5,1710,2120,512.5,787.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.152796,0.0836456,-7.14244e-05,2.91097e-08,-4.68948e-12,-8603.08,35.8436], Tmin=(100,'K'), Tmax=(1477.2,'K')), NASAPolynomial(coeffs=[21.364,0.0253823,-1.22621e-05,2.40968e-09,-1.70825e-13,-14960,-76.3573], Tmin=(1477.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-72.8236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(CCJC(C)=C=O) + radical(C=COJ)"""),
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
    label = '[O]CC=[C]C(=O)CC=C=O(7190)',
    structure = SMILES('[O]CC=C=C([O])CC=C=O'),
    E0 = (102.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,350,440,435,1725,195.54,195.546,195.55,195.553,3108.09],'cm^-1')),
        HinderedRotor(inertia=(0.360876,'amu*angstrom^2'), symmetry=1, barrier=(9.7933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.360899,'amu*angstrom^2'), symmetry=1, barrier=(9.79338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.360894,'amu*angstrom^2'), symmetry=1, barrier=(9.79321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.341449,0.105773,-0.000166753,1.46504e-07,-5.05263e-11,12490.5,36.5076], Tmin=(100,'K'), Tmax=(842.425,'K')), NASAPolynomial(coeffs=[9.0219,0.0423766,-2.01527e-05,3.80537e-09,-2.59844e-13,11584.9,-3.07134], Tmin=(842.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(CCOJ)"""),
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
    label = 'O=C=[C]C[C](O)C=CC=O(7191)',
    structure = SMILES('[O]C=CC=C(O)C[C]=C=O'),
    E0 = (-25.5699,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05852,'amu*angstrom^2'), symmetry=1, barrier=(24.3375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05948,'amu*angstrom^2'), symmetry=1, barrier=(24.3596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05917,'amu*angstrom^2'), symmetry=1, barrier=(24.3524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05896,'amu*angstrom^2'), symmetry=1, barrier=(24.3476,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.03824,0.111244,-0.000126036,6.66927e-08,-1.31383e-11,-2839.46,37.4315], Tmin=(100,'K'), Tmax=(1402.53,'K')), NASAPolynomial(coeffs=[29.9916,0.00641359,4.98582e-07,-3.06733e-10,2.60143e-14,-10498,-123.202], Tmin=(1402.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.5699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(CCCJ=C=O) + radical(C=COJ)"""),
)

species(
    label = 'O=C=C[CH]C(=O)[CH]CC=O(7192)',
    structure = SMILES('[O]C([CH]C=C=O)=CCC=O'),
    E0 = (-92.8592,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,350,440,435,1725,2120,512.5,787.5,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,364.36,364.36,364.36,364.36],'cm^-1')),
        HinderedRotor(inertia=(0.214141,'amu*angstrom^2'), symmetry=1, barrier=(20.1739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214142,'amu*angstrom^2'), symmetry=1, barrier=(20.1739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214142,'amu*angstrom^2'), symmetry=1, barrier=(20.1739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214142,'amu*angstrom^2'), symmetry=1, barrier=(20.1739,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0308801,0.0892708,-8.72542e-05,4.30646e-08,-8.62543e-12,-11027.2,31.4002], Tmin=(100,'K'), Tmax=(1184.66,'K')), NASAPolynomial(coeffs=[16.5491,0.0334983,-1.66373e-05,3.32569e-09,-2.3946e-13,-14940.9,-51.0899], Tmin=(1184.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-92.8592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = 'O=C=CC[C](O)C=[C]C=O(7193)',
    structure = SMILES('[O]C=[C]C=C(O)CC=C=O'),
    E0 = (-28.8552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04564,'amu*angstrom^2'), symmetry=1, barrier=(24.0413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04645,'amu*angstrom^2'), symmetry=1, barrier=(24.06,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04742,'amu*angstrom^2'), symmetry=1, barrier=(24.0823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04431,'amu*angstrom^2'), symmetry=1, barrier=(24.0107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78216,0.113091,-0.000135593,7.74361e-08,-1.66333e-11,-3250.15,36.3925], Tmin=(100,'K'), Tmax=(1250.71,'K')), NASAPolynomial(coeffs=[27.618,0.010009,-1.10505e-06,-3.87715e-11,9.96208e-15,-9896.15,-109.192], Tmin=(1250.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.8552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = 'O=C=CCC(=O)[C]=C[CH]O(7194)',
    structure = SMILES('[O]C(=[C]C=CO)CC=C=O'),
    E0 = (-32.513,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.935761,'amu*angstrom^2'), symmetry=1, barrier=(21.515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.935945,'amu*angstrom^2'), symmetry=1, barrier=(21.5192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.936211,'amu*angstrom^2'), symmetry=1, barrier=(21.5253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937766,'amu*angstrom^2'), symmetry=1, barrier=(21.5611,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33262,0.110109,-0.000124898,5.91778e-08,-6.95323e-12,-3711.48,34.8804], Tmin=(100,'K'), Tmax=(875.404,'K')), NASAPolynomial(coeffs=[26.076,0.0115129,-1.60775e-06,5.3118e-11,2.8598e-15,-9531.09,-99.5341], Tmin=(875.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'O=CC=CC(=O)[CH]C=[C]O(7195)',
    structure = SMILES('[O]C(C=CC=O)=CC=[C]O'),
    E0 = (-22.9314,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2782.5,750,1395,475,1775,1000,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.882573,'amu*angstrom^2'), symmetry=1, barrier=(20.2921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.878579,'amu*angstrom^2'), symmetry=1, barrier=(20.2003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.882401,'amu*angstrom^2'), symmetry=1, barrier=(20.2881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.878293,'amu*angstrom^2'), symmetry=1, barrier=(20.1937,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3609,0.107644,-0.000125731,7.03995e-08,-1.50569e-11,-2556.17,36.6807], Tmin=(100,'K'), Tmax=(1160.28,'K')), NASAPolynomial(coeffs=[25.8376,0.0138795,-4.51281e-06,7.51377e-10,-5.02249e-14,-8867.78,-98.5799], Tmin=(1160.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.9314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJO) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=C=[C]CC(=O)C[CH]C=O(7196)',
    structure = SMILES('[O]C=CCC(=O)C[C]=C=O'),
    E0 = (-17.8449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,375,552.5,462.5,1710,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.456965,0.0851462,-7.30407e-05,2.9744e-08,-4.73421e-12,-1975.12,35.9207], Tmin=(100,'K'), Tmax=(1513.07,'K')), NASAPolynomial(coeffs=[23.6586,0.0213937,-9.839e-06,1.89707e-09,-1.33147e-13,-9272.82,-90.41], Tmin=(1513.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.8449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(CCCJ=C=O) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]CCC(=O)[C]=CC=O(7197)',
    structure = SMILES('[O]C(=C=CC=O)CC[C]=O'),
    E0 = (-22.8039,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,350,440,435,1725,540,610,2055,1855,455,950,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.410635,'amu*angstrom^2'), symmetry=1, barrier=(9.4413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.410703,'amu*angstrom^2'), symmetry=1, barrier=(9.44288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.410848,'amu*angstrom^2'), symmetry=1, barrier=(9.4462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.41085,'amu*angstrom^2'), symmetry=1, barrier=(9.44626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.130569,0.100718,-0.000134589,9.22395e-08,-2.18155e-11,-2603.18,35.5417], Tmin=(100,'K'), Tmax=(620.313,'K')), NASAPolynomial(coeffs=[11.382,0.0398993,-1.99688e-05,3.92608e-09,-2.77448e-13,-4289.62,-16.5831], Tmin=(620.313,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.8039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(CCCJ=O) + radical(C=C(C)OJ)"""),
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
    label = 'O=C=[C]CC(=O)[CH]CC=O(7198)',
    structure = SMILES('[O]C(=CCC=O)C[C]=C=O'),
    E0 = (-7.49503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,350,440,435,1725,2120,512.5,787.5,1685,370,3010,987.5,1337.5,450,1655,314.975,315.207,315.363,315.694],'cm^-1')),
        HinderedRotor(inertia=(0.00167466,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227886,'amu*angstrom^2'), symmetry=1, barrier=(16.097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228939,'amu*angstrom^2'), symmetry=1, barrier=(16.1177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228324,'amu*angstrom^2'), symmetry=1, barrier=(16.1008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.273368,0.0858536,-8.4849e-05,4.35282e-08,-9.18283e-12,-770.491,33.044], Tmin=(100,'K'), Tmax=(1120.3,'K')), NASAPolynomial(coeffs=[14.3409,0.0356261,-1.75982e-05,3.50876e-09,-2.52338e-13,-3922.46,-36.4216], Tmin=(1120.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.49503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(CCCJ=C=O)"""),
)

species(
    label = 'O=C[CH]CC(=O)[C]=CC=O(7199)',
    structure = SMILES('[O]C=CCC([O])=C=CC=O'),
    E0 = (-38.2434,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.835144,'amu*angstrom^2'), symmetry=1, barrier=(19.2016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.836087,'amu*angstrom^2'), symmetry=1, barrier=(19.2233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.836112,'amu*angstrom^2'), symmetry=1, barrier=(19.2239,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.661376,0.0941895,-9.5637e-05,4.74274e-08,-9.1763e-12,-4424.79,36.3234], Tmin=(100,'K'), Tmax=(1261.47,'K')), NASAPolynomial(coeffs=[22.2628,0.0214987,-9.20085e-06,1.74712e-09,-1.23273e-13,-10208.4,-79.5973], Tmin=(1261.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-38.2434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]C=C[C](O)CC=C=O(7200)',
    structure = SMILES('O=[C][CH]C=C(O)CC=C=O'),
    E0 = (-25.4913,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1855,455,950,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.884629,0.0943237,-9.36048e-05,4.42253e-08,-8.06739e-12,-2879.03,37.1969], Tmin=(100,'K'), Tmax=(1343.93,'K')), NASAPolynomial(coeffs=[25.1133,0.0169458,-7.24219e-06,1.38499e-09,-9.82793e-14,-9867,-95.9134], Tmin=(1343.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.4913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CCC(=O)C=[C]C=O(7201)',
    structure = SMILES('[O]C=C=CC(=O)CC[C]=O'),
    E0 = (-25.2809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,375,552.5,462.5,1710,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,284.086,284.087,284.087,284.087],'cm^-1')),
        HinderedRotor(inertia=(0.199527,'amu*angstrom^2'), symmetry=1, barrier=(11.4269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199526,'amu*angstrom^2'), symmetry=1, barrier=(11.4269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199527,'amu*angstrom^2'), symmetry=1, barrier=(11.4269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199527,'amu*angstrom^2'), symmetry=1, barrier=(11.4269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0691943,0.099584,-0.00012879,8.30408e-08,-1.63709e-11,-2903.54,34.6267], Tmin=(100,'K'), Tmax=(608.413,'K')), NASAPolynomial(coeffs=[11.1294,0.0407409,-2.01606e-05,3.94294e-09,-2.77824e-13,-4539.8,-16.0838], Tmin=(608.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.2809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]CC=CC(=O)[CH]C=C=O(7202)',
    structure = SMILES('[O]CC=CC([O])=CC=C=O'),
    E0 = (40.6416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2120,512.5,787.5,350,440,435,1725,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.243801,'amu*angstrom^2'), symmetry=1, barrier=(5.60547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243036,'amu*angstrom^2'), symmetry=1, barrier=(5.58787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24336,'amu*angstrom^2'), symmetry=1, barrier=(5.59533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.229961,0.101138,-0.000150242,1.26697e-07,-4.26604e-11,5032.59,35.2939], Tmin=(100,'K'), Tmax=(828.707,'K')), NASAPolynomial(coeffs=[9.83945,0.0405446,-1.88632e-05,3.54829e-09,-2.42685e-13,3775.39,-8.9092], Tmin=(828.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.6416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + radical(CCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=CC=[C]C(=O)CC=[C]O(7203)',
    structure = SMILES('[O]C(=C=CC=O)CC=[C]O'),
    E0 = (60.0382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.687773,'amu*angstrom^2'), symmetry=1, barrier=(15.8133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.689749,'amu*angstrom^2'), symmetry=1, barrier=(15.8587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.688782,'amu*angstrom^2'), symmetry=1, barrier=(15.8365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.689196,'amu*angstrom^2'), symmetry=1, barrier=(15.846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.59487,0.0997639,-0.000115017,6.61698e-08,-1.49413e-11,7387.62,37.8042], Tmin=(100,'K'), Tmax=(1082.76,'K')), NASAPolynomial(coeffs=[19.7379,0.0246504,-1.09598e-05,2.10187e-09,-1.48774e-13,2984.45,-61.9067], Tmin=(1082.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.0382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=CJO)"""),
)

species(
    label = 'O=C[C]=CC(=O)C[CH]C=O(7204)',
    structure = SMILES('[O]C=C=CC(=O)CC=C[O]'),
    E0 = (-35.8075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,540,610,2055,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.12699,'amu*angstrom^2'), symmetry=1, barrier=(25.9117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12441,'amu*angstrom^2'), symmetry=1, barrier=(25.8524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12774,'amu*angstrom^2'), symmetry=1, barrier=(25.9289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.547211,0.0841719,-5.85936e-05,8.25312e-09,3.8547e-12,-4129.41,34.7418], Tmin=(100,'K'), Tmax=(1086.57,'K')), NASAPolynomial(coeffs=[24.1608,0.0210957,-1.00081e-05,2.05947e-09,-1.53933e-13,-11144.7,-94.0852], Tmin=(1086.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.8075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'O=C=C[CH]C(=O)C=C[CH]O(7205)',
    structure = SMILES('[O]C([CH]C=C=O)=CC=CO'),
    E0 = (-114.592,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2120,512.5,787.5,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16713,'amu*angstrom^2'), symmetry=1, barrier=(26.8346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17114,'amu*angstrom^2'), symmetry=1, barrier=(26.9268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17029,'amu*angstrom^2'), symmetry=1, barrier=(26.9072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16816,'amu*angstrom^2'), symmetry=1, barrier=(26.8582,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.26805,0.117005,-0.000137298,7.4874e-08,-1.516e-11,-13538.6,35.832], Tmin=(100,'K'), Tmax=(1367.65,'K')), NASAPolynomial(coeffs=[30.993,0.0051624,1.34093e-06,-4.92367e-10,3.98486e-14,-21274.5,-130.068], Tmin=(1367.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(C=CCJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]CC=CC(=O)C[C]=C=O(7206)',
    structure = SMILES('[O]CC=CC(=O)C[C]=C=O'),
    E0 = (107.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,375,552.5,462.5,1710,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.313301,0.0930251,-0.000144928,1.40192e-07,-5.40796e-11,13055.8,36.4303], Tmin=(100,'K'), Tmax=(800.202,'K')), NASAPolynomial(coeffs=[3.29989,0.0543149,-2.77867e-05,5.46003e-09,-3.83286e-13,13339.2,27.4449], Tmin=(800.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + radical(CCCJ=C=O) + radical(CCOJ)"""),
)

species(
    label = 'O=[C]C=CC(=O)CC[C]=O(7207)',
    structure = SMILES('O=[C]CCC(=O)[CH]C=C=O'),
    E0 = (-62.297,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,2120,512.5,787.5,3025,407.5,1350,352.5,1855,455,950,3010,987.5,1337.5,450,1655,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.313262,0.107688,-0.000178355,1.64562e-07,-5.90365e-11,-7349.65,37.9575], Tmin=(100,'K'), Tmax=(837.548,'K')), NASAPolynomial(coeffs=[7.04391,0.0468307,-2.32999e-05,4.47624e-09,-3.08307e-13,-7679.93,9.15315], Tmin=(837.548,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJC(C)=C=O)"""),
)

species(
    label = 'O=C[C]=CC(=O)CC=[C]O(7208)',
    structure = SMILES('[O]C=C=CC(=O)CC=[C]O'),
    E0 = (62.4741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,540,610,2055,3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.931033,'amu*angstrom^2'), symmetry=1, barrier=(21.4063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.932939,'amu*angstrom^2'), symmetry=1, barrier=(21.4501,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.932197,'amu*angstrom^2'), symmetry=1, barrier=(21.433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.931508,'amu*angstrom^2'), symmetry=1, barrier=(21.4172,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.774497,0.0931107,-8.91415e-05,4.04102e-08,-7.10832e-12,7695.67,37.2815], Tmin=(100,'K'), Tmax=(1384.16,'K')), NASAPolynomial(coeffs=[24.8585,0.0190357,-8.86722e-06,1.74693e-09,-1.25163e-13,599.645,-94.7158], Tmin=(1384.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.4741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'O=C=[C]CC(=O)C=C[CH]O(7209)',
    structure = SMILES('[O]C(=CC=CO)C[C]=C=O'),
    E0 = (-29.2277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.962127,'amu*angstrom^2'), symmetry=1, barrier=(22.1212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.960572,'amu*angstrom^2'), symmetry=1, barrier=(22.0854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.959271,'amu*angstrom^2'), symmetry=1, barrier=(22.0555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.962759,'amu*angstrom^2'), symmetry=1, barrier=(22.1357,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.97803,0.113113,-0.000133569,7.39849e-08,-1.5258e-11,-3284.23,37.2988], Tmin=(100,'K'), Tmax=(1337.3,'K')), NASAPolynomial(coeffs=[29.1887,0.00662339,7.58199e-07,-3.97748e-10,3.42616e-14,-10433.8,-117.686], Tmin=(1337.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.2277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(CCCJ=C=O) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=C=CCC(=O)[C]CC=O(7214)',
    structure = SMILES('O=C=CCC(=O)[C]CC=O'),
    E0 = (62.6849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.921657,0.0804236,-5.76346e-05,-3.88584e-08,6.79862e-11,7637.17,43.2131], Tmin=(100,'K'), Tmax=(488.17,'K')), NASAPolynomial(coeffs=[7.50242,0.0479293,-2.36297e-05,4.61766e-09,-3.25431e-13,6739.34,13.5687], Tmin=(488.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.6849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsCsHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'O=C=CCC(=O)C[C]C=O(7215)',
    structure = SMILES('O=C=CCC(=O)C[C]C=O'),
    E0 = (62.6849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.921657,0.0804236,-5.76346e-05,-3.88584e-08,6.79862e-11,7637.17,43.2131], Tmin=(100,'K'), Tmax=(488.17,'K')), NASAPolynomial(coeffs=[7.50242,0.0479293,-2.36297e-05,4.61766e-09,-3.25431e-13,6739.34,13.5687], Tmin=(488.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.6849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsCsHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'O=C[C]CC(=O)C=CC=O(7216)',
    structure = SMILES('O=C[C]CC(=O)C=CC=O'),
    E0 = (50.2663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0779576,0.0984421,-0.000153459,1.49985e-07,-5.94431e-11,6174.99,45.216], Tmin=(100,'K'), Tmax=(767.476,'K')), NASAPolynomial(coeffs=[3.65807,0.0577764,-3.09684e-05,6.22228e-09,-4.43254e-13,6273.57,33.1139], Tmin=(767.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.2663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'O=CC=C[C]1C[CH]C(=O)O1(7217)',
    structure = SMILES('[O]C=CC=C1C[CH]C(=O)O1'),
    E0 = (-131.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.498976,0.0500049,4.67367e-05,-1.10592e-07,4.96013e-11,-15678,32.6164], Tmin=(100,'K'), Tmax=(943.75,'K')), NASAPolynomial(coeffs=[25.3239,0.0128191,-2.29127e-06,4.25191e-10,-4.15348e-14,-23393.4,-101.764], Tmin=(943.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(CCJCO) + radical(C=COJ)"""),
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
    label = 'C=C(C=CC=O)OC=C=O(7211)',
    structure = SMILES('C=C(C=CC=O)OC=C=O'),
    E0 = (-134.319,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.805487,'amu*angstrom^2'), symmetry=1, barrier=(18.5197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805542,'amu*angstrom^2'), symmetry=1, barrier=(18.521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805661,'amu*angstrom^2'), symmetry=1, barrier=(18.5237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805476,'amu*angstrom^2'), symmetry=1, barrier=(18.5195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.148678,0.0948826,-0.000106066,6.14051e-08,-1.43487e-11,-16008.4,31.7042], Tmin=(100,'K'), Tmax=(1031.01,'K')), NASAPolynomial(coeffs=[15.8089,0.0329732,-1.59967e-05,3.16568e-09,-2.26957e-13,-19299,-45.7693], Tmin=(1031.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-(Cdd-O2d)OsH) + group(Cds-O2d(Cds-Cds)H)"""),
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
    label = 'O=C=CCC(O)=C=CC=O(7212)',
    structure = SMILES('O=C=CCC(O)=C=CC=O'),
    E0 = (-195.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.765941,'amu*angstrom^2'), symmetry=1, barrier=(17.6105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.76593,'amu*angstrom^2'), symmetry=1, barrier=(17.6102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.762809,'amu*angstrom^2'), symmetry=1, barrier=(17.5385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.763316,'amu*angstrom^2'), symmetry=1, barrier=(17.5501,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.447547,0.104084,-0.000133452,8.933e-08,-2.40067e-11,-23322.6,32.0647], Tmin=(100,'K'), Tmax=(904.619,'K')), NASAPolynomial(coeffs=[15.4686,0.0337054,-1.67512e-05,3.32439e-09,-2.3782e-13,-26202.1,-43.1259], Tmin=(904.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsOs) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
)

species(
    label = 'O=C=CCC(=O)C=C=CO(7213)',
    structure = SMILES('O=C=CCC(=O)C=C=CO'),
    E0 = (-199.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,2120,512.5,787.5,540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.820018,'amu*angstrom^2'), symmetry=1, barrier=(18.8538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.8201,'amu*angstrom^2'), symmetry=1, barrier=(18.8557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.819645,'amu*angstrom^2'), symmetry=1, barrier=(18.8453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.819645,'amu*angstrom^2'), symmetry=1, barrier=(18.8453,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.335487,0.0936421,-0.000100097,5.40493e-08,-1.15719e-11,-23835.9,33.3572], Tmin=(100,'K'), Tmax=(1133.28,'K')), NASAPolynomial(coeffs=[18.4784,0.0272369,-1.22031e-05,2.34459e-09,-1.65912e-13,-28100.2,-59.7627], Tmin=(1133.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-199.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

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
    E0 = (62.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (92.4019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (149.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (158.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (148.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (157.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (159.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (116.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (5.60925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-58.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (173.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (165.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (27.5985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-21.9624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (203.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (20.2704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (198.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (191.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-32.2444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-7.78788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-53.6342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (31.8128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (6.71195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (16.2936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-0.0628953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (16.4211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (69.3768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (31.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (0.981629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (9.02667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (13.9441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (79.8666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (99.2632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (3.41749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-75.3668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (146.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-37.3237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (101.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (9.99731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (99.7229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (99.7229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (137.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (18.7043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (-175.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (54.8076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (179.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (-117.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (-52.0883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (-55.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (-98.4333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction223',
    reactants = ['S(780)(779)', 'CH2CHCO(2811)'],
    products = ['S(558)(557)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.12e+08,'m^3/(mol*s)'), n=-0.5, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_sec_rad;C_pri_rad] for rate rule [CO_rad/OneDe;C_rad/H2/Cd]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction224',
    reactants = ['HCCO(47)(47)', '[CH2]C(=O)C=CC=O(6095)'],
    products = ['S(558)(557)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;Cd_allenic] for rate rule [C_rad/H2/CO;Cd_allenic]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction225',
    reactants = ['H(3)(3)', 'O=C=C[CH]C(=O)C=CC=O(7144)'],
    products = ['S(558)(557)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(78817,'m^3/(mol*s)'), n=0.288419, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H/TwoDe;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction226',
    reactants = ['CHCHCHO(1679)', 'O=[C]CC=C=O(7179)'],
    products = ['S(558)(557)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.6647e+07,'m^3/(mol*s)'), n=-0.0666667, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;CO_rad/NonDe] + [Cd_pri_rad;CO_rad] for rate rule [Cd_pri_rad;CO_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction227',
    reactants = ['H(3)(3)', 'O=C=CCC(=O)[C]=CC=O(7180)'],
    products = ['S(558)(557)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction228',
    reactants = ['HCO(14)(15)', '[CH]=CC(=O)CC=C=O(7181)'],
    products = ['S(558)(557)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.81e+13,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 91 used for Cd_pri_rad;CO_pri_rad
Exact match found for rate rule [CO_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction229',
    reactants = ['H(3)(3)', 'O=C=CCC(=O)C=[C]C=O(7182)'],
    products = ['S(558)(557)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction230',
    reactants = ['H(3)(3)', 'O=C=[C]CC(=O)C=CC=O(7183)'],
    products = ['S(558)(557)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction231',
    reactants = ['H(3)(3)', 'O=[C]C=CC(=O)CC=C=O(7184)'],
    products = ['S(558)(557)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;CO_sec_rad] for rate rule [H_rad;CO_rad/OneDe]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction202',
    reactants = ['S(557)(556)'],
    products = ['S(558)(557)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction232',
    reactants = ['O=[C]C[CH]C(=O)C=CC=O(7185)'],
    products = ['S(558)(557)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction233',
    reactants = ['[O]C([C]=CC=O)CC=C=O(7151)'],
    products = ['S(558)(557)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction234',
    reactants = ['[O]C[C]=CC(=O)CC=C=O(7186)'],
    products = ['S(558)(557)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction235',
    reactants = ['[O]C=[C]CC(=O)C=CC=O(7187)'],
    products = ['S(558)(557)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction236',
    reactants = ['O=[C]C[CH]C(=O)CC=C=O(7188)'],
    products = ['S(558)(557)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction237',
    reactants = ['[O]C(C=CC=O)C[C]=C=O(7149)'],
    products = ['S(558)(557)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=C=C[CH]C(=O)C[CH]C=O(7189)'],
    products = ['S(558)(557)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction239',
    reactants = ['[O]C(C=[C]C=O)CC=C=O(7153)'],
    products = ['S(558)(557)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction240',
    reactants = ['[O]CC=[C]C(=O)CC=C=O(7190)'],
    products = ['S(558)(557)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction241',
    reactants = ['[O][C](C=CC=O)C=CC=O(7174)'],
    products = ['S(558)(557)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction242',
    reactants = ['O=C=[C]C[C](O)C=CC=O(7191)'],
    products = ['S(558)(557)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_De] for rate rule [R4radEndo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction243',
    reactants = ['O=C=C[CH]C(=O)[CH]CC=O(7192)'],
    products = ['S(558)(557)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction244',
    reactants = ['O=C=CC[C](O)C=[C]C=O(7193)'],
    products = ['S(558)(557)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(60.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction245',
    reactants = ['O=C=CCC(=O)[C]=C[CH]O(7194)'],
    products = ['S(558)(557)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction246',
    reactants = ['O=CC=CC(=O)[CH]C=[C]O(7195)'],
    products = ['S(558)(557)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction247',
    reactants = ['O=C=[C]CC(=O)C[CH]C=O(7196)'],
    products = ['S(558)(557)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_De] for rate rule [R4radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction248',
    reactants = ['O=[C]CCC(=O)[C]=CC=O(7197)'],
    products = ['S(558)(557)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction249',
    reactants = ['[O]C(C=C[C]=O)CC=C=O(7156)'],
    products = ['S(558)(557)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction250',
    reactants = ['O=C=[C]CC(=O)[CH]CC=O(7198)'],
    products = ['S(558)(557)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction251',
    reactants = ['O=C[CH]CC(=O)[C]=CC=O(7199)'],
    products = ['S(558)(557)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction252',
    reactants = ['O=[C]C=C[C](O)CC=C=O(7200)'],
    products = ['S(558)(557)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction253',
    reactants = ['O=[C]CCC(=O)C=[C]C=O(7201)'],
    products = ['S(558)(557)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction254',
    reactants = ['[O]CC=CC(=O)[CH]C=C=O(7202)'],
    products = ['S(558)(557)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction255',
    reactants = ['O=CC=[C]C(=O)CC=[C]O(7203)'],
    products = ['S(558)(557)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction256',
    reactants = ['O=C[C]=CC(=O)C[CH]C=O(7204)'],
    products = ['S(558)(557)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction257',
    reactants = ['O=C=C[CH]C(=O)C=C[CH]O(7205)'],
    products = ['S(558)(557)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction258',
    reactants = ['[O]CC=CC(=O)C[C]=C=O(7206)'],
    products = ['S(558)(557)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction259',
    reactants = ['O=[C]C=CC(=O)CC[C]=O(7207)'],
    products = ['S(558)(557)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction260',
    reactants = ['O=C[C]=CC(=O)CC=[C]O(7208)'],
    products = ['S(558)(557)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad_De;XH_Rrad] for rate rule [R7radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction261',
    reactants = ['O=C=[C]CC(=O)C=C[CH]O(7209)'],
    products = ['S(558)(557)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad_De;XH_Rrad] for rate rule [R7radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction267',
    reactants = ['O=C=CCC(=O)[C]CC=O(7214)'],
    products = ['S(558)(557)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2C;CH2(C)] + [CsJ2-C;CsJ2C;CH] for rate rule [CsJ2-C;CsJ2C;CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction268',
    reactants = ['O=C=CCC(=O)C[C]C=O(7215)'],
    products = ['S(558)(557)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2C;CH2(C)] + [CsJ2-C;CsJ2C;CH] for rate rule [CsJ2-C;CsJ2C;CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction269',
    reactants = ['O=C[C]CC(=O)C=CC=O(7216)'],
    products = ['S(558)(557)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.61832e+16,'s^-1'), n=-0.885455, Ea=(87.4392,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [CsJ2-C;CsJ2(CsC);CH]
Euclidian distance = 0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction271',
    reactants = ['O=CC=C[C]1C[CH]C(=O)O1(7217)'],
    products = ['S(558)(557)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction272',
    reactants = ['S(559)(558)'],
    products = ['S(558)(557)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction262',
    reactants = ['CO(10)(11)', 'C=CC(=O)CC=C=O(7210)'],
    products = ['S(558)(557)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction263',
    reactants = ['C=C(C=CC=O)OC=C=O(7211)'],
    products = ['S(558)(557)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction264',
    reactants = ['O=C=CC=C(O)C=CC=O(7162)'],
    products = ['S(558)(557)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction265',
    reactants = ['O=C=CCC(O)=C=CC=O(7212)'],
    products = ['S(558)(557)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction266',
    reactants = ['O=C=CCC(=O)C=C=CO(7213)'],
    products = ['S(558)(557)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction270',
    reactants = ['S(560)(559)'],
    products = ['S(558)(557)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

network(
    label = '425',
    isomers = [
        'S(558)(557)',
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
    label = '425',
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

