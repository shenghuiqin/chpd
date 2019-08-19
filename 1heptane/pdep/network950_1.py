species(
    label = '[CH]=C(O)C([CH2])C(1364)',
    structure = SMILES('[CH]=C(O)C([CH2])C'),
    E0 = (189.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.454116,'amu*angstrom^2'), symmetry=1, barrier=(16.341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.463541,'amu*angstrom^2'), symmetry=1, barrier=(16.3221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00326786,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.44269,'amu*angstrom^2'), symmetry=1, barrier=(16.3611,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.939304,0.0564189,-2.38083e-05,-2.00181e-08,1.53985e-11,22893.4,25.3067], Tmin=(100,'K'), Tmax=(908.958,'K')), NASAPolynomial(coeffs=[17.2077,0.0143073,-2.96265e-06,3.73868e-10,-2.38553e-14,18718.2,-58.3249], Tmin=(908.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
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
    label = 'C3H6(72)',
    structure = SMILES('C=CC'),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,749.325,9.54025], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,-487.138,-4.54468], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=C(O)C(=C)C(7317)',
    structure = SMILES('[CH]=C(O)C(=C)C'),
    E0 = (89.5339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.909063,'amu*angstrom^2'), symmetry=1, barrier=(20.9012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.907868,'amu*angstrom^2'), symmetry=1, barrier=(20.8737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.907136,'amu*angstrom^2'), symmetry=1, barrier=(20.8569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.14479,0.0699057,-7.11153e-05,3.54022e-08,-6.65235e-12,10920,22.0848], Tmin=(100,'K'), Tmax=(1466.75,'K')), NASAPolynomial(coeffs=[19.0669,0.0102314,-1.83383e-06,1.60606e-10,-6.13673e-15,6237.39,-73.4919], Tmin=(1466.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.5339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
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
    label = '[CH]=C(O)C=C(6273)',
    structure = SMILES('[CH]=C(O)C=C'),
    E0 = (127.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.06241,'amu*angstrom^2'), symmetry=1, barrier=(24.427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06165,'amu*angstrom^2'), symmetry=1, barrier=(24.4094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55553,0.0432866,-1.7698e-05,-2.104e-08,1.52122e-11,15387.4,16.4164], Tmin=(100,'K'), Tmax=(907.962,'K')), NASAPolynomial(coeffs=[16.9959,0.00364253,9.14579e-07,-2.83671e-10,1.91369e-14,11413.8,-63.025], Tmin=(907.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
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
    label = 'C3H6(T)(143)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238389,'amu*angstrom^2'), symmetry=1, barrier=(5.48103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00909639,'amu*angstrom^2'), symmetry=1, barrier=(22.1005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93778,0.0190991,4.26842e-06,-1.44873e-08,5.74941e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.81,'K')), NASAPolynomial(coeffs=[5.93909,0.0171892,-6.69152e-06,1.21546e-09,-8.39795e-14,33151.2,-4.14888], Tmin=(1046.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C#CC([CH2])C(2946)',
    structure = SMILES('C#CC([CH2])C'),
    E0 = (321.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.45876,'amu*angstrom^2'), symmetry=1, barrier=(10.5478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0674332,'amu*angstrom^2'), symmetry=1, barrier=(83.1351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.60896,'amu*angstrom^2'), symmetry=1, barrier=(82.9771,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05266,0.0371709,-7.10649e-06,-1.96893e-08,1.19932e-11,38774.1,18.6599], Tmin=(100,'K'), Tmax=(877.4,'K')), NASAPolynomial(coeffs=[9.62985,0.0193968,-5.38942e-06,7.89676e-10,-4.88604e-14,36799,-20.583], Tmin=(877.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C(O)[C](C)C(7318)',
    structure = SMILES('[CH]C(O)=C(C)C'),
    E0 = (87.169,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,228.03,228.033,228.035],'cm^-1')),
        HinderedRotor(inertia=(1.36215,'amu*angstrom^2'), symmetry=1, barrier=(50.2629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36215,'amu*angstrom^2'), symmetry=1, barrier=(50.2629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36216,'amu*angstrom^2'), symmetry=1, barrier=(50.2629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36215,'amu*angstrom^2'), symmetry=1, barrier=(50.2629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770123,0.063083,-4.15552e-05,8.41729e-09,1.64833e-12,10607.2,22.4873], Tmin=(100,'K'), Tmax=(1021.34,'K')), NASAPolynomial(coeffs=[13.6079,0.02726,-1.01732e-05,1.79067e-09,-1.21629e-13,7230.91,-43.4097], Tmin=(1021.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(C)C(=C)[O](1151)',
    structure = SMILES('[CH2]C(C)C(=C)[O]'),
    E0 = (80.0546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,348.303,348.547],'cm^-1')),
        HinderedRotor(inertia=(0.00138873,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00138742,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127038,'amu*angstrom^2'), symmetry=1, barrier=(10.9498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3675.75,'J/mol'), sigma=(6.3037,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=574.14 K, Pc=33.3 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36686,0.0490894,-1.53711e-05,-1.91732e-08,1.26765e-11,9731.3,24.8689], Tmin=(100,'K'), Tmax=(922.496,'K')), NASAPolynomial(coeffs=[13.2401,0.0204788,-6.04096e-06,9.61536e-10,-6.39148e-14,6567.47,-36.7297], Tmin=(922.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.0546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2][C](C)C(=C)O(2587)',
    structure = SMILES('[CH2]C(C)=C([CH2])O'),
    E0 = (26.9004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.777219,0.0606937,-3.47937e-05,-7.51733e-09,1.02387e-11,3360.78,22.3236], Tmin=(100,'K'), Tmax=(929.747,'K')), NASAPolynomial(coeffs=[17.2507,0.0158046,-4.29313e-06,6.71495e-10,-4.57631e-14,-825.517,-61.9911], Tmin=(929.747,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.9004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C([CH2])C(=C)O(2588)',
    structure = SMILES('[CH2]C([CH2])C(=C)O'),
    E0 = (147.332,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,448.72],'cm^-1')),
        HinderedRotor(inertia=(0.00855371,'amu*angstrom^2'), symmetry=1, barrier=(1.3355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0803487,'amu*angstrom^2'), symmetry=1, barrier=(11.2014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.084487,'amu*angstrom^2'), symmetry=1, barrier=(11.2402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.579082,'amu*angstrom^2'), symmetry=1, barrier=(82.3034,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0328214,0.0669216,-6.20229e-05,2.92278e-08,-5.14499e-12,17883.9,29.4954], Tmin=(100,'K'), Tmax=(1658.31,'K')), NASAPolynomial(coeffs=[16.3309,0.0138188,-1.65886e-06,3.41639e-12,8.08117e-15,14331.1,-52.0751], Tmin=(1658.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C([O])C(C)C(2589)',
    structure = SMILES('[CH]=C([O])C(C)C'),
    E0 = (122.068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3120,650,792.5,1650,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,410.28],'cm^-1')),
        HinderedRotor(inertia=(0.0836663,'amu*angstrom^2'), symmetry=1, barrier=(9.99705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0836538,'amu*angstrom^2'), symmetry=1, barrier=(9.99689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0836214,'amu*angstrom^2'), symmetry=1, barrier=(9.99667,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25055,0.0512362,-2.07222e-05,-1.25271e-08,9.50099e-12,14788.7,23.2136], Tmin=(100,'K'), Tmax=(963.766,'K')), NASAPolynomial(coeffs=[13.7471,0.0209185,-7.07298e-06,1.22984e-09,-8.51941e-14,11379.3,-41.8056], Tmin=(963.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C(C)C(=C)O(2592)',
    structure = SMILES('C=C(C)C(=C)O'),
    E0 = (-157.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.788196,0.0597584,-3.16941e-05,-1.09827e-08,1.14352e-11,-18824.7,19.3668], Tmin=(100,'K'), Tmax=(935.057,'K')), NASAPolynomial(coeffs=[17.6031,0.0151963,-4.11324e-06,6.55661e-10,-4.56823e-14,-23165.8,-67.0241], Tmin=(935.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    label = '[CH]=C(O)C[CH2](1041)',
    structure = SMILES('[CH]=C(O)C[CH2]'),
    E0 = (221.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,3615,1277.5,1000,3000,3100,440,815,1455,1000,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(0.115587,'amu*angstrom^2'), symmetry=1, barrier=(13.8736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115587,'amu*angstrom^2'), symmetry=1, barrier=(13.8736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115587,'amu*angstrom^2'), symmetry=1, barrier=(13.8735,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3511,'J/mol'), sigma=(5.95716,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=548.41 K, Pc=37.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62102,0.0431667,-1.84841e-05,-1.48324e-08,1.1209e-11,26705.8,20.864], Tmin=(100,'K'), Tmax=(932.653,'K')), NASAPolynomial(coeffs=[14.7957,0.00990318,-2.36432e-06,3.63363e-10,-2.61429e-14,23237.5,-47.1966], Tmin=(932.653,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)[CH]CC(7319)',
    structure = SMILES('[CH]C(O)=CCC'),
    E0 = (103.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.982613,0.0547487,-1.27467e-05,-2.50605e-08,1.44077e-11,12576.7,23.8801], Tmin=(100,'K'), Tmax=(964.609,'K')), NASAPolynomial(coeffs=[14.8205,0.0247386,-8.64489e-06,1.52233e-09,-1.06136e-13,8633.61,-48.982], Tmin=(964.609,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(O)C[CH]C(1363)',
    structure = SMILES('[CH]=C(O)C[CH]C'),
    E0 = (186.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,392.826],'cm^-1')),
        HinderedRotor(inertia=(0.00109619,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0968749,'amu*angstrom^2'), symmetry=1, barrier=(10.5542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0970612,'amu*angstrom^2'), symmetry=1, barrier=(10.5617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0967014,'amu*angstrom^2'), symmetry=1, barrier=(10.5683,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05165,0.0584476,-4.50713e-05,1.47087e-08,-6.08598e-13,22564.4,26.1003], Tmin=(100,'K'), Tmax=(1008.7,'K')), NASAPolynomial(coeffs=[13.153,0.0219256,-7.81154e-06,1.35219e-09,-9.12556e-14,19539.7,-35.2781], Tmin=(1008.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(RCCJC)"""),
)

species(
    label = 'CC1CC=C1O(7320)',
    structure = SMILES('CC1CC=C1O'),
    E0 = (-101.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71639,0.0278113,6.41956e-05,-1.12469e-07,4.77959e-11,-12146,18.8604], Tmin=(100,'K'), Tmax=(933.01,'K')), NASAPolynomial(coeffs=[18.9134,0.0113533,-1.41627e-06,2.01058e-10,-2.19528e-14,-17847.7,-76.271], Tmin=(933.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
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
    label = '[CH]=C(O)[CH]C(6274)',
    structure = SMILES('[CH]C(O)=CC'),
    E0 = (126.224,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.14535,'amu*angstrom^2'), symmetry=1, barrier=(49.3258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14495,'amu*angstrom^2'), symmetry=1, barrier=(49.3165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14518,'amu*angstrom^2'), symmetry=1, barrier=(49.3218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65216,0.0426314,-9.45765e-06,-2.05987e-08,1.18491e-11,15273.9,18.5647], Tmin=(100,'K'), Tmax=(954.088,'K')), NASAPolynomial(coeffs=[12.4354,0.0190833,-6.48986e-06,1.12283e-09,-7.76725e-14,12230.4,-38.1179], Tmin=(954.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
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
    E0 = (189.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (312.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (287.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (373.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (387.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (350.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (290.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (381.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (381.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (233.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (273.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (633.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (252.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (641.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (384.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (349.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (197.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (507.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(O)C([CH2])C(1364)'],
    products = ['HCCOH(50)', 'C3H6(72)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', '[CH]=C(O)C(=C)C(7317)'],
    products = ['[CH]=C(O)C([CH2])C(1364)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(72.3521,'m^3/(mol*s)'), n=1.66655, Ea=(10.8198,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeCs_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH3(17)', '[CH]=C(O)C=C(6273)'],
    products = ['[CH]=C(O)C([CH2])C(1364)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.0143836,'m^3/(mol*s)'), n=2.40621, Ea=(24.5235,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=[C]O(172)', 'C3H6(72)'],
    products = ['[CH]=C(O)C([CH2])C(1364)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['HCCOH(50)', 'C3H6(T)(143)'],
    products = ['[CH]=C(O)C([CH2])C(1364)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;CJ] for rate rule [Ct-O_Ct;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OH(5)', 'C#CC([CH2])C(2946)'],
    products = ['[CH]=C(O)C([CH2])C(1364)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.148e+06,'cm^3/(mol*s)'), n=1.876, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 212 used for Ct-Cs_Ct-H;OJ_pri
Exact match found for rate rule [Ct-Cs_Ct-H;OJ_pri]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -1.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C(O)[C](C)C(7318)'],
    products = ['[CH]=C(O)C([CH2])C(1364)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.614e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(O)C([CH2])C(1364)'],
    products = ['[CH2]C(C)C(=C)[O](1151)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(O)C([CH2])C(1364)'],
    products = ['[CH2][C](C)C(=C)O(2587)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(O)C([CH2])C(1364)'],
    products = ['[CH2]C([CH2])C(=C)O(2588)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(O)C([CH2])C(1364)'],
    products = ['[CH]=C([O])C(C)C(2589)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]O(172)', 'C3H6(T)(143)'],
    products = ['[CH]=C(O)C([CH2])C(1364)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(O)C([CH2])C(1364)'],
    products = ['C=C(C)C(=C)O(2592)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CH2(S)(23)', '[CH]=C(O)C[CH2](1041)'],
    products = ['[CH]=C(O)C([CH2])C(1364)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(O)C([CH2])C(1364)'],
    products = ['[CH]=C(O)[CH]CC(7319)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(O)C([CH2])C(1364)'],
    products = ['[CH]=C(O)C[CH]C(1363)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(O)C([CH2])C(1364)'],
    products = ['CC1CC=C1O(7320)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH2(19)', '[CH]=C(O)[CH]C(6274)'],
    products = ['[CH]=C(O)C([CH2])C(1364)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '950',
    isomers = [
        '[CH]=C(O)C([CH2])C(1364)',
    ],
    reactants = [
        ('HCCOH(50)', 'C3H6(72)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '950',
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

