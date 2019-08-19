species(
    label = '[CH2]C(C)C=[C]C=C(26550)',
    structure = SMILES('[CH2]C(C)C=[C]C=C'),
    E0 = (407.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,302.368,302.374],'cm^-1')),
        HinderedRotor(inertia=(0.183703,'amu*angstrom^2'), symmetry=1, barrier=(11.9185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183705,'amu*angstrom^2'), symmetry=1, barrier=(11.9186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.338998,'amu*angstrom^2'), symmetry=1, barrier=(21.9946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1268,'amu*angstrom^2'), symmetry=1, barrier=(73.1071,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621362,0.0654423,-3.25224e-05,-6.9307e-09,8.77755e-12,49149.8,28.2179], Tmin=(100,'K'), Tmax=(938.04,'K')), NASAPolynomial(coeffs=[14.2172,0.0292934,-9.62003e-06,1.59778e-09,-1.06331e-13,45638.8,-41.623], Tmin=(938.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Isobutyl)"""),
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
    label = 'CH2CHCCH(26391)',
    structure = SMILES('C#CC=C'),
    E0 = (274.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,3010,987.5,1337.5,450,1655,2175,525],'cm^-1')),
        HinderedRotor(inertia=(1.46338,'amu*angstrom^2'), symmetry=1, barrier=(33.6459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.87083,0.0182042,1.06711e-05,-2.72492e-08,1.19478e-11,33023.8,11.2934], Tmin=(100,'K'), Tmax=(955.249,'K')), NASAPolynomial(coeffs=[8.52653,0.0108962,-3.56564e-06,6.31243e-10,-4.51891e-14,31196.2,-19.6435], Tmin=(955.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""CH2CHCCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C1[C]=CC(C)C1(26745)',
    structure = SMILES('[CH2]C1[C]=CC(C)C1'),
    E0 = (421.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5922,0.0344533,5.45872e-05,-9.39606e-08,3.82981e-11,50747.9,24.1745], Tmin=(100,'K'), Tmax=(953.512,'K')), NASAPolynomial(coeffs=[13.9578,0.0288541,-9.40051e-06,1.67448e-09,-1.20862e-13,46286.1,-45.9245], Tmin=(953.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(cyclopentene-vinyl) + radical(Isobutyl)"""),
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
    label = 'C=C[C]=CC(=C)C(26746)',
    structure = SMILES('C=C[C]=CC(=C)C'),
    E0 = (307.098,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.09163,'amu*angstrom^2'), symmetry=1, barrier=(25.0988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09186,'amu*angstrom^2'), symmetry=1, barrier=(25.1039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0922,'amu*angstrom^2'), symmetry=1, barrier=(25.1119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.462,0.0680021,-4.03173e-05,-1.5034e-09,7.33656e-12,37071.6,23.5262], Tmin=(100,'K'), Tmax=(952.865,'K')), NASAPolynomial(coeffs=[16.3586,0.0243804,-8.02797e-06,1.35854e-09,-9.23327e-14,32993,-57.9035], Tmin=(952.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C(C)C#CC=C(26747)',
    structure = SMILES('[CH2]C(C)C#CC=C'),
    E0 = (378.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,339.697,340.037],'cm^-1')),
        HinderedRotor(inertia=(0.00146258,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265808,'amu*angstrom^2'), symmetry=1, barrier=(21.6855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925765,'amu*angstrom^2'), symmetry=1, barrier=(75.888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.931403,'amu*angstrom^2'), symmetry=1, barrier=(75.8743,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15095,0.0510854,-8.2234e-07,-3.69592e-08,1.88964e-11,45601.4,27.4208], Tmin=(100,'K'), Tmax=(938.01,'K')), NASAPolynomial(coeffs=[13.5755,0.0270247,-8.59594e-06,1.43656e-09,-9.7738e-14,41998.2,-38.5075], Tmin=(938.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C=C=C=C(26748)',
    structure = SMILES('[CH2]C(C)C=C=C=C'),
    E0 = (437.966,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.622965,'amu*angstrom^2'), symmetry=1, barrier=(14.3232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.622865,'amu*angstrom^2'), symmetry=1, barrier=(14.3209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4212,'amu*angstrom^2'), symmetry=1, barrier=(32.6761,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.574934,0.0671709,-5.23442e-05,2.19204e-08,-3.74055e-12,52805.2,27.8779], Tmin=(100,'K'), Tmax=(1390.52,'K')), NASAPolynomial(coeffs=[14.1286,0.028182,-1.02855e-05,1.75586e-09,-1.1519e-13,49035.9,-41.9791], Tmin=(1390.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Isobutyl)"""),
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
    label = 'C=C[C]=CC=C(26601)',
    structure = SMILES('C=C[C]=CC=C'),
    E0 = (344.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.32866,'amu*angstrom^2'), symmetry=1, barrier=(30.5484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33122,'amu*angstrom^2'), symmetry=1, barrier=(30.6074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31057,0.0479268,-9.42475e-06,-2.91657e-08,1.70739e-11,41563.7,19.8807], Tmin=(100,'K'), Tmax=(926.369,'K')), NASAPolynomial(coeffs=[15.4385,0.0157821,-4.10458e-06,6.34548e-10,-4.37525e-14,37707.8,-53.8816], Tmin=(926.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=[C]C=C(4699)',
    structure = SMILES('[CH]=C=C[CH2]'),
    E0 = (451.584,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3000,3100,440,815,1455,1000,180,1024.85,1025.53,1026.61],'cm^-1')),
        HinderedRotor(inertia=(0.00938781,'amu*angstrom^2'), symmetry=1, barrier=(7.01846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76805,0.020302,8.75519e-06,-2.87666e-08,1.37354e-11,54363.7,13.5565], Tmin=(100,'K'), Tmax=(915.031,'K')), NASAPolynomial(coeffs=[9.46747,0.00887314,-1.78262e-06,2.38534e-10,-1.6263e-14,52390.1,-22.2544], Tmin=(915.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C=CJ) + radical(Allyl_P)"""),
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
    label = 'C=C[C]=C[C](C)C(26749)',
    structure = SMILES('[CH2]C=[C]C=C(C)C'),
    E0 = (299.866,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.536166,'amu*angstrom^2'), symmetry=1, barrier=(12.3275,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.491938,0.0686237,-4.97023e-05,1.96852e-08,-3.22309e-12,36198.9,26.1167], Tmin=(100,'K'), Tmax=(1431.13,'K')), NASAPolynomial(coeffs=[13.4767,0.0323305,-1.16618e-05,1.96434e-09,-1.2741e-13,32482.4,-41.1817], Tmin=(1431.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C(C)[C]=CC=C(26750)',
    structure = SMILES('[CH2]C(C)[C]=CC=C'),
    E0 = (446.422,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,251.008,253.479],'cm^-1')),
        HinderedRotor(inertia=(0.00269849,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361289,'amu*angstrom^2'), symmetry=1, barrier=(16.2261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36046,'amu*angstrom^2'), symmetry=1, barrier=(16.2141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.366915,'amu*angstrom^2'), symmetry=1, barrier=(16.2108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.653782,0.0636565,-2.88187e-05,-9.6763e-09,9.17109e-12,53821.6,28.9779], Tmin=(100,'K'), Tmax=(965.146,'K')), NASAPolynomial(coeffs=[14.6041,0.0286432,-9.84212e-06,1.69564e-09,-1.15534e-13,50066.7,-43.3315], Tmin=(965.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C=C[C]=C(26751)',
    structure = SMILES('[CH2]C(C)C=C[C]=C'),
    E0 = (407.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,301.719,303.233],'cm^-1')),
        HinderedRotor(inertia=(0.338666,'amu*angstrom^2'), symmetry=1, barrier=(21.9993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182812,'amu*angstrom^2'), symmetry=1, barrier=(11.9139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184336,'amu*angstrom^2'), symmetry=1, barrier=(11.9216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13211,'amu*angstrom^2'), symmetry=1, barrier=(73.1247,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621371,0.0654421,-3.2522e-05,-6.93131e-09,8.77783e-12,49149.8,28.2179], Tmin=(100,'K'), Tmax=(938.037,'K')), NASAPolynomial(coeffs=[14.2172,0.0292934,-9.62006e-06,1.59779e-09,-1.06331e-13,45638.8,-41.6229], Tmin=(938.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=[C]C(C)C(26752)',
    structure = SMILES('[CH2][CH]C#CC(C)C'),
    E0 = (403.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0857,0.0567873,-1.70907e-05,-1.86449e-08,1.26297e-11,48636.7,28.5319], Tmin=(100,'K'), Tmax=(897.209,'K')), NASAPolynomial(coeffs=[11.2114,0.0321265,-1.01048e-05,1.60863e-09,-1.03667e-13,45995.4,-23.8148], Tmin=(897.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(RCCJ) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH2][C](C)C=CC=C(26753)',
    structure = SMILES('[CH2]C=CC=C([CH2])C'),
    E0 = (218.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06944,0.0507094,1.06835e-05,-4.86626e-08,2.22209e-11,26448.6,25.9979], Tmin=(100,'K'), Tmax=(964.242,'K')), NASAPolynomial(coeffs=[13.6838,0.0315054,-1.0972e-05,1.93687e-09,-1.35159e-13,22476.1,-42.385], Tmin=(964.242,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC=CC([CH2])C(26754)',
    structure = SMILES('[CH]=CC=CC([CH2])C'),
    E0 = (455.676,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0480298,'amu*angstrom^2'), symmetry=1, barrier=(15.4215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246431,'amu*angstrom^2'), symmetry=1, barrier=(15.5103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.663358,'amu*angstrom^2'), symmetry=1, barrier=(15.2519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664969,'amu*angstrom^2'), symmetry=1, barrier=(15.289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621252,0.0624768,-2.02993e-05,-2.21192e-08,1.44847e-11,54937.6,29.0245], Tmin=(100,'K'), Tmax=(949.448,'K')), NASAPolynomial(coeffs=[15.9989,0.0263704,-8.5654e-06,1.45617e-09,-1.00031e-13,50724.9,-51.1733], Tmin=(949.448,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(455.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([CH2])C=CC=C(26755)',
    structure = SMILES('[CH2]C([CH2])C=CC=C'),
    E0 = (413.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,358.694,358.732],'cm^-1')),
        HinderedRotor(inertia=(0.00131016,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144242,'amu*angstrom^2'), symmetry=1, barrier=(13.1659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144017,'amu*angstrom^2'), symmetry=1, barrier=(13.1656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.839162,'amu*angstrom^2'), symmetry=1, barrier=(76.625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736503,0.0603367,-1.49392e-05,-2.88351e-08,1.77213e-11,49880.2,29.2977], Tmin=(100,'K'), Tmax=(918.924,'K')), NASAPolynomial(coeffs=[15.5205,0.0258832,-7.50646e-06,1.1816e-09,-7.82364e-14,45900.7,-47.6457], Tmin=(918.924,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C=[C][C]=CC(C)C(26756)',
    structure = SMILES('[CH2]C#C[CH]C(C)C'),
    E0 = (353.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.91625,0.0601257,-2.6834e-05,-5.64255e-09,6.46678e-12,42624,26.8196], Tmin=(100,'K'), Tmax=(971.238,'K')), NASAPolynomial(coeffs=[11.5237,0.0334092,-1.17808e-05,2.01441e-09,-1.34731e-13,39763.1,-28.1658], Tmin=(971.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Propargyl) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=C[C]=CC(C)C(26757)',
    structure = SMILES('[CH]C=C=CC(C)C'),
    E0 = (426.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,275.513,275.514,275.515,275.515],'cm^-1')),
        HinderedRotor(inertia=(0.938297,'amu*angstrom^2'), symmetry=1, barrier=(50.5425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.938299,'amu*angstrom^2'), symmetry=1, barrier=(50.5425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.938301,'amu*angstrom^2'), symmetry=1, barrier=(50.5425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.938299,'amu*angstrom^2'), symmetry=1, barrier=(50.5425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.594112,0.0648618,-2.72301e-05,-3.74499e-09,4.25163e-12,51482.9,27.3896], Tmin=(100,'K'), Tmax=(1123.7,'K')), NASAPolynomial(coeffs=[12.4773,0.0394493,-1.58502e-05,2.87758e-09,-1.97155e-13,47746,-36.0698], Tmin=(1123.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(C)C=C1[CH]C1(26758)',
    structure = SMILES('[CH2]C(C)C=C1[CH]C1'),
    E0 = (442.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26037,0.0442325,2.83417e-05,-6.78852e-08,2.93151e-11,53359.6,25.402], Tmin=(100,'K'), Tmax=(957.172,'K')), NASAPolynomial(coeffs=[14.348,0.0290276,-9.71276e-06,1.72029e-09,-1.22214e-13,49045.3,-46.6145], Tmin=(957.172,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(Isobutyl)"""),
)

species(
    label = 'CC1C=[C][CH]CC1(26759)',
    structure = SMILES('CC1[CH][C]=CCC1'),
    E0 = (322.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87678,0.0257975,7.87027e-05,-1.14147e-07,4.34775e-11,38865.3,19.4695], Tmin=(100,'K'), Tmax=(979.642,'K')), NASAPolynomial(coeffs=[12.9631,0.0328861,-1.23162e-05,2.34734e-09,-1.73065e-13,34180.9,-46.6097], Tmin=(979.642,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(Cds_S) + radical(cyclohexene-allyl)"""),
)

species(
    label = 'C=CC=CC(=C)C(26760)',
    structure = SMILES('C=CC=CC(=C)C'),
    E0 = (108.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.540181,0.0614551,-1.12031e-05,-3.44223e-08,1.91302e-11,13139.5,23.6489], Tmin=(100,'K'), Tmax=(960.578,'K')), NASAPolynomial(coeffs=[17.8944,0.0243013,-8.01425e-06,1.41741e-09,-1.01049e-13,8185.61,-67.8089], Tmin=(960.578,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C=CC(C)C(26761)',
    structure = SMILES('C=C=C=CC(C)C'),
    E0 = (232.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.630536,0.0646193,-3.68253e-05,5.53222e-09,1.52132e-12,28138.6,25.2036], Tmin=(100,'K'), Tmax=(1137.21,'K')), NASAPolynomial(coeffs=[13.7534,0.0324759,-1.29132e-05,2.35106e-09,-1.61691e-13,24247.7,-43.7783], Tmin=(1137.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    label = '[CH2]CC=[C]C=C(26493)',
    structure = SMILES('[CH2]CC=[C]C=C'),
    E0 = (439.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,269.117,270.082],'cm^-1')),
        HinderedRotor(inertia=(0.299206,'amu*angstrom^2'), symmetry=1, barrier=(15.3925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297554,'amu*angstrom^2'), symmetry=1, barrier=(15.3949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97541,'amu*angstrom^2'), symmetry=1, barrier=(102.1,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3289.69,'J/mol'), sigma=(5.72141,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=513.84 K, Pc=39.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28443,0.0523914,-2.77929e-05,-1.18158e-09,4.46299e-12,52962.9,23.8431], Tmin=(100,'K'), Tmax=(993.884,'K')), NASAPolynomial(coeffs=[11.9769,0.0246018,-8.85774e-06,1.54887e-09,-1.05452e-13,50084.6,-31.4637], Tmin=(993.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=C[CH]CC(26762)',
    structure = SMILES('[CH2]C=[C]C=CCC'),
    E0 = (316.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865944,0.0585056,-1.52563e-05,-2.01906e-08,1.18576e-11,38161.3,26.9225], Tmin=(100,'K'), Tmax=(978.016,'K')), NASAPolynomial(coeffs=[13.0042,0.0324952,-1.16115e-05,2.03325e-09,-1.39164e-13,34656.7,-37.1461], Tmin=(978.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=CC[CH]C(26546)',
    structure = SMILES('C=C[C]=CC[CH]C'),
    E0 = (404.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.867493,0.0658641,-4.79219e-05,1.96362e-08,-3.42366e-12,48814.8,28.5322], Tmin=(100,'K'), Tmax=(1316.96,'K')), NASAPolynomial(coeffs=[10.6687,0.0360945,-1.40143e-05,2.47144e-09,-1.65212e-13,46233.2,-21.4515], Tmin=(1316.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJC) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC1=CC(C)C1(26763)',
    structure = SMILES('C=CC1=CC(C)C1'),
    E0 = (162.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25068,0.0376863,5.99861e-05,-1.09104e-07,4.56988e-11,19714.8,21.9043], Tmin=(100,'K'), Tmax=(949.122,'K')), NASAPolynomial(coeffs=[18.429,0.0225385,-6.55119e-06,1.18339e-09,-9.05871e-14,13875.3,-73.6584], Tmin=(949.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
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
    label = 'C=C[C]=C[CH]C(26604)',
    structure = SMILES('[CH2]C=[C]C=CC'),
    E0 = (338.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180],'cm^-1')),
        HinderedRotor(inertia=(2.72871,'amu*angstrom^2'), symmetry=1, barrier=(62.7384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.968464,'amu*angstrom^2'), symmetry=1, barrier=(22.2669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.7301,'amu*angstrom^2'), symmetry=1, barrier=(62.7705,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53745,0.0463667,-1.1901e-05,-1.57978e-08,9.31935e-12,40858.4,21.6], Tmin=(100,'K'), Tmax=(968.725,'K')), NASAPolynomial(coeffs=[10.6019,0.0268685,-9.47277e-06,1.63757e-09,-1.11014e-13,38260.9,-26.1849], Tmin=(968.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=[C]C1CC1C(26764)',
    structure = SMILES('[CH2]C=[C]C1CC1C'),
    E0 = (434.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24404,0.0420355,3.86642e-05,-8.04136e-08,3.38016e-11,52412.3,25.3543], Tmin=(100,'K'), Tmax=(966.187,'K')), NASAPolynomial(coeffs=[15.8517,0.027281,-9.41111e-06,1.73528e-09,-1.27113e-13,47455.4,-55.6604], Tmin=(966.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(434.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C)C=C=[C]C(26765)',
    structure = SMILES('[CH2]C(C)[CH]C#CC'),
    E0 = (402.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.404054,0.0679435,-5.37436e-05,2.50147e-08,-4.7118e-12,48499.3,30.0754], Tmin=(100,'K'), Tmax=(1427.12,'K')), NASAPolynomial(coeffs=[11.9579,0.0301453,-8.32394e-06,1.13879e-09,-6.35343e-14,45753,-27.8425], Tmin=(1427.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Isobutyl) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH2]C(C)[C]=C=CC(26766)',
    structure = SMILES('[CH2]C(C)C#C[CH]C'),
    E0 = (403.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.300935,0.0667188,-5.0413e-05,2.21346e-08,-3.89592e-12,48650.9,31.3466], Tmin=(100,'K'), Tmax=(1564.55,'K')), NASAPolynomial(coeffs=[12.1931,0.0292007,-7.62221e-06,9.9481e-10,-5.35926e-14,45800.5,-28.5665], Tmin=(1564.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C)C=C=CC(26767)',
    structure = SMILES('[CH2]C(C)=C[C]=CC'),
    E0 = (299.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.491938,0.0686237,-4.97023e-05,1.96852e-08,-3.22309e-12,36198.9,26.1167], Tmin=(100,'K'), Tmax=(1431.13,'K')), NASAPolynomial(coeffs=[13.4767,0.0323305,-1.16618e-05,1.96434e-09,-1.2741e-13,32482.4,-41.1817], Tmin=(1431.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C([CH2])C=C=CC(26768)',
    structure = SMILES('[CH2]C([CH2])C=C=CC'),
    E0 = (466.444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,243.742],'cm^-1')),
        HinderedRotor(inertia=(0.124466,'amu*angstrom^2'), symmetry=1, barrier=(5.23741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12811,'amu*angstrom^2'), symmetry=1, barrier=(5.23842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124293,'amu*angstrom^2'), symmetry=1, barrier=(5.24083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67045,'amu*angstrom^2'), symmetry=1, barrier=(70.2488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621613,0.0685063,-5.37835e-05,2.40238e-08,-4.46473e-12,56226.7,29.7392], Tmin=(100,'K'), Tmax=(1271.31,'K')), NASAPolynomial(coeffs=[11.9504,0.0328623,-1.17281e-05,1.97048e-09,-1.28043e-13,53346.2,-27.635], Tmin=(1271.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=C1[CH]C(C)C1(26769)',
    structure = SMILES('[CH2]C=C1[CH]C(C)C1'),
    E0 = (324.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.2416,-0.0148959,0.00014612,-1.39043e-07,3.23432e-11,38650.5,-16.9855], Tmin=(100,'K'), Tmax=(1679.84,'K')), NASAPolynomial(coeffs=[69.0783,0.0398293,-7.67146e-05,1.84325e-08,-1.36791e-12,-8605.45,-413.177], Tmin=(1679.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(methylenecyclobutane) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C)C1[C]=CC1(26770)',
    structure = SMILES('[CH2]C(C)C1[C]=CC1'),
    E0 = (513.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30838,0.0437836,2.68903e-05,-6.50025e-08,2.79975e-11,61892.4,26.7315], Tmin=(100,'K'), Tmax=(959.388,'K')), NASAPolynomial(coeffs=[13.8199,0.0294939,-9.98523e-06,1.77127e-09,-1.25447e-13,57748.7,-42.1951], Tmin=(959.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Isobutyl) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'C=C(C)C=C=CC(26771)',
    structure = SMILES('C=C(C)C=C=CC'),
    E0 = (160.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.297461,0.0708707,-5.31742e-05,2.07207e-08,-3.26002e-12,19491.9,24.5665], Tmin=(100,'K'), Tmax=(1504.99,'K')), NASAPolynomial(coeffs=[16.4078,0.0280524,-1.04981e-05,1.81645e-09,-1.19775e-13,14642.7,-59.7423], Tmin=(1504.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH]=C=CC([CH2])C(26772)',
    structure = SMILES('C#C[CH]C([CH2])C'),
    E0 = (444.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.00233445,'amu*angstrom^2'), symmetry=1, barrier=(9.42173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.15434,'amu*angstrom^2'), symmetry=1, barrier=(72.5244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.15391,'amu*angstrom^2'), symmetry=1, barrier=(72.5145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.085203,'amu*angstrom^2'), symmetry=1, barrier=(72.5369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.909241,0.0587591,-5.09716e-05,2.51885e-08,-4.93469e-12,53552.7,24.8399], Tmin=(100,'K'), Tmax=(1390.71,'K')), NASAPolynomial(coeffs=[11.6928,0.0223336,-5.84896e-06,7.61115e-10,-4.07369e-14,51076.5,-28.8606], Tmin=(1390.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Sec_Propargyl)"""),
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
    E0 = (407.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (527.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (529.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (602.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (666.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (505.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (477.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (568.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (503.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (642.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (622.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (525.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (600.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (678.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (833.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (595.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (559.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (736.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (545.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (495.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (496.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (432.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (859.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (602.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (567.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (415.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (720.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (446.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (615.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (645.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (483.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (519.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (533.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (526.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (432.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (825.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['C3H6(72)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['[CH2]C1[C]=CC(C)C1(26745)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(9.40382e+09,'s^-1'), n=0.352, Ea=(120.148,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C[C]=CC(=C)C(26746)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(72.3521,'m^3/(mol*s)'), n=1.66655, Ea=(10.8198,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeCs_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(C)C#CC=C(26747)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.66e+08,'cm^3/(mol*s)'), n=1.64, Ea=(12.2591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2703 used for Ct-Cs_Ct-Cd;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C(C)C=C=C=C(26748)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH3(17)', 'C=C[C]=CC=C(26601)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0143836,'m^3/(mol*s)'), n=2.40621, Ea=(24.5235,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C3H6(72)', '[CH]=[C]C=C(4699)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C3H6(T)(143)', 'CH2CHCCH(26391)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0117197,'m^3/(mol*s)'), n=2.33233, Ea=(9.74423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C[C]=C[C](C)C(26749)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.614e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C)[C]=CC=C(26750)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C)C=C[C]=C(26751)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['C=C[C]=[C]C(C)C(26752)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['[CH2][C](C)C=CC=C(26753)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC=CC([CH2])C(26754)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH2])C=CC=C(26755)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.4207e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['C=[C][C]=CC(C)C(26756)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.11562e+07,'s^-1'), n=1.49763, Ea=(188.308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_2H;Cd_H_out_doubleC] for rate rule [R5HJ_3;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C[C]=CC(C)C(26757)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(45.4385,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C3H6(T)(143)', '[CH]=[C]C=C(4699)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['[CH2]C(C)C=C1[CH]C1(26758)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['CC1C=[C][CH]CC1(26759)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.63396e+08,'s^-1'), n=0.483333, Ea=(87.4777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['C=CC=CC(=C)C(26760)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['C=C=C=CC(C)C(26761)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(S)(23)', '[CH2]CC=[C]C=C(26493)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['C=C[C]=C[CH]CC(26762)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['C=C[C]=CC[CH]C(26546)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['C=CC1=CC(C)C1(26763)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(19)', 'C=C[C]=C[CH]C(26604)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['[CH2]C=[C]C1CC1C(26764)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.50867e+09,'s^-1'), n=0.788889, Ea=(38.4463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['[CH2]C(C)C=C=[C]C(26765)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['[CH2]C(C)[C]=C=CC(26766)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.0029655,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['[CH2][C](C)C=C=CC(26767)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.21176e+06,'s^-1'), n=1.41298, Ea=(75.8094,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;C_rad_out_2H;XH_out] for rate rule [R5H_SMMS;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C([CH2])C=C=CC(26768)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(10873.3,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6H;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['[CH2]C=C1[CH]C(C)C1(26769)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(125.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['[CH2]C(C)C1[C]=CC1(26770)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(118.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(C)C=[C]C=C(26550)'],
    products = ['C=C(C)C=C=CC(26771)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(19)', '[CH]=C=CC([CH2])C(26772)'],
    products = ['[CH2]C(C)C=[C]C=C(26550)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4520',
    isomers = [
        '[CH2]C(C)C=[C]C=C(26550)',
    ],
    reactants = [
        ('C3H6(72)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4520',
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

