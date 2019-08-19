species(
    label = 'C7H9O(9)',
    structure = SMILES('[O]C1CC=CC=CC1'),
    E0 = (130.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4342.65,'J/mol'), sigma=(7.12519,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=678.31 K, Pc=27.24 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.232152,0.0680108,-2.23082e-05,-2.33238e-08,1.45863e-11,15873.3,11.8777], Tmin=(100,'K'), Tmax=(995.653,'K')), NASAPolynomial(coeffs=[18.7753,0.0266676,-9.96939e-06,1.85763e-09,-1.33863e-13,10537.5,-85.7536], Tmin=(995.653,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C=CC=CCC=O(956)',
    structure = SMILES('[CH2]C=CC=CCC=O'),
    E0 = (40.6633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,286.105,645.958],'cm^-1')),
        HinderedRotor(inertia=(0.0787087,'amu*angstrom^2'), symmetry=1, barrier=(22.8412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.994264,'amu*angstrom^2'), symmetry=1, barrier=(22.8601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0789567,'amu*angstrom^2'), symmetry=1, barrier=(22.8332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.994741,'amu*angstrom^2'), symmetry=1, barrier=(22.8711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3882.24,'J/mol'), sigma=(6.33784,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=606.40 K, Pc=34.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.868578,0.0528089,9.16168e-06,-4.58481e-08,1.96156e-11,5017.29,28.9089], Tmin=(100,'K'), Tmax=(1039.72,'K')), NASAPolynomial(coeffs=[15.7717,0.0320753,-1.37307e-05,2.68843e-09,-1.96067e-13,-60.0845,-53.0844], Tmin=(1039.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.6633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C7H9(3)',
    structure = SMILES('[CH]1CC=CC=CC1'),
    E0 = (267.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3862.82,'J/mol'), sigma=(6.55596,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=603.36 K, Pc=31.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05029,0.0548564,-1.71402e-05,-1.22379e-08,7.31431e-12,32264.5,9.34222], Tmin=(100,'K'), Tmax=(1067.28,'K')), NASAPolynomial(coeffs=[12.4249,0.0319753,-1.27385e-05,2.35045e-09,-1.64068e-13,28711.7,-51.544], Tmin=(1067.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(RCCJCC)"""),
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
    label = 'C=C=CC=CCC=O(976)',
    structure = SMILES('C=C=CC=CCC=O'),
    E0 = (99.2119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,540,610,2055,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.977632,'amu*angstrom^2'), symmetry=1, barrier=(22.4777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978595,'amu*angstrom^2'), symmetry=1, barrier=(22.4998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976925,'amu*angstrom^2'), symmetry=1, barrier=(22.4614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61164,0.0614653,-2.31198e-05,-1.21261e-08,7.8506e-12,12065.4,27.3318], Tmin=(100,'K'), Tmax=(1105.39,'K')), NASAPolynomial(coeffs=[16.7958,0.0283936,-1.28352e-05,2.53455e-09,-1.83671e-13,6930.01,-59.4137], Tmin=(1105.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.2119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC=C=CCC=O(1001)',
    structure = SMILES('C=CC=C=CCC=O'),
    E0 = (99.2119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,540,610,2055,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.977632,'amu*angstrom^2'), symmetry=1, barrier=(22.4777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978595,'amu*angstrom^2'), symmetry=1, barrier=(22.4998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976925,'amu*angstrom^2'), symmetry=1, barrier=(22.4614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61164,0.0614653,-2.31198e-05,-1.21261e-08,7.8506e-12,12065.4,27.3318], Tmin=(100,'K'), Tmax=(1105.39,'K')), NASAPolynomial(coeffs=[16.7958,0.0283936,-1.28352e-05,2.53455e-09,-1.83671e-13,6930.01,-59.4137], Tmin=(1105.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.2119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC=CC=CC=O(1013)',
    structure = SMILES('C=CC=CC=CC=O'),
    E0 = (22.2051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.976937,'amu*angstrom^2'), symmetry=1, barrier=(22.4617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976304,'amu*angstrom^2'), symmetry=1, barrier=(22.4472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978369,'amu*angstrom^2'), symmetry=1, barrier=(22.4946,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279345,0.0696767,-4.23342e-05,2.69834e-09,4.17388e-12,2815.02,27.0664], Tmin=(100,'K'), Tmax=(1070.38,'K')), NASAPolynomial(coeffs=[18.1419,0.0249069,-1.04006e-05,1.99581e-09,-1.43366e-13,-2268.19,-66.2073], Tmin=(1070.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.2051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'HCO(30)',
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
    label = 'C=CC=CC=C(1014)',
    structure = SMILES('C=CC=CC=C'),
    E0 = (145.693,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08797,'amu*angstrom^2'), symmetry=1, barrier=(25.0147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0889,'amu*angstrom^2'), symmetry=1, barrier=(25.0359,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38428,0.0414412,1.94287e-05,-6.16689e-08,2.86504e-11,17631.7,19.3258], Tmin=(100,'K'), Tmax=(941.215,'K')), NASAPolynomial(coeffs=[16.962,0.0157222,-4.10118e-06,6.95736e-10,-5.26521e-14,12906.2,-64.4098], Tmin=(941.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CH2CHO(51)',
    structure = SMILES('[CH2]C=O'),
    E0 = (1.22925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,526.75,532.597,975.94,1639.13,1641.45],'cm^-1')),
        HinderedRotor(inertia=(0.00114821,'amu*angstrom^2'), symmetry=1, barrier=(2.1986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66874,0.0096233,1.60617e-05,-2.87682e-08,1.2503e-11,219.438,12.5694], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91637,0.0088465,-3.14955e-06,5.05413e-10,-3.01305e-14,-1047.8,-6.1065], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(1.22925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2CHO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=CC=C[CH2](324)',
    structure = SMILES('[CH]=CC=C[CH2]'),
    E0 = (423.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.52804,'amu*angstrom^2'), symmetry=1, barrier=(35.1326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53309,'amu*angstrom^2'), symmetry=1, barrier=(35.2488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22719,0.0272982,2.2821e-05,-5.20811e-08,2.2993e-11,50955.4,18.1254], Tmin=(100,'K'), Tmax=(937.462,'K')), NASAPolynomial(coeffs=[12.1491,0.0145308,-4.06041e-06,6.79572e-10,-4.92004e-14,47795.9,-36.031], Tmin=(937.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=CC=C[CH2](983)',
    structure = SMILES('[CH2]C=CC=C[CH2]'),
    E0 = (257.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,467.896],'cm^-1')),
        HinderedRotor(inertia=(0.359675,'amu*angstrom^2'), symmetry=1, barrier=(55.8581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359783,'amu*angstrom^2'), symmetry=1, barrier=(55.858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359749,'amu*angstrom^2'), symmetry=1, barrier=(55.8569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91527,0.0306757,4.13798e-05,-7.59835e-08,3.17674e-11,31116.9,21.5125], Tmin=(100,'K'), Tmax=(942.95,'K')), NASAPolynomial(coeffs=[12.7392,0.0229469,-7.07071e-06,1.21797e-09,-8.6992e-14,27378,-39.0729], Tmin=(942.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=CC=C[CH]C=O(984)',
    structure = SMILES('[CH2]C=CC=CC=C[O]'),
    E0 = (160.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3001,3007,3013,3019,3025,975,980,985,990,995,1000,1300,1315,1330,1345,1360,1375,400,420,440,460,480,500,1630,1640,1650,1660,1670,1680,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.38115,'amu*angstrom^2'), symmetry=1, barrier=(31.7554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3807,'amu*angstrom^2'), symmetry=1, barrier=(31.745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37859,'amu*angstrom^2'), symmetry=1, barrier=(31.6965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.377903,0.0570618,2.05373e-05,-8.36393e-08,4.10182e-11,19442.2,28.2652], Tmin=(100,'K'), Tmax=(925.175,'K')), NASAPolynomial(coeffs=[24.2658,0.0120925,-1.09251e-06,7.03681e-11,-1.00541e-14,12526.6,-98.6092], Tmin=(925.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=CC=[C]CC=O(985)',
    structure = SMILES('[CH2]C=CC=[C]CC=O'),
    E0 = (278.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,939.653,941.219],'cm^-1')),
        HinderedRotor(inertia=(0.985009,'amu*angstrom^2'), symmetry=1, barrier=(22.6473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.749895,'amu*angstrom^2'), symmetry=1, barrier=(22.6448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.985694,'amu*angstrom^2'), symmetry=1, barrier=(22.6631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.985635,'amu*angstrom^2'), symmetry=1, barrier=(22.6617,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.815756,0.0575814,-1.59213e-05,-1.66823e-08,8.92157e-12,33621.6,29.576], Tmin=(100,'K'), Tmax=(1099.34,'K')), NASAPolynomial(coeffs=[15.0696,0.0307889,-1.35718e-05,2.63719e-09,-1.89258e-13,28972.6,-47.4307], Tmin=(1099.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH2](321)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,229.711,230.18,230.787],'cm^-1')),
        HinderedRotor(inertia=(1.33306,'amu*angstrom^2'), symmetry=1, barrier=(50.5153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,44095,-3.44606], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CCC=O(986)',
    structure = SMILES('[CH]=CCC=O'),
    E0 = (153.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.12336,'amu*angstrom^2'), symmetry=1, barrier=(18.0046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122606,'amu*angstrom^2'), symmetry=1, barrier=(18.0324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51056,0.0244663,1.05068e-05,-2.60457e-08,9.69719e-12,18577.3,17.4244], Tmin=(100,'K'), Tmax=(1119.61,'K')), NASAPolynomial(coeffs=[9.62674,0.019207,-9.46251e-06,1.93136e-09,-1.41855e-13,15720,-23.3551], Tmin=(1119.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C[C]=CCC=O(987)',
    structure = SMILES('[CH2]C=C[C]=CCC=O'),
    E0 = (239.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,273.974,671.565],'cm^-1')),
        HinderedRotor(inertia=(1.8213,'amu*angstrom^2'), symmetry=1, barrier=(41.8753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1465,'amu*angstrom^2'), symmetry=1, barrier=(26.3603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.082324,'amu*angstrom^2'), symmetry=1, barrier=(26.3624,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14661,'amu*angstrom^2'), symmetry=1, barrier=(26.3628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.751815,0.0597784,-2.12599e-05,-1.15117e-08,7.3498e-12,28951.1,28.9266], Tmin=(100,'K'), Tmax=(1092.92,'K')), NASAPolynomial(coeffs=[14.6319,0.0315089,-1.33834e-05,2.54615e-09,-1.80551e-13,24571.5,-45.4259], Tmin=(1092.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=[C]C=CCC=O(988)',
    structure = SMILES('[CH2]C=[C]C=CCC=O'),
    E0 = (239.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,273.974,671.565],'cm^-1')),
        HinderedRotor(inertia=(1.8213,'amu*angstrom^2'), symmetry=1, barrier=(41.8753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1465,'amu*angstrom^2'), symmetry=1, barrier=(26.3603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.082324,'amu*angstrom^2'), symmetry=1, barrier=(26.3624,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14661,'amu*angstrom^2'), symmetry=1, barrier=(26.3628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.751815,0.0597784,-2.12599e-05,-1.15117e-08,7.3498e-12,28951.1,28.9266], Tmin=(100,'K'), Tmax=(1092.92,'K')), NASAPolynomial(coeffs=[14.6319,0.0315089,-1.33834e-05,2.54615e-09,-1.80551e-13,24571.5,-45.4259], Tmin=(1092.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=CC=CCC=O(989)',
    structure = SMILES('[CH2][C]=CC=CCC=O'),
    E0 = (278.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,939.653,941.219],'cm^-1')),
        HinderedRotor(inertia=(0.985009,'amu*angstrom^2'), symmetry=1, barrier=(22.6473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.749895,'amu*angstrom^2'), symmetry=1, barrier=(22.6448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.985694,'amu*angstrom^2'), symmetry=1, barrier=(22.6631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.985635,'amu*angstrom^2'), symmetry=1, barrier=(22.6617,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.815756,0.0575814,-1.59213e-05,-1.66823e-08,8.92157e-12,33621.6,29.576], Tmin=(100,'K'), Tmax=(1099.34,'K')), NASAPolynomial(coeffs=[15.0696,0.0307889,-1.35718e-05,2.63719e-09,-1.89258e-13,28972.6,-47.4307], Tmin=(1099.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC=CC[C]=O(990)',
    structure = SMILES('[CH2]C=CC=CC[C]=O'),
    E0 = (200.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,1855,455,950,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,398.202,398.862],'cm^-1')),
        HinderedRotor(inertia=(0.168801,'amu*angstrom^2'), symmetry=1, barrier=(18.9466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168322,'amu*angstrom^2'), symmetry=1, barrier=(18.9461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06321,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168646,'amu*angstrom^2'), symmetry=1, barrier=(18.9521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.807102,0.0556641,-5.62122e-06,-3.1726e-08,1.53577e-11,24257.1,29.9468], Tmin=(100,'K'), Tmax=(1032.71,'K')), NASAPolynomial(coeffs=[16.3527,0.0275993,-1.15522e-05,2.24654e-09,-1.63698e-13,19332,-53.852], Tmin=(1032.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]C=CC=CCC=O(999)',
    structure = SMILES('[CH]C=CC=CCC=O'),
    E0 = (293.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.568808,0.0617787,-1.37506e-05,-2.11433e-08,1.04844e-11,35409.9,29.5412], Tmin=(100,'K'), Tmax=(1102.37,'K')), NASAPolynomial(coeffs=[15.5755,0.0352773,-1.57226e-05,3.04971e-09,-2.18215e-13,30403,-52.0231], Tmin=(1102.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C2H3(44)',
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
    label = '[CH]C=CCC=O(1010)',
    structure = SMILES('[CH]C=CCC=O'),
    E0 = (241.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,399.149,399.149,399.149,399.15],'cm^-1')),
        HinderedRotor(inertia=(0.433694,'amu*angstrom^2'), symmetry=1, barrier=(49.0322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.433694,'amu*angstrom^2'), symmetry=1, barrier=(49.0322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.433693,'amu*angstrom^2'), symmetry=1, barrier=(49.0322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90949,0.0366928,6.54401e-06,-2.56812e-08,9.48674e-12,29131,22.0986], Tmin=(100,'K'), Tmax=(1159.72,'K')), NASAPolynomial(coeffs=[10.0066,0.0308595,-1.44886e-05,2.83732e-09,-2.02355e-13,25767.1,-24.5711], Tmin=(1159.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]CC=CC=CC=O(994)',
    structure = SMILES('[CH2]CC=CC=CC=O'),
    E0 = (117.008,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,261.532,261.782],'cm^-1')),
        HinderedRotor(inertia=(0.296558,'amu*angstrom^2'), symmetry=1, barrier=(14.4648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297856,'amu*angstrom^2'), symmetry=1, barrier=(14.4587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297434,'amu*angstrom^2'), symmetry=1, barrier=(14.4586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296809,'amu*angstrom^2'), symmetry=1, barrier=(14.4588,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.404746,0.0721851,-5.29331e-05,1.90749e-08,-2.76523e-12,14207.7,30.4947], Tmin=(100,'K'), Tmax=(1604.27,'K')), NASAPolynomial(coeffs=[17.6442,0.0292007,-1.2742e-05,2.37308e-09,-1.62484e-13,8676.36,-60.8241], Tmin=(1604.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJ)"""),
)

species(
    label = 'O=CC[CH]C1C=CC1(975)',
    structure = SMILES('O=CC[CH]C1C=CC1'),
    E0 = (182.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18782,0.0470578,1.66891e-05,-5.08096e-08,2.13332e-11,22034.7,29.7193], Tmin=(100,'K'), Tmax=(1011.52,'K')), NASAPolynomial(coeffs=[13.6433,0.0323759,-1.28072e-05,2.42022e-09,-1.73905e-13,17746.3,-39.2565], Tmin=(1011.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(CCJCC=O)"""),
)

species(
    label = 'C[C]=CC=CCC=O(977)',
    structure = SMILES('C[C]=CC=CCC=O'),
    E0 = (160.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.476274,0.0669099,-3.93847e-05,7.72978e-09,3.19428e-13,19433.1,29.255], Tmin=(100,'K'), Tmax=(1324.51,'K')), NASAPolynomial(coeffs=[16.9266,0.030823,-1.39103e-05,2.65587e-09,-1.85343e-13,13883.1,-59.2323], Tmin=(1324.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(Cds_S)"""),
)

species(
    label = 'CC=[C]C=CCC=O(978)',
    structure = SMILES('CC=[C]C=CCC=O'),
    E0 = (121.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.8422,'amu*angstrom^2'), symmetry=1, barrier=(19.3638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84219,'amu*angstrom^2'), symmetry=1, barrier=(19.3636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84224,'amu*angstrom^2'), symmetry=1, barrier=(19.3647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842149,'amu*angstrom^2'), symmetry=1, barrier=(19.3627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.310614,0.0701047,-4.73904e-05,1.53952e-08,-1.98703e-12,14767.6,28.9856], Tmin=(100,'K'), Tmax=(1806.44,'K')), NASAPolynomial(coeffs=[20.0558,0.0263827,-1.1085e-05,1.99658e-09,-1.32732e-13,7633.94,-77.9499], Tmin=(1806.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CJC=C)"""),
)

species(
    label = 'CC=C[C]=CCC=O(979)',
    structure = SMILES('CC=C[C]=CCC=O'),
    E0 = (121.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.8422,'amu*angstrom^2'), symmetry=1, barrier=(19.3638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84219,'amu*angstrom^2'), symmetry=1, barrier=(19.3636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84224,'amu*angstrom^2'), symmetry=1, barrier=(19.3647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842149,'amu*angstrom^2'), symmetry=1, barrier=(19.3627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.310614,0.0701047,-4.73904e-05,1.53952e-08,-1.98703e-12,14767.6,28.9856], Tmin=(100,'K'), Tmax=(1806.44,'K')), NASAPolynomial(coeffs=[20.0558,0.0263827,-1.1085e-05,1.99658e-09,-1.32732e-13,7633.94,-77.9499], Tmin=(1806.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CJC=C)"""),
)

species(
    label = 'CC=CC=[C]CC=O(980)',
    structure = SMILES('CC=CC=[C]CC=O'),
    E0 = (160.449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,252.151,252.154],'cm^-1')),
        HinderedRotor(inertia=(0.342913,'amu*angstrom^2'), symmetry=1, barrier=(15.4635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.342876,'amu*angstrom^2'), symmetry=1, barrier=(15.4633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.342718,'amu*angstrom^2'), symmetry=1, barrier=(15.4632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.342603,'amu*angstrom^2'), symmetry=1, barrier=(15.464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.476274,0.0669099,-3.93847e-05,7.72978e-09,3.19428e-13,19433.1,29.255], Tmin=(100,'K'), Tmax=(1324.51,'K')), NASAPolynomial(coeffs=[16.9266,0.030823,-1.39103e-05,2.65587e-09,-1.85343e-13,13883.1,-59.2323], Tmin=(1324.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(Cds_S)"""),
)

species(
    label = 'CC=CC=C[CH]C=O(981)',
    structure = SMILES('CC=C[CH]C=CC=O'),
    E0 = (30.9611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.785265,0.0596031,-2.07619e-05,-9.46316e-09,5.90608e-12,3848.7,28.9271], Tmin=(100,'K'), Tmax=(1146.43,'K')), NASAPolynomial(coeffs=[13.8929,0.0345032,-1.49184e-05,2.8382e-09,-2.0001e-13,-512.637,-42.0144], Tmin=(1146.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.9611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC=CC=CC[C]=O(982)',
    structure = SMILES('CC=CC=CC[C]=O'),
    E0 = (82.5682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,181.694,183.186],'cm^-1')),
        HinderedRotor(inertia=(0.646968,'amu*angstrom^2'), symmetry=1, barrier=(15.1685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.651044,'amu*angstrom^2'), symmetry=1, barrier=(15.1677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.653886,'amu*angstrom^2'), symmetry=1, barrier=(15.1711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.636691,'amu*angstrom^2'), symmetry=1, barrier=(15.1733,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.486477,0.0649411,-2.96674e-05,-5.54895e-09,5.61235e-12,10067.5,29.5463], Tmin=(100,'K'), Tmax=(1121.12,'K')), NASAPolynomial(coeffs=[16.5176,0.0302112,-1.32602e-05,2.56928e-09,-1.83826e-13,5060.98,-55.9246], Tmin=(1121.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.5682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CCC=CC1[CH]C1(991)',
    structure = SMILES('O=CCC=CC1[CH]C1'),
    E0 = (204.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37409,0.0406664,3.45858e-05,-6.59909e-08,2.50784e-11,24758.7,29.7165], Tmin=(100,'K'), Tmax=(1051.57,'K')), NASAPolynomial(coeffs=[13.8803,0.0344141,-1.54343e-05,3.08572e-09,-2.26972e-13,19843.9,-42.1101], Tmin=(1051.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(cyclopropane)"""),
)

species(
    label = 'O=CCC1[CH]C=CC1(992)',
    structure = SMILES('O=CCC1[CH]C=CC1'),
    E0 = (9.90119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49319,0.0380587,4.25197e-05,-7.54543e-08,2.93049e-11,1296.08,24.1063], Tmin=(100,'K'), Tmax=(1008.07,'K')), NASAPolynomial(coeffs=[12.8624,0.0344409,-1.38423e-05,2.65377e-09,-1.92713e-13,-3104.48,-41.2925], Tmin=(1008.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.90119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopentene) + radical(cyclopentene-allyl)"""),
)

species(
    label = '[CH]1CC=CC=CCO1(993)',
    structure = SMILES('[CH]1CC=CC=CCO1'),
    E0 = (153.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63134,0.0275449,8.72609e-05,-1.31975e-07,5.20142e-11,18521.2,21.6968], Tmin=(100,'K'), Tmax=(961.243,'K')), NASAPolynomial(coeffs=[16.4202,0.0277052,-9.27243e-06,1.75252e-09,-1.32938e-13,12827.5,-63.8939], Tmin=(961.243,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclooctane) + radical(CCsJOCs)"""),
)

species(
    label = 'CO(28)',
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
    label = '[CH2]C=CC=CC(995)',
    structure = SMILES('[CH2]C=CC=CC'),
    E0 = (139.926,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,208.671],'cm^-1')),
        HinderedRotor(inertia=(0.617553,'amu*angstrom^2'), symmetry=1, barrier=(19.0743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618352,'amu*angstrom^2'), symmetry=1, barrier=(19.0735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.84468,'amu*angstrom^2'), symmetry=1, barrier=(87.8036,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62189,0.0397423,1.75012e-05,-4.91152e-08,2.12944e-11,16926.1,21.7005], Tmin=(100,'K'), Tmax=(970.029,'K')), NASAPolynomial(coeffs=[12.1207,0.026819,-9.4764e-06,1.7006e-09,-1.20078e-13,13460.4,-35.9954], Tmin=(970.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=CC=COC=C(996)',
    structure = SMILES('C=CC=C[CH]OC=C'),
    E0 = (86.7631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.135973,0.0665063,-8.28904e-06,-4.70386e-08,2.55813e-11,10591.1,29.2256], Tmin=(100,'K'), Tmax=(953.19,'K')), NASAPolynomial(coeffs=[22.031,0.0198182,-5.9365e-06,1.05698e-09,-7.89965e-14,4364.01,-86.1252], Tmin=(953.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.7631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH2]C=CC=CC=CO(997)',
    structure = SMILES('[CH2]C=CC=CC=CO'),
    E0 = (18.9307,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3001,3007,3013,3019,3025,975,980,985,990,995,1000,1300,1315,1330,1345,1360,1375,400,420,440,460,480,500,1630,1640,1650,1660,1670,1680,3615,1277.5,1000,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23287,'amu*angstrom^2'), symmetry=1, barrier=(28.3461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23237,'amu*angstrom^2'), symmetry=1, barrier=(28.3347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23235,'amu*angstrom^2'), symmetry=1, barrier=(28.3342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23257,'amu*angstrom^2'), symmetry=1, barrier=(28.3391,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.03794,0.0640684,1.70946e-05,-8.99754e-08,4.59105e-11,2445.3,28.3485], Tmin=(100,'K'), Tmax=(913.516,'K')), NASAPolynomial(coeffs=[27.607,0.0082577,1.61541e-06,-5.04222e-10,3.11802e-14,-5327.6,-117.421], Tmin=(913.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.9307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'CH2(T)(33)',
    structure = SMILES('[CH2]'),
    E0 = (381.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([971.045,2816.03,3444.23],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71758,0.00127391,2.17347e-06,-3.48858e-09,1.65209e-12,45872.4,1.75298], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.14632,0.00303671,-9.96474e-07,1.50484e-10,-8.57336e-15,46041.3,4.72342], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(381.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=CC=CCC=O(998)',
    structure = SMILES('[CH]=CC=CCC=O'),
    E0 = (205.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.88079,'amu*angstrom^2'), symmetry=1, barrier=(20.2511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.880746,'amu*angstrom^2'), symmetry=1, barrier=(20.2501,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.880591,'amu*angstrom^2'), symmetry=1, barrier=(20.2465,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17864,0.0494311,-9.2841e-06,-2.22687e-08,1.10615e-11,24855.9,24.8368], Tmin=(100,'K'), Tmax=(1068.93,'K')), NASAPolynomial(coeffs=[15.325,0.0234295,-1.05938e-05,2.12109e-09,-1.55935e-13,20292.8,-51.5523], Tmin=(1068.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1C=CCC1[O](1000)',
    structure = SMILES('C=CC1C=CCC1[O]'),
    E0 = (151.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0311,0.0490888,1.81025e-05,-5.62172e-08,2.41154e-11,18353.2,24.4341], Tmin=(100,'K'), Tmax=(997.852,'K')), NASAPolynomial(coeffs=[15.159,0.0307352,-1.18511e-05,2.23973e-09,-1.622e-13,13627.9,-53.244], Tmin=(997.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopentene) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC[C]=CCC=O(1002)',
    structure = SMILES('C=CC[C]=CCC=O'),
    E0 = (191.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,180,521.634,521.781],'cm^-1')),
        HinderedRotor(inertia=(0.0810792,'amu*angstrom^2'), symmetry=1, barrier=(1.86417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664181,'amu*angstrom^2'), symmetry=1, barrier=(15.2708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0792741,'amu*angstrom^2'), symmetry=1, barrier=(15.2691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0790749,'amu*angstrom^2'), symmetry=1, barrier=(15.27,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.755749,0.0614967,-2.8708e-05,1.51168e-10,2.14473e-12,23133.3,31.2771], Tmin=(100,'K'), Tmax=(1315.48,'K')), NASAPolynomial(coeffs=[15.4816,0.0332695,-1.53928e-05,2.96693e-09,-2.07918e-13,17827,-49.2474], Tmin=(1315.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC=CCC=O(1003)',
    structure = SMILES('C=[C]CC=CCC=O'),
    E0 = (191.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,180,521.634,521.781],'cm^-1')),
        HinderedRotor(inertia=(0.0810792,'amu*angstrom^2'), symmetry=1, barrier=(1.86417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664181,'amu*angstrom^2'), symmetry=1, barrier=(15.2708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0792741,'amu*angstrom^2'), symmetry=1, barrier=(15.2691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0790749,'amu*angstrom^2'), symmetry=1, barrier=(15.27,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.755749,0.0614967,-2.8708e-05,1.51168e-10,2.14473e-12,23133.3,31.2771], Tmin=(100,'K'), Tmax=(1315.48,'K')), NASAPolynomial(coeffs=[15.4816,0.0332695,-1.53928e-05,2.96693e-09,-2.07918e-13,17827,-49.2474], Tmin=(1315.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CCC=[C]CC=O(1004)',
    structure = SMILES('C=CCC=[C]CC=O'),
    E0 = (191.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,180,521.634,521.781],'cm^-1')),
        HinderedRotor(inertia=(0.0810792,'amu*angstrom^2'), symmetry=1, barrier=(1.86417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664181,'amu*angstrom^2'), symmetry=1, barrier=(15.2708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0792741,'amu*angstrom^2'), symmetry=1, barrier=(15.2691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0790749,'amu*angstrom^2'), symmetry=1, barrier=(15.27,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.755749,0.0614967,-2.8708e-05,1.51168e-10,2.14473e-12,23133.3,31.2771], Tmin=(100,'K'), Tmax=(1315.48,'K')), NASAPolynomial(coeffs=[15.4816,0.0332695,-1.53928e-05,2.96693e-09,-2.07918e-13,17827,-49.2474], Tmin=(1315.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CCC=CCC=O(1005)',
    structure = SMILES('[CH]=CCC=CCC=O'),
    E0 = (200.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,322.682,323.798],'cm^-1')),
        HinderedRotor(inertia=(0.00160135,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225942,'amu*angstrom^2'), symmetry=1, barrier=(16.6429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223747,'amu*angstrom^2'), symmetry=1, barrier=(16.6465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226579,'amu*angstrom^2'), symmetry=1, barrier=(16.6449,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.711407,0.0606206,-2.20042e-05,-8.89489e-09,5.58199e-12,24249.5,31.3549], Tmin=(100,'K'), Tmax=(1187.53,'K')), NASAPolynomial(coeffs=[15.375,0.0332829,-1.53303e-05,2.99698e-09,-2.13722e-13,19211.8,-48.4563], Tmin=(1187.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC=C[CH]C=O(1006)',
    structure = SMILES('C=CCC=CC=C[O]'),
    E0 = (73.1917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.96808,'amu*angstrom^2'), symmetry=1, barrier=(22.2581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.968374,'amu*angstrom^2'), symmetry=1, barrier=(22.2648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.96485,'amu*angstrom^2'), symmetry=1, barrier=(22.1838,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.314498,0.0613528,4.8243e-06,-6.06285e-08,3.05529e-11,8953.59,29.9564], Tmin=(100,'K'), Tmax=(947.388,'K')), NASAPolynomial(coeffs=[21.999,0.0187248,-5.14903e-06,9.01695e-10,-6.87071e-14,2649.15,-85.0752], Tmin=(947.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.1917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'C=CCC=CC[C]=O(1007)',
    structure = SMILES('C=CCC=CC[C]=O'),
    E0 = (113.422,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1855,455,950,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,344,344.004,344.005],'cm^-1')),
        HinderedRotor(inertia=(0.0014245,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174252,'amu*angstrom^2'), symmetry=1, barrier=(14.6323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174243,'amu*angstrom^2'), symmetry=1, barrier=(14.6323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174242,'amu*angstrom^2'), symmetry=1, barrier=(14.6323,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.724019,0.0600026,-2.05543e-05,-1.12438e-08,6.69853e-12,13769.5,31.7201], Tmin=(100,'K'), Tmax=(1146.51,'K')), NASAPolynomial(coeffs=[15.2778,0.0323215,-1.4554e-05,2.83666e-09,-2.02834e-13,8914.38,-47.1031], Tmin=(1146.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC1[CH]C1CC=O(1008)',
    structure = SMILES('C=CC1[CH]C1CC=O'),
    E0 = (200.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08594,0.0523977,-3.97677e-06,-2.58721e-08,1.17321e-11,24269.6,29.9899], Tmin=(100,'K'), Tmax=(1062.68,'K')), NASAPolynomial(coeffs=[12.4466,0.0346306,-1.41798e-05,2.66263e-09,-1.87975e-13,20443.7,-32.1497], Tmin=(1062.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(cyclopropane)"""),
)

species(
    label = 'C=CC1C=CC[CH]O1(1009)',
    structure = SMILES('C=CC1C=CC[CH]O1'),
    E0 = (118.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924222,0.0476811,3.58813e-05,-8.64214e-08,3.85385e-11,14337.8,24.2325], Tmin=(100,'K'), Tmax=(937.921,'K')), NASAPolynomial(coeffs=[18.3732,0.0237327,-6.52958e-06,1.09292e-09,-7.96746e-14,8844.9,-70.664], Tmin=(937.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(36dihydro2hpyran) + radical(CCsJOCs)"""),
)

species(
    label = 'C=CC=CC1CC1[O](1011)',
    structure = SMILES('C=CC=CC1CC1[O]'),
    E0 = (223.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.741229,0.049624,3.46778e-05,-8.82748e-08,3.94592e-11,27076.6,27.4356], Tmin=(100,'K'), Tmax=(950.07,'K')), NASAPolynomial(coeffs=[20.7926,0.0202163,-5.74869e-06,1.03987e-09,-8.04453e-14,20783.7,-81.3405], Tmin=(950.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1C=CC1CC=O(1012)',
    structure = SMILES('[CH2]C1C=CC1CC=O'),
    E0 = (179.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06831,0.0506478,8.8427e-06,-4.55277e-08,2.05938e-11,21702.6,28.4934], Tmin=(100,'K'), Tmax=(981.794,'K')), NASAPolynomial(coeffs=[13.8394,0.0312511,-1.13831e-05,2.06282e-09,-1.45632e-13,17622,-40.8955], Tmin=(981.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC=CC[CH]C=O(1015)',
    structure = SMILES('C=CC=CCC=C[O]'),
    E0 = (73.1917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.96808,'amu*angstrom^2'), symmetry=1, barrier=(22.2581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.968374,'amu*angstrom^2'), symmetry=1, barrier=(22.2648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.96485,'amu*angstrom^2'), symmetry=1, barrier=(22.1838,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.314498,0.0613528,4.8243e-06,-6.06285e-08,3.05529e-11,8953.59,29.9564], Tmin=(100,'K'), Tmax=(947.388,'K')), NASAPolynomial(coeffs=[21.999,0.0187248,-5.14903e-06,9.01695e-10,-6.87071e-14,2649.15,-85.0752], Tmin=(947.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.1917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'C=CC=[C]CCC=O(1016)',
    structure = SMILES('C=CC=[C]CCC=O'),
    E0 = (166.512,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,302.978,303.158,303.305],'cm^-1')),
        HinderedRotor(inertia=(0.188191,'amu*angstrom^2'), symmetry=1, barrier=(12.2542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18787,'amu*angstrom^2'), symmetry=1, barrier=(12.2524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188139,'amu*angstrom^2'), symmetry=1, barrier=(12.252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1877,'amu*angstrom^2'), symmetry=1, barrier=(12.2497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.675351,0.0726889,-5.87869e-05,2.56885e-08,-4.70419e-12,20146.7,29.3962], Tmin=(100,'K'), Tmax=(1264.94,'K')), NASAPolynomial(coeffs=[12.2051,0.0362299,-1.55531e-05,2.90289e-09,-2.00935e-13,17229.8,-28.9378], Tmin=(1264.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CC=CCC[C]=O(1017)',
    structure = SMILES('C=CC=CCC[C]=O'),
    E0 = (88.6311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1855,455,950,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.809227,'amu*angstrom^2'), symmetry=1, barrier=(18.6057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.80933,'amu*angstrom^2'), symmetry=1, barrier=(18.6081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.809322,'amu*angstrom^2'), symmetry=1, barrier=(18.6079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.80921,'amu*angstrom^2'), symmetry=1, barrier=(18.6053,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.312095,0.0748064,-6.18868e-05,2.6974e-08,-4.77202e-12,10797.8,31.0486], Tmin=(100,'K'), Tmax=(1342.19,'K')), NASAPolynomial(coeffs=[15.3321,0.0300443,-1.18624e-05,2.12718e-09,-1.44036e-13,6765.8,-45.8348], Tmin=(1342.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.6311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C[C]=CCCC=O(1018)',
    structure = SMILES('C=C[C]=CCCC=O'),
    E0 = (127.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,231.405,231.41,232.311],'cm^-1')),
        HinderedRotor(inertia=(0.341281,'amu*angstrom^2'), symmetry=1, barrier=(12.9732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339947,'amu*angstrom^2'), symmetry=1, barrier=(12.979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343093,'amu*angstrom^2'), symmetry=1, barrier=(12.9765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.73995,'amu*angstrom^2'), symmetry=1, barrier=(104.787,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.593571,0.0750347,-6.44205e-05,3.10292e-08,-6.29377e-12,15477.2,28.8157], Tmin=(100,'K'), Tmax=(1154.05,'K')), NASAPolynomial(coeffs=[11.4035,0.0375665,-1.57198e-05,2.89577e-09,-1.9919e-13,12982.2,-24.8846], Tmin=(1154.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=CCCC=O(1019)',
    structure = SMILES('C=[C]C=CCCC=O'),
    E0 = (127.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,231.405,231.41,232.311],'cm^-1')),
        HinderedRotor(inertia=(0.341281,'amu*angstrom^2'), symmetry=1, barrier=(12.9732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339947,'amu*angstrom^2'), symmetry=1, barrier=(12.979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343093,'amu*angstrom^2'), symmetry=1, barrier=(12.9765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.73995,'amu*angstrom^2'), symmetry=1, barrier=(104.787,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.593571,0.0750347,-6.44205e-05,3.10292e-08,-6.29377e-12,15477.2,28.8157], Tmin=(100,'K'), Tmax=(1154.05,'K')), NASAPolynomial(coeffs=[11.4035,0.0375665,-1.57198e-05,2.89577e-09,-1.9919e-13,12982.2,-24.8846], Tmin=(1154.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CC=CCCC=O(1020)',
    structure = SMILES('[CH]=CC=CCCC=O'),
    E0 = (175.767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,225.781,226.054],'cm^-1')),
        HinderedRotor(inertia=(0.375479,'amu*angstrom^2'), symmetry=1, barrier=(13.5696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374879,'amu*angstrom^2'), symmetry=1, barrier=(13.5711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374441,'amu*angstrom^2'), symmetry=1, barrier=(13.5681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.37533,'amu*angstrom^2'), symmetry=1, barrier=(13.5677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.404869,0.0742727,-5.97156e-05,2.5173e-08,-4.33988e-12,21273.1,30.2989], Tmin=(100,'K'), Tmax=(1363.08,'K')), NASAPolynomial(coeffs=[14.8176,0.0319778,-1.31717e-05,2.40873e-09,-1.64691e-13,17344,-43.6982], Tmin=(1363.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CC=CC1C[CH]O1(1021)',
    structure = SMILES('C=CC=CC1C[CH]O1'),
    E0 = (203.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.649518,0.0525061,3.36808e-05,-9.77319e-08,4.74779e-11,24574.1,27.1103], Tmin=(100,'K'), Tmax=(893.49,'K')), NASAPolynomial(coeffs=[21.8988,0.0150062,-1.12352e-07,-3.29749e-10,2.63313e-14,18476.5,-85.8847], Tmin=(893.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Oxetane) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C(C=O)C=CC=C(1022)',
    structure = SMILES('[CH2]C(C=O)C=CC=C'),
    E0 = (132.286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.607137,'amu*angstrom^2'), symmetry=1, barrier=(13.9593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.607146,'amu*angstrom^2'), symmetry=1, barrier=(13.9595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60716,'amu*angstrom^2'), symmetry=1, barrier=(13.9598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.607146,'amu*angstrom^2'), symmetry=1, barrier=(13.9595,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.374593,0.0796253,-7.32954e-05,3.66513e-08,-7.56143e-12,16041,29.5512], Tmin=(100,'K'), Tmax=(1148.92,'K')), NASAPolynomial(coeffs=[13.2116,0.0349325,-1.49452e-05,2.79303e-09,-1.93966e-13,13091.2,-34.162], Tmin=(1148.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O)"""),
)

species(
    label = 'CH2CHCHCH(170)',
    structure = SMILES('[CH]=CC=C'),
    E0 = (346.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.31937,'amu*angstrom^2'), symmetry=1, barrier=(30.3349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64255,0.0163337,3.86225e-05,-6.71377e-08,2.83603e-11,41729.6,13.282], Tmin=(100,'K'), Tmax=(937.724,'K')), NASAPolynomial(coeffs=[12.9705,0.00669127,-1.0007e-06,1.67602e-10,-1.71436e-14,38279.7,-43.9476], Tmin=(937.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH2CHCHCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]CC=O(1023)',
    structure = SMILES('[CH]CC=O'),
    E0 = (243.878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,216.297,1936.8,1938.25],'cm^-1')),
        HinderedRotor(inertia=(0.00364564,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254291,'amu*angstrom^2'), symmetry=1, barrier=(8.73847,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.86518,0.0272923,-2.85547e-05,2.34994e-08,-9.02537e-12,29370.4,14.7057], Tmin=(100,'K'), Tmax=(711.957,'K')), NASAPolynomial(coeffs=[3.83606,0.0199904,-9.27842e-06,1.80503e-09,-1.27782e-13,29279,10.6805], Tmin=(711.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet)"""),
)

species(
    label = 'O=C1CC=CC=CC1(229)',
    structure = SMILES('O=C1CC=CC=CC1'),
    E0 = (-22.6359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85754,-0.0195366,0.000292265,-4.0156e-07,1.61897e-10,-2583.52,21.0166], Tmin=(100,'K'), Tmax=(925.678,'K')), NASAPolynomial(coeffs=[42.0753,-0.0216079,1.73679e-05,-3.18444e-09,1.85544e-13,-17386.2,-209.643], Tmin=(925.678,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.6359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cycloheptane)"""),
)

species(
    label = 'O[C]1CC=CC=CC1(957)',
    structure = SMILES('O[C]1CC=CC=CC1'),
    E0 = (77.0078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.298958,0.0831094,-7.2786e-05,3.25107e-08,-5.74692e-12,9426.25,13.4096], Tmin=(100,'K'), Tmax=(1368.93,'K')), NASAPolynomial(coeffs=[19.4959,0.025269,-9.40753e-06,1.64549e-09,-1.10181e-13,4006.69,-88.3054], Tmin=(1368.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.0078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C2CsJOH)"""),
)

species(
    label = 'OC1[CH]C=CC=CC1(958)',
    structure = SMILES('OC1C=C[CH]C=CC1'),
    E0 = (19.2665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8195,-0.0205219,0.000307529,-4.31225e-07,1.78374e-10,2460.02,27.104], Tmin=(100,'K'), Tmax=(904.714,'K')), NASAPolynomial(coeffs=[44.0968,-0.028869,2.52974e-05,-5.08494e-09,3.31997e-13,-12497.9,-213.016], Tmin=(904.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.2665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C)"""),
)

species(
    label = 'OC1C[C]=CC=CC1(959)',
    structure = SMILES('OC1C[C]=CC=CC1'),
    E0 = (138.222,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.000533947,0.0770197,-5.43466e-05,1.25355e-08,1.44743e-12,16777.7,13.3795], Tmin=(100,'K'), Tmax=(1048.49,'K')), NASAPolynomial(coeffs=[17.9647,0.0274399,-1.05379e-05,1.92535e-09,-1.34527e-13,11968.4,-79.1126], Tmin=(1048.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(138.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Cds_S)"""),
)

species(
    label = 'OC1CC=[C]C=CC1(960)',
    structure = SMILES('OC1CC=[C]C=CC1'),
    E0 = (99.3755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.043201,0.0789598,-5.87619e-05,1.64797e-08,4.08669e-13,12106.3,12.6542], Tmin=(100,'K'), Tmax=(1024.74,'K')), NASAPolynomial(coeffs=[17.4739,0.0282536,-1.04048e-05,1.84761e-09,-1.26935e-13,7588.45,-76.811], Tmin=(1024.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.3755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C)"""),
)

species(
    label = '[O][C]1CC=CC=CC1(230)',
    structure = SMILES('[O][C]1CC=CC=CC1'),
    E0 = (307.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00465397,0.0775611,-6.19516e-05,2.18255e-08,-2.04961e-12,37120.5,12.1214], Tmin=(100,'K'), Tmax=(1097.33,'K')), NASAPolynomial(coeffs=[18.0631,0.0256507,-1.0016e-05,1.83007e-09,-1.2715e-13,32319.4,-80.4955], Tmin=(1097.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1[CH]C=CC=CC1(952)',
    structure = SMILES('[O]C1C=C[CH]C=CC1'),
    E0 = (249.627,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97778,-0.0244582,0.000313212,-4.3593e-07,1.79818e-10,30160.7,26.3437], Tmin=(100,'K'), Tmax=(905.282,'K')), NASAPolynomial(coeffs=[44.0319,-0.0306622,2.58834e-05,-5.17246e-09,3.36985e-13,15186.6,-213.008], Tmin=(905.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1C[C]=CC=CC1(953)',
    structure = SMILES('[O]C1C[C]=CC=CC1'),
    E0 = (368.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.171546,0.0729177,-4.80703e-05,7.03871e-09,3.23952e-12,44477.9,12.57], Tmin=(100,'K'), Tmax=(1040.09,'K')), NASAPolynomial(coeffs=[17.8407,0.0257457,-1.00084e-05,1.85106e-09,-1.30629e-13,39678.3,-78.7713], Tmin=(1040.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1CC=[C]C=CC1(954)',
    structure = SMILES('[O]C1CC=[C]C=CC1'),
    E0 = (329.736,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.126375,0.0748854,-5.2572e-05,1.10786e-08,2.16772e-12,39806.6,11.8537], Tmin=(100,'K'), Tmax=(1018.77,'K')), NASAPolynomial(coeffs=[17.3721,0.0265232,-9.85499e-06,1.76862e-09,-1.22655e-13,35288.6,-76.5956], Tmin=(1018.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]1CC2CC=CC1O2(961)',
    structure = SMILES('[CH]1CC2CC=CC1O2'),
    E0 = (110.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01737,0.0183651,0.000104923,-1.42653e-07,5.33021e-11,13409.2,20.6713], Tmin=(100,'K'), Tmax=(985.73,'K')), NASAPolynomial(coeffs=[14.5163,0.0315587,-1.24114e-05,2.47897e-09,-1.88531e-13,7839.99,-55.1997], Tmin=(985.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_ene_1) + radical(CCJCO)"""),
)

species(
    label = '[CH]1C=CCC2CC1O2(962)',
    structure = SMILES('[CH]1C=CCC2CC1O2'),
    E0 = (-44.1223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.363393,0.0491554,0.000104695,-2.06771e-07,9.17334e-11,-5106.39,-3.85289], Tmin=(100,'K'), Tmax=(934.066,'K')), NASAPolynomial(coeffs=[41.8244,-0.0116606,9.89825e-06,-1.74874e-09,9.49267e-14,-18215.8,-232.493], Tmin=(934.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.1223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(C=CCJCO)"""),
)

species(
    label = '[O]C1C=CC=CCC1(219)',
    structure = SMILES('[O]C1C=CC=CCC1'),
    E0 = (130.404,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281652,0.0647965,-9.96522e-06,-3.82168e-08,2.02738e-11,15833,12.6088], Tmin=(100,'K'), Tmax=(985.759,'K')), NASAPolynomial(coeffs=[19.8086,0.024911,-9.15138e-06,1.72864e-09,-1.27033e-13,10071.4,-91.0149], Tmin=(985.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1CC2C=CC2C1(955)',
    structure = SMILES('[O]C1CC2C=CC2C1'),
    E0 = (211.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66306,0.0299414,7.18699e-05,-1.09685e-07,4.2294e-11,25554.9,19.7883], Tmin=(100,'K'), Tmax=(984.223,'K')), NASAPolynomial(coeffs=[14.4163,0.0315574,-1.20477e-05,2.33106e-09,-1.73408e-13,20455.8,-54.6864], Tmin=(984.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + radical(CC(C)OJ)"""),
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
    E0 = (322.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (327.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (239.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (183.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (425.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (290.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (372.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (490.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (530.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (455.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (451.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (490.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (412.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (505.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (527.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (162.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (236.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (133.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (248.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (314.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (225.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (193.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (141.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (163.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (266.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (199.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (153.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (378.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (400.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (162.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (586.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (151.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (348.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (348.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (340.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (345.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (183.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (175.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (285.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (118.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (223.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (236.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (155.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (324.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (213.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (297.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (231.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (208.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (203.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (255.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (590.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (260.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (296.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (272.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (239.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (242.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (519.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (467.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (580.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (545.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (261.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (157.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (246.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (253.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (510.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction17',
    reactants = ['H(19)', 'C=C=CC=CCC=O(976)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(19)', 'C=CC=C=CCC=O(1001)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction55',
    reactants = ['H(19)', 'C=CC=CC=CC=O(1013)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.165e+08,'cm^3/(mol*s)'), n=1.533, Ea=(5.18398,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 2853 used for Cds-COH_Cds;HJ
Exact match found for rate rule [Cds-COH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction56',
    reactants = ['HCO(30)', 'C=CC=CC=C(1014)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0534234,'m^3/(mol*s)'), n=2.459, Ea=(4.91982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-CdH;CJ] for rate rule [Cds-HH_Cds-CdH;CO_pri_rad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2CHO(51)', '[CH]=CC=C[CH2](324)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;C_pri_rad] for rate rule [Cd_pri_rad;C_rad/H2/CO]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HCO(30)', '[CH2]C=CC=C[CH2](983)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.62e+13,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [C_pri_rad;CO_pri_rad] for rate rule [C_rad/H2/Cd;CO_pri_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(19)', '[CH2]C=CC=C[CH]C=O(984)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(78817,'m^3/(mol*s)'), n=0.288419, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H/TwoDe;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(19)', '[CH2]C=CC=[C]CC=O(985)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C[CH2](321)', '[CH]=CCC=O(986)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using an average for rate rule [Cd_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(19)', '[CH2]C=C[C]=CCC=O(987)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(19)', '[CH2]C=[C]C=CCC=O(988)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(19)', '[CH2][C]=CC=CCC=O(989)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(19)', '[CH2]C=CC=CC[C]=O(990)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_rad/NonDe;Y_rad] for rate rule [CO_rad/NonDe;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(19)', '[CH]C=CC=CCC=O(999)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C2H3(44)', '[CH]C=CCC=O(1010)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using an average for rate rule [Cd_pri_rad;Cd_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]CC=CC=CC=O(994)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.25e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;CdHCO_2] for rate rule [1_3_pentadiene;CH(CJ)_1;CdHCO_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['O=CC[CH]C1C=CC1(975)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra_HNd_pri;radadd_intra_cs] for rate rule [R5_SD_D;doublebond_intra_HNd_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['C7H9O(9)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.19e+11,'s^-1'), n=0.08, Ea=(92.3387,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8;multiplebond_intra;radadd_intra_cs2H] for rate rule [R8;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['C[C]=CC=CCC=O(977)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC=[C]C=CCC=O(978)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 288 used for R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CC=C[C]=CCC=O(979)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.33e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 284 used for R4H_SDS;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC=CC=[C]CC=O(980)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSMS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['CC=CC=C[CH]C=O(981)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(500000,'s^-1'), n=1.95, Ea=(100.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SMSMS;C_rad_out_2H;Cs_H_out_H/OneDe] for rate rule [R6H_SMSMS;C_rad_out_2H;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC=CC=CC[C]=O(982)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(244756,'s^-1'), n=1.235, Ea=(80.8557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;Y_rad_out;Cs_H_out_2H] for rate rule [R7H;CO_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['O=CCC=CC1[CH]C1(991)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 142 used for R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['O=CCC1[CH]C=CC1(992)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.09944e+09,'s^-1'), n=0.4512, Ea=(158.872,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 95 used for R5_SD_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H
Exact match found for rate rule [R5_SD_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['[CH]1CC=CC=CCO1(993)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.18057e+11,'s^-1'), n=0.420859, Ea=(112.437,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;multiplebond_intra;radadd_intra_cs2H] + [R6plus;carbonyl_intra_H;radadd_intra] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic
Ea raised from 109.3 to 112.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction37',
    reactants = ['CO(28)', '[CH2]C=CC=CC(995)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO;C_pri] for rate rule [CO;C_pri/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C=CC=COC=C(996)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C=CC=CC=CO(997)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CH2(T)(33)', '[CH]=CC=CCC=O(998)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['C=CC1C=CCC1[O](1000)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(110.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra_cs] for rate rule [R6_SMS_CO;carbonylbond_intra_H;radadd_intra_csHCd]
Euclidian distance = 3.0
family: Intra_R_Add_Exocyclic
Ea raised from 110.4 to 110.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=CC[C]=CCC=O(1002)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=[C]CC=CCC=O(1003)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=CCC=[C]CC=O(1004)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=CCC=CCC=O(1005)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=CCC=C[CH]C=O(1006)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.86616e+06,'s^-1'), n=1.7775, Ea=(110.667,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H;C_rad_out_H/OneDe;Cs_H_out_H/Cd] + [R4H_SDS;C_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_SDS;C_rad_out_H/OneDe;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=CCC=CC[C]=O(1007)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(60205.5,'s^-1'), n=1.86417, Ea=(61.9987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5H_SSMS;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['C=CC1[CH]C1CC=O(1008)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;doublebond_intra_pri_HNd_Cs;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['C=CC1C=CC[CH]O1(1009)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.00901e+09,'s^-1'), n=0.463766, Ea=(77.4753,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHCd] for rate rule [R6_SMS_CO;carbonyl_intra_H;radadd_intra_csHCd]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic
Ea raised from 77.1 to 77.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['C=CC=CC1CC1[O](1011)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(183.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCd] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 183.2 to 183.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['[CH2]C1C=CC1CC=O(1012)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra_csHNd] for rate rule [R5_SD_D;doublebond_intra_2H_pri;radadd_intra_csHNd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction57',
    reactants = ['C=CC=CC[CH]C=O(1015)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4.02296e+08,'s^-1'), n=1.43567, Ea=(82.7609,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/OneDe;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['C=CC=[C]CCC=O(1016)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['C=CC=CCC[C]=O(1017)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;CO_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['C=C[C]=CCCC=O(1018)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(6.1583e+09,'s^-1'), n=0.92705, Ea=(170.178,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] + [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['C=[C]C=CCCC=O(1019)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(2.22e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_Cd;Cs_H_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH]=CC=CCCC=O(1020)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R5H_DSMS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2]C=CC=CCC=O(956)'],
    products = ['C=CC=CC1C[CH]O1(1021)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(6.89861e+07,'s^-1'), n=1.13751, Ea=(163.006,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_csH(CdCdCd)]
Euclidian distance = 3.74165738677
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2]C(C=O)C=CC=C(1022)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO] for rate rule [cCs(-HC)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction65',
    reactants = ['CH2CHCHCH(170)', '[CH]CC=O(1023)'],
    products = ['[CH2]C=CC=CCC=O(956)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction1',
    reactants = ['H(19)', 'O=C1CC=CC=CC1(229)'],
    products = ['C7H9O(9)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(0.0366254,'m^3/(mol*s)'), n=1.743, Ea=(71.4418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsCs_O;YJ] for rate rule [CO-CsCs_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C7H9O(9)'],
    products = ['O[C]1CC=CC=CC1(957)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C7H9O(9)'],
    products = ['OC1[CH]C=CC=CC1(958)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(4.15038e+07,'s^-1'), n=1.60667, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OC1C[C]=CC=CC1(959)'],
    products = ['C7H9O(9)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C7H9O(9)'],
    products = ['OC1CC=[C]C=CC1(960)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5H_SSSD;O_rad_out;Cd_H_out_single] for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(19)', '[O][C]1CC=CC=CC1(230)'],
    products = ['C7H9O(9)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', '[O]C1[CH]C=CC=CC1(952)'],
    products = ['C7H9O(9)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(2.71464e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', '[O]C1C[C]=CC=CC1(953)'],
    products = ['C7H9O(9)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(19)', '[O]C1CC=[C]C=CC1(954)'],
    products = ['C7H9O(9)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C7H9O(9)'],
    products = ['[CH]1CC2CC=CC1O2(961)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(7.92069e+12,'s^-1'), n=-0.233729, Ea=(130.948,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn1cx_gamma;doublebond_intra_pri_HCd;radadd_intra] for rate rule [Rn1c7_gamma_short;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C7H9O(9)'],
    products = ['[CH]1C=CCC2CC1O2(962)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(220000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn1c7_beta_long;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C7H9O(9)'],
    products = ['[O]C1C=CC=CCC1(219)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(4.05241e+09,'s^-1'), n=1.0765, Ea=(115.616,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_pentadiene;CH_end;CdHC_2] + [1_3_pentadiene;CH2(C)_1;unsaturated_end] for rate rule [1_3_pentadiene;CH2(C)_1;CdHC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C7H9O(9)'],
    products = ['[O]C1CC2C=CC2C1(955)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C7H9(3)', 'O(20)'],
    products = ['C7H9O(9)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '174',
    isomers = [
        'C7H9O(9)',
        '[CH2]C=CC=CCC=O(956)',
    ],
    reactants = [
        ('C7H9(3)', 'O(20)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '174',
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

