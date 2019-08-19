species(
    label = '[CH2]OC1[CH]CC=CC=C1(1357)',
    structure = SMILES('[CH2]OC1[CH]CC=CC=C1'),
    E0 = (304.276,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.60699,0.0853836,-4.6141e-05,-7.84275e-09,1.06243e-11,36775.9,19.2417], Tmin=(100,'K'), Tmax=(1004.75,'K')), NASAPolynomial(coeffs=[22.3626,0.027721,-1.04884e-05,1.96358e-09,-1.41649e-13,30455,-100.168], Tmin=(1004.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CCJCO) + radical(CsJOCC2)"""),
)

species(
    label = 'C7H8(11)',
    structure = SMILES('C1=CC=CCC=C1'),
    E0 = (165.503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3872.77,'J/mol'), sigma=(6.33517,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=604.92 K, Pc=34.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.976321,0.0470646,2.43918e-05,-7.34711e-08,3.37246e-11,20032,17.1406], Tmin=(100,'K'), Tmax=(948.557,'K')), NASAPolynomial(coeffs=[19.6256,0.016921,-4.63468e-06,8.3158e-10,-6.48659e-14,14312.2,-83.3476], Tmin=(948.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
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
    label = '[CH2]OC1C=CC=CCC=1(1423)',
    structure = SMILES('[CH2]OC1C=CC=CCC=1'),
    E0 = (162.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.156,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.519087,0.0800559,-2.73497e-05,-3.51889e-08,2.28222e-11,19700.7,26.5641], Tmin=(100,'K'), Tmax=(955.208,'K')), NASAPolynomial(coeffs=[25.1052,0.0195369,-5.7819e-06,1.03322e-09,-7.78318e-14,12671,-107.056], Tmin=(955.208,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=COCJ)"""),
)

species(
    label = '[CH2]OC1C=CC=CC=C1(1412)',
    structure = SMILES('[CH2]OC1C=CC=CC=C1'),
    E0 = (192.001,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.156,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.647238,0.0859127,-4.91154e-05,-7.66558e-09,1.1506e-11,23274.2,26.8], Tmin=(100,'K'), Tmax=(988.1,'K')), NASAPolynomial(coeffs=[23.6724,0.0232861,-8.42668e-06,1.57352e-09,-1.14938e-13,16719.3,-99.0864], Tmin=(988.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(CsJOCC2)"""),
)

species(
    label = '[CH]1C=C[CH]CC=C1(125)',
    structure = SMILES('[CH]1C=C[CH]CC=C1'),
    E0 = (333.186,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.70798,-0.040239,0.000339088,-4.54205e-07,1.84685e-10,40184,19.5412], Tmin=(100,'K'), Tmax=(904.636,'K')), NASAPolynomial(coeffs=[40.627,-0.029508,2.54911e-05,-5.11056e-09,3.33198e-13,26023.7,-199.942], Tmin=(904.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH2]O[C]1C=CC=CCC1(1414)',
    structure = SMILES('[CH2]OC1=C[CH]C=CCC1'),
    E0 = (188.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0887,-0.00641741,0.000295392,-4.28613e-07,1.79247e-10,22806.9,31.3308], Tmin=(100,'K'), Tmax=(905.468,'K')), NASAPolynomial(coeffs=[47.4356,-0.0271498,2.49066e-05,-5.02579e-09,3.27602e-13,6870.56,-229.318], Tmin=(905.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(C=COCJ)"""),
)

species(
    label = '[CH2]OC1[CH]C=CC=CC1(1358)',
    structure = SMILES('[CH2]OC1C=C[CH]C=CC1'),
    E0 = (223.597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4098.12,'J/mol'), sigma=(6.95722,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=640.12 K, Pc=27.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.932879,0.000905919,0.00026984,-3.97719e-07,1.66909e-10,27065.4,31.0403], Tmin=(100,'K'), Tmax=(906.969,'K')), NASAPolynomial(coeffs=[45.7935,-0.0231727,2.22712e-05,-4.49746e-09,2.9164e-13,11780.8,-220.407], Tmin=(906.969,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(CsJOCC2)"""),
)

species(
    label = 'CO[C]1[CH]CC=CC=C1(1424)',
    structure = SMILES('COC1[CH]C=C[CH]CC=1'),
    E0 = (152.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.641986,0.00671005,0.000260894,-3.91449e-07,1.64974e-10,18525.7,26.87], Tmin=(100,'K'), Tmax=(909.19,'K')), NASAPolynomial(coeffs=[47.3631,-0.023678,2.20422e-05,-4.41007e-09,2.83684e-13,2790.3,-233.898], Tmin=(909.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]OC1[C]=CC=CCC1(1416)',
    structure = SMILES('[CH2]OC1[C]=CC=CCC1'),
    E0 = (342.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0511,0.0975468,-8.6755e-05,3.86948e-08,-6.79828e-12,41352.3,18.8254], Tmin=(100,'K'), Tmax=(1379.24,'K')), NASAPolynomial(coeffs=[23.083,0.0275539,-1.06332e-05,1.90043e-09,-1.28918e-13,34695,-105.367], Tmin=(1379.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CsJOCC2) + radical(Cds_S)"""),
)

species(
    label = '[CH2]OC1C=CC=[C]CC1(1413)',
    structure = SMILES('[CH2]OC1C=CC=[C]CC1'),
    E0 = (342.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0511,0.0975468,-8.6755e-05,3.86948e-08,-6.79828e-12,41352.3,18.8254], Tmin=(100,'K'), Tmax=(1379.24,'K')), NASAPolynomial(coeffs=[23.083,0.0275539,-1.06332e-05,1.90043e-09,-1.28918e-13,34695,-105.367], Tmin=(1379.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Cds_S) + radical(CsJOCC2)"""),
)

species(
    label = 'COC1[C]=CC=CC[CH]1(1425)',
    structure = SMILES('COC1[C]=CC=CC[CH]1'),
    E0 = (358.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.395939,0.0840288,-5.66488e-05,1.30701e-08,8.96762e-13,43290.6,19.5466], Tmin=(100,'K'), Tmax=(1112.68,'K')), NASAPolynomial(coeffs=[19.1043,0.0329278,-1.33751e-05,2.4901e-09,-1.74432e-13,37774.8,-81.8985], Tmin=(1112.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Cds_S) + radical(CCJCO)"""),
)

species(
    label = '[CH2]OC1C=[C]C=CCC1(1417)',
    structure = SMILES('[CH2]OC1C=[C]C=CCC1'),
    E0 = (303.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06343,0.0991899,-9.04042e-05,4.20098e-08,-7.71293e-12,36679.5,17.9871], Tmin=(100,'K'), Tmax=(1321.59,'K')), NASAPolynomial(coeffs=[22.0728,0.0291647,-1.09259e-05,1.9176e-09,-1.28855e-13,30564.1,-100.083], Tmin=(1321.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CsJOCC2) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]OC1C=C[C]=CCC1(1415)',
    structure = SMILES('[CH2]OC1C=C[C]=CCC1'),
    E0 = (303.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06343,0.0991899,-9.04042e-05,4.20098e-08,-7.71292e-12,36679.5,17.9871], Tmin=(100,'K'), Tmax=(1321.59,'K')), NASAPolynomial(coeffs=[22.0728,0.0291647,-1.09259e-05,1.91759e-09,-1.28855e-13,30564.1,-100.083], Tmin=(1321.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(CsJOCC2)"""),
)

species(
    label = 'COC1[CH]C=C[CH]C=C1(1388)',
    structure = SMILES('COC1[CH]C=C[CH]C=C1'),
    E0 = (110.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49885,-0.0198966,0.000335849,-4.6769e-07,1.91691e-10,13408.8,30.3893], Tmin=(100,'K'), Tmax=(910.909,'K')), NASAPolynomial(coeffs=[46.8329,-0.0245376,2.33211e-05,-4.63703e-09,2.95023e-13,-2916.72,-228.368], Tmin=(910.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC(O)C=C) + radical(C=CCJC=C)"""),
)

species(
    label = 'COC1[CH]CC=C[C]=C1(1426)',
    structure = SMILES('COC1[CH]CC=C[C]=C1'),
    E0 = (319.687,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.457197,0.086196,-6.18928e-05,1.81321e-08,-6.34895e-13,38620,18.8874], Tmin=(100,'K'), Tmax=(1104.52,'K')), NASAPolynomial(coeffs=[18.6432,0.0336855,-1.32076e-05,2.40385e-09,-1.66114e-13,33384.3,-79.7607], Tmin=(1104.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CCJCO) + radical(C=CJC=C)"""),
)

species(
    label = 'COC1[CH]C[C]=CC=C1(1427)',
    structure = SMILES('COC1[CH]C[C]=CC=C1'),
    E0 = (358.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.395939,0.0840288,-5.66488e-05,1.30701e-08,8.96762e-13,43290.6,19.5466], Tmin=(100,'K'), Tmax=(1112.68,'K')), NASAPolynomial(coeffs=[19.1043,0.0329278,-1.33751e-05,2.4901e-09,-1.74432e-13,37774.8,-81.8985], Tmin=(1112.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Cds_S) + radical(CCJCO)"""),
)

species(
    label = 'COC1[CH]CC=[C]C=C1(1428)',
    structure = SMILES('COC1[CH]CC=[C]C=C1'),
    E0 = (319.687,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.457197,0.086196,-6.18928e-05,1.81321e-08,-6.34895e-13,38620,18.8874], Tmin=(100,'K'), Tmax=(1104.52,'K')), NASAPolynomial(coeffs=[18.6432,0.0336855,-1.32076e-05,2.40385e-09,-1.66114e-13,33384.3,-79.7607], Tmin=(1104.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(CCJCO)"""),
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
    label = '[CH2]OC1[CH]C2C=CCC21(1429)',
    structure = SMILES('[CH2]OC1[CH]C2C=CCC21'),
    E0 = (356.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.764083,0.0503289,3.80387e-05,-8.24144e-08,3.38225e-11,43021,27.2963], Tmin=(100,'K'), Tmax=(996.348,'K')), NASAPolynomial(coeffs=[17.2396,0.0343631,-1.34678e-05,2.59562e-09,-1.90894e-13,37247.3,-64.6272], Tmin=(996.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(CCJCO) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]OC1C=CC2[CH]CC21(1430)',
    structure = SMILES('[CH2]OC1C=CC2[CH]CC21'),
    E0 = (344.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.774017,0.0529779,2.35504e-05,-6.19778e-08,2.51348e-11,41526.4,28.0234], Tmin=(100,'K'), Tmax=(1031.38,'K')), NASAPolynomial(coeffs=[15.4766,0.0380335,-1.591e-05,3.0843e-09,-2.23934e-13,36255.7,-54.2116], Tmin=(1031.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(cyclobutane) + radical(CsJOCC2)"""),
)

species(
    label = '[CH]1CC=CC2[CH]C1OC2(1431)',
    structure = SMILES('[CH]1CC=CC2[CH]C1OC2'),
    E0 = (212.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.297165,0.0792006,-3.46384e-05,-1.20648e-08,9.9827e-12,25778.3,2.34717], Tmin=(100,'K'), Tmax=(1047.26,'K')), NASAPolynomial(coeffs=[19.5245,0.0341015,-1.38846e-05,2.63271e-09,-1.88196e-13,19948.1,-102.21], Tmin=(1047.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Tetrahydrofuran) + ring(Cycloheptene) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(CCJCO)"""),
)

species(
    label = '[CH2]OC1C=C[CH]C2CC21(1432)',
    structure = SMILES('[CH2]OC1[CH]C=CC2CC21'),
    E0 = (254.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.246402,0.0626308,1.30371e-05,-6.10882e-08,2.71385e-11,30790.2,23.42], Tmin=(100,'K'), Tmax=(1004.03,'K')), NASAPolynomial(coeffs=[19.1033,0.0339478,-1.34942e-05,2.59777e-09,-1.90172e-13,24662.8,-79.287], Tmin=(1004.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + radical(C=CCJCO) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]OC1C2[CH]C=CCC21(1433)',
    structure = SMILES('[CH2]OC1C2C=C[CH]CC21'),
    E0 = (276.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.668853,0.0533831,2.99609e-05,-7.36704e-08,3.05921e-11,33383.1,23.5979], Tmin=(100,'K'), Tmax=(1002.5,'K')), NASAPolynomial(coeffs=[17.0183,0.0354786,-1.40671e-05,2.7022e-09,-1.97439e-13,27726.7,-67.1817], Tmin=(1002.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + radical(cyclohexene-allyl) + radical(CsJOCC2)"""),
)

species(
    label = '[CH]1C[CH]C2C=CC1CO2(1434)',
    structure = SMILES('[CH]1C[CH]C2C=CC1CO2'),
    E0 = (213.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.538163,0.0896152,-6.82215e-05,2.46126e-08,-3.51225e-12,25802.9,3.22084], Tmin=(100,'K'), Tmax=(1647.82,'K')), NASAPolynomial(coeffs=[23.9763,0.0301079,-1.40526e-05,2.69736e-09,-1.87368e-13,17723.8,-127.291], Tmin=(1647.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(oxepane) + ring(36dihydro2hpyran) + radical(Cs_S) + radical(CCJCO)"""),
)

species(
    label = '[CH]1C=CC2[CH]CC1CO2(1435)',
    structure = SMILES('[CH]1C=CC2C[CH]C1OC2'),
    E0 = (136.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.843968,0.0586371,9.81535e-05,-2.05418e-07,9.2033e-11,16656.2,2.71169], Tmin=(100,'K'), Tmax=(935.063,'K')), NASAPolynomial(coeffs=[43.3399,-0.00794112,8.55619e-06,-1.51148e-09,7.86095e-14,3040.9,-236.104], Tmin=(935.063,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(Cycloheptane) + ring(Oxane) + radical(C=CCJCO) + radical(CCJCO)"""),
)

species(
    label = '[CH]1C=CC[CH]C2OCC12(1436)',
    structure = SMILES('[CH]1C=CC[CH]C2OCC12'),
    E0 = (153.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.295196,0.0794532,-3.55599e-05,-9.48913e-09,8.55282e-12,18578.5,-1.29534], Tmin=(100,'K'), Tmax=(1070.63,'K')), NASAPolynomial(coeffs=[19.0779,0.0359134,-1.49651e-05,2.84714e-09,-2.02905e-13,12777.3,-103.801], Tmin=(1070.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(Cycloheptene) + ring(Oxetane) + radical(CCJCO) + radical(Allyl_S)"""),
)

species(
    label = '[CH]=CC=CC(C=C)O[CH2](1437)',
    structure = SMILES('[CH]=CC=CC(C=C)O[CH2]'),
    E0 = (414.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3120,650,792.5,1650,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.854801,0.0959475,-8.70481e-05,3.9902e-08,-7.23409e-12,49990.4,37.1268], Tmin=(100,'K'), Tmax=(1334.27,'K')), NASAPolynomial(coeffs=[21.6305,0.0285378,-1.12641e-05,2.03595e-09,-1.39069e-13,43990.2,-77.8359], Tmin=(1334.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(414.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CsJOCC2) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC=CCC=CO[CH2](1438)',
    structure = SMILES('[CH]=CC=CCC=CO[CH2]'),
    E0 = (389.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.606225,0.0819032,-2.81202e-05,-3.61668e-08,2.36854e-11,47031.6,36.3739], Tmin=(100,'K'), Tmax=(948.132,'K')), NASAPolynomial(coeffs=[25.3428,0.0200462,-5.59234e-06,9.62881e-10,-7.17547e-14,39970.7,-98.7201], Tmin=(948.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COCJ) + radical(Cds_P)"""),
)

species(
    label = 'COC1C=CC=CCC=1(1439)',
    structure = SMILES('COC1C=CC=CCC=1'),
    E0 = (-29.7054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.637492,0.0821783,-2.57634e-05,-3.77776e-08,2.37279e-11,-3387.8,24.5036], Tmin=(100,'K'), Tmax=(957.016,'K')), NASAPolynomial(coeffs=[25.1168,0.0223226,-6.8497e-06,1.22471e-09,-9.09333e-14,-10505.7,-110.048], Tmin=(957.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.7054,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'COC1C=CC=CC=C1(1419)',
    structure = SMILES('COC1C=CC=CC=C1'),
    E0 = (8.41607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.345171,0.0793504,-3.31025e-05,-1.76431e-08,1.31396e-11,1182.9,26.2999], Tmin=(100,'K'), Tmax=(1011.42,'K')), NASAPolynomial(coeffs=[20.8685,0.0301871,-1.1702e-05,2.20459e-09,-1.58893e-13,-4884.84,-85.0671], Tmin=(1011.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.41607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = '[CH2]OC1[CH]CCC=CC=1(1440)',
    structure = SMILES('[CH2]OC1[CH]C=CCCC=1'),
    E0 = (203.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.558087,0.00985901,0.000250203,-3.80201e-07,1.61035e-10,24648.7,30.7633], Tmin=(100,'K'), Tmax=(908.489,'K')), NASAPolynomial(coeffs=[47.3036,-0.0239698,2.20889e-05,-4.42549e-09,2.85625e-13,9057.65,-229.333], Tmin=(908.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=COCJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]O[CH]C1C=CC=CC1(1411)',
    structure = SMILES('[CH2]O[CH]C1C=CC=CC1'),
    E0 = (285.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.238873,0.0641657,4.87151e-06,-5.5911e-08,2.70154e-11,34521.8,29.9445], Tmin=(100,'K'), Tmax=(968.622,'K')), NASAPolynomial(coeffs=[19.4473,0.0293461,-1.01245e-05,1.84346e-09,-1.33331e-13,28713,-72.8895], Tmin=(968.622,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(CsJOCC) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]OC1C=CC=CC1[CH2](1421)',
    structure = SMILES('[CH2]OC1C=CC=CC1[CH2]'),
    E0 = (293.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.226445,0.0633952,1.01583e-05,-6.3097e-08,2.99696e-11,35420.6,29.9238], Tmin=(100,'K'), Tmax=(963.571,'K')), NASAPolynomial(coeffs=[20.0501,0.0286016,-9.61988e-06,1.74513e-09,-1.26938e-13,29395.2,-76.4203], Tmin=(963.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Isobutyl) + radical(CsJOCC2)"""),
)

species(
    label = 'C1C=CC2OCC2CC=1(1362)',
    structure = SMILES('C1C=CC2OCC2CC=1'),
    E0 = (-92.4314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31041,0.14906,-0.000297552,3.18134e-07,-1.25222e-10,-10957,-17.4167], Tmin=(100,'K'), Tmax=(841.396,'K')), NASAPolynomial(coeffs=[-4.66842,0.0868291,-4.72083e-05,9.3251e-09,-6.49218e-13,-7624.01,14.652], Tmin=(841.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-92.4314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(1,3-Cycloheptadiene) + ring(Oxetane)"""),
)

species(
    label = '[CH2]OC1[CH]CC2C=CC21(1441)',
    structure = SMILES('[CH2]OC1[CH]CC2C=CC21'),
    E0 = (385.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.784065,0.0504148,3.60892e-05,-7.98136e-08,3.28529e-11,46497.3,26.3866], Tmin=(100,'K'), Tmax=(995.959,'K')), NASAPolynomial(coeffs=[16.9211,0.03445,-1.3432e-05,2.57711e-09,-1.88942e-13,40860.4,-63.5621], Tmin=(995.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + radical(CCJCO) + radical(CsJOCC2)"""),
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
    label = '[O]C1[CH]CC=CC=C1(1442)',
    structure = SMILES('[O]C1[CH]CC=CC=C1'),
    E0 = (330.307,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.480261,0.0595434,-1.25097e-06,-4.77339e-08,2.40858e-11,39869.5,14.392], Tmin=(100,'K'), Tmax=(972.786,'K')), NASAPolynomial(coeffs=[20.2652,0.020795,-7.19811e-06,1.3641e-09,-1.02536e-13,34004.3,-90.8745], Tmin=(972.786,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(CCJCO)"""),
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
    E0 = (304.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (383.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (424.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (311.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (465.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (423.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (421.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (483.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (483.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (402.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (347.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (347.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (363.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (352.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (501.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (462.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (526.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (435.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (435.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (385.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (305.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (305.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (400.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (331.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (371.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (455.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (423.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (367.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (329.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (350.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (464.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (549.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (312.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (426.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (711.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['C7H8(11)', 'CH2O(35)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(19)', '[CH2]OC1C=CC=CCC=1(1423)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(87.5179,'m^3/(mol*s)'), n=1.66467, Ea=(9.4546,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;HJ] for rate rule [Cds-CdOs_Cds;HJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', '[CH2]OC1C=CC=CC=C1(1412)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.0545849,'m^3/(mol*s)'), n=2.81111, Ea=(21.1569,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 26 used for Cds-CdH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CsH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2O(35)', '[CH]1C=C[CH]CC=C1(125)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4660,'cm^3/(mol*s)'), n=3.17, Ea=(97.3791,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-HH;YJ] for rate rule [Od_CO-HH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH2]O[C]1C=CC=CCC1(1414)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.41e+08,'s^-1'), n=1.52, Ea=(161.502,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_(CdCdCd)]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]OC1[CH]C=CC=CC1(1358)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.022e+10,'s^-1'), n=1.34, Ea=(199.577,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for R2H_S;C_rad_out_H/(Cd-Cd-Cd);Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/(Cd-Cd-Cd);Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['CO[C]1[CH]CC=CC=C1(1424)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.30996,'s^-1'), n=3.5644, Ea=(117.059,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_2H;XH_out] for rate rule [R3H_SS_O;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]OC1[C]=CC=CCC1(1416)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 206 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]OC1C=CC=[C]CC1(1413)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 206 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['COC1[C]=CC=CC[CH]1(1425)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]OC1C=[C]C=CCC1(1417)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]OC1C=C[C]=CCC1(1415)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['COC1[CH]C=C[CH]C=C1(1388)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(62296.1,'s^-1'), n=1.86, Ea=(59.4128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_H/Cd] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['COC1[CH]CC=C[C]=C1(1426)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['COC1[CH]C[C]=CC=C1(1427)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['COC1[CH]CC=[C]C=C1(1428)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6H_SMSSR;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', '[CH2][O](354)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH2]OC1[CH]C2C=CCC21(1429)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 116 used for Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH2]OC1C=CC2[CH]CC21(1430)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 116 used for Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH]1CC=CC2[CH]C1OC2(1431)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.9e+10,'s^-1'), n=0.06, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri_HCd;radadd_intra_cs2H] for rate rule [Rn2c7_beta_short;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH2]OC1C=C[CH]C2CC21(1432)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(1.2552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH2]OC1C2[CH]C=CCC21(1433)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(1.2552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH]1C[CH]C2C=CC1CO2(1434)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.54269e+17,'s^-1'), n=-0.500262, Ea=(96.3509,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_cyclic;doublebond_intra_pri_HCd;radadd_intra_cs] for rate rule [Rn2c7_gamma_long;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH]1C=CC2[CH]CC1CO2(1435)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R7;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H] for rate rule [Rn2c7_gamma_long;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH]1C=CC[CH]C2OCC12(1436)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.54811e+07,'s^-1'), n=0.785425, Ea=(66.8554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=CC=CC(C=C)O[CH2](1437)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=CC=CCC=CO[CH2](1438)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.25555e+08,'s^-1'), n=0.67, Ea=(34.1205,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;doublebond_intra_pri;radadd_intra_cdsingleH] + [R7_linear;doublebond_intra_pri_HNd;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_HNd_O;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['COC1C=CC=CCC=1(1439)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['COC1C=CC=CC=C1(1419)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH2]OC1[CH]CCC=CC=1(1440)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.125e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;CdHC_2] for rate rule [1_3_pentadiene;CH(CJ)_1;CdHC_2]
Euclidian distance = 1.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH2]O[CH]C1C=CC=CC1(1411)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH2]OC1C=CC=CC1[CH2](1421)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['C1C=CC2OCC2CC=1(1362)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    products = ['[CH2]OC1[CH]CC2C=CC21(1441)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CH2(T)(33)', '[O]C1[CH]CC=CC=C1(1442)'],
    products = ['[CH2]OC1[CH]CC=CC=C1(1357)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

network(
    label = '247',
    isomers = [
        '[CH2]OC1[CH]CC=CC=C1(1357)',
    ],
    reactants = [
        ('C7H8(11)', 'CH2O(35)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '247',
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

