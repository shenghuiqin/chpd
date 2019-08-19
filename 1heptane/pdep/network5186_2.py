species(
    label = 'CC[CH]CC(C)=O(2056)',
    structure = SMILES('CC[CH]CC(C)=O'),
    E0 = (-105.852,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3638.59,'J/mol'), sigma=(6.31136,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=568.34 K, Pc=32.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1706,0.0623114,-3.6678e-05,1.07041e-08,-1.29542e-12,-12630.4,28.1137], Tmin=(100,'K'), Tmax=(1782.23,'K')), NASAPolynomial(coeffs=[12.3522,0.0372159,-1.55568e-05,2.80359e-09,-1.87198e-13,-16616.1,-32.2928], Tmin=(1782.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]C(=O)CCCC(28404)',
    structure = SMILES('C=C([O])CCCC'),
    E0 = (-140.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,318.496,318.506,318.53,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00166199,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206814,'amu*angstrom^2'), symmetry=1, barrier=(14.8896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2068,'amu*angstrom^2'), symmetry=1, barrier=(14.8897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206801,'amu*angstrom^2'), symmetry=1, barrier=(14.8894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3795.8,'J/mol'), sigma=(6.58551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.90 K, Pc=30.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493565,0.0673092,-3.64475e-05,1.63285e-09,3.70933e-12,-16804,28.7123], Tmin=(100,'K'), Tmax=(1055.96,'K')), NASAPolynomial(coeffs=[14.1679,0.0330108,-1.25861e-05,2.26305e-09,-1.5565e-13,-20667.6,-42.6231], Tmin=(1055.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-140.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC[CH]CC(140)',
    structure = SMILES('CC[CH]CC'),
    E0 = (26.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,1748.69,1749.84],'cm^-1')),
        HinderedRotor(inertia=(0.100684,'amu*angstrom^2'), symmetry=1, barrier=(6.89816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104397,'amu*angstrom^2'), symmetry=1, barrier=(6.91808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100753,'amu*angstrom^2'), symmetry=1, barrier=(6.91387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103404,'amu*angstrom^2'), symmetry=1, barrier=(6.85304,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1408,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3051.59,'J/mol'), sigma=(5.73385,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=476.65 K, Pc=36.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9332,0.0431369,-1.41036e-05,-6.73618e-10,8.04699e-13,3248.31,21.1643], Tmin=(100,'K'), Tmax=(1648.16,'K')), NASAPolynomial(coeffs=[9.63706,0.0333954,-1.33881e-05,2.33367e-09,-1.51518e-13,-507.453,-23.5416], Tmin=(1648.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC)"""),
)

species(
    label = 'CO(12)',
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
    label = 'CH3CO(57)',
    structure = SMILES('C[C]=O'),
    E0 = (-22.2282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,2865.56],'cm^-1')),
        HinderedRotor(inertia=(0.0188671,'amu*angstrom^2'), symmetry=1, barrier=(18.7749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03587,0.000877295,3.071e-05,-3.92476e-08,1.52969e-11,-2682.07,7.86177], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.31372,0.00917378,-3.32204e-06,5.39475e-10,-3.24524e-14,-3645.04,-1.67576], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-22.2282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH3CO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'butene1(127)',
    structure = SMILES('C=CCC'),
    E0 = (-16.4325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,252.555],'cm^-1')),
        HinderedRotor(inertia=(0.178654,'amu*angstrom^2'), symmetry=1, barrier=(7.72883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185589,'amu*angstrom^2'), symmetry=1, barrier=(7.72103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58773,0.0232778,1.93412e-05,-3.55496e-08,1.36906e-11,-1918.73,14.5751], Tmin=(100,'K'), Tmax=(1007.28,'K')), NASAPolynomial(coeffs=[7.20517,0.0236362,-9.0315e-06,1.65393e-09,-1.16019e-13,-3797.34,-12.4426], Tmin=(1007.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.4325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene1""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2CO(28)',
    structure = SMILES('C=C=O'),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'butyl_1(82)',
    structure = SMILES('[CH2]CCC'),
    E0 = (63.0573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0977402,'amu*angstrom^2'), symmetry=1, barrier=(2.24724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0976865,'amu*angstrom^2'), symmetry=1, barrier=(2.246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0977534,'amu*angstrom^2'), symmetry=1, barrier=(2.24754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.1143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25388,0.0316763,2.89994e-06,-1.98049e-08,8.20503e-12,7652.64,17.2725], Tmin=(100,'K'), Tmax=(1050.57,'K')), NASAPolynomial(coeffs=[7.59591,0.0260842,-1.01719e-05,1.85189e-09,-1.28169e-13,5716.37,-12.6366], Tmin=(1050.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.0573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), label="""butyl_1""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C(=O)CCC(2148)',
    structure = SMILES('C=C([O])CCC'),
    E0 = (-117.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,335.188,336.29,336.687],'cm^-1')),
        HinderedRotor(inertia=(0.00148026,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0981052,'amu*angstrom^2'), symmetry=1, barrier=(7.83672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0992182,'amu*angstrom^2'), symmetry=1, barrier=(7.84609,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3654.1,'J/mol'), sigma=(6.27192,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.76 K, Pc=33.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26061,0.0511177,-1.77487e-05,-1.26003e-08,8.41823e-12,-13971.9,23.7131], Tmin=(100,'K'), Tmax=(999.148,'K')), NASAPolynomial(coeffs=[12.5223,0.0258178,-9.4698e-06,1.69475e-09,-1.17565e-13,-17209.9,-35.5505], Tmin=(999.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC1CC1(C)[O](31202)',
    structure = SMILES('CCC1CC1(C)[O]'),
    E0 = (0.152271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803968,0.0537246,1.35123e-05,-5.6067e-08,2.53307e-11,148.194,24.3316], Tmin=(100,'K'), Tmax=(976.873,'K')), NASAPolynomial(coeffs=[16.2932,0.0299024,-1.07166e-05,1.96665e-09,-1.41427e-13,-4767.54,-59.7037], Tmin=(976.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.152271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ)"""),
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
    label = 'CC=CCC(C)=O(3026)',
    structure = SMILES('CC=CCC(C)=O'),
    E0 = (-184.035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,180,543.365],'cm^-1')),
        HinderedRotor(inertia=(0.4901,'amu*angstrom^2'), symmetry=1, barrier=(11.2684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.014344,'amu*angstrom^2'), symmetry=1, barrier=(2.99928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0521676,'amu*angstrom^2'), symmetry=1, barrier=(11.2673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0540453,'amu*angstrom^2'), symmetry=1, barrier=(11.2903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16124,0.0529964,-1.29258e-05,-1.14222e-08,5.2111e-12,-22024.5,25.2775], Tmin=(100,'K'), Tmax=(1255.61,'K')), NASAPolynomial(coeffs=[12.7424,0.0360882,-1.66024e-05,3.20671e-09,-2.25647e-13,-26508.3,-39.5045], Tmin=(1255.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = 'CCC=CC(C)=O(3025)',
    structure = SMILES('CCC=CC(C)=O'),
    E0 = (-185.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,231.762,231.763],'cm^-1')),
        HinderedRotor(inertia=(0.00313875,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138919,'amu*angstrom^2'), symmetry=1, barrier=(5.29335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138918,'amu*angstrom^2'), symmetry=1, barrier=(5.2935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138934,'amu*angstrom^2'), symmetry=1, barrier=(5.29348,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25698,0.0561983,-2.8305e-05,5.63056e-09,-2.76307e-13,-22167.9,26.0136], Tmin=(100,'K'), Tmax=(1790.29,'K')), NASAPolynomial(coeffs=[16.6092,0.0297385,-1.27054e-05,2.26802e-09,-1.48379e-13,-28921.5,-60.502], Tmin=(1790.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H)"""),
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
    label = 'C=CCC(C)=O(760)',
    structure = SMILES('C=CCC(C)=O'),
    E0 = (-148.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,554.838,554.981],'cm^-1')),
        HinderedRotor(inertia=(0.102125,'amu*angstrom^2'), symmetry=1, barrier=(2.34805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0647671,'amu*angstrom^2'), symmetry=1, barrier=(14.1475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0647668,'amu*angstrom^2'), symmetry=1, barrier=(14.1476,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84389,0.0371822,7.58146e-06,-2.8484e-08,1.08178e-11,-17715.3,21.0088], Tmin=(100,'K'), Tmax=(1132.83,'K')), NASAPolynomial(coeffs=[10.9325,0.0291216,-1.35654e-05,2.68706e-09,-1.93817e-13,-21316.4,-30.7776], Tmin=(1132.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C[CH]CCC(C)=O(31054)',
    structure = SMILES('C[CH]CCC(C)=O'),
    E0 = (-111.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,180,400,3200],'cm^-1')),
        HinderedRotor(inertia=(0.0388552,'amu*angstrom^2'), symmetry=1, barrier=(3.57343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388552,'amu*angstrom^2'), symmetry=1, barrier=(3.57343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388552,'amu*angstrom^2'), symmetry=1, barrier=(3.57343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388552,'amu*angstrom^2'), symmetry=1, barrier=(3.57343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388552,'amu*angstrom^2'), symmetry=1, barrier=(3.57343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3638.59,'J/mol'), sigma=(6.31136,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=568.34 K, Pc=32.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14769,0.0679392,-6.561e-05,5.13391e-08,-1.91728e-11,-13289.5,28.4504], Tmin=(100,'K'), Tmax=(744.413,'K')), NASAPolynomial(coeffs=[3.09061,0.05233,-2.3741e-05,4.51471e-09,-3.14736e-13,-13435.5,20.6124], Tmin=(744.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(RCCJC)"""),
)

species(
    label = 'CCC[CH]C(C)=O(31053)',
    structure = SMILES('CCCC=C(C)[O]'),
    E0 = (-153.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.489125,0.0680772,-4.26908e-05,1.07671e-08,-8.57781e-14,-18277.4,28.2935], Tmin=(100,'K'), Tmax=(1149.79,'K')), NASAPolynomial(coeffs=[13.7144,0.0336781,-1.29609e-05,2.31149e-09,-1.56833e-13,-22086.1,-40.6939], Tmin=(1149.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-153.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CCCC(C)=O(22489)',
    structure = SMILES('[CH2]CCCC(C)=O'),
    E0 = (-100.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32464,0.0645169,-4.26217e-05,1.57126e-08,-2.65687e-12,-11997,26.8605], Tmin=(100,'K'), Tmax=(1234.52,'K')), NASAPolynomial(coeffs=[6.8557,0.0465955,-2.08465e-05,3.95351e-09,-2.75556e-13,-13362.7,-0.988973], Tmin=(1234.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC(C)=O(2610)',
    structure = SMILES('[CH2][CH]CC(C)=O'),
    E0 = (123.174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,180,2281.1],'cm^-1')),
        HinderedRotor(inertia=(0.103588,'amu*angstrom^2'), symmetry=1, barrier=(2.38168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103003,'amu*angstrom^2'), symmetry=1, barrier=(2.36824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102727,'amu*angstrom^2'), symmetry=1, barrier=(2.36191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00532682,'amu*angstrom^2'), symmetry=1, barrier=(2.36323,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00318,0.0465797,-2.77939e-05,8.38015e-09,-1.07604e-12,14883.2,24.6626], Tmin=(100,'K'), Tmax=(1631.1,'K')), NASAPolynomial(coeffs=[8.44585,0.03078,-1.3264e-05,2.44141e-09,-1.65801e-13,12781.5,-9.57157], Tmin=(1631.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]CC(C)=O(3022)',
    structure = SMILES('C[CH][CH]CC(C)=O'),
    E0 = (88.5939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,785.244,837.921,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0970954,'amu*angstrom^2'), symmetry=1, barrier=(2.23242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0970954,'amu*angstrom^2'), symmetry=1, barrier=(2.23242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0970954,'amu*angstrom^2'), symmetry=1, barrier=(2.23242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0970954,'amu*angstrom^2'), symmetry=1, barrier=(2.23242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0970954,'amu*angstrom^2'), symmetry=1, barrier=(2.23242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34212,0.0458583,3.8762e-05,-1.7912e-07,1.63511e-10,10705.1,26.9239], Tmin=(100,'K'), Tmax=(408.519,'K')), NASAPolynomial(coeffs=[3.53262,0.0482834,-2.18477e-05,4.16793e-09,-2.91946e-13,10490.3,20.808], Tmin=(408.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.5939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]CC(130)',
    structure = SMILES('[CH2][CH]CC'),
    E0 = (255.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2031.24],'cm^-1')),
        HinderedRotor(inertia=(0.244974,'amu*angstrom^2'), symmetry=1, barrier=(5.63244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00192352,'amu*angstrom^2'), symmetry=1, barrier=(5.63177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244928,'amu*angstrom^2'), symmetry=1, barrier=(5.63137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98997,0.0287412,-9.51469e-06,4.19232e-10,1.90526e-13,30780.1,16.8971], Tmin=(100,'K'), Tmax=(2154.56,'K')), NASAPolynomial(coeffs=[12.4231,0.0182241,-7.06316e-06,1.16769e-09,-7.11818e-14,25091.4,-39.6212], Tmin=(2154.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'CC[CH][CH]C(C)=O(3021)',
    structure = SMILES('CC[CH]C=C(C)[O]'),
    E0 = (-11.9687,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,459.681,460.446,462.732],'cm^-1')),
        HinderedRotor(inertia=(0.000794365,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0582067,'amu*angstrom^2'), symmetry=1, barrier=(8.91453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0595382,'amu*angstrom^2'), symmetry=1, barrier=(8.93158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0597371,'amu*angstrom^2'), symmetry=1, barrier=(8.8913,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.747727,0.0621889,-3.16656e-05,-1.4599e-11,3.68442e-12,-1314.48,26.2558], Tmin=(100,'K'), Tmax=(1072.96,'K')), NASAPolynomial(coeffs=[13.2184,0.0320635,-1.24285e-05,2.2478e-09,-1.54801e-13,-4932.59,-39.1757], Tmin=(1072.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-11.9687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C[CH]CC(C)=O(3023)',
    structure = SMILES('[CH2]C[CH]CC(C)=O'),
    E0 = (99.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,210.595,408.613,3281.19],'cm^-1')),
        HinderedRotor(inertia=(0.04897,'amu*angstrom^2'), symmetry=1, barrier=(2.46151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.04897,'amu*angstrom^2'), symmetry=1, barrier=(2.46151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.04897,'amu*angstrom^2'), symmetry=1, barrier=(2.46151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.04897,'amu*angstrom^2'), symmetry=1, barrier=(2.46151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.04897,'amu*angstrom^2'), symmetry=1, barrier=(2.46151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40653,0.0607418,-3.94153e-05,1.36127e-08,-2.06236e-12,12044,29.0537], Tmin=(100,'K'), Tmax=(1403.11,'K')), NASAPolynomial(coeffs=[8.51067,0.040489,-1.77639e-05,3.32533e-09,-2.29383e-13,10050.5,-7.62579], Tmin=(1403.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(RCCJ)"""),
)

species(
    label = 'CC[CH]C[C]=O(2431)',
    structure = SMILES('CC[CH]C[C]=O'),
    E0 = (108.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,231.874,3966.35],'cm^-1')),
        HinderedRotor(inertia=(0.229635,'amu*angstrom^2'), symmetry=1, barrier=(8.76113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00313542,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.501267,'amu*angstrom^2'), symmetry=1, barrier=(19.1256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00313556,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77134,0.0497357,-3.44094e-05,1.30647e-08,-2.15083e-12,13186.5,25.1913], Tmin=(100,'K'), Tmax=(1345.56,'K')), NASAPolynomial(coeffs=[8.29233,0.0303504,-1.2799e-05,2.35767e-09,-1.61497e-13,11431.6,-8.20425], Tmin=(1345.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJCC=O)"""),
)

species(
    label = 'C=C([O])C[CH]CC(1126)',
    structure = SMILES('C=C([O])C[CH]CC'),
    E0 = (53.6223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,180,180,1286.34,2707.76],'cm^-1')),
        HinderedRotor(inertia=(0.42994,'amu*angstrom^2'), symmetry=1, barrier=(13.8863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.601606,'amu*angstrom^2'), symmetry=1, barrier=(13.8321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0117751,'amu*angstrom^2'), symmetry=1, barrier=(13.9572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00266568,'amu*angstrom^2'), symmetry=1, barrier=(13.8472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3795.8,'J/mol'), sigma=(6.58551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.90 K, Pc=30.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1185,0.0634505,-4.53405e-05,1.82342e-08,-3.18071e-12,6552.72,29.0296], Tmin=(100,'K'), Tmax=(1288.34,'K')), NASAPolynomial(coeffs=[9.34835,0.0378989,-1.55913e-05,2.84021e-09,-1.9357e-13,4432.14,-12.7598], Tmin=(1288.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.6223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJCC)"""),
)

species(
    label = 'CCC1C[C](C)O1(31203)',
    structure = SMILES('CCC1C[C](C)O1'),
    E0 = (-12.2422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09235,0.0559511,-7.37549e-06,-2.90303e-08,1.63049e-11,-1360.22,23.9005], Tmin=(100,'K'), Tmax=(894.307,'K')), NASAPolynomial(coeffs=[10.556,0.0359496,-1.12757e-05,1.79323e-09,-1.15506e-13,-3945.72,-25.6906], Tmin=(894.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.2422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOCs)"""),
)

species(
    label = 'C[CH]CC(C)=O(1920)',
    structure = SMILES('C[CH]CC(C)=O'),
    E0 = (-82.0721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,182.282,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00278205,'amu*angstrom^2'), symmetry=1, barrier=(31.5875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.253136,'amu*angstrom^2'), symmetry=1, barrier=(5.9765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33618,'amu*angstrom^2'), symmetry=1, barrier=(31.5887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00507255,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66914,0.049028,-2.71333e-05,7.19381e-09,-7.67636e-13,-9785.77,24.0978], Tmin=(100,'K'), Tmax=(2066.03,'K')), NASAPolynomial(coeffs=[13.509,0.0261052,-1.04907e-05,1.82356e-09,-1.17809e-13,-14678.1,-41.6139], Tmin=(2066.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-82.0721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]C(C)CC(C)=O(31204)',
    structure = SMILES('[CH2]C(C)CC(C)=O'),
    E0 = (-109.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,306.175,1376.02],'cm^-1')),
        HinderedRotor(inertia=(0.00179824,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0913426,'amu*angstrom^2'), symmetry=1, barrier=(6.0765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0913403,'amu*angstrom^2'), symmetry=1, barrier=(6.0765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.091343,'amu*angstrom^2'), symmetry=1, barrier=(6.0765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0913429,'amu*angstrom^2'), symmetry=1, barrier=(6.0765,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04634,0.0661533,-4.58101e-05,1.81312e-08,-3.17524e-12,-13090.4,27.5735], Tmin=(100,'K'), Tmax=(1262.19,'K')), NASAPolynomial(coeffs=[8.55011,0.0423729,-1.75489e-05,3.20397e-09,-2.18598e-13,-14984.6,-10.3751], Tmin=(1262.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(CC)C(C)=O(2028)',
    structure = SMILES('[CH2]C(CC)C(C)=O'),
    E0 = (-101.004,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360.849,3625.56],'cm^-1')),
        HinderedRotor(inertia=(0.0892061,'amu*angstrom^2'), symmetry=1, barrier=(8.24176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0891953,'amu*angstrom^2'), symmetry=1, barrier=(8.24159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237333,'amu*angstrom^2'), symmetry=1, barrier=(21.9289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0891997,'amu*angstrom^2'), symmetry=1, barrier=(8.24155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0891979,'amu*angstrom^2'), symmetry=1, barrier=(8.24169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.881823,0.0730601,-6.16129e-05,3.1203e-08,-7.03934e-12,-12039.5,26.8912], Tmin=(100,'K'), Tmax=(1011.05,'K')), NASAPolynomial(coeffs=[7.95836,0.0450632,-2.00762e-05,3.81445e-09,-2.6699e-13,-13470.4,-7.32674], Tmin=(1011.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CJC(C)C=O)"""),
)

species(
    label = 'CC[CH]C=C(C)O(31205)',
    structure = SMILES('CC[CH]C=C(C)O'),
    E0 = (-149.774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,393.728,393.752],'cm^-1')),
        HinderedRotor(inertia=(0.00108772,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129214,'amu*angstrom^2'), symmetry=1, barrier=(14.2116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129199,'amu*angstrom^2'), symmetry=1, barrier=(14.2115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12923,'amu*angstrom^2'), symmetry=1, barrier=(14.2109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129163,'amu*angstrom^2'), symmetry=1, barrier=(14.2108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.433755,0.0655857,-2.23356e-05,-1.91241e-08,1.24142e-11,-17873.8,25.8794], Tmin=(100,'K'), Tmax=(993.635,'K')), NASAPolynomial(coeffs=[16.453,0.0295345,-1.08397e-05,1.96453e-09,-1.38249e-13,-22461.1,-58.3657], Tmin=(993.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_S)"""),
)

species(
    label = 'C=C(O)C[CH]CC(3448)',
    structure = SMILES('C=C(O)C[CH]CC'),
    E0 = (-84.1825,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.311615,0.0724501,-5.46265e-05,2.18887e-08,-3.57965e-12,-9984.9,30.4355], Tmin=(100,'K'), Tmax=(1441.05,'K')), NASAPolynomial(coeffs=[15.1689,0.0312093,-1.1698e-05,2.02865e-09,-1.34187e-13,-14266.9,-46.6706], Tmin=(1441.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.1825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJCC)"""),
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
    label = '[CH]CC(C)=O(31206)',
    structure = SMILES('[CH]CC(C)=O'),
    E0 = (190.021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,392.831,392.834,1941.81,1941.81],'cm^-1')),
        HinderedRotor(inertia=(0.00109242,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737644,'amu*angstrom^2'), symmetry=1, barrier=(8.07777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737636,'amu*angstrom^2'), symmetry=1, barrier=(8.07778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37006,0.0389834,-2.93791e-05,1.33512e-08,-2.82496e-12,22910.3,17.8458], Tmin=(100,'K'), Tmax=(1033.13,'K')), NASAPolynomial(coeffs=[5.39921,0.0272553,-1.23508e-05,2.36297e-09,-1.65954e-13,22284.4,3.13325], Tmin=(1033.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]CC(144)',
    structure = SMILES('[CH]CC'),
    E0 = (328.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,328.86,330.153,2894.73],'cm^-1')),
        HinderedRotor(inertia=(0.00155435,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00154734,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00371,0.0161454,1.44268e-05,-2.62511e-08,1.01943e-11,39516.8,12.6846], Tmin=(100,'K'), Tmax=(1003.27,'K')), NASAPolynomial(coeffs=[6.63826,0.0156478,-5.75067e-06,1.05874e-09,-7.51289e-14,38083.3,-8.37163], Tmin=(1003.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet)"""),
)

species(
    label = 'propen2oxy(1009)',
    structure = SMILES('C=C(C)[O]'),
    E0 = (-43.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,459.471],'cm^-1')),
        HinderedRotor(inertia=(0.0528875,'amu*angstrom^2'), symmetry=1, barrier=(7.86742,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3365.98,'J/mol'), sigma=(5.64088,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.76 K, Pc=42.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84406,0.0189152,1.20804e-05,-2.62117e-08,1.05279e-11,-5215.42,14.5362], Tmin=(100,'K'), Tmax=(1005.62,'K')), NASAPolynomial(coeffs=[7.59011,0.0155447,-6.02335e-06,1.12468e-09,-8.0242e-14,-6954.09,-12.286], Tmin=(1005.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""propen2oxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CC[C]CC(C)=O(31207)',
    structure = SMILES('CC[C]CC(C)=O'),
    E0 = (142.461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13646,0.0678214,-5.44299e-05,2.59856e-08,-5.60703e-12,17233.1,25.6684], Tmin=(100,'K'), Tmax=(1038.59,'K')), NASAPolynomial(coeffs=[7.34992,0.0438912,-1.98686e-05,3.80102e-09,-2.66993e-13,15942.5,-4.54326], Tmin=(1038.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=C(O)[CH]CCC(31048)',
    structure = SMILES('[CH2]C(O)=CCCC'),
    E0 = (-131.97,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,316.139,316.267],'cm^-1')),
        HinderedRotor(inertia=(0.00168652,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190321,'amu*angstrom^2'), symmetry=1, barrier=(13.4749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190071,'amu*angstrom^2'), symmetry=1, barrier=(13.4745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190099,'amu*angstrom^2'), symmetry=1, barrier=(13.4692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189851,'amu*angstrom^2'), symmetry=1, barrier=(13.4712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.201522,0.0713412,-3.45101e-05,-1.13445e-08,1.11962e-11,-15724.5,28.5085], Tmin=(100,'K'), Tmax=(962.577,'K')), NASAPolynomial(coeffs=[17.8216,0.0265066,-8.87786e-06,1.53891e-09,-1.06649e-13,-20431.7,-62.6571], Tmin=(962.577,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH]=C(O)CCCC(31049)',
    structure = SMILES('[CH]=C(O)CCCC'),
    E0 = (-31.5447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,350,440,435,1725,3120,650,792.5,1650,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,273.307,273.314],'cm^-1')),
        HinderedRotor(inertia=(0.00225672,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25092,'amu*angstrom^2'), symmetry=1, barrier=(13.3009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250918,'amu*angstrom^2'), symmetry=1, barrier=(13.3009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250932,'amu*angstrom^2'), symmetry=1, barrier=(13.3009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250918,'amu*angstrom^2'), symmetry=1, barrier=(13.3009,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0996828,0.074292,-4.39491e-05,5.69428e-11,6.4938e-12,-3643.36,29.0263], Tmin=(100,'K'), Tmax=(993.019,'K')), NASAPolynomial(coeffs=[17.7305,0.0275066,-9.88385e-06,1.76273e-09,-1.22742e-13,-8339.73,-61.9247], Tmin=(993.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.5447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)CC[CH]C(31050)',
    structure = SMILES('C=C(O)CC[CH]C'),
    E0 = (-84.1946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.354972,0.0711555,-4.72451e-05,1.17244e-08,5.31001e-13,-9987.34,30.4438], Tmin=(100,'K'), Tmax=(1065.11,'K')), NASAPolynomial(coeffs=[14.5699,0.031933,-1.19516e-05,2.11693e-09,-1.4402e-13,-13818.7,-42.8027], Tmin=(1065.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.1946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CCCC(=C)O(31051)',
    structure = SMILES('[CH2]CCCC(=C)O'),
    E0 = (-73.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2950,3100,1380,975,1025,1650,3615,1277.5,1000,3000,3100,440,815,1455,1000,350,440,435,1725,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.214045,0.0714612,-3.75888e-05,-5.2864e-09,8.07015e-12,-8680.57,29.9994], Tmin=(100,'K'), Tmax=(994.317,'K')), NASAPolynomial(coeffs=[17.2628,0.0280813,-1.01703e-05,1.8238e-09,-1.27408e-13,-13316.9,-58.4192], Tmin=(994.317,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.3945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC(=C)[O](1018)',
    structure = SMILES('[CH2]CC(=C)[O]'),
    E0 = (111.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,350,440,435,1725,460.539,4000],'cm^-1')),
        HinderedRotor(inertia=(0.765416,'amu*angstrom^2'), symmetry=1, barrier=(17.5984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0401086,'amu*angstrom^2'), symmetry=1, barrier=(17.5479,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3511,'J/mol'), sigma=(5.95716,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=548.41 K, Pc=37.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0447,0.0358733,-1.01186e-05,-1.3997e-08,8.54798e-12,13543.8,20.4407], Tmin=(100,'K'), Tmax=(962.15,'K')), NASAPolynomial(coeffs=[10.8859,0.015978,-5.38751e-06,9.38125e-10,-6.51388e-14,11062.1,-25.9275], Tmin=(962.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = 'npropyl(83)',
    structure = SMILES('[CH2]CC'),
    E0 = (87.0621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0928812,'amu*angstrom^2'), symmetry=1, barrier=(2.13552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.092914,'amu*angstrom^2'), symmetry=1, barrier=(2.13628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0877,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.02815,0.0147023,2.4051e-05,-3.66738e-08,1.38611e-11,10512.1,12.4699], Tmin=(100,'K'), Tmax=(984.464,'K')), NASAPolynomial(coeffs=[6.16543,0.0184495,-6.79029e-06,1.23049e-09,-8.63866e-14,9095.06,-6.67607], Tmin=(984.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.0621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""npropyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(=C)[O](2376)',
    structure = SMILES('[CH2]C(=C)[O]'),
    E0 = (110.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1477.82],'cm^-1')),
        HinderedRotor(inertia=(0.530916,'amu*angstrom^2'), symmetry=1, barrier=(12.2068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.69662,0.0242907,-7.63674e-06,-9.31844e-09,6.21033e-12,13330.3,14.3522], Tmin=(100,'K'), Tmax=(924.136,'K')), NASAPolynomial(coeffs=[8.66414,0.00984176,-2.65664e-06,4.14902e-10,-2.77495e-14,11741.3,-16.5962], Tmin=(924.136,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CCC(=C)[O](2378)',
    structure = SMILES('[CH2]CCC(=C)[O]'),
    E0 = (88.1906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,350,440,435,1725,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00963226,'amu*angstrom^2'), symmetry=1, barrier=(1.72298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0972939,'amu*angstrom^2'), symmetry=1, barrier=(17.4018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0973056,'amu*angstrom^2'), symmetry=1, barrier=(17.4022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29222,0.0519151,-2.84102e-05,-1.00054e-10,3.88868e-12,10711,25.3862], Tmin=(100,'K'), Tmax=(1016.22,'K')), NASAPolynomial(coeffs=[12.3476,0.0234711,-8.67163e-06,1.54519e-09,-1.06384e-13,7685.87,-31.9563], Tmin=(1016.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.1906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])CC[CH]C(3009)',
    structure = SMILES('C=C([O])CC[CH]C'),
    E0 = (53.6102,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,292.399,292.4,1711.38,3959.71],'cm^-1')),
        HinderedRotor(inertia=(0.314661,'amu*angstrom^2'), symmetry=1, barrier=(19.0907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175916,'amu*angstrom^2'), symmetry=1, barrier=(10.6729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.314658,'amu*angstrom^2'), symmetry=1, barrier=(19.0907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00810084,'amu*angstrom^2'), symmetry=1, barrier=(90.133,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.909678,0.0649853,-4.71311e-05,1.88688e-08,-3.19336e-12,6561.35,29.9516], Tmin=(100,'K'), Tmax=(1355.13,'K')), NASAPolynomial(coeffs=[11.1222,0.0348408,-1.37641e-05,2.45373e-09,-1.65049e-13,3793.49,-22.4214], Tmin=(1355.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.6102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C][O](173)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97566e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02476e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CsCJ=O) + radical(CJC=O)"""),
)

species(
    label = 'C=C([O])[CH]CCC(2856)',
    structure = SMILES('[CH2]C([O])=CCCC'),
    E0 = (5.83526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,373.018,373.028,373.332],'cm^-1')),
        HinderedRotor(inertia=(0.00121155,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0813913,'amu*angstrom^2'), symmetry=1, barrier=(8.04061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0814427,'amu*angstrom^2'), symmetry=1, barrier=(8.04097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0813616,'amu*angstrom^2'), symmetry=1, barrier=(8.04126,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.525374,0.0678657,-4.37604e-05,7.97548e-09,2.2261e-12,834.377,28.8471], Tmin=(100,'K'), Tmax=(1013.82,'K')), NASAPolynomial(coeffs=[14.3435,0.0294284,-1.06846e-05,1.87219e-09,-1.27258e-13,-2793.92,-42.0832], Tmin=(1013.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.83526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CCCC(=C)[O](3010)',
    structure = SMILES('[CH2]CCCC(=C)[O]'),
    E0 = (64.4103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,350,440,435,1725,354.403,354.423,354.423,4000],'cm^-1')),
        HinderedRotor(inertia=(0.196046,'amu*angstrom^2'), symmetry=1, barrier=(17.4719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242259,'amu*angstrom^2'), symmetry=1, barrier=(2.15946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.024219,'amu*angstrom^2'), symmetry=1, barrier=(2.15946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196001,'amu*angstrom^2'), symmetry=1, barrier=(17.4719,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479778,0.0686037,-4.86553e-05,1.58413e-08,-1.4134e-12,7880.91,30.5507], Tmin=(100,'K'), Tmax=(1137.91,'K')), NASAPolynomial(coeffs=[14.4173,0.0299786,-1.14069e-05,2.02594e-09,-1.37357e-13,4037.72,-41.4399], Tmin=(1137.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.4103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C([O])CCCC(3012)',
    structure = SMILES('[CH]=C([O])CCCC'),
    E0 = (106.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,244.272,244.296,4000],'cm^-1')),
        HinderedRotor(inertia=(0.365703,'amu*angstrom^2'), symmetry=1, barrier=(15.4857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36551,'amu*angstrom^2'), symmetry=1, barrier=(15.4843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36548,'amu*angstrom^2'), symmetry=1, barrier=(15.4853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00282267,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281285,0.0723536,-5.78814e-05,2.43982e-08,-4.1484e-12,12921.9,29.8842], Tmin=(100,'K'), Tmax=(1403.53,'K')), NASAPolynomial(coeffs=[15.9272,0.0277633,-1.02261e-05,1.76224e-09,-1.16411e-13,8529.98,-50.902], Tmin=(1403.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCCC[C]1CO1(31052)',
    structure = SMILES('CCCC[C]1CO1'),
    E0 = (5.33346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.419173,0.0709109,-5.47799e-05,2.54535e-08,-4.94666e-12,776.982,27.4557], Tmin=(100,'K'), Tmax=(1263.21,'K')), NASAPolynomial(coeffs=[11.2879,0.0353372,-1.11633e-05,1.70914e-09,-1.03874e-13,-1876.54,-27.1534], Tmin=(1263.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.33346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(C2CsJO)"""),
)

species(
    label = 'C=[C]CCCC(373)',
    structure = SMILES('C=[C]CCCC'),
    E0 = (173.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,314.678,314.68,4000],'cm^-1')),
        HinderedRotor(inertia=(0.117801,'amu*angstrom^2'), symmetry=1, barrier=(8.27753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228969,'amu*angstrom^2'), symmetry=1, barrier=(16.0885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117808,'amu*angstrom^2'), symmetry=1, barrier=(8.27753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398702,'amu*angstrom^2'), symmetry=1, barrier=(28.0163,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.927626,0.0586121,-3.39837e-05,9.70484e-09,-1.11655e-12,21016.8,26.1062], Tmin=(100,'K'), Tmax=(1977.18,'K')), NASAPolynomial(coeffs=[16.2843,0.027544,-1.04134e-05,1.75728e-09,-1.11626e-13,14944.3,-58.4487], Tmin=(1977.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = 'CCCCC[C]=O(22490)',
    structure = SMILES('CCCCC[C]=O'),
    E0 = (-114.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09011,0.0677712,-4.97544e-05,2.11382e-08,-4.03303e-12,-13694.1,27.392], Tmin=(100,'K'), Tmax=(1156.73,'K')), NASAPolynomial(coeffs=[7.85583,0.0443751,-1.94152e-05,3.65251e-09,-2.53898e-13,-15259.3,-6.23367], Tmin=(1156.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(=C)OCCC(31055)',
    structure = SMILES('[CH2]C(=C)OCCC'),
    E0 = (-92.1621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,350,440,435,1725,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.128093,0.0730398,-3.83072e-05,-6.34755e-09,8.87267e-12,-10934.4,26.8259], Tmin=(100,'K'), Tmax=(984.806,'K')), NASAPolynomial(coeffs=[17.7611,0.0279169,-9.93768e-06,1.76876e-09,-1.23377e-13,-15692.3,-64.4967], Tmin=(984.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-92.1621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(O)CJ)"""),
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
    label = 'CCCC[C]=O(2388)',
    structure = SMILES('CCCC[C]=O'),
    E0 = (-90.9231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,301.503,3786.82],'cm^-1')),
        HinderedRotor(inertia=(0.144354,'amu*angstrom^2'), symmetry=1, barrier=(9.31236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144568,'amu*angstrom^2'), symmetry=1, barrier=(9.31535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267294,'amu*angstrom^2'), symmetry=1, barrier=(17.2386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144364,'amu*angstrom^2'), symmetry=1, barrier=(9.3122,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7077,0.053355,-3.72805e-05,1.49242e-08,-2.69409e-12,-10855.6,22.9276], Tmin=(100,'K'), Tmax=(1204.56,'K')), NASAPolynomial(coeffs=[6.95453,0.0359319,-1.55844e-05,2.91653e-09,-2.01974e-13,-12119.6,-3.36187], Tmin=(1204.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.9231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
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
    E0 = (344.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (0.152271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (33.4891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (31.8114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (11.2528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-24.0256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (53.1396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (19.2496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (46.3499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-65.9424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (259.455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (300.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (233.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (199.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (311.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (245.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (265.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (19.4952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (337.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (139.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (21.7459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-67.7253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (119.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (297.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (284.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (354.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (56.4218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (25.4115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (160.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-0.084133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-55.1576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-28.6257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (219.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (197.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (265.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (224.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (265.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (223.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (223.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (276.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (318.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (90.3802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (302.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (416.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (0.583161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (-82.2712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (-50.7185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (130.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (221.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (290.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction43',
    reactants = ['CC[CH]CC(140)', 'CO(12)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1076,'cm^3/(mol*s)'), n=3.29, Ea=(437.228,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [CO;Cs_Cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CC[CH]CC(C)=O(2056)'],
    products = ['CCC1CC1(C)[O](31202)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.44049e+09,'s^-1'), n=0.679905, Ea=(106.005,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHNd] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 104.9 to 106.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', 'CC=CCC(C)=O(3026)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', 'CCC=CC(C)=O(3025)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.165e+08,'cm^3/(mol*s)'), n=1.533, Ea=(5.18398,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 2853 used for Cds-COH_Cds;HJ
Exact match found for rate rule [Cds-COH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH3(17)', 'C=CCC(C)=O(760)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.105698,'m^3/(mol*s)'), n=2.13, Ea=(23.0748,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH3CO(57)', 'butene1(127)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0172287,'m^3/(mol*s)'), n=2.32603, Ea=(14.6351,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-CsH;CJ] for rate rule [Cds-HH_Cds-CsH;CO_rad/NonDe]
Euclidian distance = 3.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CC[CH]CC(C)=O(2056)'],
    products = ['C[CH]CCC(C)=O(31054)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CC[CH]CC(C)=O(2056)'],
    products = ['CCC[CH]C(C)=O(31053)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]CCCC(C)=O(22489)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(=O)CCCC(28404)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(754000,'s^-1'), n=1.63, Ea=(74.8936,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 110 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CH3(17)', '[CH2][CH]CC(C)=O(2610)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.34468e+07,'m^3/(mol*s)'), n=0.0549843, Ea=(0.0928142,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', 'C[CH][CH]CC(C)=O(3022)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH]CC(130)', 'CH3CO(57)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;CO_rad/NonDe]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(3)', 'CC[CH][CH]C(C)=O(3021)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(3)', '[CH2]C[CH]CC(C)=O(3023)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['CH3(17)', 'CC[CH]C[C]=O(2431)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.2e+13,'cm^3/(mol*s)','+|-',8.4e+12), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 72 used for C_methyl;CO_rad/NonDe
Exact match found for rate rule [CO_rad/NonDe;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(3)', 'C=C([O])C[CH]CC(1126)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.83648e+06,'m^3/(mol*s)'), n=0.3308, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;H_rad] for rate rule [C_rad/H2/CO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -1.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction41',
    reactants = ['CC[CH]CC(C)=O(2056)'],
    products = ['CCC1C[C](C)O1(31203)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.7342e+08,'s^-1'), n=0.920995, Ea=(125.348,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCs] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CH2(S)(23)', 'C[CH]CC(C)=O(1920)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(215646,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction44',
    reactants = ['CC[CH]CC(C)=O(2056)'],
    products = ['[CH2]C(C)CC(C)=O(31204)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;CH3] for rate rule [cCs(-HH)CJ;CsJ-CsH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C(CC)C(C)=O(2028)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO] for rate rule [cCs(-HC)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction46',
    reactants = ['CC[CH]C=C(C)O(31205)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(104,'s^-1'), n=3.21, Ea=(82.0482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C(O)C[CH]CC(3448)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_Cs;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C2H5(29)', '[CH]CC(C)=O(31206)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]CC(144)', 'propen2oxy(1009)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/CO;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction50',
    reactants = ['H(3)', 'CC[C]CC(C)=O(31207)'],
    products = ['CC[CH]CC(C)=O(2056)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction1',
    reactants = ['CH2CO(28)', 'butyl_1(82)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(8.04,'m^3/(mol*s)'), n=1.68, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Od;CsJ] for rate rule [Ck_O;CsJ-CsHH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(O)[CH]CCC(31048)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(187863,'s^-1'), n=2.03383, Ea=(157.381,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeC;O_H_out] + [R3H_SS_2Cd;C_rad_out_H/NonDeC;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(O)CCCC(31049)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C(O)C[CH]CC(3448)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;O_H_out] for rate rule [R4H_SS(Cd)S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C(O)CC[CH]C(31050)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.0508,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_H/NonDeC;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CCCC(=C)O(31051)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(104,'s^-1'), n=2.1, Ea=(44.7688,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSS(Cd)S;C_rad_out_2H;XH_out] for rate rule [R6H_SSSS(Cd)S;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C2H5(29)', '[CH2]CC(=C)[O](1018)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(8.73e+14,'cm^3/(mol*s)'), n=-0.699, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 11 used for C_rad/H2/Cs;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['npropyl(83)', '[CH2]C(=C)[O](2376)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 126 used for C_rad/H2/Cd;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cd;C_rad/H2/Cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C([O])C[CH]CC(1126)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH3(17)', '[CH2]CCC(=C)[O](2378)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=C([O])CC[CH]C(3009)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C][O](173)', 'butyl_1(82)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(5.25602e+06,'m^3/(mol*s)'), n=0.135617, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'C=C([O])[CH]CCC(2856)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.71464e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2]CCCC(=C)[O](3010)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]=C([O])CCCC(3012)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(=O)CCCC(28404)'],
    products = ['CCCC[C]1CO1(31052)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['CH2(S)(23)', '[CH2]C(=O)CCC(2148)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C]CCCC(373)', 'O(4)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(=O)CCCC(28404)'],
    products = ['CCC[CH]C(C)=O(31053)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C[CH]CCC(C)=O(31054)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(0.1524,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R5H_CC(O2d)CC;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]CCCC(C)=O(22489)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(3690,'s^-1'), n=1.79, Ea=(49.7896,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 113 used for R6H_SSSSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R6H_SSSSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CCCCC[C]=O(22490)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(=C)OCCC(31055)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(14080,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(19)', 'CCCC[C]=O(2388)'],
    products = ['[CH2]C(=O)CCCC(28404)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '5186',
    isomers = [
        'CC[CH]CC(C)=O(2056)',
        '[CH2]C(=O)CCCC(28404)',
    ],
    reactants = [
        ('CC[CH]CC(140)', 'CO(12)'),
        ('CH3CO(57)', 'butene1(127)'),
        ('CH2CO(28)', 'butyl_1(82)'),
        ('CH2(S)(23)', '[CH2]C(=O)CCC(2148)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '5186',
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

