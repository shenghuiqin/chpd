species(
    label = 'CC[CH]CC([O])=O(872)',
    structure = SMILES('CC[CH]CC([O])=O'),
    E0 = (-90.9629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2775,0.0477396,-2.22726e-05,3.94155e-09,-1.96743e-13,-10887.8,24.0178], Tmin=(100,'K'), Tmax=(2323.93,'K')), NASAPolynomial(coeffs=[24.2284,0.0180429,-8.32369e-06,1.43721e-09,-8.83974e-14,-23273.8,-105.091], Tmin=(2323.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.9629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCJCC=O) + radical(CCOJ)"""),
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
    label = 'CO2(13)',
    structure = SMILES('O=C=O'),
    E0 = (-403.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([459.166,1086.67,1086.68,2300.05],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2779,0.00275783,7.12787e-06,-1.07855e-08,4.14228e-12,-48475.6,5.97856], Tmin=(100,'K'), Tmax=(988.185,'K')), NASAPolynomial(coeffs=[4.55071,0.00290728,-1.14643e-06,2.25798e-10,-1.69526e-14,-48986,-1.45662], Tmin=(988.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.131,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'CCC1CC1([O])[O](2669)',
    structure = SMILES('CCC1CC1([O])[O]'),
    E0 = (64.8136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46579,0.047136,-1.07989e-05,-1.13476e-08,5.35354e-12,7893.54,24.2627], Tmin=(100,'K'), Tmax=(1179.3,'K')), NASAPolynomial(coeffs=[10.6293,0.0329101,-1.41434e-05,2.66264e-09,-1.85749e-13,4560.17,-26.4266], Tmin=(1179.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.8136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
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
    label = 'CC=CCC([O])=O(2681)',
    structure = SMILES('CC=CCC([O])=O'),
    E0 = (-169.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99522,0.041277,-6.92608e-06,-9.28011e-09,3.27172e-12,-20269.3,22.187], Tmin=(100,'K'), Tmax=(1481.02,'K')), NASAPolynomial(coeffs=[11.9654,0.032422,-1.62621e-05,3.16194e-09,-2.1938e-13,-25204.6,-36.5205], Tmin=(1481.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + radical(CCOJ)"""),
)

species(
    label = 'CCC=CC([O])=O(2682)',
    structure = SMILES('CCC=CC([O])=O'),
    E0 = (-86.8105,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,973.192,1600,2464.9,3200],'cm^-1')),
        HinderedRotor(inertia=(0.126308,'amu*angstrom^2'), symmetry=1, barrier=(2.90407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126308,'amu*angstrom^2'), symmetry=1, barrier=(2.90407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126308,'amu*angstrom^2'), symmetry=1, barrier=(2.90407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36599,0.0403472,-1.76014e-05,2.4552e-09,-3.70405e-14,-10444.7,19.0109], Tmin=(100,'K'), Tmax=(2577.02,'K')), NASAPolynomial(coeffs=[34.3218,0.0043,-3.60531e-06,6.41678e-10,-3.64286e-14,-30384.7,-167.369], Tmin=(2577.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.8105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + radical(CCOJ)"""),
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
    label = 'C=CCC([O])=O(2504)',
    structure = SMILES('C=CCC([O])=O'),
    E0 = (-133.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,542.206,542.229,542.364,542.403,542.505,542.64,542.664],'cm^-1')),
        HinderedRotor(inertia=(0.000572427,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000573398,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55128,0.0269555,8.43028e-06,-1.98824e-08,6.26926e-12,-15954.8,18.3716], Tmin=(100,'K'), Tmax=(1311.84,'K')), NASAPolynomial(coeffs=[8.76298,0.0273589,-1.41494e-05,2.83278e-09,-2.01709e-13,-19249,-19.6263], Tmin=(1311.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-133.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = '[O][C]=O(669)',
    structure = SMILES('[O][C]=O'),
    E0 = (31.9507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90478,-0.000175995,8.15126e-06,-1.13656e-08,4.4768e-12,3848.25,8.04855], Tmin=(100,'K'), Tmax=(975.388,'K')), NASAPolynomial(coeffs=[5.59398,-0.00122084,7.11747e-07,-9.7712e-11,3.97995e-15,3238.91,-1.49318], Tmin=(975.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.9507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(OJC=O) + radical((O)CJOH)"""),
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
    label = 'C[CH]CCC([O])=O(2683)',
    structure = SMILES('C[CH]CCC([O])=O'),
    E0 = (-96.4186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63189,0.0607512,-7.73694e-05,7.92687e-08,-3.33176e-11,-11519.7,26.5874], Tmin=(100,'K'), Tmax=(800.761,'K')), NASAPolynomial(coeffs=[-0.885262,0.053252,-2.57208e-05,4.96463e-09,-3.46308e-13,-10473,42.1903], Tmin=(800.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.4186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(RCCJC) + radical(CCOJ)"""),
)

species(
    label = 'CCC[CH]C([O])=O(2567)',
    structure = SMILES('CCCC=C([O])[O]'),
    E0 = (-123.306,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,374.346,374.346,374.35,374.35],'cm^-1')),
        HinderedRotor(inertia=(0.00120296,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12703,'amu*angstrom^2'), symmetry=1, barrier=(12.6323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127029,'amu*angstrom^2'), symmetry=1, barrier=(12.6323,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.720821,0.0645773,-5.02095e-05,1.9969e-08,-3.20328e-12,-14706.1,26.5685], Tmin=(100,'K'), Tmax=(1473.1,'K')), NASAPolynomial(coeffs=[15.2262,0.0251895,-1.01021e-05,1.81791e-09,-1.22821e-13,-18979.6,-49.0305], Tmin=(1473.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-123.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]CCCC([O])=O(2684)',
    structure = SMILES('[CH2]CCCC([O])=O'),
    E0 = (-85.6186,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.73248,0.0472999,-2.20852e-05,3.91914e-09,-2.16378e-13,-10271.6,21.6077], Tmin=(100,'K'), Tmax=(2587.01,'K')), NASAPolynomial(coeffs=[33.208,0.00909484,-5.10271e-06,8.74953e-10,-5.09347e-14,-29023.1,-160.152], Tmin=(2587.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-85.6186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = 'CC[CH][CH]C(=O)O(2685)',
    structure = SMILES('CC[CH]C=C([O])O'),
    E0 = (-123.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.555939,0.065916,-4.40933e-05,5.90404e-09,3.5986e-12,-14739.8,24.6339], Tmin=(100,'K'), Tmax=(1007.14,'K')), NASAPolynomial(coeffs=[16.5011,0.022202,-8.20026e-06,1.48227e-09,-1.03839e-13,-18946.3,-57.3443], Tmin=(1007.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-123.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(Allyl_S) + radical(C=COJ)"""),
)

species(
    label = 'C[CH][CH]CC(=O)O(2686)',
    structure = SMILES('C[CH][CH]CC(=O)O'),
    E0 = (-122.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17012,0.0474662,-2.22047e-06,-6.36635e-08,5.79051e-11,-14641.3,27.12], Tmin=(100,'K'), Tmax=(471.07,'K')), NASAPolynomial(coeffs=[3.70061,0.0438014,-2.02627e-05,3.91947e-09,-2.77489e-13,-14889,19.7895], Tmin=(471.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCJCC=O) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C[CH]CC(=O)O(2687)',
    structure = SMILES('[CH2]C[CH]CC(=O)O'),
    E0 = (-111.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55985,0.0558,-3.64256e-05,1.21343e-08,-1.71749e-12,-13315.6,28.2184], Tmin=(100,'K'), Tmax=(1528.68,'K')), NASAPolynomial(coeffs=[9.88296,0.0340215,-1.50557e-05,2.8148e-09,-1.93373e-13,-15860.3,-15.4682], Tmin=(1528.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(RCCJ) + radical(CCJCC=O)"""),
)

species(
    label = 'CC[CH]C[C]1OO1(2688)',
    structure = SMILES('CC[CH]C[C]1OO1'),
    E0 = (265.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20588,0.0539316,-2.53979e-05,-3.85724e-09,5.27207e-12,32072.3,28.5402], Tmin=(100,'K'), Tmax=(994.847,'K')), NASAPolynomial(coeffs=[11.5057,0.0281814,-1.01882e-05,1.77576e-09,-1.20301e-13,29247.9,-24.993], Tmin=(994.847,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(dioxirane) + radical(RCCJCC) + radical(Cs_P)"""),
)

species(
    label = 'CCC1C[C]([O])O1(2689)',
    structure = SMILES('CCC1C[C]([O])O1'),
    E0 = (74.1976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60626,0.0450723,-1.94563e-06,-2.67019e-08,1.35955e-11,9017.03,24.9164], Tmin=(100,'K'), Tmax=(926.56,'K')), NASAPolynomial(coeffs=[9.31835,0.0310068,-1.0303e-05,1.70806e-09,-1.12886e-13,6762.52,-16.1557], Tmin=(926.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.1976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Oxetane) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = 'CCC=CC(=O)O(884)',
    structure = SMILES('CCC=CC(=O)O'),
    E0 = (-312.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62724,0.0575218,-4.91779e-05,2.76684e-08,-7.45401e-12,-37506.1,25.2291], Tmin=(100,'K'), Tmax=(830.598,'K')), NASAPolynomial(coeffs=[5.04631,0.0410563,-1.94425e-05,3.80188e-09,-2.70506e-13,-38074.1,9.36865], Tmin=(830.598,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s)"""),
)

species(
    label = 'CC=CCC(=O)O(2690)',
    structure = SMILES('CC=CCC(=O)O'),
    E0 = (-394.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4171,0.0469088,-6.28232e-06,-1.70677e-08,7.07774e-12,-47388.6,24.071], Tmin=(100,'K'), Tmax=(1205.94,'K')), NASAPolynomial(coeffs=[12.4713,0.0321994,-1.52966e-05,3.01337e-09,-2.15088e-13,-51651.3,-37.9489], Tmin=(1205.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-394.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs)"""),
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
    label = 'C[CH]CC([O])=O(928)',
    structure = SMILES('C[CH]CC([O])=O'),
    E0 = (-67.1826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4032.06,'J/mol'), sigma=(6.27406,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=629.80 K, Pc=37.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76998,0.0344919,-1.27638e-05,4.31799e-10,3.35012e-13,-8042.72,20.0272], Tmin=(100,'K'), Tmax=(2034.82,'K')), NASAPolynomial(coeffs=[15.1596,0.0199348,-9.25563e-06,1.64888e-09,-1.05264e-13,-15113.3,-53.5318], Tmin=(2034.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.1826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]C(C)CC([O])=O(2691)',
    structure = SMILES('[CH2]C(C)CC([O])=O'),
    E0 = (-94.8241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53957,0.0586487,-5.57417e-05,4.33021e-08,-1.63059e-11,-11320.4,25.6945], Tmin=(100,'K'), Tmax=(734.139,'K')), NASAPolynomial(coeffs=[3.11797,0.0458581,-2.10454e-05,4.01923e-09,-2.80924e-13,-11439.2,19.3366], Tmin=(734.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.8241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(Isobutyl) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(CC)C([O])=O(873)',
    structure = SMILES('[CH2]C(CC)C([O])=O'),
    E0 = (-86.1147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4190.81,'J/mol'), sigma=(6.62078,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=654.59 K, Pc=32.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18859,0.0680442,-8.15552e-05,7.11743e-08,-2.73015e-11,-10261.9,25.6609], Tmin=(100,'K'), Tmax=(751.672,'K')), NASAPolynomial(coeffs=[3.71807,0.0464834,-2.2365e-05,4.34135e-09,-3.05204e-13,-10413.4,15.7019], Tmin=(751.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.1147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CJC(C)C=O)"""),
)

species(
    label = 'CCC1CC(=O)O1(2644)',
    structure = SMILES('CCC1CC(=O)O1'),
    E0 = (-356.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98376,0.0277322,5.75071e-05,-9.39017e-08,3.84791e-11,-42809.4,22.1927], Tmin=(100,'K'), Tmax=(932.065,'K')), NASAPolynomial(coeffs=[12.2631,0.0258749,-7.50904e-06,1.24277e-09,-8.73182e-14,-46561.2,-36.5225], Tmin=(932.065,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-356.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone)"""),
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
    label = '[CH2]C([O])=O(1440)',
    structure = SMILES('[CH2]C([O])=O'),
    E0 = (-8.20105,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,841.694,844.034,844.61,844.694,845.141],'cm^-1')),
        HinderedRotor(inertia=(0.00480388,'amu*angstrom^2'), symmetry=1, barrier=(2.45209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.60873,0.0103835,4.164e-06,-6.72683e-09,1.73621e-12,-973.928,12.8424], Tmin=(100,'K'), Tmax=(1564.29,'K')), NASAPolynomial(coeffs=[5.70762,0.0133434,-6.65897e-06,1.2886e-09,-8.86347e-14,-2649.38,-1.47903], Tmin=(1564.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.20105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJCO) + radical(CCOJ)"""),
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
    E0 = (-90.9629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (64.8136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (48.3786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (130.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (26.1423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (30.1533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-90.9629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (68.0291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (34.1391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (61.2394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (20.0817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-40.0854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (2.54954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (287.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (265.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (74.1976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-27.5627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-65.9896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (352.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (154.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (36.6355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-82.6785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (320.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (351.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['butene1(127)', 'CO2(13)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['CCC1CC1([O])[O](2669)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.44049e+09,'s^-1'), n=0.679905, Ea=(155.776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHNd] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_csHNd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic
Ea raised from 154.5 to 155.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'CC=CCC([O])=O(2681)'],
    products = ['CC[CH]CC([O])=O(872)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'CCC=CC([O])=O(2682)'],
    products = ['CC[CH]CC([O])=O(872)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.165e+08,'cm^3/(mol*s)'), n=1.533, Ea=(5.18398,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 2853 used for Cds-COH_Cds;HJ
Exact match found for rate rule [Cds-COH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH3(17)', 'C=CCC([O])=O(2504)'],
    products = ['CC[CH]CC([O])=O(872)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.105698,'m^3/(mol*s)'), n=2.13, Ea=(23.0748,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['butene1(127)', '[O][C]=O(669)'],
    products = ['CC[CH]CC([O])=O(872)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0172287,'m^3/(mol*s)'), n=2.32603, Ea=(14.6351,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH]CC(130)', 'CO2(13)'],
    products = ['CC[CH]CC([O])=O(872)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.04,'m^3/(mol*s)'), n=1.68, Ea=(56.4997,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Od;CJ] for rate rule [CO2;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 54.2 to 56.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['C[CH]CCC([O])=O(2683)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['CCC[CH]C([O])=O(2567)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]CCCC([O])=O(2684)'],
    products = ['CC[CH]CC([O])=O(872)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['CC[CH][CH]C(=O)O(2685)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.19965e+08,'s^-1'), n=1.57622, Ea=(111.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['C[CH][CH]CC(=O)O(2686)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.39199e+09,'s^-1'), n=0.63, Ea=(50.8774,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['[CH2]C[CH]CC(=O)O(2687)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.36e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6Hall;O_rad_out;Cs_H_out_2H] for rate rule [R6HJ_3;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]CC(130)', '[O][C]=O(669)'],
    products = ['CC[CH]CC([O])=O(872)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['CC[CH]C[C]1OO1(2688)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(356.736,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 356.5 to 356.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['CCC1C[C]([O])O1(2689)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.46839e+08,'s^-1'), n=0.920995, Ea=(165.16,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCs] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHCs]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 163.9 to 165.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['CCC=CC(=O)O(884)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['CC=CCC(=O)O(2690)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.50442e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH2(S)(23)', 'C[CH]CC([O])=O(928)'],
    products = ['CC[CH]CC([O])=O(872)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(215646,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['[CH2]C(C)CC([O])=O(2691)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;CH3] for rate rule [cCs(-HH)CJ;CsJ-CsH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(CC)C([O])=O(873)'],
    products = ['CC[CH]CC([O])=O(872)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO] for rate rule [cCs(-HC)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC[CH]CC([O])=O(872)'],
    products = ['CCC1CC(=O)O1(2644)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]CC(144)', '[CH2]C([O])=O(1440)'],
    products = ['CC[CH]CC([O])=O(872)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/CO;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(4)', 'CC[CH]C[C]=O(2431)'],
    products = ['CC[CH]CC([O])=O(872)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '416',
    isomers = [
        'CC[CH]CC([O])=O(872)',
    ],
    reactants = [
        ('butene1(127)', 'CO2(13)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '416',
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

