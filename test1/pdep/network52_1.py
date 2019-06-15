species(
    label = 'C=CC=CC=CC(99)',
    structure = SMILES('C=CC=CC=CC'),
    E0 = (109.668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.907165,'amu*angstrom^2'), symmetry=1, barrier=(20.8575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.901128,'amu*angstrom^2'), symmetry=1, barrier=(20.7187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.902516,'amu*angstrom^2'), symmetry=1, barrier=(20.7506,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3523.79,'J/mol'), sigma=(5.9032,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=550.41 K, Pc=38.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.707351,0.0573834,-2.46786e-06,-4.14155e-08,2.1071e-11,13322.1,24.2542], Tmin=(100,'K'), Tmax=(964.439,'K')), NASAPolynomial(coeffs=[17.1633,0.0251949,-8.49298e-06,1.52028e-09,-1.0885e-13,8470.73,-63.2356], Tmin=(964.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC1C=CC=CC1(203)',
    structure = SMILES('CC1C=CC=CC1'),
    E0 = (43.3206,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3742.82,'J/mol'), sigma=(6.34669,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=584.62 K, Pc=33.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86632,0.0262021,7.67233e-05,-1.12919e-07,4.34216e-11,5305.66,19.7608], Tmin=(100,'K'), Tmax=(974.026,'K')), NASAPolynomial(coeffs=[13.1932,0.0315886,-1.15014e-05,2.17362e-09,-1.60274e-13,637.056,-47.2255], Tmin=(974.026,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.3206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]=CC=CC=C(342)',
    structure = SMILES('[CH]=CC=CC=C'),
    E0 = (392.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.10935,'amu*angstrom^2'), symmetry=1, barrier=(25.5061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10912,'amu*angstrom^2'), symmetry=1, barrier=(25.5008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31283,0.044935,2.87911e-06,-4.44372e-08,2.2804e-11,47351.4,20.6786], Tmin=(100,'K'), Tmax=(937.297,'K')), NASAPolynomial(coeffs=[17.2005,0.0128923,-3.06895e-06,4.97402e-10,-3.78212e-14,42802.3,-63.3209], Tmin=(937.297,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'CH3(23)',
    structure = SMILES('[CH3]'),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'CH2CHCHCH(119)',
    structure = SMILES('[CH]=CC=C'),
    E0 = (346.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.31937,'amu*angstrom^2'), symmetry=1, barrier=(30.3349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64258,0.0163334,3.86236e-05,-6.71392e-08,2.8361e-11,41729.6,13.2819], Tmin=(100,'K'), Tmax=(937.72,'K')), NASAPolynomial(coeffs=[12.9704,0.0066914,-1.00078e-06,1.6762e-10,-1.71452e-14,38279.7,-43.9471], Tmin=(937.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH2CHCHCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'propen1yl(311)',
    structure = SMILES('[CH]=CC'),
    E0 = (254.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,180,694.529,1531.96,4000],'cm^-1')),
        HinderedRotor(inertia=(0.102118,'amu*angstrom^2'), symmetry=1, barrier=(2.34789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22801,0.0118522,1.72183e-05,-2.70814e-08,1.02844e-11,30640.6,10.1787], Tmin=(100,'K'), Tmax=(988.714,'K')), NASAPolynomial(coeffs=[5.81113,0.014019,-5.21094e-06,9.49007e-10,-6.67726e-14,29513.2,-5.37277], Tmin=(988.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""propen1yl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=CC=CC(344)',
    structure = SMILES('[CH]=CC=CC'),
    E0 = (304.992,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680],'cm^-1')),
        HinderedRotor(inertia=(0.867881,'amu*angstrom^2'), symmetry=1, barrier=(19.9543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.872416,'amu*angstrom^2'), symmetry=1, barrier=(20.0586,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93302,0.0363757,-1.10365e-06,-2.51393e-08,1.24817e-11,36764.6,17.623], Tmin=(100,'K'), Tmax=(973.632,'K')), NASAPolynomial(coeffs=[11.5283,0.0184064,-6.46807e-06,1.16265e-09,-8.23217e-14,33879.4,-33.6337], Tmin=(973.632,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C2H3(32)',
    structure = SMILES('[CH]=C'),
    E0 = (286.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,677.08,1086.68,3788.01],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]C=CC[C]=CC(207)',
    structure = SMILES('[CH2]C=CC[C]=CC'),
    E0 = (406.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,192.05,1830.36],'cm^-1')),
        HinderedRotor(inertia=(0.42076,'amu*angstrom^2'), symmetry=1, barrier=(11.0206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423412,'amu*angstrom^2'), symmetry=1, barrier=(11.0321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.787584,'amu*angstrom^2'), symmetry=1, barrier=(20.5268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784888,'amu*angstrom^2'), symmetry=1, barrier=(20.5111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.634809,0.0641336,-4.13358e-05,1.33798e-08,-1.75261e-12,48964.7,29.2977], Tmin=(100,'K'), Tmax=(1760.9,'K')), NASAPolynomial(coeffs=[16.2867,0.0285794,-1.10495e-05,1.91361e-09,-1.24723e-13,43452.4,-55.0697], Tmin=(1760.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C[C]=CC=CC(90)',
    structure = SMILES('[CH2]C[C]=CC=CC'),
    E0 = (442.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,320.475,320.478],'cm^-1')),
        HinderedRotor(inertia=(0.158764,'amu*angstrom^2'), symmetry=1, barrier=(11.5711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158764,'amu*angstrom^2'), symmetry=1, barrier=(11.5711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158765,'amu*angstrom^2'), symmetry=1, barrier=(11.5711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158764,'amu*angstrom^2'), symmetry=1, barrier=(11.5711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.520891,0.0677225,-4.89514e-05,1.84228e-08,-2.82295e-12,53330.4,29.2799], Tmin=(100,'K'), Tmax=(1528.04,'K')), NASAPolynomial(coeffs=[15.204,0.0292875,-1.12231e-05,1.96295e-09,-1.30076e-13,48843,-47.7835], Tmin=(1528.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = 'C=CC=[C]C[CH]C(227)',
    structure = SMILES('C=CC=[C]C[CH]C'),
    E0 = (443.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,358.343,358.344,358.395],'cm^-1')),
        HinderedRotor(inertia=(0.111783,'amu*angstrom^2'), symmetry=1, barrier=(10.1879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111806,'amu*angstrom^2'), symmetry=1, barrier=(10.1879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0013127,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111785,'amu*angstrom^2'), symmetry=1, barrier=(10.188,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.88849,0.0641557,-4.41666e-05,1.62849e-08,-2.51557e-12,53487,29.3363], Tmin=(100,'K'), Tmax=(1475.31,'K')), NASAPolynomial(coeffs=[12.0514,0.0338898,-1.33942e-05,2.37939e-09,-1.59213e-13,50193.3,-28.8594], Tmin=(1475.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(RCCJC)"""),
)

species(
    label = 'C=[C]CC=C[CH]C(218)',
    structure = SMILES('C=[C]C[CH]C=CC'),
    E0 = (406.629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,358.727,359.512,359.667],'cm^-1')),
        HinderedRotor(inertia=(0.0719273,'amu*angstrom^2'), symmetry=1, barrier=(6.60742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0715721,'amu*angstrom^2'), symmetry=1, barrier=(6.5473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00131797,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281906,'amu*angstrom^2'), symmetry=1, barrier=(26.1891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.851833,0.0605612,-3.44136e-05,8.58458e-09,-6.27474e-13,49026.7,27.5401], Tmin=(100,'K'), Tmax=(1461.78,'K')), NASAPolynomial(coeffs=[13.9434,0.0322167,-1.30028e-05,2.32011e-09,-1.54722e-13,44400.2,-43.3229], Tmin=(1461.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C=CCC=[C]C(209)',
    structure = SMILES('[CH2]C=CCC=[C]C'),
    E0 = (406.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,192.05,1830.36],'cm^-1')),
        HinderedRotor(inertia=(0.42076,'amu*angstrom^2'), symmetry=1, barrier=(11.0206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423412,'amu*angstrom^2'), symmetry=1, barrier=(11.0321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.787584,'amu*angstrom^2'), symmetry=1, barrier=(20.5268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784888,'amu*angstrom^2'), symmetry=1, barrier=(20.5111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.634809,0.0641336,-4.13358e-05,1.33798e-08,-1.75261e-12,48964.7,29.2977], Tmin=(100,'K'), Tmax=(1760.9,'K')), NASAPolynomial(coeffs=[16.2867,0.0285794,-1.10495e-05,1.91361e-09,-1.24723e-13,43452.4,-55.0697], Tmin=(1760.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC=[C]C=CC(88)',
    structure = SMILES('[CH2]CC=[C]C=CC'),
    E0 = (403.466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,245.938,246.102],'cm^-1')),
        HinderedRotor(inertia=(0.284455,'amu*angstrom^2'), symmetry=1, barrier=(12.2116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284338,'amu*angstrom^2'), symmetry=1, barrier=(12.2037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284454,'amu*angstrom^2'), symmetry=1, barrier=(12.2056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.78521,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.532199,0.0691234,-5.19103e-05,2.10377e-08,-3.51018e-12,48656.5,28.3541], Tmin=(100,'K'), Tmax=(1406.22,'K')), NASAPolynomial(coeffs=[13.7251,0.0315959,-1.188e-05,2.05991e-09,-1.36257e-13,44946.1,-39.7915], Tmin=(1406.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(RCCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=CC[CH]C(228)',
    structure = SMILES('C=C[C]=CC[CH]C'),
    E0 = (404.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,371.557,371.706,372.077],'cm^-1')),
        HinderedRotor(inertia=(0.0863704,'amu*angstrom^2'), symmetry=1, barrier=(8.4807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0864237,'amu*angstrom^2'), symmetry=1, barrier=(8.47718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0863683,'amu*angstrom^2'), symmetry=1, barrier=(8.47631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.729437,'amu*angstrom^2'), symmetry=1, barrier=(71.6353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.86749,0.0658642,-4.7922e-05,1.96363e-08,-3.42368e-12,48814.8,28.5322], Tmin=(100,'K'), Tmax=(1316.94,'K')), NASAPolynomial(coeffs=[10.6687,0.0360946,-1.40143e-05,2.47144e-09,-1.65213e-13,46233.2,-21.4513], Tmin=(1316.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJC) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CCC=C[CH]C(220)',
    structure = SMILES('[CH]=CC[CH]C=CC'),
    E0 = (415.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,180,856.696],'cm^-1')),
        HinderedRotor(inertia=(0.136629,'amu*angstrom^2'), symmetry=1, barrier=(3.14136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136592,'amu*angstrom^2'), symmetry=1, barrier=(3.14051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.852823,'amu*angstrom^2'), symmetry=1, barrier=(19.6081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.852881,'amu*angstrom^2'), symmetry=1, barrier=(19.6094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.892077,0.0587911,-2.50278e-05,-3.37475e-09,3.83915e-12,50139,27.3069], Tmin=(100,'K'), Tmax=(1142.37,'K')), NASAPolynomial(coeffs=[12.5451,0.0341783,-1.39685e-05,2.57742e-09,-1.78439e-13,46420.2,-35.0871], Tmin=(1142.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_S) + radical(Cds_P)"""),
)

species(
    label = 'C=CC[CH]C=[C]C(224)',
    structure = SMILES('C=CC[CH]C=[C]C'),
    E0 = (406.629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,358.727,359.512,359.667],'cm^-1')),
        HinderedRotor(inertia=(0.0719273,'amu*angstrom^2'), symmetry=1, barrier=(6.60742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0715721,'amu*angstrom^2'), symmetry=1, barrier=(6.5473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00131797,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281906,'amu*angstrom^2'), symmetry=1, barrier=(26.1891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.851833,0.0605612,-3.44136e-05,8.58458e-09,-6.27474e-13,49026.7,27.5401], Tmin=(100,'K'), Tmax=(1461.78,'K')), NASAPolynomial(coeffs=[13.9434,0.0322167,-1.30028e-05,2.32011e-09,-1.54722e-13,44400.2,-43.3229], Tmin=(1461.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = 'C[CH]C=[C]C=CC(198)',
    structure = SMILES('C[CH]C=[C]C=CC'),
    E0 = (339.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,289.784,290.785],'cm^-1')),
        HinderedRotor(inertia=(0.220659,'amu*angstrom^2'), symmetry=1, barrier=(13.4783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236812,'amu*angstrom^2'), symmetry=1, barrier=(13.5071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343135,'amu*angstrom^2'), symmetry=1, barrier=(19.8796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29984,'amu*angstrom^2'), symmetry=1, barrier=(75.2831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736739,0.0627147,-3.12792e-05,-7.09226e-10,4.05418e-12,40937.5,24.7234], Tmin=(100,'K'), Tmax=(1053.11,'K')), NASAPolynomial(coeffs=[12.798,0.0332713,-1.26562e-05,2.26078e-09,-1.5461e-13,37489.5,-38.3988], Tmin=(1053.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(Allyl_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=[C]C=CCC(89)',
    structure = SMILES('[CH2]C=[C]C=CCC'),
    E0 = (316.276,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,306.362,307.249],'cm^-1')),
        HinderedRotor(inertia=(0.305649,'amu*angstrom^2'), symmetry=1, barrier=(20.9932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08553,'amu*angstrom^2'), symmetry=1, barrier=(74.0827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185113,'amu*angstrom^2'), symmetry=1, barrier=(12.353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10565,'amu*angstrom^2'), symmetry=1, barrier=(74.0671,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865961,0.0585054,-1.52556e-05,-2.01916e-08,1.18581e-11,38161.3,26.9224], Tmin=(100,'K'), Tmax=(978.011,'K')), NASAPolynomial(coeffs=[13.0041,0.0324953,-1.16116e-05,2.03327e-09,-1.39165e-13,34656.7,-37.1458], Tmin=(978.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C[CH]CC=CC(225)',
    structure = SMILES('[CH]C=CCC=CC'),
    E0 = (387.383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.809834,0.0605439,-2.09133e-05,-6.63373e-09,4.45419e-12,46714.2,29.0052], Tmin=(100,'K'), Tmax=(1167.28,'K')), NASAPolynomial(coeffs=[11.384,0.0407837,-1.66917e-05,3.04644e-09,-2.08639e-13,43123.2,-28.4525], Tmin=(1167.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C=CCC=C[CH2](211)',
    structure = SMILES('[CH2]C=CCC=C[CH2]'),
    E0 = (319.697,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,819.749,1702.88],'cm^-1')),
        HinderedRotor(inertia=(0.939523,'amu*angstrom^2'), symmetry=1, barrier=(21.6015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.939492,'amu*angstrom^2'), symmetry=1, barrier=(21.6008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.939531,'amu*angstrom^2'), symmetry=1, barrier=(21.6017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.9395,'amu*angstrom^2'), symmetry=1, barrier=(21.601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3519.58,'J/mol'), sigma=(6.10993,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.75 K, Pc=35.01 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.89194,0.0548011,-2.01022e-06,-3.44964e-08,1.64354e-11,38574.5,27.4606], Tmin=(100,'K'), Tmax=(1007.62,'K')), NASAPolynomial(coeffs=[14.652,0.0308082,-1.18922e-05,2.21116e-09,-1.57403e-13,34246.5,-46.7446], Tmin=(1007.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.697,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C=CC[CH]C(229)',
    structure = SMILES('C=[C]C=CC[CH]C'),
    E0 = (404.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,371.557,371.706,372.077],'cm^-1')),
        HinderedRotor(inertia=(0.0863704,'amu*angstrom^2'), symmetry=1, barrier=(8.4807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0864237,'amu*angstrom^2'), symmetry=1, barrier=(8.47718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0863683,'amu*angstrom^2'), symmetry=1, barrier=(8.47631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.729437,'amu*angstrom^2'), symmetry=1, barrier=(71.6353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.86749,0.0658642,-4.7922e-05,1.96363e-08,-3.42368e-12,48814.8,28.5322], Tmin=(100,'K'), Tmax=(1316.94,'K')), NASAPolynomial(coeffs=[10.6687,0.0360946,-1.40143e-05,2.47144e-09,-1.65213e-13,46233.2,-21.4513], Tmin=(1316.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJC) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]CC=C[C]=CC(86)',
    structure = SMILES('[CH2]CC=C[C]=CC'),
    E0 = (403.466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,245.938,246.102],'cm^-1')),
        HinderedRotor(inertia=(0.284455,'amu*angstrom^2'), symmetry=1, barrier=(12.2116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284338,'amu*angstrom^2'), symmetry=1, barrier=(12.2037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284454,'amu*angstrom^2'), symmetry=1, barrier=(12.2056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.78521,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.532199,0.0691234,-5.19103e-05,2.10377e-08,-3.51018e-12,48656.5,28.3541], Tmin=(100,'K'), Tmax=(1406.22,'K')), NASAPolynomial(coeffs=[13.7251,0.0315959,-1.188e-05,2.05991e-09,-1.36257e-13,44946.1,-39.7915], Tmin=(1406.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(RCCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][CH]C=CCC=C(108)',
    structure = SMILES('[CH2]C=C[CH]CC=C'),
    E0 = (320.286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,382.277,382.289,382.292],'cm^-1')),
        HinderedRotor(inertia=(0.00115341,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00115229,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.299501,'amu*angstrom^2'), symmetry=1, barrier=(31.0695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.299303,'amu*angstrom^2'), symmetry=1, barrier=(31.0696,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03012,0.0521083,2.03792e-06,-3.57957e-08,1.61466e-11,38640,26.6828], Tmin=(100,'K'), Tmax=(1021.44,'K')), NASAPolynomial(coeffs=[13.6024,0.0325733,-1.28876e-05,2.41084e-09,-1.71393e-13,34522.3,-41.822], Tmin=(1021.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CC=CCC(91)',
    structure = SMILES('C=[C]C=C[CH]CC'),
    E0 = (351.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,180,446.295,446.301],'cm^-1')),
        HinderedRotor(inertia=(2.36815,'amu*angstrom^2'), symmetry=1, barrier=(77.4776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0962383,'amu*angstrom^2'), symmetry=1, barrier=(13.6037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0963124,'amu*angstrom^2'), symmetry=1, barrier=(13.6037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167911,'amu*angstrom^2'), symmetry=1, barrier=(23.7276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.708447,0.0622886,-2.60167e-05,-8.89312e-09,7.58607e-12,42412.4,25.2624], Tmin=(100,'K'), Tmax=(1012.75,'K')), NASAPolynomial(coeffs=[13.6629,0.0319467,-1.19192e-05,2.12966e-09,-1.46746e-13,38720.6,-42.672], Tmin=(1012.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Allyl_S)"""),
)

species(
    label = 'C[CH]C=C[C]=CC(197)',
    structure = SMILES('C[CH]C=C[C]=CC'),
    E0 = (339.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,289.784,290.785],'cm^-1')),
        HinderedRotor(inertia=(0.220659,'amu*angstrom^2'), symmetry=1, barrier=(13.4783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236812,'amu*angstrom^2'), symmetry=1, barrier=(13.5071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343135,'amu*angstrom^2'), symmetry=1, barrier=(19.8796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29984,'amu*angstrom^2'), symmetry=1, barrier=(75.2831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736739,0.0627147,-3.12792e-05,-7.09226e-10,4.05418e-12,40937.5,24.7234], Tmin=(100,'K'), Tmax=(1053.11,'K')), NASAPolynomial(coeffs=[12.798,0.0332713,-1.26562e-05,2.26078e-09,-1.5461e-13,37489.5,-38.3988], Tmin=(1053.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(Allyl_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CC=CC[CH]C(230)',
    structure = SMILES('[CH]=CC=CC[CH]C'),
    E0 = (453.012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,264.711,945.887],'cm^-1')),
        HinderedRotor(inertia=(0.150498,'amu*angstrom^2'), symmetry=1, barrier=(7.51289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286962,'amu*angstrom^2'), symmetry=1, barrier=(14.3356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150988,'amu*angstrom^2'), symmetry=1, barrier=(7.5151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286525,'amu*angstrom^2'), symmetry=1, barrier=(14.3354,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.552667,0.0663651,-4.67292e-05,1.72963e-08,-2.615e-12,54616.7,30.4848], Tmin=(100,'K'), Tmax=(1546.29,'K')), NASAPolynomial(coeffs=[14.7971,0.029517,-1.09842e-05,1.88521e-09,-1.23373e-13,50211.5,-44.4451], Tmin=(1546.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC=CC=[C]C(84)',
    structure = SMILES('[CH2]CC=CC=[C]C'),
    E0 = (442.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,320.475,320.478],'cm^-1')),
        HinderedRotor(inertia=(0.158764,'amu*angstrom^2'), symmetry=1, barrier=(11.5711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158764,'amu*angstrom^2'), symmetry=1, barrier=(11.5711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158765,'amu*angstrom^2'), symmetry=1, barrier=(11.5711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158764,'amu*angstrom^2'), symmetry=1, barrier=(11.5711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.520891,0.0677225,-4.89514e-05,1.84228e-08,-2.82295e-12,53330.4,29.2799], Tmin=(100,'K'), Tmax=(1528.04,'K')), NASAPolynomial(coeffs=[15.204,0.0292875,-1.12231e-05,1.96295e-09,-1.30076e-13,48843,-47.7835], Tmin=(1528.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C[CH]C=CCC(109)',
    structure = SMILES('[CH]C=CC=CCC'),
    E0 = (369.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.656789,0.0608466,-9.09463e-06,-2.78421e-08,1.40435e-11,44621.2,27.6288], Tmin=(100,'K'), Tmax=(1011.62,'K')), NASAPolynomial(coeffs=[13.8761,0.0363665,-1.40025e-05,2.54775e-09,-1.77656e-13,40524.7,-43.3277], Tmin=(1011.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C[C]=CC=C[CH]C(196)',
    structure = SMILES('C[C]=C[CH]C=CC'),
    E0 = (356.265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,367.319,374.44],'cm^-1')),
        HinderedRotor(inertia=(0.111713,'amu*angstrom^2'), symmetry=1, barrier=(10.7695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111173,'amu*angstrom^2'), symmetry=1, barrier=(10.8078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0012105,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.390715,'amu*angstrom^2'), symmetry=1, barrier=(36.5662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15472,0.0522902,-7.39826e-06,-2.1515e-08,1.03951e-11,42960.1,26.7936], Tmin=(100,'K'), Tmax=(1041.48,'K')), NASAPolynomial(coeffs=[11.3464,0.0348705,-1.35968e-05,2.48019e-09,-1.7224e-13,39659,-28.446], Tmin=(1041.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC=C[CH]C(83)',
    structure = SMILES('[CH2]C=CC=C[CH]C'),
    E0 = (258.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,390.692,391.064],'cm^-1')),
        HinderedRotor(inertia=(0.00110601,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257141,'amu*angstrom^2'), symmetry=1, barrier=(27.8738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257475,'amu*angstrom^2'), symmetry=1, barrier=(27.8724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.687607,'amu*angstrom^2'), symmetry=1, barrier=(74.5388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3519.58,'J/mol'), sigma=(6.10993,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.75 K, Pc=35.01 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16033,0.0465168,2.35882e-05,-6.26226e-08,2.70681e-11,31194,25.163], Tmin=(100,'K'), Tmax=(970.337,'K')), NASAPolynomial(coeffs=[14.5417,0.0300013,-1.06227e-05,1.92701e-09,-1.37629e-13,26777.8,-48.3663], Tmin=(970.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Allyl_S)"""),
)

species(
    label = 'CH2(S)(28)',
    structure = SMILES('[CH2]'),
    E0 = (418.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1358.21,2621.43,3089.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19331,-0.00233105,8.15676e-06,-6.62986e-09,1.93233e-12,50366.2,-0.746734], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.13502,0.00289594,-8.16668e-07,1.13573e-10,-6.36263e-15,50504.1,4.06031], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(418.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=CC=CC=C(114)',
    structure = SMILES('C=CC=CC=C'),
    E0 = (145.693,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.09105,'amu*angstrom^2'), symmetry=1, barrier=(25.0855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08569,'amu*angstrom^2'), symmetry=1, barrier=(24.9622,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38428,0.0414412,1.94287e-05,-6.16688e-08,2.86503e-11,17631.7,19.3258], Tmin=(100,'K'), Tmax=(941.216,'K')), NASAPolynomial(coeffs=[16.962,0.0157222,-4.10118e-06,6.95735e-10,-5.26521e-14,12906.2,-64.4098], Tmin=(941.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC=CCC=C(111)',
    structure = SMILES('C=CC=CCC=C'),
    E0 = (140.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.766652,'amu*angstrom^2'), symmetry=1, barrier=(17.6268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.76638,'amu*angstrom^2'), symmetry=1, barrier=(17.6206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766347,'amu*angstrom^2'), symmetry=1, barrier=(17.6198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3438.54,'J/mol'), sigma=(5.83085,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=537.09 K, Pc=39.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.928268,0.0526924,5.52854e-06,-4.52831e-08,2.12066e-11,17024.7,26.4844], Tmin=(100,'K'), Tmax=(980.344,'K')), NASAPolynomial(coeffs=[15.6612,0.0277121,-1.00058e-05,1.83677e-09,-1.3177e-13,12447.8,-52.9117], Tmin=(980.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC=C=CC=CC(202)',
    structure = SMILES('CC=C=CC=CC'),
    E0 = (162.449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.812922,'amu*angstrom^2'), symmetry=1, barrier=(18.6907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.811973,'amu*angstrom^2'), symmetry=1, barrier=(18.6689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.811315,'amu*angstrom^2'), symmetry=1, barrier=(18.6537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.616316,0.0651576,-3.9368e-05,8.05849e-09,7.36507e-13,19667.6,24.6176], Tmin=(100,'K'), Tmax=(1142.69,'K')), NASAPolynomial(coeffs=[13.8171,0.0318313,-1.25326e-05,2.26888e-09,-1.55483e-13,15809.6,-44.5099], Tmin=(1142.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C=CC=CCC(100)',
    structure = SMILES('C=C=CC=CCC'),
    E0 = (175.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.84029,'amu*angstrom^2'), symmetry=1, barrier=(19.3199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.839309,'amu*angstrom^2'), symmetry=1, barrier=(19.2974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.841543,'amu*angstrom^2'), symmetry=1, barrier=(19.3487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3533.3,'J/mol'), sigma=(6.00633,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=551.89 K, Pc=37 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.701239,0.0605125,-1.83838e-05,-1.8941e-08,1.14648e-11,21276.7,25.4137], Tmin=(100,'K'), Tmax=(1007.05,'K')), NASAPolynomial(coeffs=[15.0975,0.0294814,-1.11145e-05,2.03249e-09,-1.43106e-13,17051.1,-50.725], Tmin=(1007.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC=CC[C]C(346)',
    structure = SMILES('C=CC=CC[C]C'),
    E0 = (438.999,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04649,'amu*angstrom^2'), symmetry=1, barrier=(24.0608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04961,'amu*angstrom^2'), symmetry=1, barrier=(24.1326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.30951,'amu*angstrom^2'), symmetry=1, barrier=(53.1003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04808,'amu*angstrom^2'), symmetry=1, barrier=(24.0974,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.606412,0.0708854,-5.20749e-05,2.0297e-08,-3.29473e-12,52924.4,35.7688], Tmin=(100,'K'), Tmax=(1416.96,'K')), NASAPolynomial(coeffs=[13.1676,0.0354261,-1.45378e-05,2.63626e-09,-1.78811e-13,49364.6,-29.2096], Tmin=(1416.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=CC=C[C]CC(347)',
    structure = SMILES('C=CC=C[C]CC'),
    E0 = (366.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,253.383,253.386,253.387,253.389],'cm^-1')),
        HinderedRotor(inertia=(0.332939,'amu*angstrom^2'), symmetry=1, barrier=(15.1691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332926,'amu*angstrom^2'), symmetry=1, barrier=(15.1691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332947,'amu*angstrom^2'), symmetry=1, barrier=(15.1691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332927,'amu*angstrom^2'), symmetry=1, barrier=(15.169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.399047,0.0695791,-4.76447e-05,1.36745e-08,-6.20531e-13,44257.9,28.7816], Tmin=(100,'K'), Tmax=(1143.17,'K')), NASAPolynomial(coeffs=[15.2107,0.0300008,-1.17839e-05,2.13398e-09,-1.46423e-13,40071.1,-48.1584], Tmin=(1143.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CsJ2_singlet-(Cds-Cds-Cds-Cds)Cs_6_ring)"""),
)

species(
    label = 'C=C[C]CC=CC(348)',
    structure = SMILES('C=C[C]CC=CC'),
    E0 = (454.069,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34067,'amu*angstrom^2'), symmetry=1, barrier=(7.83268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39306,0.06489,-3.76961e-05,-9.59337e-09,2.08474e-11,54698.6,34.2517], Tmin=(100,'K'), Tmax=(541.777,'K')), NASAPolynomial(coeffs=[4.62093,0.0503822,-2.33434e-05,4.51086e-09,-3.18863e-13,54212,19.3946], Tmin=(541.777,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=CC[C]C=CC(349)',
    structure = SMILES('C=CC[C]C=CC'),
    E0 = (454.069,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34067,'amu*angstrom^2'), symmetry=1, barrier=(7.83268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39306,0.06489,-3.76961e-05,-9.59337e-09,2.08474e-11,54698.6,34.2517], Tmin=(100,'K'), Tmax=(541.777,'K')), NASAPolynomial(coeffs=[4.62093,0.0503822,-2.33434e-05,4.51086e-09,-3.18863e-13,54212,19.3946], Tmin=(541.777,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C[C]C=CC=CC(350)',
    structure = SMILES('C[C]C=CC=CC'),
    E0 = (354.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.658752,'amu*angstrom^2'), symmetry=1, barrier=(15.146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.657318,'amu*angstrom^2'), symmetry=1, barrier=(15.113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.656346,'amu*angstrom^2'), symmetry=1, barrier=(15.0907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658443,'amu*angstrom^2'), symmetry=1, barrier=(15.1389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.439105,0.0697918,-5.17993e-05,1.99011e-08,-3.10452e-12,42782.6,28.2051], Tmin=(100,'K'), Tmax=(1505.1,'K')), NASAPolynomial(coeffs=[15.6755,0.0292987,-1.14429e-05,2.02545e-09,-1.35313e-13,38196.2,-51.531], Tmin=(1505.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(CsJ2_singlet-(Cds-Cds-Cds-Cds)Cs_6_ring)"""),
)

species(
    label = '[CH]CC=CC=CC(351)',
    structure = SMILES('[CH]CC=CC=CC'),
    E0 = (449.19,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.828894,0.0681591,-5.07171e-05,2.00354e-08,-3.31458e-12,54140.2,25.3726], Tmin=(100,'K'), Tmax=(1382.66,'K')), NASAPolynomial(coeffs=[12.2539,0.0351068,-1.48597e-05,2.74629e-09,-1.88519e-13,50980.8,-33.4481], Tmin=(1382.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=CC1C=CC1C(223)',
    structure = SMILES('C=CC1C=CC1C'),
    E0 = (188.222,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3542.51,'J/mol'), sigma=(6.04563,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.33 K, Pc=36.38 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50606,0.0365431,4.94138e-05,-8.65193e-08,3.45852e-11,22744,21.4164], Tmin=(100,'K'), Tmax=(977.442,'K')), NASAPolynomial(coeffs=[13.8527,0.0310902,-1.13885e-05,2.12877e-09,-1.5497e-13,18177.2,-48.8813], Tmin=(977.442,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'CC=CC1C=CC1(215)',
    structure = SMILES('CC=CC1C=CC1'),
    E0 = (173.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3664.15,'J/mol'), sigma=(6.2193,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=572.33 K, Pc=34.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60062,0.0347553,5.19335e-05,-8.71662e-08,3.43053e-11,20941.7,21.9782], Tmin=(100,'K'), Tmax=(982.245,'K')), NASAPolynomial(coeffs=[13.1106,0.0322671,-1.20461e-05,2.26085e-09,-1.64258e-13,16539.5,-44.2436], Tmin=(982.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: FFCM1(-)"""),
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
    E0 = (273.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (528.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (600.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (591.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (433.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (465.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (466.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (434.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (469.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (492.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (504.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (494.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (414.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (378.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (355.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (412.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (337.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (465.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (442.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (345.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (390.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (378.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (487.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (450.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (394.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (364.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (283.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (564.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (325.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (326.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (361.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (456.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (425.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (471.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (471.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (454.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (466.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (232.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (232.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction283',
    reactants = ['C=CC=CC=CC(99)'],
    products = ['CC1C=CC=CC1(203)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.67327e+12,'s^-1'), n=-0.101958, Ea=(164.285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_5_unsaturated_hexane] for rate rule [linear_1_3_5_hexatriene]
Euclidian distance = 1.0
family: Intra_Diels_alder_monocyclic"""),
)

reaction(
    label = 'reaction359',
    reactants = ['[CH]=CC=CC=C(342)', 'CH3(23)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.81e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 69 used for C_methyl;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction362',
    reactants = ['CH2CHCHCH(119)', 'propen1yl(311)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 89 used for Cd_pri_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction364',
    reactants = ['[CH]=CC=CC(344)', 'C2H3(32)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 89 used for Cd_pri_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction369',
    reactants = ['[CH2]C=CC[C]=CC(207)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction370',
    reactants = ['[CH2]C[C]=CC=CC(90)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction371',
    reactants = ['C=CC=[C]C[CH]C(227)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction372',
    reactants = ['C=[C]CC=C[CH]C(218)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction373',
    reactants = ['[CH2]C=CCC=[C]C(209)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction374',
    reactants = ['[CH2]CC=[C]C=CC(88)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction375',
    reactants = ['C=C[C]=CC[CH]C(228)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.94203e+09,'s^-1'), n=0.37, Ea=(99.6901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad_NDe] + [R3radExo;Y_rad;XH_Rrad_NDe] for rate rule [R3radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction376',
    reactants = ['[CH]=CCC=C[CH]C(220)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction377',
    reactants = ['C=CC[CH]C=[C]C(224)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction378',
    reactants = ['C[CH]C=[C]C=CC(198)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.68987e+09,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction379',
    reactants = ['[CH2]C=[C]C=CCC(89)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction380',
    reactants = ['[CH]=C[CH]CC=CC(225)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction381',
    reactants = ['[CH2]C=CCC=C[CH2](211)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.036e+09,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction382',
    reactants = ['C=[C]C=CC[CH]C(229)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(60.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction383',
    reactants = ['[CH2]CC=C[C]=CC(86)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction384',
    reactants = ['[CH2][CH]C=CCC=C(108)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction385',
    reactants = ['[CH2][C]=CC=CCC(91)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction386',
    reactants = ['C[CH]C=C[C]=CC(197)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction387',
    reactants = ['[CH]=CC=CC[CH]C(230)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction388',
    reactants = ['[CH2]CC=CC=[C]C(84)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction389',
    reactants = ['[CH]=C[CH]C=CCC(109)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction390',
    reactants = ['C[C]=CC=C[CH]C(196)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad_NDe;XH_Rrad] for rate rule [R6radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction235',
    reactants = ['[CH2]C=CC=C[CH]C(83)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction391',
    reactants = ['CH2(S)(28)', 'C=CC=CC=C(114)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction392',
    reactants = ['C=CC=CCC=C(111)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(9.66e+08,'s^-1'), n=1.162, Ea=(185.28,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH2(C)_1;unsaturated_end] for rate rule [1_3_pentadiene;CH2(C)_1;CdH2_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction393',
    reactants = ['CC=C=CC=CC(202)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.53605e+09,'s^-1'), n=1.02346, Ea=(163.761,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH_end;CddC_2] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH3_1;CddC_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction394',
    reactants = ['C=C=CC=CCC(100)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(9.66e+08,'s^-1'), n=1.162, Ea=(185.28,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH2(C)_1;unsaturated_end] for rate rule [1_3_pentadiene;CH2(C)_1;CddC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction395',
    reactants = ['C=CC=CC[C]C(346)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;CsJ2C;CH2(C=C)] for rate rule [CsJ2-C;CsJ2C;CH2(C=C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction396',
    reactants = ['C=CC=C[C]CC(347)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.22315e+12,'s^-1'), n=0.271316, Ea=(58.4033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2(C=C);CH2(C)] + [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction397',
    reactants = ['C=C[C]CC=CC(348)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;CsJ2(C=C);CH2(C=C)] for rate rule [CsJ2-C;CsJ2(C=C);CH2(C=C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction398',
    reactants = ['C=CC[C]C=CC(349)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;CsJ2(C=C);CH2(C=C)] for rate rule [CsJ2-C;CsJ2(C=C);CH2(C=C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction399',
    reactants = ['C[C]C=CC=CC(350)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.08234e+13,'s^-1'), n=0.129132, Ea=(99.579,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction400',
    reactants = ['[CH]CC=CC=CC(351)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;singletcarbene;CH2(C=C)] for rate rule [CsJ2-C;CsJ2H;CH2(C=C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction401',
    reactants = ['C=CC=CC=CC(99)'],
    products = ['C=CC1C=CC1C(223)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction356',
    reactants = ['C=CC=CC=CC(99)'],
    products = ['CC=CC1C=CC1(215)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH2_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

network(
    label = '52',
    isomers = [
        'C=CC=CC=CC(99)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '52',
    Tmin = (300,'K'),
    Tmax = (3000,'K'),
    Tcount = 8,
    Tlist = ([302.617,324.619,374.997,470.374,649.057,1000.02,1706.11,2761.25],'K'),
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

