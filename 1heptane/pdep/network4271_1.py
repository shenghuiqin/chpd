species(
    label = 'C=[C]C(C)C[C]=CC(24186)',
    structure = SMILES('C=[C]C(C)C[C]=CC'),
    E0 = (471.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.502649,0.0767529,-5.27539e-05,1.9321e-08,-3.01417e-12,56846.6,32.34], Tmin=(100,'K'), Tmax=(1427.7,'K')), NASAPolynomial(coeffs=[12.0444,0.0444165,-1.87802e-05,3.45701e-09,-2.36293e-13,53550.9,-27.4519], Tmin=(1427.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'CH3CHCCH2(18175)',
    structure = SMILES('C=C=CC'),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.759584,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74635,0.0218189,8.22353e-06,-2.14768e-08,8.55624e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.6,'K')), NASAPolynomial(coeffs=[6.82078,0.0192338,-7.45622e-06,1.36536e-09,-9.53195e-14,16028,-10.4333], Tmin=(1025.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=C=C(C)C[C]=CC(25256)',
    structure = SMILES('C=C=C(C)C[C]=CC'),
    E0 = (392.09,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1685,370,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.507862,'amu*angstrom^2'), symmetry=1, barrier=(11.6767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507735,'amu*angstrom^2'), symmetry=1, barrier=(11.6738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507483,'amu*angstrom^2'), symmetry=1, barrier=(11.668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507383,'amu*angstrom^2'), symmetry=1, barrier=(11.6657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635202,0.0771949,-5.95459e-05,2.57891e-08,-4.86746e-12,47275.8,28.799], Tmin=(100,'K'), Tmax=(1194.08,'K')), NASAPolynomial(coeffs=[9.99112,0.0458534,-2.01742e-05,3.80724e-09,-2.65133e-13,45041.5,-17.9973], Tmin=(1194.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C)C=C=CC(25257)',
    structure = SMILES('C=[C]C(C)C=C=CC'),
    E0 = (397.575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.565879,'amu*angstrom^2'), symmetry=1, barrier=(13.0107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.566296,'amu*angstrom^2'), symmetry=1, barrier=(13.0203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.565576,'amu*angstrom^2'), symmetry=1, barrier=(13.0037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.565655,'amu*angstrom^2'), symmetry=1, barrier=(13.0055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.552205,0.0812487,-7.28532e-05,4.00401e-08,-9.79801e-12,47936.8,27.9436], Tmin=(100,'K'), Tmax=(941.688,'K')), NASAPolynomial(coeffs=[8.0225,0.0495176,-2.23098e-05,4.25848e-09,-2.98801e-13,46529.8,-7.6475], Tmin=(941.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(397.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C)CC#CC(25258)',
    structure = SMILES('C=[C]C(C)CC#CC'),
    E0 = (393.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.123978,0.0753364,-5.47889e-05,2.22056e-08,-3.72716e-12,47500.4,32.6961], Tmin=(100,'K'), Tmax=(1401.96,'K')), NASAPolynomial(coeffs=[14.0211,0.035686,-1.23658e-05,2.03248e-09,-1.29854e-13,43603.8,-39.0448], Tmin=(1401.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(C)C[C]=CC(25259)',
    structure = SMILES('C#CC(C)C[C]=CC'),
    E0 = (398.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,349.895,349.898],'cm^-1')),
        HinderedRotor(inertia=(0.126716,'amu*angstrom^2'), symmetry=1, barrier=(11.0088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126715,'amu*angstrom^2'), symmetry=1, barrier=(11.0088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126717,'amu*angstrom^2'), symmetry=1, barrier=(11.0088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126717,'amu*angstrom^2'), symmetry=1, barrier=(11.0089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.926067,'amu*angstrom^2'), symmetry=1, barrier=(80.4543,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.313192,0.0782057,-6.11937e-05,2.6813e-08,-4.94762e-12,48089.2,30.3601], Tmin=(100,'K'), Tmax=(1260.49,'K')), NASAPolynomial(coeffs=[12.3223,0.0400972,-1.58451e-05,2.82884e-09,-1.90816e-13,45061.7,-30.3573], Tmin=(1260.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH]C(18176)',
    structure = SMILES('[CH2][C]=CC'),
    E0 = (361.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.352622,'amu*angstrom^2'), symmetry=1, barrier=(8.10748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828631,'amu*angstrom^2'), symmetry=1, barrier=(19.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12013e-13,43485.7,14.8304], Tmin=(100,'K'), Tmax=(2065.83,'K')), NASAPolynomial(coeffs=[10.7464,0.014324,-5.20136e-06,8.69079e-10,-5.48385e-14,40045.6,-31.3799], Tmin=(2065.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Allyl_P)"""),
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
    label = 'C=C=CC[C]=CC(25260)',
    structure = SMILES('C=C=CC[C]=CC'),
    E0 = (431.145,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.588706,'amu*angstrom^2'), symmetry=1, barrier=(13.5355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587887,'amu*angstrom^2'), symmetry=1, barrier=(13.5167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.588356,'amu*angstrom^2'), symmetry=1, barrier=(13.5275,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18516,0.0607069,-4.13777e-05,1.47067e-08,-2.19637e-12,51956.5,26.0623], Tmin=(100,'K'), Tmax=(1498.96,'K')), NASAPolynomial(coeffs=[11.3785,0.0335056,-1.41571e-05,2.60013e-09,-1.7718e-13,48900.7,-27.2404], Tmin=(1498.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S)"""),
)

species(
    label = 'C#CCC(C)[C]=C(25261)',
    structure = SMILES('C#CCC(C)[C]=C'),
    E0 = (435.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,1685,370,366.019,368.689],'cm^-1')),
        HinderedRotor(inertia=(0.105158,'amu*angstrom^2'), symmetry=1, barrier=(9.94134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203971,'amu*angstrom^2'), symmetry=1, barrier=(19.3556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102652,'amu*angstrom^2'), symmetry=1, barrier=(9.92666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.770037,'amu*angstrom^2'), symmetry=1, barrier=(74.248,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.638197,0.0660572,-5.17386e-05,2.209e-08,-3.85427e-12,52553.4,27.4274], Tmin=(100,'K'), Tmax=(1361.7,'K')), NASAPolynomial(coeffs=[13.51,0.0282463,-1.00876e-05,1.69839e-09,-1.10497e-13,49047.9,-38.6457], Tmin=(1361.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S)"""),
)

species(
    label = 'C=C[C](C)C[C]=CC(24593)',
    structure = SMILES('C=C[C](C)C[C]=CC'),
    E0 = (365.745,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,360,370,350,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.162713,0.0738924,-4.52069e-05,1.39696e-08,-1.75962e-12,44135.8,31.3778], Tmin=(100,'K'), Tmax=(1807.37,'K')), NASAPolynomial(coeffs=[17.3243,0.0359112,-1.36849e-05,2.34244e-09,-1.51314e-13,37932.4,-61.5741], Tmin=(1807.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C)[CH]C=CC(24759)',
    structure = SMILES('C=[C]C(C)[CH]C=CC'),
    E0 = (374.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.19058,0.0730289,-3.90832e-05,6.06405e-09,1.03838e-12,45233.2,31.9015], Tmin=(100,'K'), Tmax=(1232.37,'K')), NASAPolynomial(coeffs=[14.5694,0.0403487,-1.63344e-05,2.96959e-09,-2.02583e-13,40626.8,-44.7828], Tmin=(1232.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C)CC=[C]C(24758)',
    structure = SMILES('C=[C]C(C)CC=[C]C'),
    E0 = (471.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.502649,0.0767529,-5.27539e-05,1.9321e-08,-3.01417e-12,56846.6,32.34], Tmin=(100,'K'), Tmax=(1427.7,'K')), NASAPolynomial(coeffs=[12.0444,0.0444165,-1.87802e-05,3.45701e-09,-2.36293e-13,53550.9,-27.4519], Tmin=(1427.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C)C[C]=CC(24595)',
    structure = SMILES('[CH]=CC(C)C[C]=CC'),
    E0 = (480.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,272.433,272.435],'cm^-1')),
        HinderedRotor(inertia=(0.219067,'amu*angstrom^2'), symmetry=1, barrier=(11.538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00227136,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21907,'amu*angstrom^2'), symmetry=1, barrier=(11.538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219067,'amu*angstrom^2'), symmetry=1, barrier=(11.5379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219068,'amu*angstrom^2'), symmetry=1, barrier=(11.5379,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.186983,0.0787809,-5.4884e-05,1.99679e-08,-3.01418e-12,57975.2,33.4115], Tmin=(100,'K'), Tmax=(1512.36,'K')), NASAPolynomial(coeffs=[14.9216,0.0398099,-1.62316e-05,2.92948e-09,-1.97653e-13,53518.4,-43.7695], Tmin=(1512.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(480.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(C)[CH][C]=CC(24594)',
    structure = SMILES('C=CC(C)[CH][C]=CC'),
    E0 = (374.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.19058,0.0730289,-3.90832e-05,6.06405e-09,1.03838e-12,45233.2,31.9015], Tmin=(100,'K'), Tmax=(1232.37,'K')), NASAPolynomial(coeffs=[14.5694,0.0403487,-1.63344e-05,2.96959e-09,-2.02583e-13,40626.8,-44.7828], Tmin=(1232.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C=C)C[C]=CC(19769)',
    structure = SMILES('[CH2]C(C=C)C[C]=CC'),
    E0 = (438.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3602.68,'J/mol'), sigma=(6.38391,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.73 K, Pc=31.42 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.174654,0.0781796,-5.51456e-05,2.10209e-08,-3.34942e-12,52923.4,34.8343], Tmin=(100,'K'), Tmax=(1442.65,'K')), NASAPolynomial(coeffs=[13.8865,0.04016,-1.56135e-05,2.75215e-09,-1.83509e-13,48967.2,-36.3424], Tmin=(1442.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][C](C)CC=CC(24760)',
    structure = SMILES('C=[C][C](C)CC=CC'),
    E0 = (365.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.162713,0.0738924,-4.52069e-05,1.39696e-08,-1.75962e-12,44135.8,31.3778], Tmin=(100,'K'), Tmax=(1807.37,'K')), NASAPolynomial(coeffs=[17.3243,0.0359112,-1.36849e-05,2.34244e-09,-1.51314e-13,37932.4,-61.5741], Tmin=(1807.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_T)"""),
)

species(
    label = '[CH2]C=CCC(C)[C]=C(19768)',
    structure = SMILES('[CH2]C=CCC(C)[C]=C'),
    E0 = (385.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3602.68,'J/mol'), sigma=(6.38391,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.73 K, Pc=31.42 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0487822,0.0752307,-3.82888e-05,5.74021e-10,3.94178e-12,46488.5,33.7883], Tmin=(100,'K'), Tmax=(1104.33,'K')), NASAPolynomial(coeffs=[15.4613,0.0389347,-1.55154e-05,2.84005e-09,-1.96475e-13,41893.6,-47.4894], Tmin=(1104.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([C]=C)CC=CC(19775)',
    structure = SMILES('[CH2]C([C]=C)CC=CC'),
    E0 = (438.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.174654,0.0781796,-5.51456e-05,2.10209e-08,-3.34942e-12,52923.4,34.8343], Tmin=(100,'K'), Tmax=(1442.65,'K')), NASAPolynomial(coeffs=[13.8865,0.04016,-1.56135e-05,2.75215e-09,-1.83509e-13,48967.2,-36.3424], Tmin=(1442.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(C)CC=CC(24761)',
    structure = SMILES('[CH]=[C]C(C)CC=CC'),
    E0 = (480.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,272.433,272.435],'cm^-1')),
        HinderedRotor(inertia=(0.219067,'amu*angstrom^2'), symmetry=1, barrier=(11.538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00227136,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21907,'amu*angstrom^2'), symmetry=1, barrier=(11.538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219067,'amu*angstrom^2'), symmetry=1, barrier=(11.5379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219068,'amu*angstrom^2'), symmetry=1, barrier=(11.5379,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.186983,0.0787809,-5.4884e-05,1.99679e-08,-3.01418e-12,57975.2,33.4115], Tmin=(100,'K'), Tmax=(1512.36,'K')), NASAPolynomial(coeffs=[14.9216,0.0398099,-1.62316e-05,2.92948e-09,-1.97653e-13,53518.4,-43.7695], Tmin=(1512.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(480.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C=CC(C)C[C]=[C]C(24596)',
    structure = SMILES('C=CC(C)C[C]=[C]C'),
    E0 = (471.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.502649,0.0767529,-5.27539e-05,1.9321e-08,-3.01417e-12,56846.6,32.34], Tmin=(100,'K'), Tmax=(1427.7,'K')), NASAPolynomial(coeffs=[12.0444,0.0444165,-1.87802e-05,3.45701e-09,-2.36293e-13,53550.9,-27.4519], Tmin=(1427.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]CC(C)C=C(19770)',
    structure = SMILES('[CH2]C=[C]CC(C)C=C'),
    E0 = (385.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0487822,0.0752307,-3.82888e-05,5.74021e-10,3.94178e-12,46488.5,33.7883], Tmin=(100,'K'), Tmax=(1104.33,'K')), NASAPolynomial(coeffs=[15.4613,0.0389347,-1.55154e-05,2.84005e-09,-1.96475e-13,41893.6,-47.4894], Tmin=(1104.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=C=C(C)CC=CC(24763)',
    structure = SMILES('C=C=C(C)CC=CC'),
    E0 = (154.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.23043,0.0776752,-5.19848e-05,1.78821e-08,-2.54405e-12,18691.3,29.7805], Tmin=(100,'K'), Tmax=(1592.49,'K')), NASAPolynomial(coeffs=[15.3251,0.0397606,-1.62723e-05,2.93164e-09,-1.97035e-13,13883.6,-50.0659], Tmin=(1592.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC(C)C=C=CC(24598)',
    structure = SMILES('C=CC(C)C=C=CC'),
    E0 = (159.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.481233,0.078198,-5.4689e-05,2.0522e-08,-3.29718e-12,19336.9,27.699], Tmin=(100,'K'), Tmax=(1385.61,'K')), NASAPolynomial(coeffs=[11.7158,0.0457658,-1.95792e-05,3.62946e-09,-2.49313e-13,16223.6,-30.1654], Tmin=(1385.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = 'C=[C]CC[C]=CC(25262)',
    structure = SMILES('C=[C]CC[C]=CC'),
    E0 = (503.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,258.632,258.634,258.638],'cm^-1')),
        HinderedRotor(inertia=(0.0025201,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193607,'amu*angstrom^2'), symmetry=1, barrier=(9.19049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193604,'amu*angstrom^2'), symmetry=1, barrier=(9.19049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193613,'amu*angstrom^2'), symmetry=1, barrier=(9.19049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31451,0.0626486,-4.30471e-05,1.63142e-08,-2.7426e-12,60633.3,27.4296], Tmin=(100,'K'), Tmax=(1287.49,'K')), NASAPolynomial(coeffs=[8.20129,0.0412524,-1.81191e-05,3.40633e-09,-2.36169e-13,58860,-7.53537], Tmin=(1287.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(503.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC(C)[C]=C(25263)',
    structure = SMILES('C=[C]CC(C)[C]=C'),
    E0 = (507.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,359.963,359.999,360.096],'cm^-1')),
        HinderedRotor(inertia=(0.129103,'amu*angstrom^2'), symmetry=1, barrier=(11.8725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00130112,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129091,'amu*angstrom^2'), symmetry=1, barrier=(11.8723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12911,'amu*angstrom^2'), symmetry=1, barrier=(11.8717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.828695,0.0647983,-4.41964e-05,1.56832e-08,-2.30418e-12,61171.9,29.3746], Tmin=(100,'K'), Tmax=(1552.34,'K')), NASAPolynomial(coeffs=[13.2905,0.0326872,-1.31679e-05,2.35769e-09,-1.58129e-13,57302.9,-36.2265], Tmin=(1552.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=C([CH]C)C[C]=CC(24184)',
    structure = SMILES('[CH2]C(=CC)C[C]=CC'),
    E0 = (366.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1685,370,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,579.702],'cm^-1')),
        HinderedRotor(inertia=(0.147406,'amu*angstrom^2'), symmetry=1, barrier=(3.38916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64226,'amu*angstrom^2'), symmetry=1, barrier=(14.7668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64164,'amu*angstrom^2'), symmetry=1, barrier=(14.7526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.643937,'amu*angstrom^2'), symmetry=1, barrier=(14.8054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145327,'amu*angstrom^2'), symmetry=1, barrier=(3.34136,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.29648,0.0786067,-5.42868e-05,1.96375e-08,-2.97459e-12,44273.2,31.2372], Tmin=(100,'K'), Tmax=(1490.43,'K')), NASAPolynomial(coeffs=[13.9025,0.0420909,-1.75363e-05,3.199e-09,-2.17227e-13,40217.5,-39.8334], Tmin=(1490.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(C)[CH]C[C]=CC(25264)',
    structure = SMILES('[CH2]C(C)=CC[C]=CC'),
    E0 = (366.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.29648,0.0786067,-5.42868e-05,1.96375e-08,-2.97459e-12,44273.2,31.2372], Tmin=(100,'K'), Tmax=(1490.43,'K')), NASAPolynomial(coeffs=[13.9025,0.0420909,-1.75363e-05,3.199e-09,-2.17227e-13,40217.5,-39.8334], Tmin=(1490.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C)C(=C)[CH]C(24183)',
    structure = SMILES('[CH2]C(=CC)C(C)[C]=C'),
    E0 = (369.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,345.333,347.343],'cm^-1')),
        HinderedRotor(inertia=(0.119405,'amu*angstrom^2'), symmetry=1, barrier=(9.93037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281457,'amu*angstrom^2'), symmetry=1, barrier=(24.022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116909,'amu*angstrom^2'), symmetry=1, barrier=(9.94809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117447,'amu*angstrom^2'), symmetry=1, barrier=(9.9744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116555,'amu*angstrom^2'), symmetry=1, barrier=(9.93684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.299693,0.0839308,-6.74533e-05,3.06742e-08,-6.02582e-12,44564.4,29.0122], Tmin=(100,'K'), Tmax=(1163.73,'K')), NASAPolynomial(coeffs=[10.857,0.0476425,-2.06788e-05,3.8782e-09,-2.69295e-13,42107.3,-23.5217], Tmin=(1163.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C1C(=CC)CC1C(25265)',
    structure = SMILES('C=C1C(=CC)CC1C'),
    E0 = (118.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.689924,0.0550304,2.3689e-05,-6.56265e-08,2.77602e-11,14372.8,24.9628], Tmin=(100,'K'), Tmax=(993.204,'K')), NASAPolynomial(coeffs=[15.3775,0.0380508,-1.43595e-05,2.66472e-09,-1.90565e-13,9375.16,-56.2678], Tmin=(993.204,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane)"""),
)

species(
    label = 'H2CC(41)',
    structure = SMILES('[C]=C'),
    E0 = (401.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2480.69,'J/mol'), sigma=(4.48499,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=387.48 K, Pc=62.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28155,0.00697643,-2.38528e-06,-1.21078e-09,9.82042e-13,48319.2,5.92036], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.27807,0.00475623,-1.63007e-06,2.54623e-10,-1.4886e-14,48014,0.639979], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(401.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C[CH]C[C]=CC(24192)',
    structure = SMILES('C[CH]C[C]=CC'),
    E0 = (355.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,272.37,2221.18],'cm^-1')),
        HinderedRotor(inertia=(0.00227236,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148422,'amu*angstrom^2'), symmetry=1, barrier=(7.81357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148424,'amu*angstrom^2'), symmetry=1, barrier=(7.81357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148422,'amu*angstrom^2'), symmetry=1, barrier=(7.81357,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25446,0.047533,-2.24803e-05,4.72442e-09,-3.81653e-13,42864.6,23.5426], Tmin=(100,'K'), Tmax=(2781.8,'K')), NASAPolynomial(coeffs=[19.681,0.022476,-8.96952e-06,1.48664e-09,-9.0685e-14,33168.8,-78.3599], Tmin=(2781.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(C)[C]=C(2594)',
    structure = SMILES('[CH2]C(C)[C]=C'),
    E0 = (394.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0392304,'amu*angstrom^2'), symmetry=1, barrier=(13.128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570934,'amu*angstrom^2'), symmetry=1, barrier=(13.1269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0125848,'amu*angstrom^2'), symmetry=1, barrier=(79.8348,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94624,0.0392084,-1.10826e-05,-1.04207e-08,6.34909e-12,47544.7,21.7038], Tmin=(100,'K'), Tmax=(985.332,'K')), NASAPolynomial(coeffs=[8.76025,0.0246346,-8.82083e-06,1.52962e-09,-1.03281e-13,45566.5,-14.2933], Tmin=(985.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[C]=CC(24199)',
    structure = SMILES('[C]=CC'),
    E0 = (564.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.40488,'amu*angstrom^2'), symmetry=1, barrier=(9.30899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28451,0.0144722,-2.9151e-06,-2.04635e-09,7.91551e-13,67868.8,10.0416], Tmin=(100,'K'), Tmax=(1455.86,'K')), NASAPolynomial(coeffs=[5.67821,0.0121697,-4.94658e-06,9.00486e-10,-6.07653e-14,66718.9,-3.96135], Tmin=(1455.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CdCdJ2_triplet)"""),
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
    E0 = (471.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (615.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (617.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (623.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (623.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (526.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (599.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (602.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (654.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (671.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (705.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (586.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (613.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (623.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (613.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (633.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (515.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (789.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (714.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (614.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (722.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (564.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (564.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (923.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (927.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (566.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (616.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (641.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (479.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (757.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (958.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['CH3CHCCH2(18175)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'C=C=C(C)C[C]=CC(25256)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.51e+07,'cm^3/(mol*s)'), n=1.64, Ea=(11.8407,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2579 used for Cds-CsCs_Ca;HJ
Exact match found for rate rule [Cds-CsCs_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=[C]C(C)C=C=CC(25257)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.23e+08,'cm^3/(mol*s)'), n=1.64, Ea=(8.49352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2566 used for Cds-CsH_Ca;HJ
Exact match found for rate rule [Cds-CsH_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=[C]C(C)CC#CC(25258)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.02e+09,'cm^3/(mol*s)'), n=1.64, Ea=(18.4933,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2702 used for Ct-Cs_Ct-Cs;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#CC(C)C[C]=CC(25259)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH3CHCCH2(18175)', 'C=[C][CH]C(18176)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00472174,'m^3/(mol*s)'), n=2.41, Ea=(20.1294,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3(17)', 'C=C=CC[C]=CC(25260)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(10800,'cm^3/(mol*s)'), n=2.41, Ea=(32.1331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 597 used for Cds-CsH_Ca;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Ca;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3(17)', 'C#CCC(C)[C]=C(25261)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(178000,'cm^3/(mol*s)'), n=2.41, Ea=(30.1248,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2271 used for Ct-H_Ct-Cs;CsJ-HHH
Exact match found for rate rule [Ct-H_Ct-Cs;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['C=C[C](C)C[C]=CC(24593)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_Cs2]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['C=[C]C(C)[CH]C=CC(24759)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]C(C)CC=[C]C(24758)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC(C)C[C]=CC(24595)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['C=CC(C)[CH][C]=CC(24594)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['[CH2]C(C=C)C[C]=CC(19769)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['C=[C][C](C)CC=CC(24760)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['[CH2]C=CCC(C)[C]=C(19768)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['[CH2]C([C]=C)CC=CC(19775)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]C(C)CC=CC(24761)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=CC(C)C[C]=[C]C(24596)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.24775e+08,'s^-1'), n=1.35208, Ea=(243.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_single;Cd_H_out_doubleC] + [RnH;Cd_rad_out_Cs;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cs;Cd_H_out_doubleC]
Euclidian distance = 2.82842712475
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['[CH2]C=[C]CC(C)C=C(19770)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6HJ_3;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C][CH]C(18176)', 'C=[C][CH]C(18176)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['C=C=C(C)CC=CC(24763)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['C=CC(C)C=C=CC(24598)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(S)(23)', 'C=[C]CC[C]=CC(25262)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(S)(23)', 'C=[C]CC(C)[C]=C(25263)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['C=C([CH]C)C[C]=CC(24184)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['C=C(C)[CH]C[C]=CC(25264)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.95888e+10,'s^-1'), n=0.7315, Ea=(144.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HC)CJ;CJ;CH3] + [cCs(-HC)CJ;CdsJ;C] for rate rule [cCs(-HC)CJ;CdsJ;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['C=[C]C(C)C(=C)[CH]C(24183)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C]C(C)C[C]=CC(24186)'],
    products = ['C=C1C(=CC)CC1C(25265)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H2CC(41)', 'C[CH]C[C]=CC(24192)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C)[C]=C(2594)', '[C]=CC(24199)'],
    products = ['C=[C]C(C)C[C]=CC(24186)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4271',
    isomers = [
        'C=[C]C(C)C[C]=CC(24186)',
    ],
    reactants = [
        ('CH3CHCCH2(18175)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4271',
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

