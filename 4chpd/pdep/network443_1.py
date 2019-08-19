species(
    label = '[CH2]OC([CH2])C=CC(2190)',
    structure = SMILES('[CH2]OC([CH2])C=CC'),
    E0 = (189.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.416414,0.0804102,-7.85683e-05,4.21963e-08,-9.36582e-12,22917.7,28.5276], Tmin=(100,'K'), Tmax=(1072.09,'K')), NASAPolynomial(coeffs=[12.5228,0.0352419,-1.53729e-05,2.89972e-09,-2.02444e-13,20321.9,-30.7216], Tmin=(1072.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CsJOCC2) + radical(CJC(C)OC)"""),
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
    label = 'C=CC=CC(381)',
    structure = SMILES('C=CC=CC'),
    E0 = (57.8956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.831076,'amu*angstrom^2'), symmetry=1, barrier=(19.1081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833175,'amu*angstrom^2'), symmetry=1, barrier=(19.1563,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3140.68,'J/mol'), sigma=(5.4037,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=490.57 K, Pc=45.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00727,0.0328459,1.55855e-05,-4.25745e-08,1.84259e-11,7044.82,16.9534], Tmin=(100,'K'), Tmax=(972.32,'K')), NASAPolynomial(coeffs=[11.2869,0.0212416,-7.50361e-06,1.3618e-09,-9.72233e-14,3984.25,-34.0139], Tmin=(972.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.8956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]OC1CC1[CH]C(3169)',
    structure = SMILES('[CH2]OC1CC1[CH]C'),
    E0 = (211.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.821486,0.0579303,-1.50416e-05,-1.9264e-08,1.06439e-11,25607.2,27.5088], Tmin=(100,'K'), Tmax=(1038.46,'K')), NASAPolynomial(coeffs=[14.3813,0.0308591,-1.22801e-05,2.29343e-09,-1.62355e-13,21434.4,-44.953], Tmin=(1038.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + ring(Cyclopropane) + radical(Cs_S) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]C1OCC1[CH]C(3170)',
    structure = SMILES('[CH2]C1OCC1[CH]C'),
    E0 = (215.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0724,0.0561098,-1.11713e-05,-2.70651e-08,1.63163e-11,26056.8,25.7071], Tmin=(100,'K'), Tmax=(889.225,'K')), NASAPolynomial(coeffs=[11.7997,0.0314563,-9.39591e-06,1.45137e-09,-9.23886e-14,23215.9,-30.0333], Tmin=(889.225,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CJC(C)OC) + radical(Cs_S)"""),
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
    label = '[CH2]OC(=C)C=CC(3171)',
    structure = SMILES('[CH2]OC(=C)C=CC'),
    E0 = (54.6944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,349.243,349.26,349.332],'cm^-1')),
        HinderedRotor(inertia=(0.170463,'amu*angstrom^2'), symmetry=1, barrier=(14.7466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170287,'amu*angstrom^2'), symmetry=1, barrier=(14.7479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170403,'amu*angstrom^2'), symmetry=1, barrier=(14.7483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17034,'amu*angstrom^2'), symmetry=1, barrier=(14.746,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.497052,0.0660087,-3.67337e-05,-3.59099e-09,7.25214e-12,6714.11,25.737], Tmin=(100,'K'), Tmax=(991.852,'K')), NASAPolynomial(coeffs=[16.8589,0.0237023,-8.56201e-06,1.54259e-09,-1.08467e-13,2303.71,-58.9367], Tmin=(991.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.6944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=COCJ)"""),
)

species(
    label = '[CH2]C=C[CH]C(377)',
    structure = SMILES('[CH2]C=C[CH]C'),
    E0 = (240.064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.0180055,'amu*angstrom^2'), symmetry=1, barrier=(19.7234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34503,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0180001,'amu*angstrom^2'), symmetry=1, barrier=(19.7225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18178,0.0283568,2.70949e-05,-5.14684e-08,2.05693e-11,28948.9,17.5848], Tmin=(100,'K'), Tmax=(990.212,'K')), NASAPolynomial(coeffs=[10.2369,0.0240425,-9.12514e-06,1.70243e-09,-1.22294e-13,25969.9,-28.1844], Tmin=(990.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_S)"""),
)

species(
    label = 'CH3CHCH(79)',
    structure = SMILES('[CH]=CC'),
    E0 = (254.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,180,694.504,1531.97,4000],'cm^-1')),
        HinderedRotor(inertia=(0.102104,'amu*angstrom^2'), symmetry=1, barrier=(2.34757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22802,0.0118521,1.72185e-05,-2.70818e-08,1.02846e-11,30640.6,10.1787], Tmin=(100,'K'), Tmax=(988.71,'K')), NASAPolynomial(coeffs=[5.81111,0.014019,-5.21097e-06,9.49013e-10,-6.67731e-14,29513.2,-5.37262], Tmin=(988.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""propen1yl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2OCHCH2(3035)',
    structure = SMILES('[CH2]OC=C'),
    E0 = (76.6924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,319.986,320.005],'cm^-1')),
        HinderedRotor(inertia=(0.00164603,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.277624,'amu*angstrom^2'), symmetry=1, barrier=(20.1806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59643,0.0244142,9.3695e-07,-1.83877e-08,8.69513e-12,9280.22,14.7231], Tmin=(100,'K'), Tmax=(988.882,'K')), NASAPolynomial(coeffs=[9.07407,0.0132332,-4.8876e-06,8.99482e-10,-6.42037e-14,7264.66,-20.1688], Tmin=(988.882,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.6924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH2OCHCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]O[C](C)C=CC(3167)',
    structure = SMILES('[CH2]OC(C)=C[CH]C'),
    E0 = (84.6984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719827,0.0599303,-1.78758e-05,-1.88477e-08,1.1185e-11,10315.6,26.7296], Tmin=(100,'K'), Tmax=(1018.38,'K')), NASAPolynomial(coeffs=[15.1645,0.0293485,-1.13544e-05,2.1019e-09,-1.48752e-13,6017.35,-49.8797], Tmin=(1018.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.6984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(C=COCJ)"""),
)

species(
    label = '[CH2][C](C=CC)OC(3172)',
    structure = SMILES('[CH2]C(=C[CH]C)OC'),
    E0 = (51.6075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.573917,0.0625207,-1.93693e-05,-2.21943e-08,1.38007e-11,6341.75,25.4559], Tmin=(100,'K'), Tmax=(979.057,'K')), NASAPolynomial(coeffs=[16.5155,0.0267377,-9.50916e-06,1.70782e-09,-1.20403e-13,1813.67,-58.2986], Tmin=(979.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.6075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]OC(C)[C]=CC(3165)',
    structure = SMILES('[CH2]OC(C)[C]=CC'),
    E0 = (216.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,1685,370,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.749135,0.0721615,-5.98504e-05,2.69731e-08,-5.11045e-12,26193.6,28.0484], Tmin=(100,'K'), Tmax=(1223.39,'K')), NASAPolynomial(coeffs=[11.7324,0.0362509,-1.58207e-05,2.97994e-09,-2.0747e-13,23506.2,-27.154], Tmin=(1223.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]C([C]=CC)OC(3173)',
    structure = SMILES('[CH2]C([C]=CC)OC'),
    E0 = (243.744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,1685,370,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.698185,0.0781875,-8.602e-05,5.9459e-08,-1.78251e-11,29429.5,28.5822], Tmin=(100,'K'), Tmax=(792.381,'K')), NASAPolynomial(coeffs=[7.61952,0.0432476,-1.9877e-05,3.80914e-09,-2.67053e-13,28332.6,-3.19847], Tmin=(792.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]OC(C)C=[C]C(3164)',
    structure = SMILES('[CH2]OC(C)C=[C]C'),
    E0 = (216.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,1685,370,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.749135,0.0721615,-5.98504e-05,2.69731e-08,-5.11045e-12,26193.6,28.0484], Tmin=(100,'K'), Tmax=(1223.39,'K')), NASAPolynomial(coeffs=[11.7324,0.0362509,-1.58207e-05,2.97994e-09,-2.0747e-13,23506.2,-27.154], Tmin=(1223.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]C(C=[C]C)OC(3174)',
    structure = SMILES('[CH2]C(C=[C]C)OC'),
    E0 = (243.744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,1685,370,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.698185,0.0781875,-8.602e-05,5.9459e-08,-1.78251e-11,29429.5,28.5822], Tmin=(100,'K'), Tmax=(792.381,'K')), NASAPolynomial(coeffs=[7.61952,0.0432476,-1.9877e-05,3.80914e-09,-2.67053e-13,28332.6,-3.19847], Tmin=(792.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]OC(C)[CH]C=C(2188)',
    structure = SMILES('[CH2]OC(C)[CH]C=C'),
    E0 = (110.567,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3377.38,'J/mol'), sigma=(6.0858,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=527.54 K, Pc=34 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0722922,0.0790799,-6.60602e-05,2.8405e-08,-4.9209e-12,13445.5,26.7373], Tmin=(100,'K'), Tmax=(1373.3,'K')), NASAPolynomial(coeffs=[17.0414,0.0296539,-1.20738e-05,2.19724e-09,-1.49941e-13,8784.79,-60.5116], Tmin=(1373.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CsJOCC2) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C([CH]C=C)OC(3153)',
    structure = SMILES('[CH2]C([CH]C=C)OC'),
    E0 = (137.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.498349,0.0796804,-7.41294e-05,3.84456e-08,-8.38877e-12,16660.4,25.5464], Tmin=(100,'K'), Tmax=(1077.19,'K')), NASAPolynomial(coeffs=[11.4816,0.0388951,-1.73348e-05,3.29541e-09,-2.3086e-13,14294.2,-28.2582], Tmin=(1077.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(C=CCJCO)"""),
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
    label = '[CH2]OC1[CH]C(C)C1(3175)',
    structure = SMILES('[CH2]OC1[CH]C(C)C1'),
    E0 = (212.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1432,0.0445969,3.21539e-05,-7.35533e-08,3.11479e-11,25714.8,26.3062], Tmin=(100,'K'), Tmax=(975.445,'K')), NASAPolynomial(coeffs=[16.1177,0.0275787,-9.93375e-06,1.86196e-09,-1.36664e-13,20681.7,-56.3894], Tmin=(975.445,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + ring(Cyclobutane) + radical(CsJOCC2) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1[CH]C(C)CO1(3176)',
    structure = SMILES('[CH2]C1[CH]C(C)CO1'),
    E0 = (142.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40159,0.0405915,4.31251e-05,-9.05795e-08,4.08234e-11,17260.6,22.2657], Tmin=(100,'K'), Tmax=(896.051,'K')), NASAPolynomial(coeffs=[14.7971,0.0251418,-5.25143e-06,6.47565e-10,-3.97318e-14,13079.6,-50.8243], Tmin=(896.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CJC(C)OC) + radical(CCJCO)"""),
)

species(
    label = 'C=C(C=CC)OC(3177)',
    structure = SMILES('C=C(C=CC)OC'),
    E0 = (-137.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.37697,0.0681495,-3.52039e-05,-6.12083e-09,8.14046e-12,-16374.3,23.6826], Tmin=(100,'K'), Tmax=(994.323,'K')), NASAPolynomial(coeffs=[16.8856,0.0264628,-9.61558e-06,1.73075e-09,-1.21296e-13,-20879.5,-62.0139], Tmin=(994.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-137.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'CH2(S)(40)',
    structure = SMILES('[CH2]'),
    E0 = (418.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1358.21,2621.43,3089.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19331,-0.00233105,8.15676e-06,-6.62986e-09,1.93233e-12,50366.2,-0.746734], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.13502,0.00289594,-8.16668e-07,1.13573e-10,-6.36263e-15,50504.1,4.06031], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(418.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]OC([CH2])C=C(3178)',
    structure = SMILES('[CH2]OC([CH2])C=C'),
    E0 = (225.512,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,349.595,349.642,349.758],'cm^-1')),
        HinderedRotor(inertia=(0.156527,'amu*angstrom^2'), symmetry=1, barrier=(13.5616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00137768,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156453,'amu*angstrom^2'), symmetry=1, barrier=(13.5648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156539,'amu*angstrom^2'), symmetry=1, barrier=(13.5669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.885471,0.06699,-6.57955e-05,3.41638e-08,-7.15714e-12,27236.3,25.0328], Tmin=(100,'K'), Tmax=(1148.24,'K')), NASAPolynomial(coeffs=[13.0988,0.0244437,-1.02153e-05,1.89403e-09,-1.31206e-13,24431.5,-35.5777], Tmin=(1148.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]O[CH]CC=CC(3179)',
    structure = SMILES('[CH2]O[CH]CC=CC'),
    E0 = (176.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.60255,0.0724938,-6.15194e-05,2.8655e-08,-5.53052e-12,21364,29.2543], Tmin=(100,'K'), Tmax=(1221.41,'K')), NASAPolynomial(coeffs=[12.539,0.0334032,-1.35127e-05,2.45216e-09,-1.67306e-13,18448.1,-30.7196], Tmin=(1221.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CsJOCC) + radical(CCsJOCs)"""),
)

species(
    label = 'CC=CC1CCO1(2193)',
    structure = SMILES('CC=CC1CCO1'),
    E0 = (-65.1025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57047,0.0369815,4.84148e-05,-9.10459e-08,3.92787e-11,-7727.25,23.5268], Tmin=(100,'K'), Tmax=(912.253,'K')), NASAPolynomial(coeffs=[13.3865,0.0279983,-7.2344e-06,1.0843e-09,-7.23277e-14,-11665.2,-42.161], Tmin=(912.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.1025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Oxetane)"""),
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
    label = '[CH2]C([O])C=CC(3180)',
    structure = SMILES('[CH2]C([O])C=CC'),
    E0 = (216.596,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,234.486,4000],'cm^-1')),
        HinderedRotor(inertia=(0.483831,'amu*angstrom^2'), symmetry=1, barrier=(18.8849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48383,'amu*angstrom^2'), symmetry=1, barrier=(18.885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.483831,'amu*angstrom^2'), symmetry=1, barrier=(18.885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23576,0.0525332,-2.93684e-05,1.95313e-09,2.57797e-12,26157,25.5176], Tmin=(100,'K'), Tmax=(1087.33,'K')), NASAPolynomial(coeffs=[12.793,0.0241196,-9.62576e-06,1.77652e-09,-1.23939e-13,22810,-35.0406], Tmin=(1087.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]O[CH]C=CC(3181)',
    structure = SMILES('[CH2]OC=C[CH]C'),
    E0 = (126.49,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,466.973,467.005,467.01],'cm^-1')),
        HinderedRotor(inertia=(0.134306,'amu*angstrom^2'), symmetry=1, barrier=(20.7943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000772565,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134418,'amu*angstrom^2'), symmetry=1, barrier=(20.7906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134313,'amu*angstrom^2'), symmetry=1, barrier=(20.7931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37754,0.040857,2.24417e-05,-6.51459e-08,2.98064e-11,15323,22.3632], Tmin=(100,'K'), Tmax=(945.736,'K')), NASAPolynomial(coeffs=[17.312,0.0155718,-4.24331e-06,7.45743e-10,-5.7268e-14,10425.8,-63.5788], Tmin=(945.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COCJ) + radical(Allyl_S)"""),
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
    E0 = (189.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (227.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (247.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (275.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (189.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (339.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (317.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (306.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (368.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (288.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (261.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (276.784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (245.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (245.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (432.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (312.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (247.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (252.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (644.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (348.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (197.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (597.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (507.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OC([CH2])C=CC(2190)'],
    products = ['CH2O(35)', 'C=CC=CC(381)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]OC([CH2])C=CC(2190)'],
    products = ['[CH2]OC1CC1[CH]C(3169)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(38.4928,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 338 used for R4_S_D;doublebond_intra_HNd_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_HNd_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]OC([CH2])C=CC(2190)'],
    products = ['[CH2]C1OCC1[CH]C(3170)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.65e+06,'s^-1'), n=1.3, Ea=(57.7392,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 346 used for R5_SS_D;doublebond_intra_HNd_pri;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_HNd_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', '[CH2]OC(=C)C=CC(3171)'],
    products = ['[CH2]OC([CH2])C=CC(2190)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(87.5179,'m^3/(mol*s)'), n=1.66467, Ea=(9.4546,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;HJ] for rate rule [Cds-CdOs_Cds;HJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2O(35)', '[CH2]C=C[CH]C(377)'],
    products = ['[CH2]OC([CH2])C=CC(2190)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2330,'cm^3/(mol*s)'), n=3.17, Ea=(68.4767,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-HH;YJ] for rate rule [Od_CO-HH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 67.5 to 68.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH3CHCH(79)', 'CH2OCHCH2(3035)'],
    products = ['[CH2]OC([CH2])C=CC(2190)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00977588,'m^3/(mol*s)'), n=2.40996, Ea=(8.29121,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CdsJ-H] for rate rule [Cds-OsH_Cds;CdsJ-H]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]OC([CH2])C=CC(2190)'],
    products = ['[CH2]O[C](C)C=CC(3167)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.08414e+09,'s^-1'), n=1.12201, Ea=(127.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_single;Cs_H_out_Cd] + [R2H_S;C_rad_out_2H;Cs_H_out_OneDe] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]OC([CH2])C=CC(2190)'],
    products = ['[CH2][C](C=CC)OC(3172)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.30996,'s^-1'), n=3.5644, Ea=(117.059,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_2H;XH_out] for rate rule [R3H_SS_O;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]OC(C)[C]=CC(3165)'],
    products = ['[CH2]OC([CH2])C=CC(2190)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([C]=CC)OC(3173)'],
    products = ['[CH2]OC([CH2])C=CC(2190)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]OC(C)C=[C]C(3164)'],
    products = ['[CH2]OC([CH2])C=CC(2190)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C=[C]C)OC(3174)'],
    products = ['[CH2]OC([CH2])C=CC(2190)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]OC([CH2])C=CC(2190)'],
    products = ['[CH2]OC(C)[CH]C=C(2188)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(121000,'s^-1'), n=1.9, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 92 used for R5H_SSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]OC([CH2])C=CC(2190)'],
    products = ['[CH2]C([CH]C=C)OC(3153)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8010,'s^-1'), n=1.94, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 93 used for R6H_RSSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R6H_RSSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C=C[CH]C(377)', '[CH2][O](354)'],
    products = ['[CH2]OC([CH2])C=CC(2190)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.35017e+06,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]OC([CH2])C=CC(2190)'],
    products = ['[CH2]OC1[CH]C(C)C1(3175)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.61e+08,'s^-1'), n=0.96, Ea=(123.01,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]OC([CH2])C=CC(2190)'],
    products = ['[CH2]C1[CH]C(C)CO1(3176)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.94e+07,'s^-1'), n=0.93, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]OC([CH2])C=CC(2190)'],
    products = ['C=C(C=CC)OC(3177)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH2(S)(40)', '[CH2]OC([CH2])C=C(3178)'],
    products = ['[CH2]OC([CH2])C=CC(2190)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.97e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]OC([CH2])C=CC(2190)'],
    products = ['[CH2]O[CH]CC=CC(3179)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]OC([CH2])C=CC(2190)'],
    products = ['CC=CC1CCO1(2193)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(T)(33)', '[CH2]C([O])C=CC(3180)'],
    products = ['[CH2]OC([CH2])C=CC(2190)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(T)(33)', '[CH2]O[CH]C=CC(3181)'],
    products = ['[CH2]OC([CH2])C=CC(2190)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

network(
    label = '443',
    isomers = [
        '[CH2]OC([CH2])C=CC(2190)',
    ],
    reactants = [
        ('CH2O(35)', 'C=CC=CC(381)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '443',
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

