species(
    label = 'CC1[C]=COO1(26431)',
    structure = SMILES('CC1[C]=COO1'),
    E0 = (201.176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3594.27,'J/mol'), sigma=(5.99159,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.42 K, Pc=37.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14932,0.0249527,3.74048e-05,-7.02522e-08,2.94238e-11,24276.8,20.3516], Tmin=(100,'K'), Tmax=(961.705,'K')), NASAPolynomial(coeffs=[14.6463,0.012195,-3.87025e-06,7.66602e-10,-6.16998e-14,20059.4,-48.8808], Tmin=(961.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(Cds_S)"""),
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
    label = 'CC1=C=COO1(26432)',
    structure = SMILES('CC1=C=COO1'),
    E0 = (393.542,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16001,0.0314729,5.40636e-07,-2.23655e-08,1.01967e-11,47406.2,16.2598], Tmin=(100,'K'), Tmax=(1034.52,'K')), NASAPolynomial(coeffs=[11.254,0.0165027,-7.0307e-06,1.38069e-09,-1.01142e-13,44444.1,-33.1444], Tmin=(1034.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene)"""),
)

species(
    label = 'CC1C#COO1(26433)',
    structure = SMILES('CC1C#COO1'),
    E0 = (468.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85658,0.0399862,-1.96103e-05,-6.97782e-10,1.85048e-12,56417.9,16.9705], Tmin=(100,'K'), Tmax=(1307.74,'K')), NASAPolynomial(coeffs=[13.3932,0.0179758,-9.59239e-06,1.9584e-09,-1.41375e-13,52265.3,-46.1227], Tmin=(1307.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtOs) + ring(Cyclopentyne)"""),
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
    label = 'C1=COOC=1(26434)',
    structure = SMILES('C1=COOC=1'),
    E0 = (435.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,750.683,751.338,751.366,751.383,751.481,751.496,751.499,751.568,751.587,751.846,751.98],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79508,0.0126821,3.97973e-05,-6.71799e-08,2.81361e-11,52414.5,11.2805], Tmin=(100,'K'), Tmax=(945.135,'K')), NASAPolynomial(coeffs=[13.4398,0.00265609,1.22671e-07,1.42216e-11,-8.78205e-15,48838.1,-47.7493], Tmin=(945.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene)"""),
)

species(
    label = '[C]1[CH]OOC=1(26438)',
    structure = SMILES('[C]1[CH]OOC=1'),
    E0 = (360.788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,819.85,819.85,819.85,819.85,819.85,819.85,819.85,819.85,819.85,819.85,819.85],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08524,0.00281293,6.91832e-05,-9.84561e-08,3.94476e-11,43441.6,15.165], Tmin=(100,'K'), Tmax=(940.003,'K')), NASAPolynomial(coeffs=[13.8419,0.00124637,1.14098e-06,-1.69689e-10,1.91103e-15,39466.3,-46.4528], Tmin=(940.003,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(360.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(Cds_S) + radical(C=CCJO)"""),
)

species(
    label = 'C[C]1[C]=COO1(26439)',
    structure = SMILES('CC1=[C][CH]OO1'),
    E0 = (318.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.46446,0.0214584,3.03075e-05,-5.39146e-08,2.1511e-11,38432.6,20.0918], Tmin=(100,'K'), Tmax=(992.782,'K')), NASAPolynomial(coeffs=[11.46,0.0154162,-6.19456e-06,1.23909e-09,-9.39135e-14,35158.2,-30.737], Tmin=(992.782,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(12dioxolene) + radical(Cds_S) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]C1[C]=COO1(26440)',
    structure = SMILES('[CH2]C1[C]=COO1'),
    E0 = (415.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07898,0.0277328,2.40244e-05,-5.73753e-08,2.54107e-11,50012,22.5992], Tmin=(100,'K'), Tmax=(956.119,'K')), NASAPolynomial(coeffs=[15.3356,0.00838229,-2.2681e-06,4.57837e-10,-3.94568e-14,45826.6,-49.393], Tmin=(956.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(Cds_S) + radical(CJCOOH)"""),
)

species(
    label = 'CC1[C]=[C]OO1(26441)',
    structure = SMILES('CC1[C]=[C]OO1'),
    E0 = (440.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35324,0.0267531,1.05207e-05,-3.17163e-08,1.34358e-11,53098.1,22.7499], Tmin=(100,'K'), Tmax=(1011.54,'K')), NASAPolynomial(coeffs=[10.7011,0.0159989,-6.53555e-06,1.27611e-09,-9.39657e-14,50270.6,-23.2477], Tmin=(1011.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(C=CJO) + radical(Cds_S)"""),
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
    label = '[C]1=COOC1(26442)',
    structure = SMILES('[C]1=COOC1'),
    E0 = (243.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,180,180,568.465,1028.48,1028.48,1028.48,1028.48,1028.49,1028.49,1028.5,1028.51,2930.72],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93792,0.0131609,3.16359e-05,-5.09516e-08,2.04045e-11,29332.9,15.7895], Tmin=(100,'K'), Tmax=(965.366,'K')), NASAPolynomial(coeffs=[9.78787,0.0102781,-3.50705e-06,6.80132e-10,-5.25119e-14,26822.1,-23.1703], Tmin=(965.366,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(Cds_S)"""),
)

species(
    label = 'C[C]1C=COO1(26435)',
    structure = SMILES('CC1=C[CH]OO1'),
    E0 = (81.1547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49683,0.0168717,5.50241e-05,-8.30537e-08,3.24042e-11,9829.3,19.5015], Tmin=(100,'K'), Tmax=(975.788,'K')), NASAPolynomial(coeffs=[12.6044,0.0159904,-5.95881e-06,1.19982e-09,-9.33835e-14,5926.12,-38.9064], Tmin=(975.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.1547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(12dioxolene) + radical(C=CCJO)"""),
)

species(
    label = 'CC1C=[C]OO1(26436)',
    structure = SMILES('CC1C=[C]OO1'),
    E0 = (203.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38922,0.0221283,3.53449e-05,-6.09514e-08,2.43478e-11,24494.6,22.1464], Tmin=(100,'K'), Tmax=(983.272,'K')), NASAPolynomial(coeffs=[11.8024,0.0166442,-6.33984e-06,1.24614e-09,-9.41975e-14,21057.4,-31.173], Tmin=(983.272,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1C=COO1(26437)',
    structure = SMILES('[CH2]C1C=COO1'),
    E0 = (177.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10907,0.023165,4.87186e-05,-8.65618e-08,3.63638e-11,21408.8,22.0175], Tmin=(100,'K'), Tmax=(950.281,'K')), NASAPolynomial(coeffs=[16.5314,0.00887181,-1.98459e-06,4.07473e-10,-3.80181e-14,16572.1,-57.8533], Tmin=(950.281,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(CJCOOH)"""),
)

species(
    label = 'C#CC(C)O[O](26146)',
    structure = SMILES('C#CC(C)O[O]'),
    E0 = (200.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,492.5,1135,1000,1380,1390,370,380,2900,435],'cm^-1')),
        HinderedRotor(inertia=(0.72109,'amu*angstrom^2'), symmetry=1, barrier=(16.5793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.722427,'amu*angstrom^2'), symmetry=1, barrier=(16.61,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.19122,'amu*angstrom^2'), symmetry=1, barrier=(50.3804,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3669.66,'J/mol'), sigma=(6.11465,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=573.19 K, Pc=36.42 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39394,0.052523,-5.20004e-05,2.72774e-08,-5.66058e-12,24168.1,21.8356], Tmin=(100,'K'), Tmax=(1177.97,'K')), NASAPolynomial(coeffs=[12.1857,0.0158771,-5.33574e-06,8.67386e-10,-5.552e-14,21625.7,-31.9959], Tmin=(1177.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ)"""),
)

species(
    label = 'CC=C=CO[O](26151)',
    structure = SMILES('CC=C=CO[O]'),
    E0 = (214.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.690031,'amu*angstrom^2'), symmetry=1, barrier=(15.8652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691142,'amu*angstrom^2'), symmetry=1, barrier=(15.8907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3805.99,'J/mol'), sigma=(6.15953,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=594.49 K, Pc=36.95 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75181,0.0496222,-4.93503e-05,2.71166e-08,-6.10875e-12,25868.3,20.7585], Tmin=(100,'K'), Tmax=(1064.11,'K')), NASAPolynomial(coeffs=[9.44376,0.0207069,-8.58861e-06,1.57823e-09,-1.08556e-13,24231.4,-16.8283], Tmin=(1064.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
)

species(
    label = 'CC1[CH]OOC=1(26443)',
    structure = SMILES('CC1[CH]OOC=1'),
    E0 = (83.8912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28116,0.0181387,6.36799e-05,-1.00975e-07,4.11328e-11,10169.6,18.3341], Tmin=(100,'K'), Tmax=(950.27,'K')), NASAPolynomial(coeffs=[15.9308,0.0103825,-2.52771e-06,5.10912e-10,-4.58198e-14,5331.46,-58.6283], Tmin=(950.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.8912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(C=CCJO)"""),
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
    E0 = (617.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (695.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (605.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (497.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (530.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (626.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (652.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (663.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (361.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (435.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (353.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (258.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (284.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (345.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', 'CC1=C=COO1(26432)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', 'CC1C#COO1(26433)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH3(17)', 'C1=COOC=1(26434)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.021171,'m^3/(mol*s)'), n=2.41, Ea=(33.7021,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;CsJ-HHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH3(17)', '[C]1[CH]OOC=1(26438)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.68936e+07,'m^3/(mol*s)'), n=0.0549843, Ea=(0.0928142,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', 'C[C]1[C]=COO1(26439)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2]C1[C]=COO1(26440)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', 'CC1[C]=[C]OO1(26441)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(S)(23)', '[C]1=COOC1(26442)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC1[C]=COO1(26431)'],
    products = ['C[C]1C=COO1(26435)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.49816e+10,'s^-1'), n=0.9445, Ea=(160.64,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S_cy5;Cd_rad_out_Cd;Cs_H_out] for rate rule [R2H_S_cy5;Cd_rad_out_Cd;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC1[C]=COO1(26431)'],
    products = ['CC1C=[C]OO1(26436)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CC1[C]=COO1(26431)'],
    products = ['[CH2]C1C=COO1(26437)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R3H_SS_12cy5;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#CC(C)O[O](26146)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra_H;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CC=C=CO[O](26151)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.01304e+12,'s^-1'), n=-0.3725, Ea=(69.6427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_CdCdd;radadd_intra] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CC1[C]=COO1(26431)'],
    products = ['CC1[CH]OOC=1(26443)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.95888e+10,'s^-1'), n=0.7315, Ea=(144.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HR!H)CJ;CJ;CH3] + [cCs(-HR!H)CJ;CdsJ;C] for rate rule [cCs(-HR!H)CJ;CdsJ;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

network(
    label = '4417',
    isomers = [
        'CC1[C]=COO1(26431)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4417',
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

