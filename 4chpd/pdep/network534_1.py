species(
    label = 'CC([O])C1[CH]CC=CC=C1(2280)',
    structure = SMILES('CC([O])C1[CH]CC=CC=C1'),
    E0 = (269.751,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04707,0.0957413,-5.7185e-05,5.53456e-09,4.56276e-12,32638.2,23.2735], Tmin=(100,'K'), Tmax=(1073.77,'K')), NASAPolynomial(coeffs=[21.1771,0.0395547,-1.58582e-05,2.94928e-09,-2.07305e-13,26331.9,-92.6686], Tmin=(1073.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(Cs_S)"""),
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
    label = 'CH3CHO(52)',
    structure = SMILES('CC=O'),
    E0 = (-178.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1305.64,1305.66,1305.67,3976.84],'cm^-1')),
        HinderedRotor(inertia=(0.136163,'amu*angstrom^2'), symmetry=1, barrier=(3.13064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-178.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'CC([O])C1C=CC=CCC=1(2855)',
    structure = SMILES('CC([O])C1C=CC=CCC=1'),
    E0 = (162.821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06887,0.0897397,-2.63711e-05,-4.10731e-08,2.49812e-11,19785,31.7624], Tmin=(100,'K'), Tmax=(969.228,'K')), NASAPolynomial(coeffs=[26.4536,0.0274107,-9.23421e-06,1.70188e-09,-1.25713e-13,12042.4,-112.577], Tmin=(969.228,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(=O)C1[CH]CC=CC=C1(2856)',
    structure = SMILES('CC(=O)C1[CH]CC=CC=C1'),
    E0 = (104.407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,375,552.5,462.5,1710,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.964595,0.0962746,-7.52212e-05,2.99087e-08,-4.77094e-12,12746.6,21.5838], Tmin=(100,'K'), Tmax=(1488.27,'K')), NASAPolynomial(coeffs=[21.6105,0.0355993,-1.40671e-05,2.51459e-09,-1.69241e-13,6027.1,-96.3041], Tmin=(1488.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CCJCC=O)"""),
)

species(
    label = 'CC([O])C1C=CC=CC=C1(2819)',
    structure = SMILES('CC([O])C1C=CC=CC=C1'),
    E0 = (166.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29987,0.101346,-6.94014e-05,1.05337e-08,5.06903e-12,20174.6,29.3078], Tmin=(100,'K'), Tmax=(1010.44,'K')), NASAPolynomial(coeffs=[23.6217,0.0327445,-1.21783e-05,2.21553e-09,-1.55921e-13,13604,-98.7753], Tmin=(1010.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(CC(C)OJ)"""),
)

species(
    label = 'C[CH][O](2420)',
    structure = SMILES('C[CH][O]'),
    E0 = (157.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1642.51],'cm^-1')),
        HinderedRotor(inertia=(0.123965,'amu*angstrom^2'), symmetry=1, barrier=(2.85019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65562,0.0114444,2.34936e-06,-4.83164e-09,1.17966e-12,18963.9,10.3625], Tmin=(100,'K'), Tmax=(1718.65,'K')), NASAPolynomial(coeffs=[6.06294,0.0136322,-6.35953e-06,1.18407e-09,-7.90642e-14,16985.9,-5.90233], Tmin=(1718.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CCsJOH) + radical(CCOJ)"""),
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
    label = 'CH3(34)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=CC1[CH]CC=CC=C1(1518)',
    structure = SMILES('O=CC1[CH]CC=CC=C1'),
    E0 = (159.278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.156,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.123818,0.0790013,-5.27292e-05,1.22085e-08,7.76167e-13,19315,16.9618], Tmin=(100,'K'), Tmax=(1113.52,'K')), NASAPolynomial(coeffs=[17.7847,0.032054,-1.29051e-05,2.38611e-09,-1.66355e-13,14249,-76.2015], Tmin=(1113.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(CCJCC=O)"""),
)

species(
    label = 'CC([O])[C]1C=CC=CCC1(2857)',
    structure = SMILES('CC([O])C1=C[CH]C=CCC1'),
    E0 = (187.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.498708,0.00727809,0.000280581,-4.14618e-07,1.73361e-10,22793.7,36.0137], Tmin=(100,'K'), Tmax=(909.852,'K')), NASAPolynomial(coeffs=[47.107,-0.0160167,1.95787e-05,-3.99547e-09,2.55499e-13,6795.21,-225.751], Tmin=(909.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(CC(C)OJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'C[C](O)C1[CH]CC=CC=C1(2858)',
    structure = SMILES('C[C](O)C1[CH]CC=CC=C1'),
    E0 = (216.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13543,0.105795,-9.08066e-05,4.0569e-08,-7.33215e-12,26171.6,23.2049], Tmin=(100,'K'), Tmax=(1316.09,'K')), NASAPolynomial(coeffs=[20.2276,0.0408671,-1.68069e-05,3.08488e-09,-2.11904e-13,20548.4,-85.7272], Tmin=(1316.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(C2CsJOH)"""),
)

species(
    label = 'CC([O])C1[CH]C=CC=CC1(2281)',
    structure = SMILES('CC([O])C1C=C[CH]C=CC1'),
    E0 = (194.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.557254,0.00663242,0.000280611,-4.13634e-07,1.72925e-10,23573.3,34.7764], Tmin=(100,'K'), Tmax=(908.77,'K')), NASAPolynomial(coeffs=[46.4148,-0.015159,1.93868e-05,-3.98434e-09,2.55989e-13,7803.55,-222.98], Tmin=(908.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(=O)C1C=C[CH][CH]CC1(2796)',
    structure = SMILES('CC(=O)C1C=C[CH][CH]CC1'),
    E0 = (144.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,375,552.5,462.5,1710,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.52199,0.028412,0.000175159,-2.73397e-07,1.15004e-10,17517.4,34.4582], Tmin=(100,'K'), Tmax=(915.565,'K')), NASAPolynomial(coeffs=[34.0182,0.00685558,6.03697e-06,-1.38925e-09,8.30154e-14,6153.69,-152.749], Tmin=(915.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(RCCJCC) + radical(Allyl_S)"""),
)

species(
    label = 'CC([O])C1[C]=CC=CCC1(2859)',
    structure = SMILES('CC([O])C1[C]=CC=CCC1'),
    E0 = (313.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28476,0.101742,-7.13092e-05,1.76974e-08,9.71955e-13,37853.9,22.0429], Tmin=(100,'K'), Tmax=(1075.6,'K')), NASAPolynomial(coeffs=[21.9879,0.0382793,-1.49996e-05,2.74986e-09,-1.91626e-13,31512.1,-98.1378], Tmin=(1075.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'CC(O)[C]1[CH]CC=CC=C1(2860)',
    structure = SMILES('CC(O)C1[CH]CC=C[CH]C=1'),
    E0 = (98.6787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.540007,0.0059734,0.000283908,-4.18471e-07,1.74922e-10,12058.5,34.951], Tmin=(100,'K'), Tmax=(909.363,'K')), NASAPolynomial(coeffs=[47.2419,-0.0167548,2.00352e-05,-4.08904e-09,2.6209e-13,-3989.34,-227.456], Tmin=(909.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.6787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C(O)C1[CH]CC=CC=C1(2861)',
    structure = SMILES('[CH2]C(O)C1[CH]CC=CC=C1'),
    E0 = (250.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60373,0.107577,-9.11894e-05,3.92272e-08,-6.68948e-12,30400.6,26.1097], Tmin=(100,'K'), Tmax=(1413.76,'K')), NASAPolynomial(coeffs=[24.514,0.03368,-1.27836e-05,2.25392e-09,-1.51275e-13,23015.9,-108.936], Tmin=(1413.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(CJCO)"""),
)

species(
    label = 'CC([O])C1C=CC=[C]CC1(2862)',
    structure = SMILES('CC([O])C1C=CC=[C]CC1'),
    E0 = (313.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28476,0.101742,-7.13092e-05,1.76974e-08,9.71955e-13,37853.9,22.0429], Tmin=(100,'K'), Tmax=(1075.6,'K')), NASAPolynomial(coeffs=[21.9879,0.0382793,-1.49996e-05,2.74986e-09,-1.91626e-13,31512.1,-98.1378], Tmin=(1075.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'CC(O)C1[C]=CC=CC[CH]1(2863)',
    structure = SMILES('CC(O)C1[C]=CC=CC[CH]1'),
    E0 = (277.232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3615,1277.5,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.21063,0.10388,-8.58637e-05,3.64997e-08,-6.23834e-12,33539.6,24.5302], Tmin=(100,'K'), Tmax=(1393.43,'K')), NASAPolynomial(coeffs=[21.5935,0.038418,-1.53948e-05,2.78465e-09,-1.89389e-13,27184.5,-93.0519], Tmin=(1393.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Cds_S) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C([O])C1C=CC=CCC1(2864)',
    structure = SMILES('[CH2]C([O])C1C=CC=CCC1'),
    E0 = (286.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35978,0.101857,-6.48629e-05,6.17465e-09,6.18283e-12,34700.7,22.4695], Tmin=(100,'K'), Tmax=(1023.29,'K')), NASAPolynomial(coeffs=[23.4889,0.0358223,-1.36502e-05,2.50808e-09,-1.76916e-13,27987.1,-105.939], Tmin=(1023.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'CC([O])C1C=[C]C=CCC1(2865)',
    structure = SMILES('CC([O])C1C=[C]C=CCC1'),
    E0 = (274.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33379,0.103764,-7.60391e-05,2.20912e-08,-2.76144e-13,33182.8,21.34], Tmin=(100,'K'), Tmax=(1060.68,'K')), NASAPolynomial(coeffs=[21.4888,0.0391029,-1.48704e-05,2.67277e-09,-1.84071e-13,27137.1,-95.7868], Tmin=(1060.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C1C=C[C]=CCC1(2866)',
    structure = SMILES('CC([O])C1C=C[C]=CCC1'),
    E0 = (274.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33379,0.103764,-7.60391e-05,2.20912e-08,-2.76144e-13,33182.8,21.34], Tmin=(100,'K'), Tmax=(1060.68,'K')), NASAPolynomial(coeffs=[21.4888,0.0391029,-1.48704e-05,2.67277e-09,-1.84071e-13,27137.1,-95.7868], Tmin=(1060.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'CC(O)C1C=C[CH][CH]C=C1(2832)',
    structure = SMILES('CC(O)C1[CH]C=C[CH]C=C1'),
    E0 = (105.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.598442,0.00532904,0.000283933,-4.17482e-07,1.74484e-10,12838.1,33.7141], Tmin=(100,'K'), Tmax=(908.288,'K')), NASAPolynomial(coeffs=[46.5503,-0.0158981,1.98439e-05,-4.07807e-09,2.62592e-13,-2981.25,-224.689], Tmin=(908.288,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC(O)C1[CH]CC=C[C]=C1(2867)',
    structure = SMILES('CC(O)C1[CH]CC=C[C]=C1'),
    E0 = (238.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.22681,0.105555,-8.9576e-05,3.98542e-08,-7.1594e-12,28867,23.7068], Tmin=(100,'K'), Tmax=(1328.54,'K')), NASAPolynomial(coeffs=[20.4824,0.0401927,-1.57791e-05,2.823e-09,-1.91055e-13,23098.7,-87.195], Tmin=(1328.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(Cs_S)"""),
)

species(
    label = 'CC(O)C1[CH]C[C]=CC=C1(2868)',
    structure = SMILES('CC(O)C1[CH]C[C]=CC=C1'),
    E0 = (277.232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3615,1277.5,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.21063,0.10388,-8.58637e-05,3.64997e-08,-6.23834e-12,33539.6,24.5302], Tmin=(100,'K'), Tmax=(1393.43,'K')), NASAPolynomial(coeffs=[21.5935,0.038418,-1.53948e-05,2.78465e-09,-1.89389e-13,27184.5,-93.0519], Tmin=(1393.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Cds_S) + radical(Cs_S)"""),
)

species(
    label = 'CC(O)C1[CH]CC=[C]C=C1(2869)',
    structure = SMILES('CC(O)C1[CH]CC=[C]C=C1'),
    E0 = (238.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.22681,0.105555,-8.9576e-05,3.98542e-08,-7.1594e-12,28867,23.7068], Tmin=(100,'K'), Tmax=(1328.54,'K')), NASAPolynomial(coeffs=[20.4824,0.0401927,-1.57791e-05,2.823e-09,-1.91055e-13,23098.7,-87.195], Tmin=(1328.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(Cs_S)"""),
)

species(
    label = 'CC([O])C1[CH]C2C=CCC21(2870)',
    structure = SMILES('CC([O])C1[CH]C2C=CCC21'),
    E0 = (313.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.389615,0.0590863,3.24858e-05,-7.52961e-08,3.00371e-11,37905.8,31.7565], Tmin=(100,'K'), Tmax=(1026.58,'K')), NASAPolynomial(coeffs=[15.9406,0.0464225,-1.90425e-05,3.6457e-09,-2.62735e-13,32187.4,-55.9766], Tmin=(1026.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(CC(C)OJ) + radical(cyclobutane)"""),
)

species(
    label = 'CC([O])C1C=CC2[CH]CC21(2871)',
    structure = SMILES('CC([O])C1C=CC2[CH]CC21'),
    E0 = (315.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.388655,0.0588232,3.38913e-05,-7.73262e-08,3.09117e-11,38034.8,31.7941], Tmin=(100,'K'), Tmax=(1022.64,'K')), NASAPolynomial(coeffs=[16.1229,0.0460013,-1.87666e-05,3.5906e-09,-2.59005e-13,32269.1,-56.9231], Tmin=(1022.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(CC(C)OJ) + radical(cyclobutane)"""),
)

species(
    label = 'CC1OC2[CH]C1[CH]CC=C2(2602)',
    structure = SMILES('CC1OC2[CH]C1[CH]CC=C2'),
    E0 = (169.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02909,0.09427,-5.10482e-05,1.02092e-09,5.25223e-12,20610.8,7.1225], Tmin=(100,'K'), Tmax=(1118.73,'K')), NASAPolynomial(coeffs=[21.4784,0.0416031,-1.77175e-05,3.37741e-09,-2.39499e-13,13834.6,-111.766], Tmin=(1118.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.748,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Tetrahydrofuran) + ring(Cycloheptene) + ring(Tetrahydrofuran) + radical(Cs_S) + radical(CCJCO)"""),
)

species(
    label = 'CC([O])C1C=C[CH]C2CC21(2872)',
    structure = SMILES('CC([O])C1[CH]C=CC2CC21'),
    E0 = (245.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.333103,0.0562871,5.12054e-05,-1.01839e-07,4.1164e-11,29722.3,28.0652], Tmin=(100,'K'), Tmax=(994.321,'K')), NASAPolynomial(coeffs=[18.5222,0.0421001,-1.63755e-05,3.1332e-09,-2.29295e-13,23189.3,-74.2458], Tmin=(994.321,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + radical(cyclohexene-allyl) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C1C2[CH]C=CCC21(2873)',
    structure = SMILES('CC([O])C1C2C=C[CH]CC21'),
    E0 = (245.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.333022,0.0562881,5.1202e-05,-1.01834e-07,4.1162e-11,29722.3,28.0655], Tmin=(100,'K'), Tmax=(994.33,'K')), NASAPolynomial(coeffs=[18.5226,0.0420996,-1.63752e-05,3.13312e-09,-2.29289e-13,23189.2,-74.2476], Tmin=(994.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + radical(CC(C)OJ) + radical(cyclohexene-allyl)"""),
)

species(
    label = 'CC1OC2[CH]C[CH]C1C=C2(2594)',
    structure = SMILES('CC1OC2[CH]C[CH]C1C=C2'),
    E0 = (176.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.44117,0.105073,-8.05177e-05,2.94667e-08,-4.25808e-12,21458.6,7.67944], Tmin=(100,'K'), Tmax=(1637.28,'K')), NASAPolynomial(coeffs=[27.729,0.0338081,-1.52284e-05,2.88234e-09,-1.98868e-13,11906.6,-147.432], Tmin=(1637.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(oxepane) + ring(36dihydro2hpyran) + radical(CCJCO) + radical(Cs_S)"""),
)

species(
    label = 'CC1OC2C=C[CH]C1[CH]C2(2529)',
    structure = SMILES('CC1OC2[CH]C=CC1[CH]C2'),
    E0 = (94.8868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56369,0.07605,7.2243e-05,-1.80203e-07,8.2252e-11,11655,6.56067], Tmin=(100,'K'), Tmax=(944.571,'K')), NASAPolynomial(coeffs=[43.38,0.00262022,3.22024e-06,-4.71581e-10,6.20152e-15,-2050.21,-235.307], Tmin=(944.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.8868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(Cycloheptane) + ring(Oxane) + radical(Cs_S) + radical(C=CCJCO)"""),
)

species(
    label = 'CC1OC2[CH]C=CC[CH]C12(2617)',
    structure = SMILES('CC1OC2[CH]C=CC[CH]C12'),
    E0 = (87.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.63897,0.108064,-8.31686e-05,3.15088e-08,-4.74922e-12,10691.8,3.12696], Tmin=(100,'K'), Tmax=(1573.11,'K')), NASAPolynomial(coeffs=[26.2675,0.0371062,-1.55092e-05,2.83586e-09,-1.9253e-13,1911.72,-144.149], Tmin=(1573.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(Cycloheptene) + ring(Oxetane) + radical(C=CCJCO) + radical(Cs_S)"""),
)

species(
    label = '[CH]=CC=CC(C=C)C(C)[O](2874)',
    structure = SMILES('[CH]=CC=CC(C=C)C(C)[O]'),
    E0 = (388.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24952,0.108431,-9.74456e-05,4.58726e-08,-8.70274e-12,46879.5,38.7031], Tmin=(100,'K'), Tmax=(1262.7,'K')), NASAPolynomial(coeffs=[20.365,0.0399609,-1.61079e-05,2.929e-09,-2.00446e-13,41420.9,-70.6162], Tmin=(1262.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=CC=CCC=CC(C)[O](2875)',
    structure = SMILES('[CH]=CC=CCC=CC(C)[O]'),
    E0 = (386.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0548,0.0951595,-5.41923e-05,-1.43158e-09,8.17937e-12,46681.3,41.8001], Tmin=(100,'K'), Tmax=(1024.54,'K')), NASAPolynomial(coeffs=[22.3003,0.0357622,-1.37669e-05,2.54519e-09,-1.80135e-13,40227.4,-79.5822], Tmin=(1024.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'CC(O)C1C=CC=CCC=1(2876)',
    structure = SMILES('CC(O)C1C=CC=CCC=1'),
    E0 = (-67.5394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.23032,0.0937183,-3.22273e-05,-3.60985e-08,2.34e-11,-7915.63,32.5337], Tmin=(100,'K'), Tmax=(968.874,'K')), NASAPolynomial(coeffs=[26.5153,0.0292086,-9.8226e-06,1.78991e-09,-1.30739e-13,-15640.6,-112.566], Tmin=(968.874,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.5394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'CC(=O)C1C=CC=CCC1(2798)',
    structure = SMILES('CC(=O)C1C=CC=CCC1'),
    E0 = (-95.4953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.959691,0.0992419,-7.6406e-05,3.02067e-08,-4.84251e-12,-11299,19.062], Tmin=(100,'K'), Tmax=(1465.75,'K')), NASAPolynomial(coeffs=[20.4474,0.0408224,-1.66212e-05,3.01472e-09,-2.04596e-13,-17574.4,-92.3999], Tmin=(1465.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.4953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = 'CC(O)C1C=CC=CC=C1(2838)',
    structure = SMILES('CC(O)C1C=CC=CC=C1'),
    E0 = (-64.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46677,0.105389,-7.54844e-05,1.58014e-08,3.36415e-12,-7525.72,30.0985], Tmin=(100,'K'), Tmax=(1013.04,'K')), NASAPolynomial(coeffs=[23.7087,0.0344998,-1.27422e-05,2.2978e-09,-1.60471e-13,-14089.7,-98.9067], Tmin=(1013.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
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
    label = '[O]CC1[CH]CC=CC=C1(1354)',
    structure = SMILES('[O]CC1[CH]CC=CC=C1'),
    E0 = (301.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4471.47,'J/mol'), sigma=(7.43333,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=698.43 K, Pc=24.7 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.189778,0.0797784,-5.59892e-05,1.98905e-08,-2.90621e-12,36404.1,17.9849], Tmin=(100,'K'), Tmax=(1558.64,'K')), NASAPolynomial(coeffs=[16.1011,0.0389452,-1.6693e-05,3.08284e-09,-2.10382e-13,31444,-65.8403], Tmin=(1558.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(CCOJ)"""),
)

species(
    label = 'CC([O])C1[CH]CCC=CC=1(2877)',
    structure = SMILES('CC([O])C1[CH]C=CCCC=1'),
    E0 = (187.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.498708,0.00727809,0.000280581,-4.14618e-07,1.73361e-10,22793.7,36.0137], Tmin=(100,'K'), Tmax=(909.852,'K')), NASAPolynomial(coeffs=[47.107,-0.0160167,1.95787e-05,-3.99547e-09,2.55499e-13,6795.21,-225.751], Tmin=(909.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])[CH]C1C=CC=CC1(2878)',
    structure = SMILES('CC([O])[CH]C1C=CC=CC1'),
    E0 = (266.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0378215,0.0656404,2.83352e-05,-8.36346e-08,3.64056e-11,32261,34.9142], Tmin=(100,'K'), Tmax=(982.445,'K')), NASAPolynomial(coeffs=[20.4535,0.0377711,-1.39438e-05,2.61886e-09,-1.91103e-13,25553.3,-77.2278], Tmin=(982.445,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1C=CC=CC1C(C)[O](2879)',
    structure = SMILES('[CH2]C1C=CC=CC1C(C)[O]'),
    E0 = (264.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.156544,0.069205,2.06622e-05,-7.87182e-08,3.58921e-11,31928.8,33.6865], Tmin=(100,'K'), Tmax=(964.043,'K')), NASAPolynomial(coeffs=[20.7266,0.0365216,-1.24504e-05,2.24552e-09,-1.61535e-13,25394.7,-79.3039], Tmin=(964.043,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Isobutyl) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC1OC2CC=CC=CC12(2290)',
    structure = SMILES('CC1OC2CC=CC=CC12'),
    E0 = (-127.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.05174,0.165119,-0.000313237,3.26967e-07,-1.27374e-10,-15141.7,-14.4235], Tmin=(100,'K'), Tmax=(838.793,'K')), NASAPolynomial(coeffs=[-3.03582,0.0938719,-5.00246e-05,9.83355e-09,-6.83921e-13,-12305.1,6.07561], Tmin=(838.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(1,3-Cycloheptadiene) + ring(Oxetane)"""),
)

species(
    label = 'CC([O])C1[CH]CC2C=CC21(2880)',
    structure = SMILES('CC([O])C1[CH]CC2C=CC21'),
    E0 = (349.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.438106,0.0573418,3.74712e-05,-8.09184e-08,3.2144e-11,42188.5,30.9532], Tmin=(100,'K'), Tmax=(1021.22,'K')), NASAPolynomial(coeffs=[16.1751,0.0456217,-1.86373e-05,3.57629e-09,-2.58632e-13,36371.3,-58.0438], Tmin=(1021.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(Cs_S)"""),
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
    label = 'C[CH]C1[CH]CC=CC=C1(553)',
    structure = SMILES('C[CH]C1[CH]CC=CC=C1'),
    E0 = (406.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.289159,0.0808398,-3.92104e-05,-1.70626e-09,4.84397e-12,49044.6,20.9484], Tmin=(100,'K'), Tmax=(1127.98,'K')), NASAPolynomial(coeffs=[17.2799,0.041386,-1.7129e-05,3.20135e-09,-2.23915e-13,43627.5,-72.3714], Tmin=(1127.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(Cs_S)"""),
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
    E0 = (269.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (385.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (387.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (398.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (339.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (276.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (327.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (380.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (435.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (369.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (413.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (454.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (345.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (408.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (447.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (378.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (373.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (318.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (318.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (330.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (381.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (428.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (348.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (490.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (400.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (400.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (350.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (271.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (271.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (366.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (296.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (319.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (429.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (413.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (333.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (333.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (294.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (720.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (315.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (429.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (515.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (278.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (392.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (649.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['C7H8(11)', 'CH3CHO(52)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(19)', 'CC([O])C1C=CC=CCC=1(2855)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.11e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.1294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2594 used for Cds-CdCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', 'CC(=O)C1[CH]CC=CC=C1(2856)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.0366254,'m^3/(mol*s)'), n=1.743, Ea=(71.4418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsCs_O;YJ] for rate rule [CO-CsCs_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'CC([O])C1C=CC=CC=C1(2819)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0545849,'m^3/(mol*s)'), n=2.81111, Ea=(21.1569,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 26 used for Cds-CdH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CsH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C7H8(11)', 'C[CH][O](2420)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0147295,'m^3/(mol*s)'), n=2.37207, Ea=(16.7338,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', 'CH3CHO(52)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0403742,'m^3/(mol*s)'), n=2.2105, Ea=(122.124,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3(34)', 'O=CC1[CH]CC=CC=C1(1518)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.61258,'m^3/(mol*s)'), n=1.485, Ea=(32.0285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-HHH] for rate rule [CO-CsH_O;CsJ-HHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CC([O])[C]1C=CC=CCC1(2857)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.26e+10,'s^-1'), n=0.77, Ea=(192.464,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 164 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['C[C](O)C1[CH]CC=CC=C1(2858)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.2802e+09,'s^-1'), n=0.983917, Ea=(99.768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_1H;Cs_H_out_H/(Cd-Cd-Cd)] + [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC(=O)C1C=C[CH][CH]CC1(2796)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(588307,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC([O])C1[C]=CC=CCC1(2859)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 206 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC(O)[C]1[CH]CC=CC=C1(2860)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['[CH2]C(O)C1[CH]CC=CC=C1(2861)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC([O])C1C=CC=[C]CC1(2862)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.37e+07,'s^-1'), n=1.41, Ea=(177.82,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 207 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cd_H_out_doubleC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CC(O)C1[C]=CC=CC[CH]1(2863)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([O])C1C=CC=CCC1(2864)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC([O])C1C=[C]C=CCC1(2865)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC([O])C1C=C[C]=CCC1(2866)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC(O)C1C=C[CH][CH]C=C1(2832)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(61085.3,'s^-1'), n=1.8624, Ea=(60.8997,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC(O)C1[CH]CC=C[C]=C1(2867)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.234e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5H_SSSD;O_rad_out;Cd_H_out_single] for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC(O)C1[CH]C[C]=CC=C1(2868)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R6HJ_2;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC(O)C1[CH]CC=[C]C=C1(2869)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(146928,'s^-1'), n=1.47286, Ea=(78.6904,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6H;O_rad_out;XH_out] + [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', 'C[CH][O](2420)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC([O])C1[CH]C2C=CCC21(2870)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 116 used for Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC([O])C1C=CC2[CH]CC21(2871)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 116 used for Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC1OC2[CH]C1[CH]CC=C2(2602)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(9.9e+10,'s^-1'), n=0.06, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri_HCd;radadd_intra] for rate rule [Rn2c7_beta_short;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC([O])C1C=C[CH]C2CC21(2872)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(1.2552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC([O])C1C2[CH]C=CCC21(2873)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(1.2552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC1OC2[CH]C[CH]C1C=C2(2594)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.54269e+17,'s^-1'), n=-0.500262, Ea=(96.3509,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_cyclic;doublebond_intra_pri_HCd;radadd_intra] for rate rule [Rn2c7_gamma_long;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC1OC2C=C[CH]C1[CH]C2(2529)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_gamma_long;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC1OC2[CH]C=CC[CH]C12(2617)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=CC=CC(C=C)C(C)[O](2874)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=CC=CCC=CC(C)[O](2875)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC(O)C1C=CC=CCC=1(2876)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC(=O)C1C=CC=CCC1(2798)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC(O)C1C=CC=CC=C1(2838)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CH2(S)(40)', '[O]CC1[CH]CC=CC=C1(1354)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC([O])C1[CH]CCC=CC=1(2877)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.125e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;CdHC_2] for rate rule [1_3_pentadiene;CH(CJ)_1;CdHC_2]
Euclidian distance = 1.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC([O])[CH]C1C=CC=CC1(2878)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['[CH2]C1C=CC=CC1C(C)[O](2879)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC1OC2CC=CC=CC12(2290)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC([O])C1[CH]CC2C=CC21(2880)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction44',
    reactants = ['O(20)', 'C[CH]C1[CH]CC=CC=C1(553)'],
    products = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '534',
    isomers = [
        'CC([O])C1[CH]CC=CC=C1(2280)',
    ],
    reactants = [
        ('C7H8(11)', 'CH3CHO(52)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '534',
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

