species(
    label = 'CC([O])C1[CH]C=CCCC1(2274)',
    structure = SMILES('CC([O])C1[CH]C=CCCC1'),
    E0 = (119.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0250325,0.0254093,0.00023072,-3.54621e-07,1.49113e-10,14562.7,34.9913], Tmin=(100,'K'), Tmax=(914.331,'K')), NASAPolynomial(coeffs=[43.2332,-0.00250275,1.21954e-05,-2.56804e-09,1.58832e-13,-73.2337,-206.421], Tmin=(914.331,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'CHPD(1)',
    structure = SMILES('C1C=CCCCC=1'),
    E0 = (72.8504,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3862.82,'J/mol'), sigma=(6.55596,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=603.36 K, Pc=31.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.815077,0.0542245,7.00763e-06,-4.80543e-08,2.21793e-11,8890.61,7.61919], Tmin=(100,'K'), Tmax=(987.481,'K')), NASAPolynomial(coeffs=[16.3155,0.0286278,-1.06044e-05,1.97595e-09,-1.42758e-13,4016.04,-76.1476], Tmin=(987.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.8504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene)"""),
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
    label = 'CC1OC2[CH]CCCC1[CH]2(2416)',
    structure = SMILES('CC1OC2[CH]CCCC1[CH]2'),
    E0 = (62.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53507,0.10528,-6.45471e-05,1.04848e-08,2.71024e-12,7784.96,7.41495], Tmin=(100,'K'), Tmax=(1122.16,'K')), NASAPolynomial(coeffs=[22.4166,0.0454216,-1.86436e-05,3.47856e-09,-2.43515e-13,802.768,-118.057], Tmin=(1122.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + polycyclic(PolycyclicRing) + radical(CCJCO) + radical(CCJCO)"""),
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
    label = 'CC([O])C1=CC=CCCC1(2417)',
    structure = SMILES('CC([O])C1=CC=CCCC1'),
    E0 = (68.704,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.199,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25501,0.0971676,-4.45487e-05,-1.48822e-08,1.32383e-11,8468.48,22.4878], Tmin=(100,'K'), Tmax=(1018.03,'K')), NASAPolynomial(coeffs=[23.3978,0.0386964,-1.49658e-05,2.7908e-09,-1.99055e-13,1459.47,-106.66], Tmin=(1018.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(=O)C1[CH]C=CCCC1(2418)',
    structure = SMILES('CC(=O)C1C=C[CH]CCC1'),
    E0 = (-50.2242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.199,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.24159,0.0282569,0.000197947,-3.07983e-07,1.29616e-10,-5854.57,32.9005], Tmin=(100,'K'), Tmax=(915.862,'K')), NASAPolynomial(coeffs=[38.4327,0.00264568,8.65687e-06,-1.87654e-09,1.13559e-13,-18771.5,-180.32], Tmin=(915.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.2242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Allyl_S)"""),
)

species(
    label = 'CC([O])C1C=C=CCCC1(2419)',
    structure = SMILES('CC([O])C1C=C=CCCC1'),
    E0 = (198.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.199,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.27727,0.0770281,-8.68167e-06,-3.57455e-08,1.67324e-11,24075.7,30.8025], Tmin=(100,'K'), Tmax=(1057.7,'K')), NASAPolynomial(coeffs=[16.8267,0.0486267,-1.98577e-05,3.73011e-09,-2.63683e-13,18428,-62.2675], Tmin=(1057.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CC(C)OJ)"""),
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
    label = '[CH]1C=C[CH]CCC1(93)',
    structure = SMILES('[CH]1C=C[CH]CCC1'),
    E0 = (259.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17636,-0.0217516,0.000290736,-3.97461e-07,1.6188e-10,31302.2,19.7886], Tmin=(100,'K'), Tmax=(909.243,'K')), NASAPolynomial(coeffs=[37.6598,-0.0173238,1.86036e-05,-3.75571e-09,2.40281e-13,18214,-184.512], Tmin=(909.243,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(Allyl_S)"""),
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
    label = 'O=CC1[CH]C=CCCC1(2350)',
    structure = SMILES('O=CC1C=C[CH]CCC1'),
    E0 = (4.64658,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.172,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.956337,0.0122609,0.000217004,-3.22788e-07,1.34625e-10,719.5,28.7439], Tmin=(100,'K'), Tmax=(912.633,'K')), NASAPolynomial(coeffs=[37.3681,-0.00521916,1.2163e-05,-2.53435e-09,1.58871e-13,-11844.8,-176.016], Tmin=(912.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.64658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptene) + radical(Allyl_S)"""),
)

species(
    label = 'CC([O])[C]1CC=CCCC1(2421)',
    structure = SMILES('CC([O])[C]1CC=CCCC1'),
    E0 = (130.871,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0961017,0.0240189,0.000232963,-3.5561e-07,1.49116e-10,15938.4,36.7309], Tmin=(100,'K'), Tmax=(914.72,'K')), NASAPolynomial(coeffs=[42.6905,-0.00166727,1.17646e-05,-2.48311e-09,1.5287e-13,1428.19,-201.686], Tmin=(914.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCJ(C)CO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C[C](O)C1[CH]C=CCCC1(2422)',
    structure = SMILES('C[C](O)C1[CH]C=CCCC1'),
    E0 = (65.6769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.279671,0.0379898,0.000188288,-3.08048e-07,1.32196e-10,8105.63,35.7015], Tmin=(100,'K'), Tmax=(912.048,'K')), NASAPolynomial(coeffs=[41.9156,-0.000622335,1.09408e-05,-2.36363e-09,1.4873e-13,-5682.06,-197.373], Tmin=(912.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.6769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(C2CsJOH) + radical(Allyl_S)"""),
)

species(
    label = 'CC([O])C1C[C]=CCCC1(2423)',
    structure = SMILES('CC([O])C1C[C]=CCCC1'),
    E0 = (216.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.200979,0.035174,0.000197175,-3.17086e-07,1.35247e-10,26200.2,37.3828], Tmin=(100,'K'), Tmax=(913.164,'K')), NASAPolynomial(coeffs=[41.9513,-0.000433311,1.08533e-05,-2.33278e-09,1.45499e-13,12288,-196.172], Tmin=(913.164,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C1[CH]CCC=CC1(2424)',
    structure = SMILES('CC([O])C1[CH]CCC=CC1'),
    E0 = (172.839,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0215685,0.0293796,0.000210431,-3.27875e-07,1.38124e-10,20985.2,38.666], Tmin=(100,'K'), Tmax=(915.937,'K')), NASAPolynomial(coeffs=[41.0796,0.000938728,9.94178e-06,-2.12135e-09,1.28854e-13,7135.55,-190.356], Tmin=(915.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cs_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C[C]([O])C1CC=CCCC1(2425)',
    structure = SMILES('C[C]([O])C1CC=CCCC1'),
    E0 = (154.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,360,370,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.320689,0.0392908,0.000184976,-3.04217e-07,1.30648e-10,18840.8,36.7632], Tmin=(100,'K'), Tmax=(912.672,'K')), NASAPolynomial(coeffs=[41.7807,0.000115746,1.04842e-05,-2.27004e-09,1.42138e-13,5102.47,-195.667], Tmin=(912.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(O)[C]1[CH]C=CCCC1(2426)',
    structure = SMILES('CC(O)[C]1[CH]C=CCCC1'),
    E0 = (41.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.137609,0.022712,0.000236296,-3.59468e-07,1.50677e-10,5203.22,35.6674], Tmin=(100,'K'), Tmax=(914.122,'K')), NASAPolynomial(coeffs=[42.8234,-0.00240188,1.22191e-05,-2.57622e-09,1.59421e-13,-9355.51,-203.38], Tmin=(914.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCJ(C)CO) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C(O)C1[CH]C=CCCC1(2427)',
    structure = SMILES('[CH2]C(O)C1[CH]C=CCCC1'),
    E0 = (100.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.284733,0.0345202,0.000205368,-3.30914e-07,1.41611e-10,12314.1,36.9293], Tmin=(100,'K'), Tmax=(911.511,'K')), NASAPolynomial(coeffs=[44.1139,-0.00449328,1.31455e-05,-2.78088e-09,1.76019e-13,-2253.16,-208.663], Tmin=(911.511,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(CJCO)"""),
)

species(
    label = 'CC([O])C1CC=[C]CCC1(2428)',
    structure = SMILES('CC([O])C1CC=[C]CCC1'),
    E0 = (216.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.200979,0.035174,0.000197175,-3.17086e-07,1.35247e-10,26200.2,37.3828], Tmin=(100,'K'), Tmax=(913.164,'K')), NASAPolynomial(coeffs=[41.9513,-0.000433311,1.08533e-05,-2.33278e-09,1.45499e-13,12288,-196.172], Tmin=(913.164,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'CC(O)C1[CH]C=CCC[CH]1(2429)',
    structure = SMILES('CC(O)C1[CH]C=CCC[CH]1'),
    E0 = (83.5913,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.06265,0.0280781,0.000213743,-3.31703e-07,1.39671e-10,10250,37.604], Tmin=(100,'K'), Tmax=(915.306,'K')), NASAPolynomial(coeffs=[41.2133,0.000202633,1.03972e-05,-2.21468e-09,1.35424e-13,-3648.47,-192.054], Tmin=(915.306,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.5913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cs_S) + radical(Allyl_S)"""),
)

species(
    label = 'CC([O])C1C[CH]CC=CC1(2430)',
    structure = SMILES('CC([O])C1C[CH]CC=CC1'),
    E0 = (172.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.106411,0.0307988,0.000198945,-3.11505e-07,1.31509e-10,20969,38.3699], Tmin=(100,'K'), Tmax=(913.884,'K')), NASAPolynomial(coeffs=[38.7459,0.00424334,8.53009e-06,-1.89898e-09,1.16653e-13,7953.15,-177.136], Tmin=(913.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(RCCJCC) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([O])C1CC=CCCC1(2431)',
    structure = SMILES('[CH2]C([O])C1CC=CCCC1'),
    E0 = (189.886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.326154,0.0358263,0.000202036,-3.27057e-07,1.4005e-10,23049.2,37.9925], Tmin=(100,'K'), Tmax=(912.107,'K')), NASAPolynomial(coeffs=[43.98,-0.0037569,1.26899e-05,-2.68754e-09,1.69447e-13,8530.98,-206.963], Tmin=(912.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'CC([O])C1CC=C[CH]CC1(2432)',
    structure = SMILES('CC([O])C1CC=C[CH]CC1'),
    E0 = (119.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0250325,0.0254093,0.00023072,-3.54621e-07,1.49113e-10,14562.7,34.9913], Tmin=(100,'K'), Tmax=(914.331,'K')), NASAPolynomial(coeffs=[43.2332,-0.00250275,1.21954e-05,-2.56804e-09,1.58832e-13,-73.2337,-206.421], Tmin=(914.331,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(O)C1[CH]C=CC[CH]C1(2433)',
    structure = SMILES('CC(O)C1[CH]C=CC[CH]C1'),
    E0 = (83.5074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.147871,0.0294924,0.000202277,-3.15363e-07,1.33069e-10,10233.9,37.3065], Tmin=(100,'K'), Tmax=(913.22,'K')), NASAPolynomial(coeffs=[38.8792,0.00350804,8.98506e-06,-1.99218e-09,1.23213e-13,-2830.72,-178.833], Tmin=(913.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.5074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(RCCJCC)"""),
)

species(
    label = 'CC(O)C1[CH][C]=CCCC1(2434)',
    structure = SMILES('CC(O)C1[CH][C]=CCCC1'),
    E0 = (126.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.159628,0.0338689,0.000200502,-3.20938e-07,1.36805e-10,15465,36.3199], Tmin=(100,'K'), Tmax=(912.539,'K')), NASAPolynomial(coeffs=[42.085,-0.00116951,1.13088e-05,-2.42611e-09,1.52069e-13,1503.93,-197.871], Tmin=(912.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'CC(O)C1[CH]C=C[CH]CC1(2435)',
    structure = SMILES('CC(O)C1[CH]C=C[CH]CC1'),
    E0 = (30.1618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.066625,0.0241013,0.000234057,-3.58486e-07,1.50677e-10,3827.48,33.9275], Tmin=(100,'K'), Tmax=(913.734,'K')), NASAPolynomial(coeffs=[43.3659,-0.00323715,1.26499e-05,-2.66112e-09,1.65381e-13,-10856.9,-208.115], Tmin=(913.734,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.1618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(Allyl_S)"""),
)

species(
    label = 'CC(O)C1[CH]C=[C]CCC1(2436)',
    structure = SMILES('CC(O)C1[CH]C=[C]CCC1'),
    E0 = (126.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.159628,0.0338689,0.000200502,-3.20938e-07,1.36805e-10,15465,36.3199], Tmin=(100,'K'), Tmax=(912.539,'K')), NASAPolynomial(coeffs=[42.085,-0.00116951,1.13088e-05,-2.42611e-09,1.52069e-13,1503.93,-197.871], Tmin=(912.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'CC([O])C1CCCC2[CH]C21(2437)',
    structure = SMILES('CC([O])C1CCCC2[CH]C21'),
    E0 = (252.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.117578,0.0667388,2.04782e-05,-6.09178e-08,2.40522e-11,30518.1,33.3618], Tmin=(100,'K'), Tmax=(1060.15,'K')), NASAPolynomial(coeffs=[14.7008,0.0544959,-2.27287e-05,4.31564e-09,-3.06708e-13,25022,-49.1839], Tmin=(1060.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + polycyclic(s2_3_6_ane) + radical(CC(C)OJ) + radical(bicyclo[4.1.0]heptane-C3-7)"""),
)

species(
    label = 'CC1OC2[CH][CH]C1CCC2(2438)',
    structure = SMILES('CC1OC2[CH][CH]C1CCC2'),
    E0 = (32.4828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26373,0.10751,-7.86657e-05,2.91011e-08,-4.38716e-12,4102.56,6.29117], Tmin=(100,'K'), Tmax=(1534.55,'K')), NASAPolynomial(coeffs=[21.6496,0.0477831,-2.02835e-05,3.73759e-09,-2.5508e-13,-2929.77,-114.065], Tmin=(1534.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.4828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + polycyclic(PolycyclicRing) + radical(cycloheptane) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CCC=CC=CC(C)[O](2439)',
    structure = SMILES('[CH2]CCC=CC=CC(C)[O]'),
    E0 = (215.599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29314,0.10072,-5.94472e-05,5.19236e-09,5.09112e-12,26134.6,43.5401], Tmin=(100,'K'), Tmax=(1062.72,'K')), NASAPolynomial(coeffs=[21.5378,0.0426392,-1.67818e-05,3.09001e-09,-2.1608e-13,19709.2,-75.3947], Tmin=(1062.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(CC(C)OJ) + radical(RCCJ)"""),
)

species(
    label = 'CC(O)C1=CC=CCCC1(2440)',
    structure = SMILES('CC(O)C1=CC=CCCC1'),
    E0 = (-161.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.42189,0.10121,-5.06281e-05,-9.62244e-09,1.15383e-11,-19231.9,23.2786], Tmin=(100,'K'), Tmax=(1020.36,'K')), NASAPolynomial(coeffs=[23.4869,0.0404484,-1.5528e-05,2.87267e-09,-2.03573e-13,-26235.2,-106.803], Tmin=(1020.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = 'CC(=O)C1CC=CCCC1(2441)',
    structure = SMILES('CC(=O)C1CC=CCCC1'),
    E0 = (-190.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0976717,0.0298817,0.000202943,-3.16647e-07,1.33433e-10,-22683.4,35.1918], Tmin=(100,'K'), Tmax=(915.304,'K')), NASAPolynomial(coeffs=[39.49,0.00299818,8.93774e-06,-1.94852e-09,1.18513e-13,-35979.6,-184.608], Tmin=(915.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-190.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene)"""),
)

species(
    label = 'CC(O)C1C=C=CCCC1(2442)',
    structure = SMILES('CC(O)C1C=C=CCCC1'),
    E0 = (-31.5858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.447708,0.0811094,-1.48788e-05,-3.03604e-08,1.49918e-11,-3624.52,31.6063], Tmin=(100,'K'), Tmax=(1062.87,'K')), NASAPolynomial(coeffs=[16.9522,0.0503194,-2.03867e-05,3.80434e-09,-2.6758e-13,-9282.9,-62.6178], Tmin=(1062.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.5858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
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
    label = '[O]CC1[CH]C=CCCC1(2163)',
    structure = SMILES('[O]CC1[CH]C=CCCC1'),
    E0 = (151.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.18,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4464.45,'J/mol'), sigma=(7.65423,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=697.34 K, Pc=22.59 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930807,0.013352,0.000218015,-3.21517e-07,1.33212e-10,18343.4,30.8934], Tmin=(100,'K'), Tmax=(915.053,'K')), NASAPolynomial(coeffs=[36.0839,7.2699e-05,9.65441e-06,-2.0528e-09,1.25379e-13,6032.6,-167.695], Tmin=(915.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CCOJ) + radical(Allyl_S)"""),
)

species(
    label = 'CC([O])[CH]C1C=CCCC1(2443)',
    structure = SMILES('CC([O])[CH]C1C=CCCC1'),
    E0 = (157.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.141774,0.0663365,4.03129e-05,-9.59858e-08,4.01173e-11,19079.2,36.5339], Tmin=(100,'K'), Tmax=(990.092,'K')), NASAPolynomial(coeffs=[19.9148,0.044727,-1.69694e-05,3.19897e-09,-2.32274e-13,12195.2,-74.7352], Tmin=(990.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = 'CC([O])C1[CH]CCCC=C1(2273)',
    structure = SMILES('CC([O])C1[CH]CCCC=C1'),
    E0 = (173.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0186916,0.0291306,0.000211833,-3.29981e-07,1.39073e-10,21114.2,38.711], Tmin=(100,'K'), Tmax=(915.528,'K')), NASAPolynomial(coeffs=[41.3152,0.000430221,1.02667e-05,-2.18782e-09,1.33512e-13,7193.8,-191.604], Tmin=(915.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cs_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC1OC2C=CCCCC12(2278)',
    structure = SMILES('CC1OC2C=CCCCC12'),
    E0 = (-224.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.40178,0.100608,-4.8216e-05,-7.37211e-09,9.00326e-12,-26814.7,3.40388], Tmin=(100,'K'), Tmax=(1079.92,'K')), NASAPolynomial(coeffs=[22.5717,0.0456698,-1.89374e-05,3.58602e-09,-2.54495e-13,-33967,-123.239], Tmin=(1079.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-224.696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(Cycloheptene) + ring(Oxetane)"""),
)

species(
    label = 'C[CH]C1[CH]C=CCCC1(481)',
    structure = SMILES('C[CH]C1[CH]C=CCCC1'),
    E0 = (256.062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.792043,0.0104933,0.000248272,-3.6056e-07,1.48478e-10,30968.6,32.6277], Tmin=(100,'K'), Tmax=(916.175,'K')), NASAPolynomial(coeffs=[38.7592,0.000250705,1.04164e-05,-2.19991e-09,1.32842e-13,17484.6,-182.839], Tmin=(916.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cs_S) + radical(Allyl_S)"""),
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
    label = 'CC([O])C1C=CC=CCC1(2444)',
    structure = SMILES('CC([O])C1C=CC=CCC1'),
    E0 = (75.2084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.199,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19768,0.0965466,-4.46619e-05,-1.36132e-08,1.2628e-11,9248.15,21.2542], Tmin=(100,'K'), Tmax=(1016.75,'K')), NASAPolynomial(coeffs=[22.6649,0.039621,-1.51953e-05,2.81063e-09,-1.99276e-13,2485.69,-103.658], Tmin=(1016.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.2084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C1C=CC[CH]CC1(2445)',
    structure = SMILES('CC([O])C1C=CC[CH]CC1'),
    E0 = (173.825,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.103615,0.0305487,0.000200352,-3.13619e-07,1.32461e-10,21098.1,38.4146], Tmin=(100,'K'), Tmax=(913.462,'K')), NASAPolynomial(coeffs=[38.9817,0.00373454,8.85517e-06,-1.96548e-09,1.21314e-13,8011.31,-178.386], Tmin=(913.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(RCCJCC) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C1C=[C]CCCC1(2446)',
    structure = SMILES('CC([O])C1C=[C]CCCC1'),
    E0 = (217.209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.204047,0.0349273,0.000198569,-3.19183e-07,1.36192e-10,26329.2,37.4285], Tmin=(100,'K'), Tmax=(912.775,'K')), NASAPolynomial(coeffs=[42.188,-0.000943679,1.11793e-05,-2.3995e-09,1.50179e-13,12345.8,-197.426], Tmin=(912.775,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'CC([O])C1C=CCC[CH]C1(2447)',
    structure = SMILES('CC([O])C1C=CCC[CH]C1'),
    E0 = (173.825,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.103615,0.0305487,0.000200352,-3.13619e-07,1.32461e-10,21098.1,38.4146], Tmin=(100,'K'), Tmax=(913.462,'K')), NASAPolynomial(coeffs=[38.9817,0.00373454,8.85517e-06,-1.96548e-09,1.21314e-13,8011.31,-178.386], Tmin=(913.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(RCCJCC) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C1[C]=CCCCC1(2448)',
    structure = SMILES('CC([O])C1[C]=CCCCC1'),
    E0 = (217.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.204047,0.0349273,0.000198569,-3.19183e-07,1.36192e-10,26329.2,37.4285], Tmin=(100,'K'), Tmax=(912.775,'K')), NASAPolynomial(coeffs=[42.188,-0.000943679,1.11793e-05,-2.3995e-09,1.50179e-13,12345.8,-197.426], Tmin=(912.775,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(O)C1[C]=C[CH]CCC1(2449)',
    structure = SMILES('CC(O)C1[C]=C[CH]CCC1'),
    E0 = (127.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.16239,0.0336183,0.000201911,-3.23055e-07,1.3776e-10,15594.1,36.3645], Tmin=(100,'K'), Tmax=(912.146,'K')), NASAPolynomial(coeffs=[42.321,-0.00167857,1.1634e-05,-2.49264e-09,1.56733e-13,1562.02,-199.121], Tmin=(912.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = 'CC([O])[C]1C=CCCCC1(2450)',
    structure = SMILES('CC([O])C1=C[CH]CCCC1'),
    E0 = (113.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0364544,0.0258059,0.000232098,-3.57728e-07,1.5051e-10,13912.1,36.2739], Tmin=(100,'K'), Tmax=(915.124,'K')), NASAPolynomial(coeffs=[44.1648,-0.00387537,1.27159e-05,-2.64649e-09,1.63071e-13,-1024.94,-210.462], Tmin=(915.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CC(C)OJ) + radical(Allyl_S)"""),
)

species(
    label = 'C[C]([O])C1C=CCCCC1(2451)',
    structure = SMILES('C[C]([O])C1C=CCCCC1'),
    E0 = (155.995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,360,370,350,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.323597,0.0390421,0.000186378,-3.06325e-07,1.31598e-10,18969.8,36.8083], Tmin=(100,'K'), Tmax=(912.271,'K')), NASAPolynomial(coeffs=[42.017,-0.000393993,1.08098e-05,-2.33667e-09,1.4681e-13,5160.4,-196.92], Tmin=(912.271,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C([O])C1C=CCCCC1(2452)',
    structure = SMILES('[CH2]C([O])C1C=CCCCC1'),
    E0 = (190.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.328618,0.0355719,0.00020346,-3.29194e-07,1.41014e-10,23178.3,38.036], Tmin=(100,'K'), Tmax=(911.717,'K')), NASAPolynomial(coeffs=[44.2152,-0.00426464,1.30144e-05,-2.75388e-09,1.74095e-13,8589.37,-208.209], Tmin=(911.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC1OC2[CH][CH]CCCC12(2453)',
    structure = SMILES('CC1OC2[CH][CH]CCCC12'),
    E0 = (32.4828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26373,0.10751,-7.86657e-05,2.91011e-08,-4.38716e-12,4102.56,6.29117], Tmin=(100,'K'), Tmax=(1534.55,'K')), NASAPolynomial(coeffs=[21.6496,0.0477831,-2.02835e-05,3.73759e-09,-2.5508e-13,-2929.77,-114.065], Tmin=(1534.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.4828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + polycyclic(PolycyclicRing) + radical(cycloheptane) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC(C=CC=C)C(C)[O](2454)',
    structure = SMILES('[CH2]CC(C=CC=C)C(C)[O]'),
    E0 = (219.074,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.41669,0.104336,-6.8893e-05,1.42549e-08,2.18707e-12,26556.2,42.5992], Tmin=(100,'K'), Tmax=(1067.13,'K')), NASAPolynomial(coeffs=[21.5811,0.0425933,-1.64883e-05,2.99651e-09,-2.07633e-13,20255.1,-76.3718], Tmin=(1067.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(=O)C1C=CCCCC1(2455)',
    structure = SMILES('CC(=O)C1C=CCCCC1'),
    E0 = (-191.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0420471,0.0334978,0.000188935,-2.99417e-07,1.26607e-10,-22820.1,34.7233], Tmin=(100,'K'), Tmax=(915.895,'K')), NASAPolynomial(coeffs=[38.3611,0.00517971,7.6127e-06,-1.69506e-09,1.01949e-13,-35670.8,-178.612], Tmin=(915.895,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-191.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene)"""),
)

species(
    label = 'CC(O)C1C=CC=CCC1(2456)',
    structure = SMILES('CC(O)C1C=CC=CCC1'),
    E0 = (-155.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.36454,0.100589,-5.07409e-05,-8.35339e-09,1.09276e-11,-18452.2,22.0449], Tmin=(100,'K'), Tmax=(1019.11,'K')), NASAPolynomial(coeffs=[22.7534,0.0413738,-1.57579e-05,2.8926e-09,-2.03802e-13,-25208.7,-103.798], Tmin=(1019.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-155.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = '[CH2]C1C=CC(CC1)C(C)[O](2457)',
    structure = SMILES('[CH2]C1C=CC(CC1)C(C)[O]'),
    E0 = (154.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.262071,0.069926,3.25202e-05,-9.086e-08,3.94872e-11,18747.1,35.3115], Tmin=(100,'K'), Tmax=(973.053,'K')), NASAPolynomial(coeffs=[20.1656,0.0435132,-1.54957e-05,2.83017e-09,-2.03072e-13,12046.7,-76.6848], Tmin=(973.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CC(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = 'CC1OC2C=CC1CCC2(2277)',
    structure = SMILES('CC1OC2C=CC1CCC2'),
    E0 = (-217.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09231,0.114437,-8.71756e-05,3.19934e-08,-4.62782e-12,-25954.2,4.97859], Tmin=(100,'K'), Tmax=(1644.94,'K')), NASAPolynomial(coeffs=[30.3979,0.0354309,-1.51314e-05,2.79524e-09,-1.90271e-13,-36643.1,-167.939], Tmin=(1644.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-217.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(oxepane) + ring(36dihydro2hpyran)"""),
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
    E0 = (119.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (137.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (285.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (233.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (426.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (232.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (136.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (172.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (256.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (285.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (416.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (300.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (279.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (194.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (258.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (377.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (167.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (214.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (246.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (537.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (145.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (278.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (173.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (304.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (416.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (364.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (212.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (242.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (197.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (208.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (144.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (570.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (317.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (333.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (127.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (499.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (288.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (298.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (417.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (301.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (377.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (229.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (228.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (199.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (233.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (168.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (246.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (127.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (137.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (314.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (126.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CHPD(1)', 'CH3CHO(52)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC1OC2[CH]CCCC1[CH]2(2416)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.2, Ea=(17.9912,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R6;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', 'CC([O])C1=CC=CCCC1(2417)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(4.93712,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2570 used for Cds-CsCs_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'CC(=O)C1[CH]C=CCCC1(2418)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0366254,'m^3/(mol*s)'), n=1.743, Ea=(71.4418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsCs_O;YJ] for rate rule [CO-CsCs_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(19)', 'CC([O])C1C=C=CCCC1(2419)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CHPD(1)', 'C[CH][O](2420)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00517591,'m^3/(mol*s)'), n=2.41374, Ea=(2.11244,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3CHO(52)', '[CH]1C=C[CH]CCC1(93)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0403742,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3(34)', 'O=CC1[CH]C=CCCC1(2350)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.61258,'m^3/(mol*s)'), n=1.485, Ea=(32.0285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-HHH] for rate rule [CO-CsH_O;CsJ-HHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC([O])[C]1CC=CCCC1(2421)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.94e+08,'s^-1'), n=1.27, Ea=(125.938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 162 used for R2H_S;C_rad_out_Cs2;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_Cs2;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['C[C](O)C1[CH]C=CCCC1(2422)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CC([O])C1C[C]=CCCC1(2423)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC([O])C1[CH]CCC=CC1(2424)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 219 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[C]([O])C1CC=CCCC1(2425)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC(O)[C]1[CH]C=CCCC1(2426)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['[CH2]C(O)C1[CH]C=CCCC1(2427)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC([O])C1CC=[C]CCC1(2428)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.96e+10,'s^-1'), n=0.83, Ea=(257.734,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 201 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC(O)C1[CH]C=CCC[CH]1(2429)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 323 used for R4H_SSS;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC([O])C1C[CH]CC=CC1(2430)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.502,'s^-1'), n=3.86, Ea=(41.6308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 332 used for R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([O])C1CC=CCCC1(2431)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(42200,'s^-1'), n=1.93, Ea=(56.484,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 86 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC([O])C1CC=C[CH]CC1(2432)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 334 used for R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC(O)C1[CH]C=CC[CH]C1(2433)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8e+10,'s^-1'), n=0, Ea=(25.7316,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""From training reaction 307 used for R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC(O)C1[CH][C]=CCCC1(2434)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC(O)C1[CH]C=C[CH]CC1(2435)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(120.747,'s^-1'), n=2.0375, Ea=(54.2874,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R6H_SSSSS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CC(O)C1[CH]C=[C]CCC1(2436)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.44313e+09,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;XH_out] for rate rule [R6HJ_2;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]1C=C[CH]CCC1(93)', 'C[CH][O](2420)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC([O])C1CCCC2[CH]C21(2437)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_cs] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC1OC2[CH][CH]C1CCC2(2438)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.42658e+11,'s^-1'), n=0.0609759, Ea=(92.8718,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]CCC=CC=CC(C)[O](2439)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 49 used for R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC(O)C1=CC=CCCC1(2440)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC(=O)C1CC=CCCC1(2441)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC(O)C1C=C=CCCC1(2442)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(S)(40)', '[O]CC1[CH]C=CCCC1(2163)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['CC([O])[CH]C1C=CCCC1(2443)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CC([O])C1[CH]CCCC=C1(2273)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC1OC2C=CCCCC12(2278)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C[CH]C1[CH]C=CCCC1(481)', 'O(20)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(19)', 'CC([O])C1C=CC=CCC1(2444)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.35e+08,'cm^3/(mol*s)'), n=1.64, Ea=(1.58992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2557 used for Cds-CsH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CC([O])C1C=CC[CH]CC1(2445)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CC([O])C1C=[C]CCCC1(2446)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CC([O])C1C=CCC[CH]C1(2447)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 219 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC([O])C1[C]=CCCCC1(2448)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.96e+10,'s^-1'), n=0.83, Ea=(257.734,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 201 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CC(O)C1[C]=C[CH]CCC1(2449)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC([O])[C]1C=CCCCC1(2450)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.86e+10,'s^-1'), n=0.58, Ea=(109.579,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;C_rad_out_single;Cs_H_out_Cs2] for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 2.82842712475
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C[C]([O])C1C=CCCCC1(2451)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(376186,'s^-1'), n=1.54152, Ea=(43.2892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_SSMS;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([O])C1C=CCCCC1(2452)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R6H_RSSMS;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC1OC2[CH][CH]CCCC12(2453)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]CC(C=CC=C)C(C)[O](2454)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.2, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 46 used for R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction48',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC(=O)C1C=CCCCC1(2455)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC(O)C1C=CC=CCC1(2456)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_De] for rate rule [R6radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C1C=CC(CC1)C(C)[O](2457)'],
    products = ['CC([O])C1[CH]C=CCCC1(2274)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction51',
    reactants = ['CC([O])C1[CH]C=CCCC1(2274)'],
    products = ['CC1OC2C=CC1CCC2(2277)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '528',
    isomers = [
        'CC([O])C1[CH]C=CCCC1(2274)',
    ],
    reactants = [
        ('CHPD(1)', 'CH3CHO(52)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '528',
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

