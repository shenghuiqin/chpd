species(
    label = 'C=CC=CC1CC1(118)',
    structure = SMILES('C=CC=CC1CC1'),
    E0 = (182.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3581.72,'J/mol'), sigma=(6.12627,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.46 K, Pc=35.35 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51944,0.0368195,4.71498e-05,-8.31273e-08,3.313e-11,22089.9,22.5016], Tmin=(100,'K'), Tmax=(980.797,'K')), NASAPolynomial(coeffs=[13.4274,0.0318047,-1.17842e-05,2.20292e-09,-1.59784e-13,17659.3,-45.395], Tmin=(980.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH2]C=CC=CC[CH2](56)',
    structure = SMILES('[CH2]C=CC=CC[CH2]'),
    E0 = (322.527,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,731.024,731.068],'cm^-1')),
        HinderedRotor(inertia=(0.129142,'amu*angstrom^2'), symmetry=1, barrier=(2.96923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.883449,'amu*angstrom^2'), symmetry=1, barrier=(20.3122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.51673,'amu*angstrom^2'), symmetry=1, barrier=(80.8566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238442,'amu*angstrom^2'), symmetry=1, barrier=(20.3117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3519.58,'J/mol'), sigma=(6.10993,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.75 K, Pc=35.01 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.992992,0.0525503,3.92274e-06,-4.15539e-08,1.95233e-11,38911.4,28.6567], Tmin=(100,'K'), Tmax=(978.364,'K')), NASAPolynomial(coeffs=[14.2853,0.0302044,-1.08772e-05,1.96082e-09,-1.38209e-13,34779,-43.0067], Tmin=(978.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(RCCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C7H10(1)',
    structure = SMILES('C1C=CCCCC=1'),
    E0 = (68.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3862.82,'J/mol'), sigma=(6.55596,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=603.36 K, Pc=31.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98716,0.0214517,9.2224e-05,-1.30105e-07,4.97519e-11,8284.97,18.5915], Tmin=(100,'K'), Tmax=(968.934,'K')), NASAPolynomial(coeffs=[13.7239,0.0305102,-1.08313e-05,2.05835e-09,-1.53811e-13,3310.9,-51.592], Tmin=(968.934,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C1=CC2CCCC12(59)',
    structure = SMILES('C1=CC2CCCC12'),
    E0 = (108.528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3725.22,'J/mol'), sigma=(6.44341,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=581.87 K, Pc=31.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30687,0.00902075,0.000132142,-1.74462e-07,6.61652e-11,13139.4,18.0527], Tmin=(100,'K'), Tmax=(958.841,'K')), NASAPolynomial(coeffs=[15.1126,0.0273996,-8.934e-06,1.72278e-09,-1.34162e-13,7383.11,-60.4009], Tmin=(958.841,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH2]C=C[CH]C1CC1(77)',
    structure = SMILES('[CH2]C=C[CH]C1CC1'),
    E0 = (347.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,600,889.952,889.952,889.952,889.952,889.952,889.952,889.952,889.952,889.952,889.952,2323.39],'cm^-1')),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3576.14,'J/mol'), sigma=(6.33006,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=558.58 K, Pc=31.99 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59431,0.0309971,7.21898e-05,-1.12462e-07,4.40737e-11,41940.3,22.5536], Tmin=(100,'K'), Tmax=(973.523,'K')), NASAPolynomial(coeffs=[15.1163,0.0301336,-1.07547e-05,2.04967e-09,-1.5301e-13,36715.6,-55.6314], Tmin=(973.523,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=CC1C=CCC1(113)',
    structure = SMILES('C=CC1C=CCC1'),
    E0 = (76.6775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3658.67,'J/mol'), sigma=(6.25269,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.48 K, Pc=33.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85063,0.02514,8.29919e-05,-1.21728e-07,4.71145e-11,9319.48,21.0388], Tmin=(100,'K'), Tmax=(968.098,'K')), NASAPolynomial(coeffs=[14.1364,0.0299506,-1.05687e-05,1.99729e-09,-1.48769e-13,4336.5,-51.2849], Tmin=(968.098,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.6775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C1=CC(C1)C1CC1(146)',
    structure = SMILES('C1=CC(C1)C1CC1'),
    E0 = (226.142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3725.22,'J/mol'), sigma=(6.44341,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=581.87 K, Pc=31.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00395,0.0190803,0.000102584,-1.43669e-07,5.52787e-11,27292.8,20.706], Tmin=(100,'K'), Tmax=(962.107,'K')), NASAPolynomial(coeffs=[14.8853,0.0282669,-9.55752e-06,1.81747e-09,-1.38108e-13,21910.4,-56.0326], Tmin=(962.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'H(10)',
    structure = SMILES('[H]'),
    E0 = (211.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(211.8,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=CC=CC1[CH]C1(122)',
    structure = SMILES('C=CC=CC1[CH]C1'),
    E0 = (408.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67388,0.0375425,2.84638e-05,-5.63439e-08,2.20483e-11,49251,24.291], Tmin=(100,'K'), Tmax=(1016.48,'K')), NASAPolynomial(coeffs=[11.0806,0.0331227,-1.31166e-05,2.47518e-09,-1.77259e-13,45654.6,-29.5285], Tmin=(1016.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane)"""),
)

species(
    label = 'cC3H5(120)',
    structure = SMILES('[CH]1CC1'),
    E0 = (280.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,894.423,894.423,894.423,894.423,894.423,894.424,894.425,3153.42],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.52414,-0.00309112,7.19158e-05,-9.14814e-08,3.49781e-11,33782.5,9.37886], Tmin=(100,'K'), Tmax=(937.435,'K')), NASAPolynomial(coeffs=[9.02695,0.0081344,-1.57947e-06,2.7853e-10,-2.51237e-14,31225.9,-24.9471], Tmin=(937.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""cC3H5""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=CC=C[C]1CC1(121)',
    structure = SMILES('[CH2]C=CC=C1CC1'),
    E0 = (301.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,382.401,963.489,963.489,963.489,963.489,963.489,963.489,963.489,963.489,963.489,2323.98],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6827,0.0341269,4.75546e-05,-8.18421e-08,3.26229e-11,36384.9,22.7207], Tmin=(100,'K'), Tmax=(974.031,'K')), NASAPolynomial(coeffs=[12.6673,0.0305764,-1.0979e-05,2.02576e-09,-1.46315e-13,32273.6,-40.1047], Tmin=(974.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 2 + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=[C]C1CC1(123)',
    structure = SMILES('C=CC=[C]C1CC1'),
    E0 = (420.63,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,2750,2850,2950,3050,3150,900,950,1000,1050,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4755,0.0415802,2.17632e-05,-5.30737e-08,2.18351e-11,50693.2,23.1284], Tmin=(100,'K'), Tmax=(1002.74,'K')), NASAPolynomial(coeffs=[12.427,0.0309912,-1.19076e-05,2.2295e-09,-1.59823e-13,46833,-38.0332], Tmin=(1002.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S)"""),
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
    label = '[CH]=CC1CC1(124)',
    structure = SMILES('[CH]=CC1CC1'),
    E0 = (374.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180,1058.12,1058.12,1058.12,1058.12,1058.12,1058.12,1058.12,1058.12],'cm^-1')),
        HinderedRotor(inertia=(0.136233,'amu*angstrom^2'), symmetry=1, barrier=(3.13225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54668,0.0169368,5.52948e-05,-8.2328e-08,3.2251e-11,45158.5,16.7784], Tmin=(100,'K'), Tmax=(961.591,'K')), NASAPolynomial(coeffs=[11.2743,0.0184604,-6.09067e-06,1.14085e-09,-8.58852e-14,41731.1,-34.0791], Tmin=(961.591,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_P)"""),
)

species(
    label = 'C=C[C]=CC1CC1(125)',
    structure = SMILES('[CH2]C=C=CC1CC1'),
    E0 = (366.465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,3000,3100,440,815,1455,1000,180,1025.72,1025.72,1025.72,1025.72,1025.72,1025.72,1025.72,1025.72,1025.72,2277.09],'cm^-1')),
        HinderedRotor(inertia=(0.0227815,'amu*angstrom^2'), symmetry=1, barrier=(0.523792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0227815,'amu*angstrom^2'), symmetry=1, barrier=(0.523792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44811,0.0380936,4.20006e-05,-8.03999e-08,3.30047e-11,44183.6,23.891], Tmin=(100,'K'), Tmax=(973.323,'K')), NASAPolynomial(coeffs=[14.8605,0.02731,-9.70767e-06,1.81694e-09,-1.33555e-13,39472.6,-51.2416], Tmin=(973.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C=CC1CC1(126)',
    structure = SMILES('C=C=C[CH]C1CC1'),
    E0 = (371.668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2750,2850,2950,3050,3150,900,950,1000,1050,1100,540,610,2055,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68893,0.0319729,5.6251e-05,-9.26901e-08,3.67102e-11,44801.4,21.1361], Tmin=(100,'K'), Tmax=(975.384,'K')), NASAPolynomial(coeffs=[13.9426,0.0284826,-1.02939e-05,1.94421e-09,-1.43441e-13,40186.6,-49.0781], Tmin=(975.384,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S)"""),
)

species(
    label = '[CH]=CC=CC1CC1(127)',
    structure = SMILES('[CH]=CC=CC1CC1'),
    E0 = (429.884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43888,0.0404631,2.99883e-05,-6.50106e-08,2.6871e-11,51809.4,23.1887], Tmin=(100,'K'), Tmax=(985.713,'K')), NASAPolynomial(coeffs=[13.7748,0.0287936,-1.06723e-05,1.99949e-09,-1.45083e-13,47512.5,-45.6076], Tmin=(985.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_P)"""),
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
    label = 'C=C[CH]C[C]1CC1(128)',
    structure = SMILES('[CH2]C=CC[C]1CC1'),
    E0 = (392.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180,951.39,951.39,951.39,951.39,951.39,951.39,951.39,951.39,951.39,2328.99],'cm^-1')),
        HinderedRotor(inertia=(0.0244741,'amu*angstrom^2'), symmetry=1, barrier=(0.562707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0244741,'amu*angstrom^2'), symmetry=1, barrier=(0.562707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0244741,'amu*angstrom^2'), symmetry=1, barrier=(0.562707,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45483,0.0438962,1.58384e-05,-4.19163e-08,1.63611e-11,47265.2,25.2997], Tmin=(100,'K'), Tmax=(1054.77,'K')), NASAPolynomial(coeffs=[10.3942,0.037909,-1.53433e-05,2.88208e-09,-2.03808e-13,43826.7,-25.6649], Tmin=(1054.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(Tertalkyl)"""),
)

species(
    label = '[CH2]C[C]=CC1CC1(129)',
    structure = SMILES('[CH2]C[C]=CC1CC1'),
    E0 = (501.582,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3000,3100,440,815,1455,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37802,0.0440767,1.96168e-05,-4.95146e-08,1.99397e-11,60432.6,26.4502], Tmin=(100,'K'), Tmax=(1024.34,'K')), NASAPolynomial(coeffs=[11.8359,0.0350798,-1.38344e-05,2.60134e-09,-1.85577e-13,56619.6,-32.4086], Tmin=(1024.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = 'C=[C]C[CH]C1CC1(130)',
    structure = SMILES('C=[C]C[CH]C1CC1'),
    E0 = (507.513,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,2950,3050,3150,900,950,1000,1050,1100,1685,370,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,180,994.478,994.478,994.478,994.478,994.478,994.478,994.478,994.478,994.478,994.478,994.478,2272.73],'cm^-1')),
        HinderedRotor(inertia=(0.0254607,'amu*angstrom^2'), symmetry=1, barrier=(0.585392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0254607,'amu*angstrom^2'), symmetry=1, barrier=(0.585392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0254607,'amu*angstrom^2'), symmetry=1, barrier=(0.585392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59818,0.038358,3.31389e-05,-6.12127e-08,2.34093e-11,61138.7,26.7101], Tmin=(100,'K'), Tmax=(1027.38,'K')), NASAPolynomial(coeffs=[11.1215,0.0361503,-1.45495e-05,2.76889e-09,-1.98857e-13,57341.6,-28.4479], Tmin=(1027.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(Cds_S)"""),
)

species(
    label = 'C=C[CH]CC1[CH]C1(131)',
    structure = SMILES('[CH2]C=CCC1[CH]C1'),
    E0 = (432.628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,260.927,918.958,918.958,918.958,918.958,918.958,918.958,918.958,918.958,918.958,2335.89],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55299,0.0368828,4.48111e-05,-7.75523e-08,3.01751e-11,52136.1,26.157], Tmin=(100,'K'), Tmax=(1000.11,'K')), NASAPolynomial(coeffs=[12.5749,0.0341876,-1.32214e-05,2.51074e-09,-1.82135e-13,47861.7,-37.3665], Tmin=(1000.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(cyclopropane)"""),
)

species(
    label = '[CH2]CC=[C]C1CC1(132)',
    structure = SMILES('[CH2]CC=[C]C1CC1'),
    E0 = (501.582,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3000,3100,440,815,1455,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37802,0.0440767,1.96168e-05,-4.95146e-08,1.99397e-11,60432.6,26.4502], Tmin=(100,'K'), Tmax=(1024.34,'K')), NASAPolynomial(coeffs=[11.8359,0.0350798,-1.38344e-05,2.60134e-09,-1.85577e-13,56619.6,-32.4086], Tmin=(1024.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CC[CH]C1CC1(133)',
    structure = SMILES('[CH]=CC[CH]C1CC1'),
    E0 = (516.767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,383.996,981.78,981.78,981.78,981.78,981.78,981.78,981.78,981.78,981.78,981.78,2287.33],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56103,0.037257,4.12569e-05,-7.29286e-08,2.83113e-11,62254.9,26.7717], Tmin=(100,'K'), Tmax=(1007.19,'K')), NASAPolynomial(coeffs=[12.4171,0.0340365,-1.33605e-05,2.54948e-09,-1.84976e-13,58044.6,-35.7257], Tmin=(1007.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(Cds_P)"""),
)

species(
    label = 'C=CC[CH]C1[CH]C1(134)',
    structure = SMILES('C=CC[CH]C1[CH]C1'),
    E0 = (495.584,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,420.241,965.238,965.238,965.238,965.238,965.238,965.238,965.238,965.238,965.238,965.238,2299.19],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79978,0.0342756,4.00296e-05,-6.47817e-08,2.3774e-11,59696.3,27.8616], Tmin=(100,'K'), Tmax=(1040.6,'K')), NASAPolynomial(coeffs=[9.79994,0.0382433,-1.57378e-05,3.00993e-09,-2.15924e-13,56151.5,-20.0852], Tmin=(1040.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]C=[C]C1CC1(135)',
    structure = SMILES('C[CH]C=[C]C1CC1'),
    E0 = (437.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1685,370,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,180,1075.69,1075.69,1075.69,1075.69,1075.69,1075.69,1075.69,1075.69,1075.69,1075.69,2260.85],'cm^-1')),
        HinderedRotor(inertia=(0.0555157,'amu*angstrom^2'), symmetry=1, barrier=(1.27641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0555157,'amu*angstrom^2'), symmetry=1, barrier=(1.27641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0555157,'amu*angstrom^2'), symmetry=1, barrier=(1.27641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54983,0.0380015,3.937e-05,-7.05994e-08,2.74427e-11,52715,22.9398], Tmin=(100,'K'), Tmax=(1006.66,'K')), NASAPolynomial(coeffs=[12.0116,0.0350086,-1.36538e-05,2.5846e-09,-1.8639e-13,48654.1,-37.3104], Tmin=(1006.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]CC1CC1(136)',
    structure = SMILES('[CH]C=CCC1CC1'),
    E0 = (425.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180,1019.76,1019.76,1019.76,1019.76,1019.76,1019.76,1019.76,1019.76,1019.76,1019.76,1019.76,1019.76,1019.76,2260.76],'cm^-1')),
        HinderedRotor(inertia=(0.0363586,'amu*angstrom^2'), symmetry=1, barrier=(0.835955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0363586,'amu*angstrom^2'), symmetry=1, barrier=(0.835955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0363586,'amu*angstrom^2'), symmetry=1, barrier=(0.835955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38746,0.0387342,5.48027e-05,-9.09516e-08,3.5301e-11,51334.6,25.2315], Tmin=(100,'K'), Tmax=(990.162,'K')), NASAPolynomial(coeffs=[12.7859,0.0386612,-1.47325e-05,2.75758e-09,-1.98487e-13,46823.7,-41.027], Tmin=(990.162,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]CC=C[C]1CC1(137)',
    structure = SMILES('[CH2]CC=C[C]1CC1'),
    E0 = (395.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,600,808.166,808.166,808.166,808.166,808.166,808.166,808.166,808.166,808.166,2400.09],'cm^-1')),
        HinderedRotor(inertia=(0.00402699,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00402699,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00402699,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70057,0.0339617,5.01788e-05,-8.1789e-08,3.15881e-11,47691.7,23.0707], Tmin=(100,'K'), Tmax=(990.897,'K')), NASAPolynomial(coeffs=[11.6108,0.0348631,-1.31091e-05,2.452e-09,-1.76645e-13,43719.5,-34.7832], Tmin=(990.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_T) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]C=C[C]1CC1(138)',
    structure = SMILES('C[CH]C=C[C]1CC1'),
    E0 = (331.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,600,781.68,781.68,781.68,781.68,781.68,781.68,781.68,781.68,781.68,2411.17],'cm^-1')),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86701,0.027942,6.97807e-05,-1.0275e-07,3.90739e-11,39974.4,19.5801], Tmin=(100,'K'), Tmax=(982.589,'K')), NASAPolynomial(coeffs=[11.857,0.0346764,-1.28636e-05,2.42023e-09,-1.76231e-13,35722.9,-40.0847], Tmin=(982.589,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_T) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]CC=CC1[CH]C1(95)',
    structure = SMILES('[CH2]CC=CC1[CH]C1'),
    E0 = (489.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180,937.018,937.018,937.018,937.018,937.018,937.018,937.018,937.018,937.018,2330.46],'cm^-1')),
        HinderedRotor(inertia=(0.0144577,'amu*angstrom^2'), symmetry=1, barrier=(0.332411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0144577,'amu*angstrom^2'), symmetry=1, barrier=(0.332411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0144577,'amu*angstrom^2'), symmetry=1, barrier=(0.332411,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57767,0.0400179,2.64226e-05,-5.29702e-08,2.02549e-11,58990.3,27.6087], Tmin=(100,'K'), Tmax=(1039.21,'K')), NASAPolynomial(coeffs=[10.5198,0.0371633,-1.50171e-05,2.84106e-09,-2.02533e-13,55427.3,-24.0766], Tmin=(1039.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]C=CC1[CH]C1(139)',
    structure = SMILES('CC=C[CH]C1[CH]C1'),
    E0 = (422.241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,600,849.306,849.306,849.306,849.306,849.306,849.306,849.306,849.306,849.306,2373.25],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74635,0.0341948,4.5141e-05,-7.26918e-08,2.7188e-11,50878.4,24.0771], Tmin=(100,'K'), Tmax=(1020.39,'K')), NASAPolynomial(coeffs=[10.5755,0.0373403,-1.49859e-05,2.85396e-09,-2.05378e-13,47111,-28.3284], Tmin=(1020.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane) + radical(Allyl_S)"""),
)

species(
    label = 'CC=CC=C1CC1(140)',
    structure = SMILES('CC=CC=C1CC1'),
    E0 = (183.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3702.39,'J/mol'), sigma=(6.24331,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=578.30 K, Pc=34.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39128,0.0431071,2.41812e-05,-5.59193e-08,2.2691e-11,22194.4,22.215], Tmin=(100,'K'), Tmax=(1004.86,'K')), NASAPolynomial(coeffs=[12.0796,0.0344,-1.33349e-05,2.48346e-09,-1.76812e-13,18337.9,-37.9025], Tmin=(1004.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 2"""),
)

species(
    label = '[CH2]C([CH2])C=CC=C(141)',
    structure = SMILES('[CH2]C([CH2])C=CC=C'),
    E0 = (413.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,358.733,358.734],'cm^-1')),
        HinderedRotor(inertia=(0.00130993,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144167,'amu*angstrom^2'), symmetry=1, barrier=(13.1657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144168,'amu*angstrom^2'), symmetry=1, barrier=(13.1657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.839026,'amu*angstrom^2'), symmetry=1, barrier=(76.6222,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736469,0.0603371,-1.49408e-05,-2.88328e-08,1.77203e-11,49880.2,29.2978], Tmin=(100,'K'), Tmax=(918.931,'K')), NASAPolynomial(coeffs=[15.5206,0.025883,-7.50636e-06,1.18158e-09,-7.82343e-14,45900.7,-47.6463], Tmin=(918.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC[C]C1CC1(142)',
    structure = SMILES('C=CC[C]C1CC1'),
    E0 = (454.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34885,0.041346,3.62282e-05,-7.274e-08,2.96097e-11,54723.7,24.0434], Tmin=(100,'K'), Tmax=(985.472,'K')), NASAPolynomial(coeffs=[13.8297,0.0318437,-1.19541e-05,2.23473e-09,-1.61476e-13,50265.3,-46.1267], Tmin=(985.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C=C[C]CC1CC1(143)',
    structure = SMILES('C=C[C]CC1CC1'),
    E0 = (454.024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35656,0.0417833,3.34874e-05,-6.8772e-08,2.79154e-11,54716.1,24.7341], Tmin=(100,'K'), Tmax=(990.883,'K')), NASAPolynomial(coeffs=[13.4049,0.0326273,-1.24185e-05,2.32423e-09,-1.67399e-13,50390.2,-43.0622], Tmin=(990.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C[C]C=CC1CC1(144)',
    structure = SMILES('C[C]C=CC1CC1'),
    E0 = (426.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30708,0.0445417,2.2869e-05,-5.59119e-08,2.29228e-11,51463.5,23.9517], Tmin=(100,'K'), Tmax=(1004.53,'K')), NASAPolynomial(coeffs=[12.69,0.0339612,-1.3216e-05,2.46963e-09,-1.76319e-13,47423.5,-39.7418], Tmin=(1004.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]CC=CC1CC1(145)',
    structure = SMILES('[CH]CC=CC1CC1'),
    E0 = (508.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23867,0.0445431,2.77777e-05,-6.42868e-08,2.66735e-11,61279.6,23.633], Tmin=(100,'K'), Tmax=(989.552,'K')), NASAPolynomial(coeffs=[13.9527,0.032002,-1.21052e-05,2.2592e-09,-1.6257e-13,56861.1,-47.1828], Tmin=(989.552,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH2]C[CH]C1C=CC1(78)',
    structure = SMILES('[CH2]C[CH]C1C=CC1'),
    E0 = (464.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3025,407.5,1350,352.5,670.258,905.37,905.37,905.37,905.37,905.37,905.37,905.37,905.37,905.37,905.37,905.37,905.37,2427.95,2427.95],'cm^-1')),
        HinderedRotor(inertia=(0.00366254,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00366254,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00366254,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77537,0.0318221,5.39347e-05,-8.39388e-08,3.16566e-11,55915.7,26.4999], Tmin=(100,'K'), Tmax=(1007.23,'K')), NASAPolynomial(coeffs=[11.6271,0.0350252,-1.38707e-05,2.66208e-09,-1.93799e-13,51784.1,-31.7583], Tmin=(1007.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(464.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(Cs_S)"""),
)

species(
    label = '[CH2][CH]C1C=CCC1(79)',
    structure = SMILES('[CH2][CH]C1C=CCC1'),
    E0 = (347.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3000,3100,440,815,1455,1000,991.107,991.107,991.107,991.107,991.107,991.107,991.107,991.107,991.107,991.107,991.107,991.107,991.107,991.107,991.107,991.107,991.107,1186.02,3709.49],'cm^-1')),
        HinderedRotor(inertia=(0.0691127,'amu*angstrom^2'), symmetry=1, barrier=(6.54241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0691127,'amu*angstrom^2'), symmetry=1, barrier=(6.54241,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9666,0.026012,7.01006e-05,-1.00233e-07,3.72857e-11,41926.4,25.393], Tmin=(100,'K'), Tmax=(999.639,'K')), NASAPolynomial(coeffs=[11.4783,0.0349056,-1.37014e-05,2.64358e-09,-1.93954e-13,37678.7,-32.2264], Tmin=(999.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC=CC=C(80)',
    structure = SMILES('[CH2]C=CC=CC=C'),
    E0 = (227.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.25076,'amu*angstrom^2'), symmetry=1, barrier=(28.7575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25232,'amu*angstrom^2'), symmetry=1, barrier=(28.7933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2496,'amu*angstrom^2'), symmetry=1, barrier=(28.7308,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0007,0.0483152,2.14276e-05,-6.83247e-08,3.15713e-11,27512.9,24.0664], Tmin=(100,'K'), Tmax=(942.08,'K')), NASAPolynomial(coeffs=[17.791,0.0213079,-6.07889e-06,1.03571e-09,-7.56056e-14,22384.3,-66.3645], Tmin=(942.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]CC=CC=C=C(81)',
    structure = SMILES('[CH2]CC=CC=C=C'),
    E0 = (381.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.786135,'amu*angstrom^2'), symmetry=1, barrier=(18.0748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786168,'amu*angstrom^2'), symmetry=1, barrier=(18.0755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786704,'amu*angstrom^2'), symmetry=1, barrier=(18.0879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.734258,0.061292,-2.89755e-05,-6.54144e-09,6.98294e-12,45959.5,27.082], Tmin=(100,'K'), Tmax=(1020.39,'K')), NASAPolynomial(coeffs=[14.9233,0.0271344,-1.03164e-05,1.88299e-09,-1.31933e-13,41946.4,-47.1337], Tmin=(1020.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CC=C[CH2](82)',
    structure = SMILES('[CH]=CC=C[CH2]'),
    E0 = (423.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(1.5274,'amu*angstrom^2'), symmetry=1, barrier=(35.118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53304,'amu*angstrom^2'), symmetry=1, barrier=(35.2477,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22721,0.027298,2.28219e-05,-5.20822e-08,2.29935e-11,50955.4,18.1253], Tmin=(100,'K'), Tmax=(937.459,'K')), NASAPolynomial(coeffs=[12.149,0.0145309,-4.06047e-06,6.79586e-10,-4.92016e-14,47795.9,-36.0307], Tmin=(937.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C2H4(33)',
    structure = SMILES('C=C'),
    E0 = (42.0619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9592,-0.00757051,5.7099e-05,-6.91588e-08,2.69884e-11,5089.78,4.0973], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.99183,0.0104834,-3.71721e-06,5.94628e-10,-3.5363e-14,4268.66,-0.269082], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(42.0619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: FFCM1(-)"""),
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
    label = '[CH2]C=CC=[C]CC(85)',
    structure = SMILES('[CH2]C=CC=[C]CC'),
    E0 = (355.122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.024078,'amu*angstrom^2'), symmetry=1, barrier=(15.5198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.673536,'amu*angstrom^2'), symmetry=1, barrier=(15.4859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.672433,'amu*angstrom^2'), symmetry=1, barrier=(15.4606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401921,'amu*angstrom^2'), symmetry=1, barrier=(83.4152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.909166,0.0565736,-1.09444e-05,-2.38821e-08,1.27322e-11,42832.6,27.6449], Tmin=(100,'K'), Tmax=(998.375,'K')), NASAPolynomial(coeffs=[13.4172,0.0318059,-1.18132e-05,2.12666e-09,-1.48021e-13,39071.9,-39.0052], Tmin=(998.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(Cds_S) + radical(C=CC=CCJ)"""),
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
    label = '[CH2]C=C[C]=CCC(87)',
    structure = SMILES('[CH2]C=C[C]=CCC'),
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
    label = 'C2H4(T)(92)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1436.16,1437.36,2688.4,2689.63],'cm^-1')),
        HinderedRotor(inertia=(0.0257474,'amu*angstrom^2'), symmetry=1, barrier=(17.2422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40737,0.0100311,6.40983e-06,-1.41299e-08,5.92706e-12,38288.2,6.11699], Tmin=(100,'K'), Tmax=(954.25,'K')), NASAPolynomial(coeffs=[5.52245,0.0085618,-2.90747e-06,5.02363e-10,-3.4458e-14,37547.8,-5.75253], Tmin=(954.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=CC[CH2](94)',
    structure = SMILES('[CH]=CC[CH2]'),
    E0 = (435.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.207738,'amu*angstrom^2'), symmetry=1, barrier=(4.77631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204754,'amu*angstrom^2'), symmetry=1, barrier=(4.7077,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52057,0.0279427,-9.18483e-06,-4.19516e-09,2.67714e-12,52484.1,16.9976], Tmin=(100,'K'), Tmax=(1112.3,'K')), NASAPolynomial(coeffs=[7.71987,0.0177227,-6.83493e-06,1.24847e-09,-8.64345e-14,50803,-10.9968], Tmin=(1112.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_P) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C[CH2](93)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,229.656,230.235,230.761],'cm^-1')),
        HinderedRotor(inertia=(1.3329,'amu*angstrom^2'), symmetry=1, barrier=(50.5171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817957,3.34737e-05,-4.36194e-08,1.58214e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.27621e-14,44095,-3.44604], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C=CC1[CH]CC1(96)',
    structure = SMILES('[CH2]C=CC1[CH]CC1'),
    E0 = (338.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77923,0.0293334,6.72036e-05,-1.01189e-07,3.86111e-11,40787.3,24.325], Tmin=(100,'K'), Tmax=(988.526,'K')), NASAPolynomial(coeffs=[12.8166,0.0333057,-1.26223e-05,2.41579e-09,-1.77687e-13,36228.9,-40.8158], Tmin=(988.526,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(cyclobutane)"""),
)

species(
    label = '[CH2]CC1[CH]C=CC1(97)',
    structure = SMILES('[CH2]CC1[CH]C=CC1'),
    E0 = (279.478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98122,0.0218598,9.08508e-05,-1.27193e-07,4.81716e-11,33706.3,21.4378], Tmin=(100,'K'), Tmax=(977.251,'K')), NASAPolynomial(coeffs=[13.4494,0.0317689,-1.16184e-05,2.23757e-09,-1.67385e-13,28750.2,-47.5148], Tmin=(977.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(cyclopentene-allyl)"""),
)

species(
    label = '[CH2]C1[CH]C=CCC1(98)',
    structure = SMILES('[CH2]C1[CH]C=CCC1'),
    E0 = (278.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98161,0.0222824,9.00515e-05,-1.27929e-07,4.92587e-11,33618,20.4276], Tmin=(100,'K'), Tmax=(962.172,'K')), NASAPolynomial(coeffs=[13.3364,0.030883,-1.03559e-05,1.92083e-09,-1.4243e-13,28849.8,-47.3386], Tmin=(962.172,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclohexene-allyl) + radical(Isobutyl)"""),
)

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
    label = 'CH2(T)(22)',
    structure = SMILES('[CH2]'),
    E0 = (381.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([971.045,2816.03,3444.23],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71758,0.00127391,2.17347e-06,-3.48858e-09,1.65209e-12,45872.4,1.75298], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.14632,0.00303671,-9.96474e-07,1.50484e-10,-8.57336e-15,46041.3,4.72342], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(381.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]C=CC=C[CH2](101)',
    structure = SMILES('[CH2]C=CC=C[CH2]'),
    E0 = (257.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,467.807],'cm^-1')),
        HinderedRotor(inertia=(0.359587,'amu*angstrom^2'), symmetry=1, barrier=(55.8572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359679,'amu*angstrom^2'), symmetry=1, barrier=(55.857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359772,'amu*angstrom^2'), symmetry=1, barrier=(55.8558,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91519,0.0306767,4.13763e-05,-7.59787e-08,3.17652e-11,31116.9,21.5128], Tmin=(100,'K'), Tmax=(942.962,'K')), NASAPolynomial(coeffs=[12.7395,0.0229465,-7.07045e-06,1.2179e-09,-8.69869e-14,27377.9,-39.0743], Tmin=(942.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC=CC[CH2](102)',
    structure = SMILES('[CH]=CC=CC[CH2]'),
    E0 = (487.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,264.94],'cm^-1')),
        HinderedRotor(inertia=(0.305291,'amu*angstrom^2'), symmetry=1, barrier=(15.284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.306246,'amu*angstrom^2'), symmetry=1, barrier=(15.2799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307463,'amu*angstrom^2'), symmetry=1, barrier=(15.2769,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30257,0.0492014,-1.47407e-05,-1.75085e-08,1.06844e-11,58750,24.5848], Tmin=(100,'K'), Tmax=(985.505,'K')), NASAPolynomial(coeffs=[13.7033,0.0217744,-7.85887e-06,1.42054e-09,-1.00261e-13,55193.4,-40.704], Tmin=(985.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC=C=CC=C(103)',
    structure = SMILES('[CH2]CC=C=CC=C'),
    E0 = (381.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.786135,'amu*angstrom^2'), symmetry=1, barrier=(18.0748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786168,'amu*angstrom^2'), symmetry=1, barrier=(18.0755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786704,'amu*angstrom^2'), symmetry=1, barrier=(18.0879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.734258,0.061292,-2.89755e-05,-6.54144e-09,6.98294e-12,45959.5,27.082], Tmin=(100,'K'), Tmax=(1020.39,'K')), NASAPolynomial(coeffs=[14.9233,0.0271344,-1.03164e-05,1.88299e-09,-1.31933e-13,41946.4,-47.1337], Tmin=(1020.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC=[C]CC=C(104)',
    structure = SMILES('[CH2]CC=[C]CC=C'),
    E0 = (473.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,467.87,470.971,471.049],'cm^-1')),
        HinderedRotor(inertia=(0.116217,'amu*angstrom^2'), symmetry=1, barrier=(2.67206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116873,'amu*angstrom^2'), symmetry=1, barrier=(2.68714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493568,'amu*angstrom^2'), symmetry=1, barrier=(11.3481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716102,'amu*angstrom^2'), symmetry=1, barrier=(11.3425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719682,0.0631224,-4.05446e-05,1.30994e-08,-1.71645e-12,57034.4,31.6017], Tmin=(100,'K'), Tmax=(1754.47,'K')), NASAPolynomial(coeffs=[15.7769,0.0287936,-1.11949e-05,1.94705e-09,-1.27306e-13,51750.9,-49.5051], Tmin=(1754.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC=CC[C]=C(105)',
    structure = SMILES('[CH2]CC=CC[C]=C'),
    E0 = (473.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,467.87,470.971,471.049],'cm^-1')),
        HinderedRotor(inertia=(0.116217,'amu*angstrom^2'), symmetry=1, barrier=(2.67206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116873,'amu*angstrom^2'), symmetry=1, barrier=(2.68714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493568,'amu*angstrom^2'), symmetry=1, barrier=(11.3481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716102,'amu*angstrom^2'), symmetry=1, barrier=(11.3425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719682,0.0631224,-4.05446e-05,1.30994e-08,-1.71645e-12,57034.4,31.6017], Tmin=(100,'K'), Tmax=(1754.47,'K')), NASAPolynomial(coeffs=[15.7769,0.0287936,-1.11949e-05,1.94705e-09,-1.27306e-13,51750.9,-49.5051], Tmin=(1754.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C[C]=CCC=C(106)',
    structure = SMILES('[CH2]C[C]=CCC=C'),
    E0 = (473.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,467.87,470.971,471.049],'cm^-1')),
        HinderedRotor(inertia=(0.116217,'amu*angstrom^2'), symmetry=1, barrier=(2.67206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116873,'amu*angstrom^2'), symmetry=1, barrier=(2.68714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493568,'amu*angstrom^2'), symmetry=1, barrier=(11.3481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716102,'amu*angstrom^2'), symmetry=1, barrier=(11.3425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719682,0.0631224,-4.05446e-05,1.30994e-08,-1.71645e-12,57034.4,31.6017], Tmin=(100,'K'), Tmax=(1754.47,'K')), NASAPolynomial(coeffs=[15.7769,0.0287936,-1.11949e-05,1.94705e-09,-1.27306e-13,51750.9,-49.5051], Tmin=(1754.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CCC=CC[CH2](107)',
    structure = SMILES('[CH]=CCC=CC[CH2]'),
    E0 = (482.421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,264.108,265.525],'cm^-1')),
        HinderedRotor(inertia=(0.00242939,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00240217,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.278996,'amu*angstrom^2'), symmetry=1, barrier=(13.5941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270166,'amu*angstrom^2'), symmetry=1, barrier=(13.5989,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.790519,0.0610148,-3.00894e-05,-9.34702e-11,3.2155e-12,58145.4,31.2573], Tmin=(100,'K'), Tmax=(1119.14,'K')), NASAPolynomial(coeffs=[13.2192,0.0324369,-1.30228e-05,2.39083e-09,-1.65477e-13,54371.3,-34.5362], Tmin=(1119.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
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
    label = '[CH2]CC1[CH]C1C=C(110)',
    structure = SMILES('[CH2]CC1[CH]C1C=C'),
    E0 = (497.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51951,0.0414836,2.29724e-05,-4.99596e-08,1.93367e-11,59910.2,27.2282], Tmin=(100,'K'), Tmax=(1040.11,'K')), NASAPolynomial(coeffs=[10.7384,0.0368769,-1.48696e-05,2.80862e-09,-1.99997e-13,56323.9,-25.6312], Tmin=(1040.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane) + radical(RCCJ)"""),
)

species(
    label = '[CH]1C=C[CH]CCC1(51)',
    structure = SMILES('[CH]1C=C[CH]CCC1'),
    E0 = (233.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33755,0.00875945,0.000132706,-1.72324e-07,6.44554e-11,28131.9,15.3206], Tmin=(100,'K'), Tmax=(967.659,'K')), NASAPolynomial(coeffs=[14.2911,0.0302458,-1.05031e-05,2.05604e-09,-1.58625e-13,22499.1,-59.1067], Tmin=(967.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Allyl_S)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.928268,0.0526924,5.52854e-06,-4.52831e-08,2.12066e-11,17024.7,26.4844], Tmin=(100,'K'), Tmax=(980.344,'K')), NASAPolynomial(coeffs=[15.6612,0.0277121,-1.00058e-05,1.83677e-09,-1.3177e-13,12447.8,-52.9117], Tmin=(980.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC=C=CCC(112)',
    structure = SMILES('C=CC=C=CCC'),
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
    label = '[CH2]CC1C=CC1[CH2](74)',
    structure = SMILES('[CH2]CC1C=CC1[CH2]'),
    E0 = (468.451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1082.06,1082.06,1082.06,1082.06,1082.06,1082.06,1082.06,1082.06,1082.06,1082.06,1082.06,1082.06,1082.06,2231.56],'cm^-1')),
        HinderedRotor(inertia=(0.0951722,'amu*angstrom^2'), symmetry=1, barrier=(2.1882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0951722,'amu*angstrom^2'), symmetry=1, barrier=(2.1882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0951722,'amu*angstrom^2'), symmetry=1, barrier=(2.1882,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44233,0.0408211,3.32396e-05,-6.77716e-08,2.77725e-11,56447.5,26.1784], Tmin=(100,'K'), Tmax=(977.732,'K')), NASAPolynomial(coeffs=[12.6765,0.0324502,-1.15859e-05,2.11338e-09,-1.50777e-13,52454,-36.9548], Tmin=(977.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC[C]=CC=C(115)',
    structure = SMILES('[CH2]CC[C]=CC=C'),
    E0 = (454.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,319.261,319.308,319.337],'cm^-1')),
        HinderedRotor(inertia=(0.00165081,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190649,'amu*angstrom^2'), symmetry=1, barrier=(13.7859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190691,'amu*angstrom^2'), symmetry=1, barrier=(13.786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190609,'amu*angstrom^2'), symmetry=1, barrier=(13.7856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.538822,0.0668724,-4.27759e-05,9.85535e-09,5.45876e-13,54803.1,29.6451], Tmin=(100,'K'), Tmax=(1112.51,'K')), NASAPolynomial(coeffs=[14.3354,0.0306229,-1.19084e-05,2.14929e-09,-1.47426e-13,50906.8,-42.1011], Tmin=(1112.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CCC=[C]C=C(116)',
    structure = SMILES('[CH2]CCC=[C]C=C'),
    E0 = (415.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,494.828,494.943],'cm^-1')),
        HinderedRotor(inertia=(0.10356,'amu*angstrom^2'), symmetry=1, barrier=(12.956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.5635,'amu*angstrom^2'), symmetry=1, barrier=(12.956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.563498,'amu*angstrom^2'), symmetry=1, barrier=(12.9559,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.599958,'amu*angstrom^2'), symmetry=1, barrier=(104.277,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.482421,0.068984,-4.78332e-05,1.46888e-08,-8.94701e-13,50132.3,28.9684], Tmin=(100,'K'), Tmax=(1096.54,'K')), NASAPolynomial(coeffs=[13.846,0.0314272,-1.1767e-05,2.06913e-09,-1.39605e-13,46528.7,-39.8026], Tmin=(1096.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCC=C[C]=C(117)',
    structure = SMILES('[CH2]CCC=C[C]=C'),
    E0 = (415.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,494.818,495.269],'cm^-1')),
        HinderedRotor(inertia=(0.103688,'amu*angstrom^2'), symmetry=1, barrier=(12.956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.5635,'amu*angstrom^2'), symmetry=1, barrier=(12.956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0745633,'amu*angstrom^2'), symmetry=1, barrier=(12.9561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.53538,'amu*angstrom^2'), symmetry=1, barrier=(104.277,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.482406,0.0689841,-4.78338e-05,1.46895e-08,-8.94981e-13,50132.3,28.9684], Tmin=(100,'K'), Tmax=(1096.55,'K')), NASAPolynomial(coeffs=[13.8461,0.0314271,-1.17669e-05,2.06911e-09,-1.39604e-13,46528.7,-39.8031], Tmin=(1096.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CC=CCC[CH2](55)',
    structure = SMILES('[CH]=CC=CCC[CH2]'),
    E0 = (463.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,381.276,381.293],'cm^-1')),
        HinderedRotor(inertia=(0.129699,'amu*angstrom^2'), symmetry=1, barrier=(13.3844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12972,'amu*angstrom^2'), symmetry=1, barrier=(13.3847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129719,'amu*angstrom^2'), symmetry=1, barrier=(13.3855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129729,'amu*angstrom^2'), symmetry=1, barrier=(13.3829,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.549812,0.0652362,-3.29546e-05,-3.79979e-09,6.14836e-12,55917.2,29.532], Tmin=(100,'K'), Tmax=(1026.67,'K')), NASAPolynomial(coeffs=[15.217,0.0291835,-1.10964e-05,2.01691e-09,-1.40639e-13,51793.9,-47.0284], Tmin=(1026.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(RCCJ)"""),
)

species(
    label = '[CH]1CC=CC=CC1(42)',
    structure = SMILES('[CH]1CC=CC=CC1'),
    E0 = (262.567,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25297,0.021798,6.88074e-05,-9.48413e-08,3.49359e-11,31657.1,20.1976], Tmin=(100,'K'), Tmax=(994.531,'K')), NASAPolynomial(coeffs=[9.55767,0.0343082,-1.32411e-05,2.51016e-09,-1.81994e-13,28132.5,-25.4184], Tmin=(994.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC)"""),
)

species(
    label = '[CH]1C=CC=CCC1(43)',
    structure = SMILES('[CH]1C=CCCC=C1'),
    E0 = (167.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2863,0.00959575,0.000128536,-1.73375e-07,6.6947e-11,20259.5,17.7361], Tmin=(100,'K'), Tmax=(947.927,'K')), NASAPolynomial(coeffs=[15.9797,0.0233649,-6.47559e-06,1.20576e-09,-9.641e-14,14448.7,-64.5505], Tmin=(947.927,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CCJC=C)"""),
)

species(
    label = '[C]1=CC=CCCC1(44)',
    structure = SMILES('[C]1=CC=CCCC1'),
    E0 = (305.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9496,0.0261433,6.70442e-05,-1.00258e-07,3.85129e-11,36888,19.8882], Tmin=(100,'K'), Tmax=(980.343,'K')), NASAPolynomial(coeffs=[12.6587,0.0298041,-1.10154e-05,2.09907e-09,-1.55011e-13,32512.7,-43.1705], Tmin=(980.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S)"""),
)

species(
    label = '[C]1=CCCCC=C1(45)',
    structure = SMILES('[C]1=CCCCC=C1'),
    E0 = (267.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89666,0.0281884,6.2347e-05,-9.60905e-08,3.74486e-11,32217.1,19.2007], Tmin=(100,'K'), Tmax=(970.591,'K')), NASAPolynomial(coeffs=[12.303,0.0303972,-1.07588e-05,1.99276e-09,-1.4509e-13,28073,-41.6352], Tmin=(970.591,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CJC=C)"""),
)

species(
    label = '[CH]1C=CCC[CH]C1(46)',
    structure = SMILES('[CH]1C=CCC[CH]C1'),
    E0 = (286.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.40795,0.0142582,0.000100665,-1.29074e-07,4.69074e-11,34538.7,19.4329], Tmin=(100,'K'), Tmax=(987.855,'K')), NASAPolynomial(coeffs=[9.96095,0.0367303,-1.40199e-05,2.69041e-09,-1.9795e-13,30457.7,-30.0165], Tmin=(987.855,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC) + radical(Allyl_S)"""),
)

species(
    label = '[C]1=CCCC[CH]C1(47)',
    structure = SMILES('[C]1=CCCC[CH]C1'),
    E0 = (383.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16764,0.0241747,6.66832e-05,-9.11431e-08,3.2962e-11,46176.9,21.8767], Tmin=(100,'K'), Tmax=(1010.06,'K')), NASAPolynomial(coeffs=[8.84398,0.0385262,-1.52071e-05,2.88957e-09,-2.08318e-13,42747.4,-20.6999], Tmin=(1010.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(RCCJCC)"""),
)

species(
    label = '[C]1=CC[CH]CCC1(48)',
    structure = SMILES('[C]1=CC[CH]CCC1'),
    E0 = (383.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16764,0.0241747,6.66832e-05,-9.11431e-08,3.2962e-11,46176.9,21.8767], Tmin=(100,'K'), Tmax=(1010.06,'K')), NASAPolynomial(coeffs=[8.84398,0.0385262,-1.52071e-05,2.88957e-09,-2.08318e-13,42747.4,-20.6999], Tmin=(1010.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(RCCJCC)"""),
)

species(
    label = '[CH]1C[CH]CC=CC1(49)',
    structure = SMILES('[CH]1C[CH]CC=CC1'),
    E0 = (339.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,301.656,802.082,802.082,802.082,802.082,802.082,802.082,802.082,802.082,802.082,802.082,802.082,802.082,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47139,0.0198093,6.85982e-05,-8.60552e-08,2.95893e-11,40945.9,22.1856], Tmin=(100,'K'), Tmax=(1032.62,'K')), NASAPolynomial(coeffs=[5.8299,0.0428903,-1.73551e-05,3.28287e-09,-2.33859e-13,38328.1,-3.44209], Tmin=(1032.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC) + radical(RCCJCC)"""),
)

species(
    label = '[C]1=C[CH]CCCC1(50)',
    structure = SMILES('[C]1=C[CH]CCCC1'),
    E0 = (329.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10295,0.0186212,9.88479e-05,-1.34437e-07,5.04689e-11,39769.8,18.4363], Tmin=(100,'K'), Tmax=(977.3,'K')), NASAPolynomial(coeffs=[13.0782,0.0321993,-1.17789e-05,2.27576e-09,-1.70675e-13,34830.9,-48.5538], Tmin=(977.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = '[C]1[CH]CCCCC=1(52)',
    structure = SMILES('[C]1[CH]CCCCC=1'),
    E0 = (329.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10295,0.0186212,9.88479e-05,-1.34437e-07,5.04689e-11,39769.8,17.7432], Tmin=(100,'K'), Tmax=(977.3,'K')), NASAPolynomial(coeffs=[13.0782,0.0321993,-1.17789e-05,2.27576e-09,-1.70675e-13,34830.9,-49.2469], Tmin=(977.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = '[CH]1[CH]CCC=CC1(53)',
    structure = SMILES('[CH]1[CH]CCC=CC1'),
    E0 = (339.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,301.656,802.082,802.082,802.082,802.082,802.082,802.082,802.082,802.082,802.082,802.082,802.082,802.082,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35,1607.35],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47139,0.0198093,6.85982e-05,-8.60552e-08,2.95893e-11,40945.9,22.8788], Tmin=(100,'K'), Tmax=(1032.62,'K')), NASAPolynomial(coeffs=[5.8299,0.0428903,-1.73551e-05,3.28287e-09,-2.33859e-13,38328.1,-2.74894], Tmin=(1032.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC) + radical(RCCJCC)"""),
)

species(
    label = '[CH]=CCCCC=[CH](54)',
    structure = SMILES('[CH]=CCCCC=[CH]'),
    E0 = (534.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,3115,3125,620,680,785,800,1600,1700,180,952.386],'cm^-1')),
        HinderedRotor(inertia=(0.0485358,'amu*angstrom^2'), symmetry=1, barrier=(12.1993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.536147,'amu*angstrom^2'), symmetry=1, barrier=(12.3271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.534169,'amu*angstrom^2'), symmetry=1, barrier=(12.2816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.533585,'amu*angstrom^2'), symmetry=1, barrier=(12.2682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.475672,0.0683398,-4.89883e-05,1.80153e-08,-2.68644e-12,64373.3,30.0646], Tmin=(100,'K'), Tmax=(1568.63,'K')), NASAPolynomial(coeffs=[15.9415,0.0289026,-1.1277e-05,1.98826e-09,-1.32168e-13,59521.2,-51.512], Tmin=(1568.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[C]1CC=CCCC1(57)',
    structure = SMILES('[C]1CC=CCCC1'),
    E0 = (347.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84216,0.0256082,8.14367e-05,-1.19204e-07,4.58793e-11,41870.3,19.8292], Tmin=(100,'K'), Tmax=(973.145,'K')), NASAPolynomial(coeffs=[13.885,0.0309105,-1.12088e-05,2.13345e-09,-1.58611e-13,36931.4,-51.2753], Tmin=(973.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[C]1C=CCCCC1(58)',
    structure = SMILES('[C]1C=CCCCC1'),
    E0 = (332.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8346,0.0259683,8.01649e-05,-1.17697e-07,4.52973e-11,40114.4,19.9179], Tmin=(100,'K'), Tmax=(973.898,'K')), NASAPolynomial(coeffs=[13.8082,0.0310745,-1.13081e-05,2.15197e-09,-1.59783e-13,35207.8,-50.7478], Tmin=(973.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C1=CC2CCC[C]12(64)',
    structure = SMILES('C1=CC2CCC[C]12'),
    E0 = (240.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.57993,0.00371806,0.000137212,-1.76672e-07,6.65775e-11,29002.2,16.7071], Tmin=(100,'K'), Tmax=(956.045,'K')), NASAPolynomial(coeffs=[13.9896,0.0262017,-8.23803e-06,1.57825e-09,-1.23493e-13,23611.4,-54.6086], Tmin=(956.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_T)"""),
)

species(
    label = '[CH]1CCC2C=CC12(65)',
    structure = SMILES('[CH]1CCC2C=CC12'),
    E0 = (303.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49742,0.00783907,0.000120566,-1.55982e-07,5.81176e-11,36527.3,20.6135], Tmin=(100,'K'), Tmax=(970.552,'K')), NASAPolynomial(coeffs=[13.2017,0.0280233,-1.00061e-05,1.96949e-09,-1.51565e-13,31421,-46.3104], Tmin=(970.552,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S)"""),
)

species(
    label = '[CH]1CC2C=CC2C1(66)',
    structure = SMILES('[CH]1CC2C=CC2C1'),
    E0 = (294.406,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.57733,0.0093218,0.000108831,-1.39241e-07,5.13193e-11,35479.7,19.6416], Tmin=(100,'K'), Tmax=(973.876,'K')), NASAPolynomial(coeffs=[10.8737,0.0313177,-1.14117e-05,2.19037e-09,-1.63639e-13,31204.8,-33.8156], Tmin=(973.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane)"""),
)

species(
    label = '[C]1=CC2CCCC12(67)',
    structure = SMILES('[C]1=CC2CCCC12'),
    E0 = (361.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27168,0.013688,0.000107027,-1.44663e-07,5.49274e-11,43503.7,19.3407], Tmin=(100,'K'), Tmax=(965.462,'K')), NASAPolynomial(coeffs=[14.0173,0.0267435,-9.1465e-06,1.77011e-09,-1.35905e-13,38359.2,-51.809], Tmin=(965.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH]1C[C]2CCCC12(68)',
    structure = SMILES('[CH]1C[C]2CCCC12'),
    E0 = (338.717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,303.443,805.836,805.836,805.836,805.836,805.836,805.836,805.836,805.836,805.836,805.836,805.836,805.836,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58897,0.0127337,9.48263e-05,-1.16717e-07,4.09598e-11,40805,22.7042], Tmin=(100,'K'), Tmax=(1010.55,'K')), NASAPolynomial(coeffs=[7.78084,0.0398249,-1.61033e-05,3.11639e-09,-2.27366e-13,37323.1,-14.4341], Tmin=(1010.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-6) + radical(bicyclo[3.2.0]heptane-tertiary)"""),
)

species(
    label = '[CH]1CC2CCC[C]12(69)',
    structure = SMILES('[CH]1CC2CCC[C]12'),
    E0 = (338.717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,303.443,805.836,805.836,805.836,805.836,805.836,805.836,805.836,805.836,805.836,805.836,805.836,805.836,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58897,0.0127337,9.48263e-05,-1.16717e-07,4.09598e-11,40805,22.7042], Tmin=(100,'K'), Tmax=(1010.55,'K')), NASAPolynomial(coeffs=[7.78084,0.0398249,-1.61033e-05,3.11639e-09,-2.27366e-13,37323.1,-14.4341], Tmin=(1010.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-6) + radical(bicyclo[3.2.0]heptane-tertiary)"""),
)

species(
    label = '[CH]1CCC2[CH]CC12(70)',
    structure = SMILES('[CH]1CCC2[CH]CC12'),
    E0 = (319.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66787,0.00592453,0.000123207,-1.51788e-07,5.46322e-11,38531.1,22.9388], Tmin=(100,'K'), Tmax=(986.711,'K')), NASAPolynomial(coeffs=[10.1908,0.0357286,-1.37713e-05,2.69645e-09,-2.01726e-13,34111.1,-28.1294], Tmin=(986.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-2) + radical(bicyclo[3.2.0]heptane-C5-6)"""),
)

species(
    label = '[CH]1CCC2C[CH]C12(71)',
    structure = SMILES('[CH]1CCC2C[CH]C12'),
    E0 = (319.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66787,0.00592453,0.000123207,-1.51788e-07,5.46322e-11,38531.1,22.9388], Tmin=(100,'K'), Tmax=(986.711,'K')), NASAPolynomial(coeffs=[10.1908,0.0357286,-1.37713e-05,2.69645e-09,-2.01726e-13,34111.1,-28.1294], Tmin=(986.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-2) + radical(bicyclo[3.2.0]heptane-C5-6)"""),
)

species(
    label = '[CH]1CC2[CH]CC2C1(72)',
    structure = SMILES('[CH]1CC2[CH]CC2C1'),
    E0 = (326.489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66787,0.00592453,0.000123207,-1.51788e-07,5.46322e-11,39336.3,22.9388], Tmin=(100,'K'), Tmax=(986.711,'K')), NASAPolynomial(coeffs=[10.1908,0.0357286,-1.37713e-05,2.69645e-09,-2.01726e-13,34916.2,-28.1294], Tmin=(986.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-3) + radical(bicyclo[3.2.0]heptane-C5-6)"""),
)

species(
    label = '[CH]=CC1[CH]CCC1(73)',
    structure = SMILES('[CH]=CC1[CH]CCC1'),
    E0 = (403.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88815,0.0268695,7.1117e-05,-1.03571e-07,3.90046e-11,48579.9,25.251], Tmin=(100,'K'), Tmax=(993.5,'K')), NASAPolynomial(coeffs=[12.4137,0.0335161,-1.29355e-05,2.49852e-09,-1.84363e-13,44069.1,-37.6363], Tmin=(993.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_P) + radical(Cs_S)"""),
)

species(
    label = '[CH2]CCC1[CH]C=C1(75)',
    structure = SMILES('[CH2]CCC1[CH]C=C1'),
    E0 = (434.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,358.221,1019.64,1019.64,1019.64,1019.64,1019.64,1019.64,1019.64,1019.64,1019.64,1019.64,1019.64,1019.64,2255.94],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78978,0.0277582,7.43533e-05,-1.10502e-07,4.23811e-11,52310.8,22.7838], Tmin=(100,'K'), Tmax=(981.019,'K')), NASAPolynomial(coeffs=[13.551,0.0319636,-1.18317e-05,2.26537e-09,-1.6792e-13,47493.2,-46.5247], Tmin=(981.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(434.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(cyclobutene-allyl)"""),
)

species(
    label = '[C]1CC2CCCC12(76)',
    structure = SMILES('[C]1CC2CCCC12'),
    E0 = (337.888,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22899,0.0110646,0.00012739,-1.69792e-07,6.44941e-11,40727.5,19.2813], Tmin=(100,'K'), Tmax=(960.675,'K')), NASAPolynomial(coeffs=[15.2641,0.0276615,-9.18408e-06,1.77644e-09,-1.37922e-13,34952.6,-60.104], Tmin=(960.675,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH2]C=[C]CC1CC1(147)',
    structure = SMILES('[CH2]C=[C]CC1CC1'),
    E0 = (444.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3000,3100,440,815,1455,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35405,0.0409231,3.81221e-05,-7.43312e-08,2.99995e-11,53578.4,24.9966], Tmin=(100,'K'), Tmax=(989.86,'K')), NASAPolynomial(coeffs=[13.9458,0.0320164,-1.19903e-05,2.25997e-09,-1.64285e-13,49029.1,-46.0106], Tmin=(989.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C[C]=C[CH]C1CC1(148)',
    structure = SMILES('C[C]=C[CH]C1CC1'),
    E0 = (434.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180,1076,1076,1076,1076,1076,1076,1076,1076,1076,1076,2261.01],'cm^-1')),
        HinderedRotor(inertia=(0.0559204,'amu*angstrom^2'), symmetry=1, barrier=(1.28572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0559204,'amu*angstrom^2'), symmetry=1, barrier=(1.28572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0559204,'amu*angstrom^2'), symmetry=1, barrier=(1.28572,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54617,0.0382553,3.8353e-05,-6.92969e-08,2.69165e-11,52320.7,22.9209], Tmin=(100,'K'), Tmax=(1008.64,'K')), NASAPolynomial(coeffs=[11.922,0.0352081,-1.37762e-05,2.60807e-09,-1.87922e-13,48289.5,-36.8331], Tmin=(1008.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(434.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC1CC1(149)',
    structure = SMILES('[CH2][C]=CCC1CC1'),
    E0 = (444.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3000,3100,440,815,1455,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35405,0.0409231,3.81221e-05,-7.43312e-08,2.99995e-11,53578.4,24.9966], Tmin=(100,'K'), Tmax=(989.86,'K')), NASAPolynomial(coeffs=[13.9458,0.0320164,-1.19903e-05,2.25997e-09,-1.64285e-13,49029.1,-46.0106], Tmin=(989.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=[C][CH]C1CC1(150)',
    structure = SMILES('CC=[C][CH]C1CC1'),
    E0 = (434.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1685,370,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,180,1076,1076,1076,1076,1076,1076,1076,1076,1076,1076,2261.01],'cm^-1')),
        HinderedRotor(inertia=(0.0559204,'amu*angstrom^2'), symmetry=1, barrier=(1.28572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0559204,'amu*angstrom^2'), symmetry=1, barrier=(1.28572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0559204,'amu*angstrom^2'), symmetry=1, barrier=(1.28572,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54617,0.0382553,3.8353e-05,-6.92969e-08,2.69165e-11,52320.7,22.9209], Tmin=(100,'K'), Tmax=(1008.64,'K')), NASAPolynomial(coeffs=[11.922,0.0352081,-1.37762e-05,2.60807e-09,-1.87922e-13,48289.5,-36.8331], Tmin=(1008.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(434.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1[CH]C1C1CC1(151)',
    structure = SMILES('[CH2]C1[CH]C1C1CC1'),
    E0 = (529.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83245,0.0312618,5.44252e-05,-8.51916e-08,3.2773e-11,63754.6,26.213], Tmin=(100,'K'), Tmax=(983.466,'K')), NASAPolynomial(coeffs=[11.0122,0.0345942,-1.26858e-05,2.34887e-09,-1.68692e-13,59982.3,-27.92], Tmin=(983.466,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane) + radical(Isobutyl)"""),
)

species(
    label = '[CH]1CC1[CH]C1CC1(152)',
    structure = SMILES('[CH]1CC1[CH]C1CC1'),
    E0 = (528.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15476,0.0221735,7.65085e-05,-1.03678e-07,3.77111e-11,63665.3,26.4853], Tmin=(100,'K'), Tmax=(1006.31,'K')), NASAPolynomial(coeffs=[10.2649,0.0366367,-1.46616e-05,2.83814e-09,-2.07777e-13,59668.5,-24.441], Tmin=(1006.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(528.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(cyclopropane)"""),
)

species(
    label = 'CC=C=CC1CC1(153)',
    structure = SMILES('CC=C=CC1CC1'),
    E0 = (214.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43934,0.0406488,3.33296e-05,-6.69733e-08,2.69945e-11,25960.9,23.6476], Tmin=(100,'K'), Tmax=(991.901,'K')), NASAPolynomial(coeffs=[12.6694,0.0331916,-1.26016e-05,2.34769e-09,-1.68332e-13,21872.1,-39.8207], Tmin=(991.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C=C=CCC1CC1(154)',
    structure = SMILES('C=C=CCC1CC1'),
    E0 = (230.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48902,0.0371768,4.74564e-05,-8.45028e-08,3.39029e-11,27836.3,22.9654], Tmin=(100,'K'), Tmax=(977.099,'K')), NASAPolynomial(coeffs=[13.8106,0.0311159,-1.13705e-05,2.11966e-09,-1.53998e-13,23309.9,-47.0346], Tmin=(977.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]C1CC1(155)',
    structure = SMILES('[CH]C1CC1'),
    E0 = (469.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,600,996.597,996.597,996.597,996.597,996.597,996.597,996.597,996.597,996.597,996.597,996.597,2278.31],'cm^-1')),
        HinderedRotor(inertia=(0.00273462,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.04416,0.00729144,6.14883e-05,-8.28423e-08,3.15551e-11,56546.5,13.8165], Tmin=(100,'K'), Tmax=(960.959,'K')), NASAPolynomial(coeffs=[9.63521,0.0147783,-4.70975e-06,9.00129e-10,-6.95554e-14,53667.3,-26.1087], Tmin=(960.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(469.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=C[CH]C1CC1(156)',
    structure = SMILES('[CH]C=CC1CC1'),
    E0 = (452.439,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,251.744,1070.93,1070.93,1070.93,1070.93,1070.93,1070.93,1070.93,1070.93,1070.93,1070.93,1070.93,1070.93,2256.49],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96642,0.0281718,5.56585e-05,-8.61724e-08,3.32159e-11,54503.8,21.1203], Tmin=(100,'K'), Tmax=(981.331,'K')), NASAPolynomial(coeffs=[11.3342,0.031072,-1.15731e-05,2.16341e-09,-1.56616e-13,50687,-33.9765], Tmin=(981.331,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]1[CH]C(C1)C1CC1(157)',
    structure = SMILES('[CH]1[CH]C(C1)C1CC1'),
    E0 = (427.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34594,0.016722,9.11e-05,-1.18079e-07,4.26676e-11,51499.1,24.4074], Tmin=(100,'K'), Tmax=(998.766,'K')), NASAPolynomial(coeffs=[9.89274,0.0367693,-1.45089e-05,2.80976e-09,-2.06585e-13,47484.2,-24.5446], Tmin=(998.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = 'C=CC[CH][C]1CC1(158)',
    structure = SMILES('C=CC[CH][C]1CC1'),
    E0 = (455.094,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180,1049.08,1049.08,1049.08,1049.08,1049.08,1049.08,1049.08,1049.08,1049.08,1049.08,2282.13],'cm^-1')),
        HinderedRotor(inertia=(0.0470076,'amu*angstrom^2'), symmetry=1, barrier=(1.0808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0470076,'amu*angstrom^2'), symmetry=1, barrier=(1.0808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0470076,'amu*angstrom^2'), symmetry=1, barrier=(1.0808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70377,0.0411851,1.18018e-05,-3.06794e-08,1.08647e-11,54825.4,27.0018], Tmin=(100,'K'), Tmax=(1148.49,'K')), NASAPolynomial(coeffs=[8.17099,0.0410978,-1.73883e-05,3.27472e-09,-2.29053e-13,51860.1,-11.5364], Tmin=(1148.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(455.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Tertalkyl) + radical(Cs_S)"""),
)

species(
    label = '[CH]1C=CCC1(159)',
    structure = SMILES('[CH]1C=CCC1'),
    E0 = (125.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22816,-0.0047173,0.000117611,-1.45673e-07,5.42071e-11,15155.9,11.3916], Tmin=(100,'K'), Tmax=(955.988,'K')), NASAPolynomial(coeffs=[11.0452,0.0181893,-5.59294e-06,1.09697e-09,-8.81254e-14,11120,-39.2607], Tmin=(955.988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-allyl)"""),
)

species(
    label = 'C=C[C]1C=CCC1(160)',
    structure = SMILES('[CH2]C=C1C=CCC1'),
    E0 = (157.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03279,0.0218341,8.50977e-05,-1.22192e-07,4.72283e-11,19066.8,21.0907], Tmin=(100,'K'), Tmax=(962.401,'K')), NASAPolynomial(coeffs=[13.3922,0.0286436,-9.71476e-06,1.81175e-09,-1.34828e-13,14378.6,-46.2742], Tmin=(962.401,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC1[CH]CC=C1(161)',
    structure = SMILES('C=CC1[CH]CC=C1'),
    E0 = (263.676,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01239,0.0257899,6.44891e-05,-9.50586e-08,3.60186e-11,31800.3,22.8012], Tmin=(100,'K'), Tmax=(988.864,'K')), NASAPolynomial(coeffs=[11.6822,0.0314452,-1.20005e-05,2.2926e-09,-1.6813e-13,27699,-34.8093], Tmin=(988.864,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-4)"""),
)

species(
    label = 'C=CC1C=C[CH]C1(162)',
    structure = SMILES('C=CC1[CH]C=CC1'),
    E0 = (201.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05627,0.0196116,9.29616e-05,-1.31303e-07,5.04559e-11,24366.8,19.0083], Tmin=(100,'K'), Tmax=(966.634,'K')), NASAPolynomial(coeffs=[14.3469,0.0271627,-9.39639e-06,1.8036e-09,-1.3703e-13,19261.7,-53.9853], Tmin=(966.634,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-allyl)"""),
)

species(
    label = 'C=CC1[C]=CCC1(163)',
    structure = SMILES('C=CC1[C]=CCC1'),
    E0 = (334.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81262,0.0298374,5.77897e-05,-9.18486e-08,3.58601e-11,40287.7,21.644], Tmin=(100,'K'), Tmax=(979.857,'K')), NASAPolynomial(coeffs=[13.0715,0.0292438,-1.07524e-05,2.0379e-09,-1.4996e-13,35903.3,-43.5582], Tmin=(979.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C=CC1C=[C]CC1(164)',
    structure = SMILES('C=CC1C=[C]CC1'),
    E0 = (334.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81262,0.0298374,5.77897e-05,-9.18486e-08,3.58601e-11,40287.7,21.644], Tmin=(100,'K'), Tmax=(979.857,'K')), NASAPolynomial(coeffs=[13.0715,0.0292438,-1.07524e-05,2.0379e-09,-1.4996e-13,35903.3,-43.5582], Tmin=(979.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C=[C]C1C=CCC1(165)',
    structure = SMILES('C=[C]C1C=CCC1'),
    E0 = (314.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,1685,370,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81262,0.0298374,5.77897e-05,-9.18486e-08,3.58601e-11,37922.6,21.644], Tmin=(100,'K'), Tmax=(979.857,'K')), NASAPolynomial(coeffs=[13.0715,0.0292438,-1.07524e-05,2.0379e-09,-1.4996e-13,33538.2,-43.5582], Tmin=(979.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC1C=CCC1(166)',
    structure = SMILES('[CH]=CC1C=CCC1'),
    E0 = (323.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77161,0.0287664,6.58843e-05,-1.0367e-07,4.08746e-11,39039,21.7204], Tmin=(100,'K'), Tmax=(970.273,'K')), NASAPolynomial(coeffs=[14.4706,0.0269615,-9.46934e-06,1.79677e-09,-1.34309e-13,34195.3,-51.4228], Tmin=(970.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_P)"""),
)

species(
    label = 'C=C[C]1C[CH]CC1(167)',
    structure = SMILES('[CH2]C=C1C[CH]CC1'),
    E0 = (260.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3771.23,'J/mol'), sigma=(6.58809,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=589.06 K, Pc=29.93 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92032,0.0283861,6.30056e-05,-9.23827e-08,3.45647e-11,31370.7,25.2891], Tmin=(100,'K'), Tmax=(996.92,'K')), NASAPolynomial(coeffs=[10.7644,0.0362078,-1.39244e-05,2.6374e-09,-1.91075e-13,27455.3,-28.1444], Tmin=(996.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(cyclopentane)"""),
)

species(
    label = '[CH2]C[C]1C=CCC1(168)',
    structure = SMILES('[CH2]CC1=C[CH]CC1'),
    E0 = (263.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92326,0.0247926,8.01745e-05,-1.14601e-07,4.33528e-11,31776,22.2283], Tmin=(100,'K'), Tmax=(983.383,'K')), NASAPolynomial(coeffs=[12.8072,0.0329833,-1.23421e-05,2.36913e-09,-1.7544e-13,27098.7,-42.9955], Tmin=(983.383,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-allyl) + radical(RCCJ)"""),
)

species(
    label = 'C=CC1[CH]C[CH]C1(169)',
    structure = SMILES('C=CC1[CH]C[CH]C1'),
    E0 = (341.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23266,0.0235713,6.49855e-05,-8.67829e-08,3.06994e-11,41201,26.178], Tmin=(100,'K'), Tmax=(1026.88,'K')), NASAPolynomial(coeffs=[8.03445,0.0401071,-1.63356e-05,3.12577e-09,-2.24969e-13,37946,-12.013], Tmin=(1026.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(Cs_S)"""),
)

species(
    label = '[CH2]CC1[CH]CC=C1(170)',
    structure = SMILES('[CH2]CC1[CH]CC=C1'),
    E0 = (340.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93987,0.0278243,6.32132e-05,-9.20454e-08,3.42118e-11,41019.3,25.2801], Tmin=(100,'K'), Tmax=(1003.69,'K')), NASAPolynomial(coeffs=[10.8994,0.0358355,-1.40944e-05,2.69993e-09,-1.96548e-13,37018.7,-28.9474], Tmin=(1003.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(340.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(cyclopentene-4)"""),
)

species(
    label = '[CH2]CC1[C]=CCC1(171)',
    structure = SMILES('[CH2]CC1[C]=CCC1'),
    E0 = (410.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73942,0.0318838,5.64504e-05,-8.87183e-08,3.39863e-11,49506.7,24.125], Tmin=(100,'K'), Tmax=(994.077,'K')), NASAPolynomial(coeffs=[12.2707,0.0336632,-1.28625e-05,2.44893e-09,-1.78677e-13,45231.2,-37.5933], Tmin=(994.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C=[C]C1C[CH]CC1(172)',
    structure = SMILES('C=[C]C1C[CH]CC1'),
    E0 = (385.225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,1685,370,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00247,0.0294825,5.1178e-05,-7.50344e-08,2.72861e-11,46416.3,24.9206], Tmin=(100,'K'), Tmax=(1022.5,'K')), NASAPolynomial(coeffs=[8.81919,0.0388759,-1.55025e-05,2.93237e-09,-2.09789e-13,43137.2,-17.3357], Tmin=(1022.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(cyclopentane)"""),
)

species(
    label = 'C=CC1[CH]CC[CH]1(173)',
    structure = SMILES('C=CC1[CH]CC[CH]1'),
    E0 = (350.586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15767,0.022042,7.68224e-05,-1.03551e-07,3.74584e-11,42248.3,26.4384], Tmin=(100,'K'), Tmax=(1011.49,'K')), NASAPolynomial(coeffs=[10.273,0.0369589,-1.50121e-05,2.92385e-09,-2.14443e-13,38201.7,-24.6937], Tmin=(1011.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]C1[CH]CC=C1(174)',
    structure = SMILES('C[CH]C1[CH]CC=C1'),
    E0 = (329.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09315,0.0258438,6.25329e-05,-8.6701e-08,3.11223e-11,39725,25.5016], Tmin=(100,'K'), Tmax=(1025.35,'K')), NASAPolynomial(coeffs=[9.27055,0.038632,-1.58446e-05,3.05498e-09,-2.21212e-13,36109,-19.7603], Tmin=(1025.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(cyclopentene-4)"""),
)

species(
    label = 'C[CH]C1[C]=CCC1(175)',
    structure = SMILES('C[CH]C1[C]=CCC1'),
    E0 = (400.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89106,0.0299284,5.5653e-05,-8.31754e-08,3.07899e-11,48212.4,24.3521], Tmin=(100,'K'), Tmax=(1014.6,'K')), NASAPolynomial(coeffs=[10.6168,0.0364993,-1.46344e-05,2.80891e-09,-2.03739e-13,44332.9,-28.2638], Tmin=(1014.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C=[C]C1[CH]CCC1(176)',
    structure = SMILES('C=[C]C1[CH]CCC1'),
    E0 = (393.886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,1685,370,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92877,0.0279363,6.30834e-05,-9.19057e-08,3.40956e-11,47463.5,25.1765], Tmin=(100,'K'), Tmax=(1006.81,'K')), NASAPolynomial(coeffs=[11.0626,0.0357204,-1.41751e-05,2.72961e-09,-1.99197e-13,43390.6,-30.044], Tmin=(1006.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(Cs_S)"""),
)

species(
    label = 'C=CC1C[CH][CH]C1(177)',
    structure = SMILES('C=CC1C[CH][CH]C1'),
    E0 = (333.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30518,0.0251211,5.31209e-05,-7.00483e-08,2.39877e-11,40153.8,25.9269], Tmin=(100,'K'), Tmax=(1051.66,'K')), NASAPolynomial(coeffs=[5.86137,0.0431497,-1.76007e-05,3.31428e-09,-2.34408e-13,37660.8,0.29507], Tmin=(1051.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(cyclopentane)"""),
)

species(
    label = 'C=C[C]1[CH]CCC1(178)',
    structure = SMILES('[CH2]C=C1[CH]CCC1'),
    E0 = (215.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85517,0.0228285,9.52346e-05,-1.35843e-07,5.21847e-11,25995.3,21.851], Tmin=(100,'K'), Tmax=(970.064,'K')), NASAPolynomial(coeffs=[15.051,0.0297956,-1.04487e-05,2.01263e-09,-1.5254e-13,20547.1,-56.2957], Tmin=(970.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]CC1C=[C]CC1(179)',
    structure = SMILES('[CH2]CC1C=[C]CC1'),
    E0 = (410.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73942,0.0318838,5.64504e-05,-8.87183e-08,3.39863e-11,49506.7,24.125], Tmin=(100,'K'), Tmax=(994.077,'K')), NASAPolynomial(coeffs=[12.2707,0.0336632,-1.28625e-05,2.44893e-09,-1.78677e-13,45231.2,-37.5933], Tmin=(994.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-vinyl) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CC1C[CH]CC1(180)',
    structure = SMILES('[CH]=CC1C[CH]CC1'),
    E0 = (394.479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96392,0.0283963,5.92531e-05,-8.67094e-08,3.21779e-11,47532.6,24.9873], Tmin=(100,'K'), Tmax=(1003.84,'K')), NASAPolynomial(coeffs=[10.1315,0.0367348,-1.42983e-05,2.70944e-09,-1.95621e-13,43832.9,-24.7077], Tmin=(1003.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(Cds_P)"""),
)

species(
    label = 'C[CH]C1C=C[CH]C1(181)',
    structure = SMILES('C[CH]C1[CH]C=CC1'),
    E0 = (268.774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13488,0.0198907,9.00477e-05,-1.21553e-07,4.48894e-11,32411.9,21.657], Tmin=(100,'K'), Tmax=(990.951,'K')), NASAPolynomial(coeffs=[11.7306,0.0347109,-1.34495e-05,2.61121e-09,-1.93559e-13,27880.7,-37.8168], Tmin=(990.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-allyl) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]C1C=[C]CC1(182)',
    structure = SMILES('C[CH]C1C=[C]CC1'),
    E0 = (400.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89106,0.0299284,5.5653e-05,-8.31754e-08,3.07899e-11,48212.4,24.3521], Tmin=(100,'K'), Tmax=(1014.6,'K')), NASAPolynomial(coeffs=[10.6168,0.0364993,-1.46344e-05,2.80891e-09,-2.03739e-13,44332.9,-28.2638], Tmin=(1014.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-vinyl) + radical(Cs_S)"""),
)

species(
    label = '[CH]=CC(C=C)C[CH2](183)',
    structure = SMILES('[CH]=CC(C=C)C[CH2]'),
    E0 = (483.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,1188.92],'cm^-1')),
        HinderedRotor(inertia=(0.588323,'amu*angstrom^2'), symmetry=1, barrier=(13.5267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127996,'amu*angstrom^2'), symmetry=1, barrier=(2.94287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587045,'amu*angstrom^2'), symmetry=1, barrier=(13.4973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587685,'amu*angstrom^2'), symmetry=1, barrier=(13.512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.887814,0.0683421,-5.20658e-05,2.16265e-08,-3.81158e-12,58292.7,28.0129], Tmin=(100,'K'), Tmax=(1295.93,'K')), NASAPolynomial(coeffs=[11.0732,0.0369036,-1.56761e-05,2.90627e-09,-2.0018e-13,55652.8,-23.7659], Tmin=(1295.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CCC[CH]C=C(184)',
    structure = SMILES('[CH]=CCCC=C[CH2]'),
    E0 = (426.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.00315625,'amu*angstrom^2'), symmetry=1, barrier=(2.24479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805887,'amu*angstrom^2'), symmetry=1, barrier=(18.5289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.8062,'amu*angstrom^2'), symmetry=1, barrier=(18.5361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805567,'amu*angstrom^2'), symmetry=1, barrier=(18.5216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.726582,0.0612088,-2.47002e-05,-8.65034e-09,6.80657e-12,51395.5,29.2831], Tmin=(100,'K'), Tmax=(1059.96,'K')), NASAPolynomial(coeffs=[14.0379,0.0318345,-1.26495e-05,2.33585e-09,-1.63408e-13,47401.8,-41.2388], Tmin=(1059.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=CC([CH2])C=C(185)',
    structure = SMILES('[CH2]C=CC([CH2])C=C'),
    E0 = (375.736,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,3851.24],'cm^-1')),
        HinderedRotor(inertia=(1.07702,'amu*angstrom^2'), symmetry=1, barrier=(24.7628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.89484,'amu*angstrom^2'), symmetry=1, barrier=(89.55,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.530729,'amu*angstrom^2'), symmetry=1, barrier=(12.2025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743502,'amu*angstrom^2'), symmetry=1, barrier=(89.5551,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.547793,0.0685847,-5.08023e-05,2.03371e-08,-3.35663e-12,45320.9,28.6747], Tmin=(100,'K'), Tmax=(1419.14,'K')), NASAPolynomial(coeffs=[13.6134,0.0317576,-1.18765e-05,2.05084e-09,-1.3524e-13,41612.6,-38.9326], Tmin=(1419.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC1[C]CCC1(186)',
    structure = SMILES('C=CC1[C]CCC1'),
    E0 = (336.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64708,0.0310185,6.76139e-05,-1.06291e-07,4.16937e-11,40553.4,22.2396], Tmin=(100,'K'), Tmax=(973.119,'K')), NASAPolynomial(coeffs=[14.3654,0.0303717,-1.09764e-05,2.07413e-09,-1.53306e-13,35633.4,-51.3336], Tmin=(973.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C=CC1C[C]CC1(187)',
    structure = SMILES('C=CC1C[C]CC1'),
    E0 = (325.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.628,0.031315,6.74292e-05,-1.06583e-07,4.19212e-11,39303.8,22.3406], Tmin=(100,'K'), Tmax=(972.075,'K')), NASAPolynomial(coeffs=[14.5506,0.0300787,-1.08097e-05,2.04096e-09,-1.51007e-13,34337.5,-52.2596], Tmin=(972.075,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C[C]C1C=CCC1(188)',
    structure = SMILES('C[C]C1C=CCC1'),
    E0 = (322.947,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64397,0.032184,6.20534e-05,-9.90897e-08,3.88124e-11,38944,21.488], Tmin=(100,'K'), Tmax=(977.051,'K')), NASAPolynomial(coeffs=[13.7363,0.031405,-1.15572e-05,2.17908e-09,-1.59818e-13,34255.2,-48.4719], Tmin=(977.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]CC1C=CCC1(189)',
    structure = SMILES('[CH]CC1C=CCC1'),
    E0 = (397.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57445,0.0330059,6.27492e-05,-1.01638e-07,4.01203e-11,47857.1,21.9008], Tmin=(100,'K'), Tmax=(974.19,'K')), NASAPolynomial(coeffs=[14.4977,0.0304445,-1.10658e-05,2.08867e-09,-1.53981e-13,42942.7,-52.4081], Tmin=(974.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(397.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]1CC2[CH]C1CC2(190)',
    structure = SMILES('[CH]1CC2[CH]C1CC2'),
    E0 = (343.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.82119,-0.000636564,0.000145277,-1.77133e-07,6.42308e-11,41415.4,21.6034], Tmin=(100,'K'), Tmax=(975.933,'K')), NASAPolynomial(coeffs=[11.2195,0.0334231,-1.23273e-05,2.42807e-09,-1.84964e-13,36514.9,-35.4173], Tmin=(975.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(2-norbornyl) + radical(7-norbornyl)"""),
)

species(
    label = '[CH]1C=CC1(258)',
    structure = SMILES('[CH]1C=CC1'),
    E0 = (309.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,718.425,1076.94,1076.94,1076.94,1076.94,1076.94,1076.94,1076.94,1076.94,1076.94,1821.09],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65584,-0.010946,0.00010676,-1.29424e-07,4.82651e-11,37295.1,8.17286], Tmin=(100,'K'), Tmax=(946.825,'K')), NASAPolynomial(coeffs=[9.86439,0.010387,-2.38638e-06,4.80812e-10,-4.35568e-14,33987.5,-32.6988], Tmin=(946.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutene-allyl)"""),
)

species(
    label = 'C1=CC(C1)[C]1CC1(259)',
    structure = SMILES('C1=CC(C1)[C]1CC1'),
    E0 = (411.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08447,0.0265787,5.5502e-05,-8.15696e-08,3.03069e-11,49581.9,21.5484], Tmin=(100,'K'), Tmax=(1004.59,'K')), NASAPolynomial(coeffs=[9.97177,0.0339373,-1.33654e-05,2.54284e-09,-1.8393e-13,46041.2,-26.2747], Tmin=(1004.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Tertalkyl)"""),
)

species(
    label = 'C1=C[C](C1)C1CC1(260)',
    structure = SMILES('C1=C[C](C1)C1CC1'),
    E0 = (358.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27728,0.0137748,0.000107661,-1.45886e-07,5.56918e-11,43155.6,18.6662], Tmin=(100,'K'), Tmax=(958.93,'K')), NASAPolynomial(coeffs=[13.7593,0.0270739,-8.86436e-06,1.67359e-09,-1.27492e-13,38139.9,-50.9167], Tmin=(958.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_T)"""),
)

species(
    label = '[CH]1CC1C1C=CC1(261)',
    structure = SMILES('[CH]1CC1C1C=CC1'),
    E0 = (452.051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16755,0.0197135,8.41131e-05,-1.16995e-07,4.41568e-11,54453.5,22.4615], Tmin=(100,'K'), Tmax=(978.183,'K')), NASAPolynomial(coeffs=[12.3972,0.0298172,-1.10207e-05,2.12008e-09,-1.58067e-13,49967.6,-39.3657], Tmin=(978.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane)"""),
)

species(
    label = '[CH]1C=CC1C1CC1(262)',
    structure = SMILES('[CH]1C=CC1C1CC1'),
    E0 = (390.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20188,0.0139001,0.000111296,-1.5175e-07,5.80415e-11,47076,18.8838], Tmin=(100,'K'), Tmax=(961.636,'K')), NASAPolynomial(coeffs=[15.0271,0.025617,-8.47137e-06,1.63977e-09,-1.27365e-13,41601,-58.1311], Tmin=(961.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutene-allyl)"""),
)

species(
    label = '[C]1=CCC1C1CC1(263)',
    structure = SMILES('[C]1=CCC1C1CC1'),
    E0 = (478.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96748,0.0237619,7.74243e-05,-1.13822e-07,4.4026e-11,57657.1,21.3055], Tmin=(100,'K'), Tmax=(970.741,'K')), NASAPolynomial(coeffs=[13.8012,0.027592,-9.75934e-06,1.86229e-09,-1.39645e-13,52881.7,-48.1971], Tmin=(970.741,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[C]1=CC(C1)C1CC1(264)',
    structure = SMILES('[C]1=CC(C1)C1CC1'),
    E0 = (478.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96748,0.0237619,7.74243e-05,-1.13822e-07,4.4026e-11,57657.1,21.3055], Tmin=(100,'K'), Tmax=(970.741,'K')), NASAPolynomial(coeffs=[13.8012,0.027592,-9.75934e-06,1.86229e-09,-1.39645e-13,52881.7,-48.1971], Tmin=(970.741,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'C=CC1C=CC1(265)',
    structure = SMILES('C=CC1C=CC1'),
    E0 = (219.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26734,0.01959,6.97114e-05,-1.02576e-07,4.00366e-11,26421.7,17.9956], Tmin=(100,'K'), Tmax=(961.591,'K')), NASAPolynomial(coeffs=[12.7616,0.022543,-7.59802e-06,1.42673e-09,-1.07239e-13,22248.7,-43.4265], Tmin=(961.591,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]1C[C](C1)C1CC1(266)',
    structure = SMILES('[CH]1C[C](C1)C1CC1'),
    E0 = (425.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26239,0.0235711,6.26587e-05,-8.3059e-08,2.908e-11,51206.8,23.4974], Tmin=(100,'K'), Tmax=(1035.64,'K')), NASAPolynomial(coeffs=[7.59327,0.040686,-1.67407e-05,3.20661e-09,-2.30345e-13,48080.6,-12.1699], Tmin=(1035.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Tertalkyl) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1C[CH]C1C1CC1(267)',
    structure = SMILES('[CH]1C[CH]C1C1CC1'),
    E0 = (427.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34594,0.016722,9.11e-05,-1.18079e-07,4.26676e-11,51499.1,24.4074], Tmin=(100,'K'), Tmax=(998.766,'K')), NASAPolynomial(coeffs=[9.89274,0.0367693,-1.45089e-05,2.80976e-09,-2.06585e-13,47484.2,-24.5446], Tmin=(998.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1CC(C1)[C]1CC1(268)',
    structure = SMILES('[CH]1CC(C1)[C]1CC1'),
    E0 = (425.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26239,0.0235711,6.26587e-05,-8.3059e-08,2.908e-11,51206.8,23.4974], Tmin=(100,'K'), Tmax=(1035.64,'K')), NASAPolynomial(coeffs=[7.59327,0.040686,-1.67407e-05,3.20661e-09,-2.30345e-13,48080.6,-12.1699], Tmin=(1035.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Tertalkyl) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1CC[C]1C1CC1(269)',
    structure = SMILES('[CH]1CC[C]1C1CC1'),
    E0 = (425.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26239,0.0235711,6.26588e-05,-8.30591e-08,2.908e-11,51206.8,24.1905], Tmin=(100,'K'), Tmax=(1035.64,'K')), NASAPolynomial(coeffs=[7.59326,0.040686,-1.67407e-05,3.20662e-09,-2.30345e-13,48080.6,-11.4767], Tmin=(1035.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutane) + radical(Tertalkyl)"""),
)

species(
    label = '[CH]1CCC1[C]1CC1(270)',
    structure = SMILES('[CH]1CCC1[C]1CC1'),
    E0 = (425.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26239,0.0235711,6.26587e-05,-8.3059e-08,2.908e-11,51206.8,23.4974], Tmin=(100,'K'), Tmax=(1035.64,'K')), NASAPolynomial(coeffs=[7.59327,0.040686,-1.67407e-05,3.20661e-09,-2.30345e-13,48080.6,-12.1699], Tmin=(1035.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Tertalkyl) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1CC(C1)C1[CH]C1(271)',
    structure = SMILES('[CH]1CC(C1)C1[CH]C1'),
    E0 = (465.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34594,0.016722,9.11e-05,-1.18079e-07,4.26676e-11,56078.4,24.4074], Tmin=(100,'K'), Tmax=(998.766,'K')), NASAPolynomial(coeffs=[9.89274,0.0367693,-1.45089e-05,2.80976e-09,-2.06585e-13,52063.5,-24.5446], Tmin=(998.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1CC1C1[CH]CC1(272)',
    structure = SMILES('[CH]1CC1C1[CH]CC1'),
    E0 = (465.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34594,0.016722,9.11e-05,-1.18079e-07,4.26676e-11,56078.4,24.4074], Tmin=(100,'K'), Tmax=(998.766,'K')), NASAPolynomial(coeffs=[9.89274,0.0367693,-1.45089e-05,2.80976e-09,-2.06585e-13,52063.5,-24.5446], Tmin=(998.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane) + radical(cyclobutane)"""),
)

species(
    label = '[CH2]C([CH2])C1C=CC1(273)',
    structure = SMILES('[CH2]C([CH2])C1C=CC1'),
    E0 = (474.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,2950,3050,3150,900,950,1000,1050,1100,1380,1390,370,380,2900,435,600,893.017,893.017,893.017,893.017,893.017,893.017,893.017,893.017,893.017,893.017,893.017,893.017,2336.88],'cm^-1')),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48771,0.0392203,3.99592e-05,-7.78546e-08,3.26208e-11,57191.1,25.5725], Tmin=(100,'K'), Tmax=(948.084,'K')), NASAPolynomial(coeffs=[13.0689,0.0304739,-9.67025e-06,1.67214e-09,-1.17752e-13,53192.2,-39.1907], Tmin=(948.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=CC([CH2])C1CC1(274)',
    structure = SMILES('[CH]=CC([CH2])C1CC1'),
    E0 = (529.38,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,3000,3100,440,815,1455,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,1380,1390,370,380,2900,435,600,831.352,831.352,831.352,831.352,831.352,831.352,831.352,831.352,831.352,2391.1],'cm^-1')),
        HinderedRotor(inertia=(0.00664932,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00664932,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00664932,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28042,0.0443502,2.66989e-05,-6.37255e-08,2.70623e-11,63781.5,25.9729], Tmin=(100,'K'), Tmax=(970.743,'K')), NASAPolynomial(coeffs=[13.6751,0.0309276,-1.07385e-05,1.93916e-09,-1.38206e-13,59601.1,-42.5937], Tmin=(970.743,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[C]1CCC1C1CC1(275)',
    structure = SMILES('[C]1CCC1C1CC1'),
    E0 = (461.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92373,0.0213567,9.68216e-05,-1.37703e-07,5.30858e-11,55570,21.5028], Tmin=(100,'K'), Tmax=(964.906,'K')), NASAPolynomial(coeffs=[14.9395,0.0287086,-9.91474e-06,1.89137e-09,-1.43174e-13,50204.2,-55.6149], Tmin=(964.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[C]1CC(C1)C1CC1(276)',
    structure = SMILES('[C]1CC(C1)C1CC1'),
    E0 = (451.551,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89631,0.0219127,9.59076e-05,-1.37259e-07,5.30599e-11,54406.8,20.4973], Tmin=(100,'K'), Tmax=(964.159,'K')), NASAPolynomial(coeffs=[15.1292,0.0284212,-9.75405e-06,1.85864e-09,-1.40836e-13,49000.8,-57.6627], Tmin=(964.159,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
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
    E0 = (620.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (627.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (514.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (632.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (661.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (582.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (587.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (641.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (571.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (417.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (524.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (530.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (496.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (564.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (580.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (503.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (445.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (450.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (404.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (339.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (498.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (430.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (347.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (421.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (330.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (471.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (512.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (526.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (525.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (305.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (347.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (518.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (518.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (449.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (604.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (472.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (414.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (530.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (507.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (596.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (741.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (507.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (411.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (475.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (427.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (741.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (812.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (548.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (462.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (470.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (394.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (347.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (347.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (330.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (672.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (902.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (609.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (630.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (630.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (622.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (627.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (434.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (502.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (569.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (349.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (361.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (340.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (329.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (518.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (560.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (612.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (585.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (519.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (496.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (474.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (385.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (517.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (483.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS79',
    E0 = (312.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS80',
    E0 = (406.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS81',
    E0 = (446.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS82',
    E0 = (348.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS83',
    E0 = (338.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS84',
    E0 = (250.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS85',
    E0 = (354.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS86',
    E0 = (364.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS87',
    E0 = (540.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS88',
    E0 = (471.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS89',
    E0 = (364.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS90',
    E0 = (391.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS91',
    E0 = (190.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS92',
    E0 = (452.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS93',
    E0 = (514.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS94',
    E0 = (506.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS95',
    E0 = (572.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS96',
    E0 = (361.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS97',
    E0 = (402.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS98',
    E0 = (383.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS99',
    E0 = (328.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS100',
    E0 = (334.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS101',
    E0 = (411.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS102',
    E0 = (241.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS103',
    E0 = (475.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS104',
    E0 = (441.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS105',
    E0 = (429.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS106',
    E0 = (518.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS107',
    E0 = (594.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS108',
    E0 = (595.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS109',
    E0 = (518.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS110',
    E0 = (644.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS111',
    E0 = (555.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS112',
    E0 = (560.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS113',
    E0 = (591.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS114',
    E0 = (590.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS115',
    E0 = (526.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS116',
    E0 = (483.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS117',
    E0 = (573.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS118',
    E0 = (573.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS119',
    E0 = (683.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS120',
    E0 = (411.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS121',
    E0 = (411.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS122',
    E0 = (372.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS123',
    E0 = (500.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS124',
    E0 = (355.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS125',
    E0 = (880.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS126',
    E0 = (867.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS127',
    E0 = (664.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS128',
    E0 = (651.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS129',
    E0 = (564.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS130',
    E0 = (599.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS131',
    E0 = (548.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS132',
    E0 = (470.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS133',
    E0 = (347.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS134',
    E0 = (670.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS135',
    E0 = (581.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS136',
    E0 = (708.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS137',
    E0 = (537.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS138',
    E0 = (571.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS139',
    E0 = (418.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS140',
    E0 = (369.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS141',
    E0 = (475.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS142',
    E0 = (420.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS143',
    E0 = (545.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS144',
    E0 = (545.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS145',
    E0 = (526.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS146',
    E0 = (535.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS147',
    E0 = (283.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS148',
    E0 = (285.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS149',
    E0 = (365.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS150',
    E0 = (403.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS151',
    E0 = (499.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS152',
    E0 = (484.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS153',
    E0 = (413.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS154',
    E0 = (337.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS155',
    E0 = (439.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS156',
    E0 = (433.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS157',
    E0 = (358.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS158',
    E0 = (318.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS159',
    E0 = (254.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS160',
    E0 = (419.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS161',
    E0 = (428.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS162',
    E0 = (307.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS163',
    E0 = (408.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS164',
    E0 = (428.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS165',
    E0 = (490.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS166',
    E0 = (433.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS167',
    E0 = (382.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS168',
    E0 = (427.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS169',
    E0 = (417.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS170',
    E0 = (414.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS171',
    E0 = (453.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS172',
    E0 = (438.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS173',
    E0 = (319.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS174',
    E0 = (590.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS175',
    E0 = (623.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS176',
    E0 = (570.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS177',
    E0 = (663.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS178',
    E0 = (608.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS179',
    E0 = (690.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS180',
    E0 = (690.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS181',
    E0 = (645.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS182',
    E0 = (448.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS183',
    E0 = (450.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS184',
    E0 = (488.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS185',
    E0 = (490.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS186',
    E0 = (488.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS187',
    E0 = (433.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS188',
    E0 = (473.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS189',
    E0 = (473.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS190',
    E0 = (482.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS191',
    E0 = (471.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS192',
    E0 = (537.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS193',
    E0 = (524.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS194',
    E0 = (552.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS195',
    E0 = (542.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction79',
    reactants = ['H(10)', 'C=CC=CC1[CH]C1(122)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction1',
    reactants = ['cC3H5(120)', 'CH2CHCHCH(119)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.79546e+07,'m^3/(mol*s)'), n=-0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/NonDeC] for rate rule [Cd_pri_rad;C_rad/H/NonDeC]
Euclidian distance = 2.0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction78',
    reactants = ['H(10)', 'C=CC=C[C]1CC1(121)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction80',
    reactants = ['H(10)', 'C=CC=[C]C1CC1(123)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction81',
    reactants = ['C2H3(32)', '[CH]=CC1CC1(124)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 89 used for Cd_pri_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction82',
    reactants = ['H(10)', 'C=C[C]=CC1CC1(125)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction83',
    reactants = ['H(10)', 'C=[C]C=CC1CC1(126)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction84',
    reactants = ['H(10)', '[CH]=CC=CC1CC1(127)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction85',
    reactants = ['CH2(S)(28)', 'C=CC=CC=C(114)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.13244e+08,'m^3/(mol*s)'), n=-0.441667, Ea=(7.08269,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;mb_db]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction86',
    reactants = ['C=C[CH]C[C]1CC1(128)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.64935e+10,'s^-1'), n=0.2551, Ea=(25.8154,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction87',
    reactants = ['[CH2]C[C]=CC1CC1(129)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction88',
    reactants = ['C=[C]C[CH]C1CC1(130)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction89',
    reactants = ['C=C[CH]CC1[CH]C1(131)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction90',
    reactants = ['[CH2]CC=[C]C1CC1(132)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction91',
    reactants = ['[CH]=CC[CH]C1CC1(133)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction92',
    reactants = ['C=CC[CH]C1[CH]C1(134)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction93',
    reactants = ['C[CH]C=[C]C1CC1(135)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction94',
    reactants = ['[CH]=C[CH]CC1CC1(136)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction95',
    reactants = ['[CH2]CC=C[C]1CC1(137)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction96',
    reactants = ['C[CH]C=C[C]1CC1(138)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction97',
    reactants = ['[CH2]CC=CC1[CH]C1(95)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction98',
    reactants = ['C[CH]C=CC1[CH]C1(139)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad_NDe;XH_Rrad] for rate rule [R6radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction99',
    reactants = ['CC=CC=C1CC1(140)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.02873e+09,'s^-1'), n=1.23767, Ea=(163.714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH3_1;Cd(C)C_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction100',
    reactants = ['[CH2]C([CH2])C=CC=C(141)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction76',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction101',
    reactants = ['C=CC[C]C1CC1(142)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;CsJ2C;CH2(C=C)] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C=C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction102',
    reactants = ['C=C[C]CC1CC1(143)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.22315e+12,'s^-1'), n=0.271316, Ea=(58.4033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2(C=C);CH2(C)] + [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction103',
    reactants = ['C[C]C=CC1CC1(144)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.08234e+13,'s^-1'), n=0.129132, Ea=(99.579,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction104',
    reactants = ['[CH]CC=CC1CC1(145)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;singletcarbene;CH2(C=C)] for rate rule [CsJ2-C;CsJ2H;CH2(C=C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction105',
    reactants = ['C=CC=CC1CC1(118)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH2_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.34e+11,'s^-1'), n=0.21, Ea=(25.3011,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 124 used for R4_S_D;doublebond_intra_HCd_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_HCd_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 23.4 to 25.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C[CH]C1C=CC1(78)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra_HNd_pri;radadd_intra_cs] for rate rule [R5_SD_D;doublebond_intra_HNd_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2][CH]C1C=CCC1(79)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R6_SSM_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(10)', '[CH2]C=CC=CC=C(80)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.24e+08,'cm^3/(mol*s)'), n=1.64, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2580 used for Cds-CdH_Cds-HH;HJ
Exact match found for rate rule [Cds-CdH_Cds-HH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(10)', '[CH2]CC=CC=C=C(81)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=CC=C[CH2](82)', 'C2H4(33)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(28600,'cm^3/(mol*s)'), n=2.41, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 234 used for Cds-HH_Cds-HH;CdsJ-H
Exact match found for rate rule [Cds-HH_Cds-HH;CdsJ-H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C=CC=C[CH]C(83)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(9.08185e+06,'s^-1'), n=1.84946, Ea=(92.0377,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_single;Cs_H_out_H/(Cd-Cd-Cd)] + [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]CC=CC=[C]C(84)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C=CC=[C]CC(85)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]CC=C[C]=CC(86)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 288 used for R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C=C[C]=CCC(87)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.71035e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]CC=[C]C=CC(88)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.33e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 284 used for R4H_SDS;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C=[C]C=CCC(89)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(484628,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;C_rad_out_2H;XH_out] for rate rule [R5H_SSMS;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C[C]=CC=CC(90)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSMS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2][C]=CC=CCC(91)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.27582e+06,'s^-1'), n=1.20683, Ea=(76.1033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_RSMSR;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC=C[CH2](82)', 'C2H4(T)(92)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(7.76856e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_pri_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=CC[CH2](94)', '[CH]=C[CH2](93)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using an average for rate rule [Cd_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]CC=CC1[CH]C1(95)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 142 used for R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C=CC1[CH]CC1(96)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(6.76e+10,'s^-1'), n=0.19, Ea=(139.746,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 232 used for R4_Cs_HH_D;doublebond_intra_pri_HCd;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_HH_D;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]CC1[CH]C=CC1(97)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(9.69966e+09,'s^-1'), n=0.160726, Ea=(148.069,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 95 used for R5_SD_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H
Exact match found for rate rule [R5_SD_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C1[CH]C=CCC1(98)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.16096e+09,'s^-1'), n=0.35461, Ea=(72.1771,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSM_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R6_SSM_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['C=CC=CC=CC(99)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['C=C=CC=CCC(100)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['C7H10(1)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R7;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction56',
    reactants = ['CH2(T)(22)', '[CH2]C=CC=C[CH2](101)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(6.30919e+07,'m^3/(mol*s)'), n=0.0885113, Ea=(33.9239,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction57',
    reactants = ['CH2(T)(22)', '[CH]=CC=CC[CH2](102)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(3.1546e+07,'m^3/(mol*s)'), n=0.0885113, Ea=(33.9239,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction58',
    reactants = ['H(10)', '[CH2]CC=C=CC=C(103)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]CC=[C]CC=C(104)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2]CC=CC[C]=C(105)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2]C[C]=CCC=C(106)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH]=CCC=CC[CH2](107)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2][CH]C=CCC=C(108)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(6.92799e+06,'s^-1'), n=1.7075, Ea=(114.432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH]=C[CH]C=CCC(109)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R7HJ_2;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]CC1[CH]C1C=C(110)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(5.18204e+14,'s^-1'), n=-0.478585, Ea=(247.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;doublebond_intra_pri_HNd_Cs;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH]1C=C[CH]CCC1(51)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.2, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 46 used for R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['C=CC=CCC=C(111)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction68',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['C=CC=C=CCC(112)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction69',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R5_SSDS;C_rad_out_2H;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction70',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]CC1C=CC1[CH2](74)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra_csHNd] for rate rule [R5_SD_D;doublebond_intra_2H_pri;radadd_intra_csHNd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction71',
    reactants = ['CH2(T)(22)', 'C=CC=CC=C(114)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(0.18583,'m^3/(mol*s)'), n=2.36967, Ea=(33.8986,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-CdH;YJ] for rate rule [Cds-HH_Cds-CdH;CH2_triplet]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH2]CC[C]=CC=C(115)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction73',
    reactants = ['[CH2]CCC=[C]C=C(116)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(6.1583e+09,'s^-1'), n=0.92705, Ea=(170.178,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] + [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction74',
    reactants = ['[CH2]CCC=C[C]=C(117)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(2.22e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_Cd;Cs_H_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction75',
    reactants = ['[CH]=CC=CCC[CH2](55)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R5H_DSMS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction1',
    reactants = ['H(10)', '[CH]1CC=CC=CC1(42)'],
    products = ['C7H10(1)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(10)', '[CH]1C=CC=CCC1(43)'],
    products = ['C7H10(1)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(5.32655e+09,'m^3/(mol*s)'), n=-0.416379, Ea=(6.36255,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(10)', '[C]1=CC=CCCC1(44)'],
    products = ['C7H10(1)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(10)', '[C]1=CCCCC=C1(45)'],
    products = ['C7H10(1)'],
    transitionState = 'TS78',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction79',
    reactants = ['[CH]1C=CCC[CH]C1(46)'],
    products = ['C7H10(1)'],
    transitionState = 'TS79',
    kinetics = Arrhenius(A=(1.64935e+10,'s^-1'), n=0.2551, Ea=(25.8154,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[C]1=CCCC[CH]C1(47)'],
    products = ['C7H10(1)'],
    transitionState = 'TS80',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]1=CC[CH]CCC1(48)'],
    products = ['C7H10(1)'],
    transitionState = 'TS81',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]1C[CH]CC=CC1(49)'],
    products = ['C7H10(1)'],
    transitionState = 'TS82',
    kinetics = Arrhenius(A=(3.104e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]1=C[CH]CCCC1(50)'],
    products = ['C7H10(1)'],
    transitionState = 'TS83',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]1C=C[CH]CCC1(51)'],
    products = ['C7H10(1)'],
    transitionState = 'TS84',
    kinetics = Arrhenius(A=(1.036e+09,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_De] for rate rule [R4radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[C]1[CH]CCCCC=1(52)'],
    products = ['C7H10(1)'],
    transitionState = 'TS85',
    kinetics = Arrhenius(A=(8.50442e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]1[CH]CCC=CC1(53)'],
    products = ['C7H10(1)'],
    transitionState = 'TS86',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CCCCC=[CH](54)'],
    products = ['C7H10(1)'],
    transitionState = 'TS87',
    kinetics = Arrhenius(A=(1.35773e+13,'s^-1'), n=0.0154583, Ea=(6.01101,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Ypri_rad_out] for rate rule [R7;CdsingleH_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=CC=CCC[CH2](55)'],
    products = ['C7H10(1)'],
    transitionState = 'TS88',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Ypri_rad_out] for rate rule [R7;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[C]1CC=CCCC1(57)'],
    products = ['C7H10(1)'],
    transitionState = 'TS89',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;CsJ2C;CH2(C=C)] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C=C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]1C=CCCCC1(58)'],
    products = ['C7H10(1)'],
    transitionState = 'TS90',
    kinetics = Arrhenius(A=(2.22315e+12,'s^-1'), n=0.271316, Ea=(58.4033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2(C=C);CH2(C)] + [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C7H10(1)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS91',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(10)', 'C1=CC2CCC[C]12(64)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS92',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(10)', '[CH]1CCC2C=CC12(65)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS93',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(10)', '[CH]1CC2C=CC2C1(66)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS94',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(10)', '[C]1=CC2CCCC12(67)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS95',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]1C[C]2CCCC12(68)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS96',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]1CC2CCC[C]12(69)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS97',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]1CCC2[CH]CC12(70)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS98',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]1CCC2C[CH]C12(71)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS99',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]1CC2[CH]CC2C1(72)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS100',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=CC1[CH]CCC1(73)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS101',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]1C=C[CH]CCC1(51)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS102',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Cpri_rad_out_single] + [R4;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SDS;C_rad_out_H/NonDeC;Cpri_rad_out_H/NonDeC]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]CC1C=CC1[CH2](74)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS103',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for R5_SSSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R5_SSSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]CCC1[CH]C=C1(75)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS104',
    kinetics = Arrhenius(A=(2.49159e+11,'s^-1'), n=0.1555, Ea=(7.322,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R5_SSSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R5_SSSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[C]1CC2CCCC12(76)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS105',
    kinetics = Arrhenius(A=(1.83662e+12,'s^-1'), n=0.345439, Ea=(91.3359,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction106',
    reactants = ['H(10)', 'C=CC=C[C]1CC1(121)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS106',
    kinetics = Arrhenius(A=(7.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(4.93712,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2570 used for Cds-CsCs_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction107',
    reactants = ['H(10)', 'C=C[C]=CC1CC1(125)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS107',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction108',
    reactants = ['H(10)', 'C=[C]C=CC1CC1(126)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS108',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction109',
    reactants = ['C=C[CH]C[C]1CC1(128)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS109',
    kinetics = Arrhenius(A=(2.94e+08,'s^-1'), n=1.27, Ea=(125.938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_Cs2;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_Cs2_cy3;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction110',
    reactants = ['[CH2]C=[C]CC1CC1(147)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS110',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction111',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['C[C]=C[CH]C1CC1(148)'],
    transitionState = 'TS111',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction112',
    reactants = ['C=C[CH]CC1[CH]C1(131)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS112',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] for rate rule [R3H_SS_12cy3;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction113',
    reactants = ['[CH2][C]=CCC1CC1(149)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS113',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction114',
    reactants = ['CC=[C][CH]C1CC1(150)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS114',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction115',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['C[CH]C=C[C]1CC1(138)'],
    transitionState = 'TS115',
    kinetics = Arrhenius(A=(1.05265e+10,'s^-1'), n=0.795, Ea=(178.656,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_2H;Cs_H_out_Cs2_cy3] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_Cs2_cy3]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction116',
    reactants = ['C[CH]C=CC1[CH]C1(139)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS116',
    kinetics = Arrhenius(A=(138.3,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R6HJ_2;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction117',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['[CH2]C1[CH]C1C1CC1(151)'],
    transitionState = 'TS117',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction118',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['[CH]1CC1[CH]C1CC1(152)'],
    transitionState = 'TS118',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction119',
    reactants = ['CH2(S)(28)', '[CH2]C=CC=C[CH2](101)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS119',
    kinetics = Arrhenius(A=(3.13244e+08,'m^3/(mol*s)'), n=-0.441667, Ea=(7.08269,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;mb_db]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction120',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['CC=C=CC1CC1(153)'],
    transitionState = 'TS120',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction121',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['C=C=CCC1CC1(154)'],
    transitionState = 'TS121',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction122',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['CC=CC=C1CC1(140)'],
    transitionState = 'TS122',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction123',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['[CH2]C=CC1[CH]CC1(96)'],
    transitionState = 'TS123',
    kinetics = Arrhenius(A=(1.53073e+11,'s^-1'), n=0.611, Ea=(152.377,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ-CdH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction124',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS124',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Cpri_rad_out_2H] + [R4;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SDS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction125',
    reactants = ['[CH]=C[CH2](93)', '[CH]C1CC1(155)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS125',
    kinetics = Arrhenius(A=(3.1546e+07,'m^3/(mol*s)'), n=0.0885113, Ea=(33.9239,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction126',
    reactants = ['CH2(T)(22)', '[CH]=C[CH]C1CC1(156)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS126',
    kinetics = Arrhenius(A=(3.1546e+07,'m^3/(mol*s)'), n=0.0885113, Ea=(33.9239,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction127',
    reactants = ['[CH2]C[C]=CC1CC1(129)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS127',
    kinetics = Arrhenius(A=(9.93038e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction128',
    reactants = ['[CH2]CC=[C]C1CC1(132)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS128',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction129',
    reactants = ['[CH2]CC=C[C]1CC1(137)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS129',
    kinetics = Arrhenius(A=(4.34621e+09,'s^-1'), n=0.843951, Ea=(168.336,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_Cs2_cy3;XH_out] + [R4H_SDS;C_rad_out_single;XH_out] for rate rule [R4H_SDS;C_rad_out_Cs2_cy3;XH_out]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction130',
    reactants = ['C[CH]C=[C]C1CC1(135)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS130',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R4HJ_2;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction131',
    reactants = ['[CH2]CC=CC1[CH]C1(95)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS131',
    kinetics = Arrhenius(A=(296.998,'s^-1'), n=2.47528, Ea=(59.2966,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;C_rad_out_H/NonDeC;XH_out] + [R5H_SSMS;C_rad_out_single;XH_out] for rate rule [R5H_SSMS;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction132',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['[CH]1[CH]C(C1)C1CC1(157)'],
    transitionState = 'TS132',
    kinetics = Arrhenius(A=(1.61e+08,'s^-1'), n=0.96, Ea=(123.01,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction133',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS133',
    kinetics = Arrhenius(A=(5.01e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 5 used for Y_12_01
Exact match found for rate rule [Y_12_01]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction134',
    reactants = ['C=[C]C[CH]C1CC1(130)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS134',
    kinetics = Arrhenius(A=(9.93038e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction135',
    reactants = ['C=CC[CH][C]1CC1(158)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS135',
    kinetics = Arrhenius(A=(1.28e+07,'s^-1'), n=1.56, Ea=(126.775,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;C_rad_out_Cs2;Cs_H_out_H/Cd] for rate rule [R3HJ;C_rad_out_Cs2_cy3;Cs_H_out_H/Cd]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction136',
    reactants = ['[CH]=CC[CH]C1CC1(133)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS136',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction137',
    reactants = ['C=CC[CH]C1[CH]C1(134)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS137',
    kinetics = Arrhenius(A=(0.502,'s^-1'), n=3.86, Ea=(41.6308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] for rate rule [R4HJ_2;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction138',
    reactants = ['[CH]=C[CH]CC1CC1(136)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS138',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R4HJ_2;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction139',
    reactants = ['C2H3(32)', '[CH]1C=CCC1(159)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS139',
    kinetics = Arrhenius(A=(5.32655e+09,'m^3/(mol*s)'), n=-0.416379, Ea=(6.36255,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;C_rad/H/CdCs] for rate rule [Cd_pri_rad;C_rad/H/CdCs]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction140',
    reactants = ['H(10)', 'C=C[C]1C=CCC1(160)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS140',
    kinetics = Arrhenius(A=(3.62e+13,'cm^3/(mol*s)'), n=0.228, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 13 used for C_rad/TwoDeCs;H_rad
Exact match found for rate rule [C_rad/TwoDeCs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction141',
    reactants = ['H(10)', 'C=CC1[CH]CC=C1(161)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS141',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction142',
    reactants = ['H(10)', 'C=CC1C=C[CH]C1(162)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS142',
    kinetics = Arrhenius(A=(2.66327e+09,'m^3/(mol*s)'), n=-0.416379, Ea=(6.36255,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction143',
    reactants = ['H(10)', 'C=CC1[C]=CCC1(163)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS143',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction144',
    reactants = ['H(10)', 'C=CC1C=[C]CC1(164)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS144',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction145',
    reactants = ['H(10)', 'C=[C]C1C=CCC1(165)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS145',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction146',
    reactants = ['H(10)', '[CH]=CC1C=CCC1(166)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS146',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction147',
    reactants = ['C=C[C]1C[CH]CC1(167)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS147',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction148',
    reactants = ['[CH2]C[C]1C=CCC1(168)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS148',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction149',
    reactants = ['C=CC1[CH]C[CH]C1(169)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS149',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction150',
    reactants = ['[CH2]CC1[CH]CC=C1(170)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS150',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction151',
    reactants = ['[CH2]CC1[C]=CCC1(171)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS151',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction152',
    reactants = ['C=[C]C1C[CH]CC1(172)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS152',
    kinetics = Arrhenius(A=(6.94203e+09,'s^-1'), n=0.37, Ea=(99.6901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad_NDe] + [R3radExo;Y_rad;XH_Rrad_NDe] for rate rule [R3radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction153',
    reactants = ['C=CC1[CH]CC[CH]1(173)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS153',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction154',
    reactants = ['C[CH]C1[CH]CC=C1(174)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS154',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction155',
    reactants = ['C[CH]C1[C]=CCC1(175)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS155',
    kinetics = Arrhenius(A=(1.34494e+09,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction156',
    reactants = ['C=[C]C1[CH]CCC1(176)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS156',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction157',
    reactants = ['C=CC1C[CH][CH]C1(177)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS157',
    kinetics = Arrhenius(A=(2.05689e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction158',
    reactants = ['[CH2]CC1[CH]C=CC1(97)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS158',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction159',
    reactants = ['C=C[C]1[CH]CCC1(178)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS159',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction160',
    reactants = ['[CH2]CC1C=[C]CC1(179)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS160',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction161',
    reactants = ['[CH]=CC1C[CH]CC1(180)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS161',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction162',
    reactants = ['C[CH]C1C=C[CH]C1(181)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS162',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction163',
    reactants = ['C[CH]C1C=[C]CC1(182)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS163',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction164',
    reactants = ['[CH]=CC1[CH]CCC1(73)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS164',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction165',
    reactants = ['[CH]=CC(C=C)C[CH2](183)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS165',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Ypri_rad_out] for rate rule [R5_SSSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction166',
    reactants = ['[CH]=CCC[CH]C=C(184)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS166',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R5_SSSD;C_rad_out_H/OneDe;CdsinglepriH_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction167',
    reactants = ['[CH2]C=CC([CH2])C=C(185)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS167',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R5_SSDS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction168',
    reactants = ['C=CC1[C]CCC1(186)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS168',
    kinetics = Arrhenius(A=(1.83662e+12,'s^-1'), n=0.345439, Ea=(91.3359,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction169',
    reactants = ['C=CC1C[C]CC1(187)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS169',
    kinetics = Arrhenius(A=(1.83662e+12,'s^-1'), n=0.345439, Ea=(91.3359,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction170',
    reactants = ['C[C]C1C=CCC1(188)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS170',
    kinetics = Arrhenius(A=(2.75493e+12,'s^-1'), n=0.345439, Ea=(91.3359,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction171',
    reactants = ['[CH]CC1C=CCC1(189)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS171',
    kinetics = Arrhenius(A=(1.57909e+12,'s^-1'), n=0.325393, Ea=(56.3425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;singletcarbene;CH2(C)] + [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH2(C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction172',
    reactants = ['[CH]1CC2[CH]C1CC2(190)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS172',
    kinetics = Arrhenius(A=(1.31474e+12,'s^-1'), n=0.239323, Ea=(94.3216,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction173',
    reactants = ['[CH]1CCC2[CH]CC12(70)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS173',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R7JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction261',
    reactants = ['cC3H5(120)', '[CH]1C=CC1(258)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS174',
    kinetics = Arrhenius(A=(6.5e+14,'cm^3/(mol*s)','*|/',2), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [C_sec_rad;C_rad/H/NonDeC] for rate rule [C_rad/H/CdCs;C_rad/H/NonDeC]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction262',
    reactants = ['H(10)', 'C1=CC(C1)[C]1CC1(259)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS175',
    kinetics = Arrhenius(A=(6.68468e+06,'m^3/(mol*s)'), n=-0.0135, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/Cs3;Y_rad] for rate rule [C_rad/Cs3;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction263',
    reactants = ['H(10)', 'C1=C[C](C1)C1CC1(260)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS176',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction264',
    reactants = ['H(10)', '[CH]1CC1C1C=CC1(261)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS177',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction265',
    reactants = ['H(10)', '[CH]1C=CC1C1CC1(262)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS178',
    kinetics = Arrhenius(A=(5.32655e+09,'m^3/(mol*s)'), n=-0.416379, Ea=(6.36255,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction266',
    reactants = ['H(10)', '[C]1=CCC1C1CC1(263)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS179',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [H_rad;Cd_rad/NonDe]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction267',
    reactants = ['H(10)', '[C]1=CC(C1)C1CC1(264)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS180',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction268',
    reactants = ['CH2(S)(28)', 'C=CC1C=CC1(265)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS181',
    kinetics = Arrhenius(A=(1.56622e+08,'m^3/(mol*s)'), n=-0.441667, Ea=(7.08269,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;mb_db]
Euclidian distance = 0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction269',
    reactants = ['[CH]1C[C](C1)C1CC1(266)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS182',
    kinetics = Arrhenius(A=(1.0206e+11,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction270',
    reactants = ['[CH]1C[CH]C1C1CC1(267)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS183',
    kinetics = Arrhenius(A=(1.0206e+11,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction271',
    reactants = ['[CH]1CC(C1)[C]1CC1(268)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS184',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction272',
    reactants = ['[CH]1[CH]C(C1)C1CC1(157)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS185',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction273',
    reactants = ['[CH]1CC[C]1C1CC1(269)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS186',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction274',
    reactants = ['[CH]1CCC1[C]1CC1(270)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS187',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction275',
    reactants = ['[CH]1CC(C1)C1[CH]C1(271)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS188',
    kinetics = Arrhenius(A=(3.104e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction276',
    reactants = ['[CH]1CC1C1[CH]CC1(272)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS189',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction277',
    reactants = ['[CH2]C([CH2])C1C=CC1(273)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS190',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction278',
    reactants = ['[CH2]C[CH]C1C=CC1(78)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS191',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction279',
    reactants = ['[CH]=CC([CH2])C1CC1(274)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS192',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction280',
    reactants = ['[CH]=CC[CH]C1CC1(133)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS193',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction281',
    reactants = ['[C]1CCC1C1CC1(275)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS194',
    kinetics = Arrhenius(A=(1.83662e+12,'s^-1'), n=0.345439, Ea=(91.3359,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction282',
    reactants = ['[C]1CC(C1)C1CC1(276)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS195',
    kinetics = Arrhenius(A=(3.67324e+12,'s^-1'), n=0.345439, Ea=(91.3359,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

network(
    label = '19',
    isomers = [
        'C=CC=CC1CC1(118)',
        '[CH2]C=CC=CC[CH2](56)',
        'C7H10(1)',
        'C1=CC2CCCC12(59)',
        '[CH2]C=C[CH]C1CC1(77)',
        'C=CC1C=CCC1(113)',
        'C1=CC(C1)C1CC1(146)',
    ],
    reactants = [
        ('H(10)', 'C=CC=CC1[CH]C1(122)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '19',
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

