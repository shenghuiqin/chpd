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
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54617,0.0382553,3.8353e-05,-6.92969e-08,2.69165e-11,52320.7,22.9209], Tmin=(100,'K'), Tmax=(1008.64,'K')), NASAPolynomial(coeffs=[11.922,0.0352081,-1.37762e-05,2.60807e-09,-1.87922e-13,48289.5,-36.8331], Tmin=(1008.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(434.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Cds_S)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39128,0.0431071,2.41812e-05,-5.59193e-08,2.2691e-11,22194.4,22.215], Tmin=(100,'K'), Tmax=(1004.86,'K')), NASAPolynomial(coeffs=[12.0796,0.0344,-1.33349e-05,2.48346e-09,-1.76812e-13,18337.9,-37.9025], Tmin=(1004.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 2"""),
)

species(
    label = '[CH2]C=CC1[CH]CC1(96)',
    structure = SMILES('[CH2]C=CC1[CH]CC1'),
    E0 = (338.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77923,0.0293334,6.72036e-05,-1.01189e-07,3.86111e-11,40787.3,24.325], Tmin=(100,'K'), Tmax=(988.526,'K')), NASAPolynomial(coeffs=[12.8166,0.0333057,-1.26223e-05,2.41579e-09,-1.77687e-13,36228.9,-40.8158], Tmin=(988.526,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(cyclobutane)"""),
)

species(
    label = 'C1=CC(C1)C1CC1(146)',
    structure = SMILES('C1=CC(C1)C1CC1'),
    E0 = (226.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00395,0.0190803,0.000102584,-1.43669e-07,5.52787e-11,27292.8,20.706], Tmin=(100,'K'), Tmax=(962.107,'K')), NASAPolynomial(coeffs=[14.8853,0.0282669,-9.55752e-06,1.81747e-09,-1.38108e-13,21910.4,-56.0326], Tmin=(962.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
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
    label = '[CH]1[CH]C(C1)C1CC1(157)',
    structure = SMILES('[CH]1[CH]C(C1)C1CC1'),
    E0 = (427.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34594,0.016722,9.11e-05,-1.18079e-07,4.26676e-11,51499.1,24.4074], Tmin=(100,'K'), Tmax=(998.766,'K')), NASAPolynomial(coeffs=[9.89274,0.0367693,-1.45089e-05,2.80976e-09,-2.06585e-13,47484.2,-24.5446], Tmin=(998.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutane) + radical(cyclobutane)"""),
)

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
    E0 = (518.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (644.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (555.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (560.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (591.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (590.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (526.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (483.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (573.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (573.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (683.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (411.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (411.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (372.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (500.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (355.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (880.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (867.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (664.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (651.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (564.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (599.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (548.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (470.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (347.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (670.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (581.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (708.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (537.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (571.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction109',
    reactants = ['C=C[CH]C[C]1CC1(128)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.94e+08,'s^-1'), n=1.27, Ea=(125.938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_Cs2;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_Cs2_cy3;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction110',
    reactants = ['[CH2]C=[C]CC1CC1(147)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction111',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['C[C]=C[CH]C1CC1(148)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction112',
    reactants = ['C=C[CH]CC1[CH]C1(131)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] for rate rule [R3H_SS_12cy3;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction113',
    reactants = ['[CH2][C]=CCC1CC1(149)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction114',
    reactants = ['CC=[C][CH]C1CC1(150)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction115',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['C[CH]C=C[C]1CC1(138)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.05265e+10,'s^-1'), n=0.795, Ea=(178.656,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_2H;Cs_H_out_Cs2_cy3] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_Cs2_cy3]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction116',
    reactants = ['C[CH]C=CC1[CH]C1(139)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(138.3,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R6HJ_2;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction117',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['[CH2]C1[CH]C1C1CC1(151)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction118',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['[CH]1CC1[CH]C1CC1(152)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction119',
    reactants = ['CH2(S)(28)', '[CH2]C=CC=C[CH2](101)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.13244e+08,'m^3/(mol*s)'), n=-0.441667, Ea=(7.08269,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;mb_db]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction120',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['CC=C=CC1CC1(153)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction121',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['C=C=CCC1CC1(154)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction122',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['CC=CC=C1CC1(140)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction123',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['[CH2]C=CC1[CH]CC1(96)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.53073e+11,'s^-1'), n=0.611, Ea=(152.377,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ-CdH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction124',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Cpri_rad_out_2H] + [R4;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SDS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction125',
    reactants = ['[CH]=C[CH2](93)', '[CH]C1CC1(155)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.1546e+07,'m^3/(mol*s)'), n=0.0885113, Ea=(33.9239,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction126',
    reactants = ['CH2(T)(22)', '[CH]=C[CH]C1CC1(156)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.1546e+07,'m^3/(mol*s)'), n=0.0885113, Ea=(33.9239,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction127',
    reactants = ['[CH2]C[C]=CC1CC1(129)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.93038e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction128',
    reactants = ['[CH2]CC=[C]C1CC1(132)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction129',
    reactants = ['[CH2]CC=C[C]1CC1(137)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34621e+09,'s^-1'), n=0.843951, Ea=(168.336,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_Cs2_cy3;XH_out] + [R4H_SDS;C_rad_out_single;XH_out] for rate rule [R4H_SDS;C_rad_out_Cs2_cy3;XH_out]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction130',
    reactants = ['C[CH]C=[C]C1CC1(135)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R4HJ_2;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction131',
    reactants = ['[CH2]CC=CC1[CH]C1(95)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(296.998,'s^-1'), n=2.47528, Ea=(59.2966,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;C_rad_out_H/NonDeC;XH_out] + [R5H_SSMS;C_rad_out_single;XH_out] for rate rule [R5H_SSMS;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction132',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['[CH]1[CH]C(C1)C1CC1(157)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.61e+08,'s^-1'), n=0.96, Ea=(123.01,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction133',
    reactants = ['[CH2]C=C[CH]C1CC1(77)'],
    products = ['C=CC=CC1CC1(118)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.01e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 5 used for Y_12_01
Exact match found for rate rule [Y_12_01]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction134',
    reactants = ['C=[C]C[CH]C1CC1(130)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(9.93038e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction135',
    reactants = ['C=CC[CH][C]1CC1(158)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.28e+07,'s^-1'), n=1.56, Ea=(126.775,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;C_rad_out_Cs2;Cs_H_out_H/Cd] for rate rule [R3HJ;C_rad_out_Cs2_cy3;Cs_H_out_H/Cd]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction136',
    reactants = ['[CH]=CC[CH]C1CC1(133)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction137',
    reactants = ['C=CC[CH]C1[CH]C1(134)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.502,'s^-1'), n=3.86, Ea=(41.6308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] for rate rule [R4HJ_2;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction138',
    reactants = ['[CH]=C[CH]CC1CC1(136)'],
    products = ['[CH2]C=C[CH]C1CC1(77)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R4HJ_2;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = '27',
    isomers = [
        '[CH2]C=C[CH]C1CC1(77)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '27',
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

