species(
    label = 'C=CC=CC1CC1(162)',
    structure = SMILES('C=CC=CC1CC1'),
    E0 = (182.792,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3581.72,'J/mol'), sigma=(6.12627,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.46 K, Pc=35.35 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51961,0.0368125,4.71719e-05,-8.31518e-08,3.3139e-11,22090,22.4964], Tmin=(100,'K'), Tmax=(980.784,'K')), NASAPolynomial(coeffs=[13.4282,0.0318026,-1.17828e-05,2.20267e-09,-1.59767e-13,17659,-45.405], Tmin=(980.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH2]C=CC=CC[CH2](56)',
    structure = SMILES('[CH2]C=CC=CC[CH2]'),
    E0 = (322.527,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,607.147,610.769],'cm^-1')),
        HinderedRotor(inertia=(0.455517,'amu*angstrom^2'), symmetry=1, barrier=(12.4008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.045695,'amu*angstrom^2'), symmetry=1, barrier=(12.4085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.443928,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3519.58,'J/mol'), sigma=(6.10993,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.75 K, Pc=35.01 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.992989,0.0525503,3.92263e-06,-4.15537e-08,1.95233e-11,38911.4,28.6567], Tmin=(100,'K'), Tmax=(978.365,'K')), NASAPolynomial(coeffs=[14.2853,0.0302044,-1.08771e-05,1.96082e-09,-1.38209e-13,34779,-43.0068], Tmin=(978.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(RCCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C7H10(1)',
    structure = SMILES('C1C=CCCCC=1'),
    E0 = (68.1123,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3862.82,'J/mol'), sigma=(6.55596,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=603.36 K, Pc=31.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98725,0.0214489,9.22319e-05,-1.30114e-07,4.97548e-11,8285.01,18.591], Tmin=(100,'K'), Tmax=(968.932,'K')), NASAPolynomial(coeffs=[13.7238,0.0305102,-1.08312e-05,2.05834e-09,-1.53811e-13,3310.88,-51.5923], Tmin=(968.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.1123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3853.12,'J/mol'), sigma=(6.77267,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=601.85 K, Pc=28.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33763,0.00875843,0.00013271,-1.72329e-07,6.44576e-11,28131.9,15.3203], Tmin=(100,'K'), Tmax=(967.651,'K')), NASAPolynomial(coeffs=[14.2908,0.0302463,-1.05034e-05,2.05611e-09,-1.58631e-13,22499.3,-59.105], Tmin=(967.651,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Allyl_S)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30701,0.00901903,0.000132148,-1.74471e-07,6.6169e-11,13139.4,18.0522], Tmin=(100,'K'), Tmax=(958.829,'K')), NASAPolynomial(coeffs=[15.1122,0.0274004,-8.93448e-06,1.7229e-09,-1.34172e-13,7383.31,-60.3982], Tmin=(958.829,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C=CC1C=CCC1(157)',
    structure = SMILES('C=CC1C=CCC1'),
    E0 = (76.6777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3658.67,'J/mol'), sigma=(6.25269,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.48 K, Pc=33.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85069,0.0251389,8.29943e-05,-1.2173e-07,4.71151e-11,9319.49,21.0387], Tmin=(100,'K'), Tmax=(968.097,'K')), NASAPolynomial(coeffs=[14.1361,0.0299511,-1.05691e-05,1.99737e-09,-1.48775e-13,4336.59,-51.2831], Tmin=(968.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.6777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'H(10)',
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
    label = 'C=CC=[C]C1CC1(167)',
    structure = SMILES('C=CC=[C]C1CC1'),
    E0 = (420.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,1685,370,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47573,0.0415725,2.17877e-05,-5.31013e-08,2.18455e-11,50693.3,23.123], Tmin=(100,'K'), Tmax=(1002.7,'K')), NASAPolynomial(coeffs=[12.4275,0.0309896,-1.19065e-05,2.22931e-09,-1.59812e-13,46832.8,-38.0416], Tmin=(1002.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S)"""),
)

species(
    label = 'CH2CHCHCH(163)',
    structure = SMILES('[CH]=CC=C'),
    E0 = (346.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.31937,'amu*angstrom^2'), symmetry=1, barrier=(30.3349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64261,0.0163331,3.86249e-05,-6.71409e-08,2.83618e-11,41729.6,13.2818], Tmin=(100,'K'), Tmax=(937.715,'K')), NASAPolynomial(coeffs=[12.9703,0.00669155,-1.00087e-06,1.67641e-10,-1.71469e-14,38279.8,-43.9466], Tmin=(937.715,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH2CHCHCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'cC3H5(164)',
    structure = SMILES('[CH]1CC1'),
    E0 = (280.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,893.61,894.287,894.502,894.515,894.566,894.634,894.855,3153.39],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.52406,-0.00309014,7.19122e-05,-9.14765e-08,3.49759e-11,33782.5,9.37914], Tmin=(100,'K'), Tmax=(937.448,'K')), NASAPolynomial(coeffs=[9.0272,0.00813397,-1.57922e-06,2.7847e-10,-2.51187e-14,31225.8,-24.9485], Tmin=(937.448,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""cC3H5""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=CC=C[C]1CC1(165)',
    structure = SMILES('[CH2]C=CC=C1CC1'),
    E0 = (301.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,382.609,963.445,963.445,963.445,963.445,963.445,963.445,963.445,963.445,963.445,2324.01],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6829,0.0341209,4.7572e-05,-8.18605e-08,3.26296e-11,36385,22.7244], Tmin=(100,'K'), Tmax=(974.021,'K')), NASAPolynomial(coeffs=[12.6673,0.030576,-1.09786e-05,2.0257e-09,-1.46313e-13,32273.5,-40.101], Tmin=(974.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 2 + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=CC1[CH]C1(166)',
    structure = SMILES('C=CC=CC1[CH]C1'),
    E0 = (408.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67409,0.0375352,2.84871e-05,-5.63699e-08,2.20579e-11,49251.1,24.2857], Tmin=(100,'K'), Tmax=(1016.45,'K')), NASAPolynomial(coeffs=[11.0812,0.033121,-1.31154e-05,2.47498e-09,-1.77247e-13,45654.4,-29.5371], Tmin=(1016.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane)"""),
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
    label = '[CH]=CC1CC1(168)',
    structure = SMILES('[CH]=CC1CC1'),
    E0 = (374.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,180,180,180,1057.71,1057.87,1057.89,1058.03,1058.22,1058.36,1058.39],'cm^-1')),
        HinderedRotor(inertia=(0.00393939,'amu*angstrom^2'), symmetry=1, barrier=(3.13224,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54668,0.0169369,5.52945e-05,-8.23277e-08,3.22509e-11,45158.5,16.7784], Tmin=(100,'K'), Tmax=(961.592,'K')), NASAPolynomial(coeffs=[11.2743,0.0184604,-6.09065e-06,1.14084e-09,-8.58848e-14,41731.1,-34.0792], Tmin=(961.592,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_P)"""),
)

species(
    label = 'C=C[C]=CC1CC1(169)',
    structure = SMILES('[CH2]C=C=CC1CC1'),
    E0 = (366.463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,540,610,2055,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,1025.8,1025.8,1025.8,1025.8,1025.8,1025.8,1025.8,1025.8,1025.8,2277.05],'cm^-1')),
        HinderedRotor(inertia=(0.0228575,'amu*angstrom^2'), symmetry=1, barrier=(0.525538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0228575,'amu*angstrom^2'), symmetry=1, barrier=(0.525538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44774,0.0381078,4.19572e-05,-8.03531e-08,3.29877e-11,44183.4,23.9022], Tmin=(100,'K'), Tmax=(973.338,'K')), NASAPolynomial(coeffs=[14.8593,0.027313,-9.70971e-06,1.81732e-09,-1.33578e-13,39473.1,-51.2233], Tmin=(973.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C=CC1CC1(170)',
    structure = SMILES('C=C=C[CH]C1CC1'),
    E0 = (371.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2950,3100,1380,975,1025,1650,540,610,2055,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68933,0.0319602,5.62899e-05,-9.27328e-08,3.6726e-11,44801.6,21.1522], Tmin=(100,'K'), Tmax=(975.36,'K')), NASAPolynomial(coeffs=[13.9431,0.0284808,-1.02926e-05,1.944e-09,-1.4343e-13,40186.3,-49.065], Tmin=(975.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S)"""),
)

species(
    label = '[CH]=CC=CC1CC1(171)',
    structure = SMILES('[CH]=CC=CC1CC1'),
    E0 = (429.885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43908,0.0404558,3.00115e-05,-6.50366e-08,2.68807e-11,51809.5,23.1834], Tmin=(100,'K'), Tmax=(985.693,'K')), NASAPolynomial(coeffs=[13.7755,0.0287917,-1.0671e-05,1.99925e-09,-1.45069e-13,47512.2,-45.617], Tmin=(985.693,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_P)"""),
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
    label = 'C=CC=CC=C(158)',
    structure = SMILES('C=CC=CC=C'),
    E0 = (145.693,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08861,'amu*angstrom^2'), symmetry=1, barrier=(25.0292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08825,'amu*angstrom^2'), symmetry=1, barrier=(25.021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38433,0.0414405,1.94312e-05,-6.16722e-08,2.86519e-11,17631.7,19.3256], Tmin=(100,'K'), Tmax=(941.207,'K')), NASAPolynomial(coeffs=[16.9618,0.0157225,-4.10136e-06,6.95777e-10,-5.26556e-14,12906.2,-64.4088], Tmin=(941.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C[CH]C[C]1CC1(172)',
    structure = SMILES('[CH2]C=CC[C]1CC1'),
    E0 = (392.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2850,1437.5,1250,1305,750,350,180,950.966,950.966,950.966,950.966,950.966,950.966,950.966,950.966,950.966,2329.19],'cm^-1')),
        HinderedRotor(inertia=(0.0242043,'amu*angstrom^2'), symmetry=1, barrier=(0.556505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242043,'amu*angstrom^2'), symmetry=1, barrier=(0.556505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242043,'amu*angstrom^2'), symmetry=1, barrier=(0.556505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45601,0.0438563,1.59603e-05,-4.20481e-08,1.64086e-11,47265.8,25.2569], Tmin=(100,'K'), Tmax=(1054.5,'K')), NASAPolynomial(coeffs=[10.3953,0.0379036,-1.53394e-05,2.88146e-09,-2.03778e-13,43826.2,-25.7143], Tmin=(1054.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(Tertalkyl)"""),
)

species(
    label = '[CH2]C[C]=CC1CC1(173)',
    structure = SMILES('[CH2]C[C]=CC1CC1'),
    E0 = (501.582,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37791,0.0440777,1.96166e-05,-4.95161e-08,1.99407e-11,60432.6,26.459], Tmin=(100,'K'), Tmax=(1024.34,'K')), NASAPolynomial(coeffs=[11.8371,0.0350777,-1.38332e-05,2.60109e-09,-1.85559e-13,56619.2,-32.4063], Tmin=(1024.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = 'C=[C]C[CH]C1CC1(174)',
    structure = SMILES('C=[C]C[CH]C1CC1'),
    E0 = (507.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2950,3100,1380,975,1025,1650,1685,370,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,180,994.326,994.326,994.326,994.326,994.326,994.326,994.326,994.326,994.326,994.326,994.326,2272.79],'cm^-1')),
        HinderedRotor(inertia=(0.025325,'amu*angstrom^2'), symmetry=1, barrier=(0.582272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.025325,'amu*angstrom^2'), symmetry=1, barrier=(0.582272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.025325,'amu*angstrom^2'), symmetry=1, barrier=(0.582272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59838,0.0383502,3.31639e-05,-6.12409e-08,2.34198e-11,61138.9,26.7213], Tmin=(100,'K'), Tmax=(1027.34,'K')), NASAPolynomial(coeffs=[11.1222,0.0361481,-1.4548e-05,2.76861e-09,-1.9884e-13,57341.4,-28.4409], Tmin=(1027.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(Cs_S)"""),
)

species(
    label = 'C=C[CH]CC1[CH]C1(175)',
    structure = SMILES('[CH2]C=CCC1[CH]C1'),
    E0 = (432.633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,263.315,918.503,918.503,918.503,918.503,918.503,918.503,918.503,918.503,918.503,2336.11],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55414,0.0368432,4.49331e-05,-7.76857e-08,3.0224e-11,52136.7,26.1142], Tmin=(100,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[12.5773,0.0341802,-1.32164e-05,2.50985e-09,-1.82084e-13,47860.6,-37.4234], Tmin=(1000,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC=[C]C1CC1(176)',
    structure = SMILES('[CH2]CC=[C]C1CC1'),
    E0 = (501.582,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37791,0.0440777,1.96166e-05,-4.95161e-08,1.99407e-11,60432.6,26.459], Tmin=(100,'K'), Tmax=(1024.34,'K')), NASAPolynomial(coeffs=[11.8371,0.0350777,-1.38332e-05,2.60109e-09,-1.85559e-13,56619.2,-32.4063], Tmin=(1024.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CC[CH]C1CC1(177)',
    structure = SMILES('[CH]=CC[CH]C1CC1'),
    E0 = (516.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,385.661,981.517,981.517,981.517,981.517,981.517,981.517,981.517,981.517,981.517,981.517,2287.42],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56122,0.0372494,4.12812e-05,-7.29561e-08,2.83216e-11,62255.1,26.783], Tmin=(100,'K'), Tmax=(1007.16,'K')), NASAPolynomial(coeffs=[12.4181,0.034034,-1.33589e-05,2.54916e-09,-1.84956e-13,58044.4,-35.7197], Tmin=(1007.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_P) + radical(Cs_S)"""),
)

species(
    label = 'C=CC[CH]C1[CH]C1(178)',
    structure = SMILES('C=CC[CH]C1[CH]C1'),
    E0 = (495.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,421.991,964.946,964.946,964.946,964.946,964.946,964.946,964.946,964.946,964.946,964.946,2299.29],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79998,0.0342678,4.00547e-05,-6.48098e-08,2.37845e-11,59696.5,27.8728], Tmin=(100,'K'), Tmax=(1040.55,'K')), NASAPolynomial(coeffs=[9.80059,0.0382412,-1.57364e-05,3.00967e-09,-2.15908e-13,56151.4,-20.0777], Tmin=(1040.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(cyclopropane)"""),
)

species(
    label = 'C[CH]C=[C]C1CC1(179)',
    structure = SMILES('C[CH]C=[C]C1CC1'),
    E0 = (437.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,180,1075.71,1075.71,1075.71,1075.71,1075.71,1075.71,1075.71,1075.71,1075.71,1075.71,2260.84],'cm^-1')),
        HinderedRotor(inertia=(0.0555318,'amu*angstrom^2'), symmetry=1, barrier=(1.27679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0555318,'amu*angstrom^2'), symmetry=1, barrier=(1.27679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0555318,'amu*angstrom^2'), symmetry=1, barrier=(1.27679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54979,0.0380016,3.93728e-05,-7.06049e-08,2.74454e-11,52715,22.9483], Tmin=(100,'K'), Tmax=(1006.65,'K')), NASAPolynomial(coeffs=[12.0125,0.035007,-1.36528e-05,2.58441e-09,-1.86378e-13,48653.9,-37.3064], Tmin=(1006.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]CC1CC1(180)',
    structure = SMILES('[CH]C=CCC1CC1'),
    E0 = (425.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,1019.59,1019.59,1019.59,1019.59,1019.59,1019.59,1019.59,1019.59,1019.59,1019.59,1019.59,1019.59,1019.59,2260.83],'cm^-1')),
        HinderedRotor(inertia=(0.0361927,'amu*angstrom^2'), symmetry=1, barrier=(0.832142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0361927,'amu*angstrom^2'), symmetry=1, barrier=(0.832142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0361927,'amu*angstrom^2'), symmetry=1, barrier=(0.832142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3886,0.0386946,5.49248e-05,-9.10853e-08,3.53502e-11,51335.2,25.1888], Tmin=(100,'K'), Tmax=(990.079,'K')), NASAPolynomial(coeffs=[12.7884,0.0386534,-1.47273e-05,2.75665e-09,-1.98433e-13,46822.5,-41.0848], Tmin=(990.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]CC=C[C]1CC1(181)',
    structure = SMILES('[CH2]CC=C[C]1CC1'),
    E0 = (395.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,600,808.221,808.221,808.221,808.221,808.221,808.221,808.221,808.221,808.221,2400.06],'cm^-1')),
        HinderedRotor(inertia=(0.0040308,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0040308,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0040308,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70054,0.0339617,5.01821e-05,-8.17952e-08,3.15911e-11,47691.7,23.0792], Tmin=(100,'K'), Tmax=(990.882,'K')), NASAPolynomial(coeffs=[11.6116,0.0348615,-1.31081e-05,2.45182e-09,-1.76633e-13,43719.3,-34.7792], Tmin=(990.882,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(Allyl_T)"""),
)

species(
    label = 'C[CH]C=C[C]1CC1(182)',
    structure = SMILES('C[CH]C=C[C]1CC1'),
    E0 = (331.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,600,781.693,781.693,781.693,781.693,781.693,781.693,781.693,781.693,781.693,2411.16],'cm^-1')),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86686,0.0279433,6.97792e-05,-1.02749e-07,3.90741e-11,39974.4,19.589], Tmin=(100,'K'), Tmax=(982.592,'K')), NASAPolynomial(coeffs=[11.8583,0.034674,-1.28622e-05,2.41995e-09,-1.7621e-13,35722.5,-40.0832], Tmin=(982.592,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_T) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]CC=CC1[CH]C1(140)',
    structure = SMILES('[CH2]CC=CC1[CH]C1'),
    E0 = (489.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,937.063,937.063,937.063,937.063,937.063,937.063,937.063,937.063,937.063,2330.43],'cm^-1')),
        HinderedRotor(inertia=(0.0144844,'amu*angstrom^2'), symmetry=1, barrier=(0.333026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0144844,'amu*angstrom^2'), symmetry=1, barrier=(0.333026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0144844,'amu*angstrom^2'), symmetry=1, barrier=(0.333026,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57761,0.0400181,2.64249e-05,-5.29749e-08,2.02571e-11,58990.3,27.6173], Tmin=(100,'K'), Tmax=(1039.19,'K')), NASAPolynomial(coeffs=[10.5207,0.0371617,-1.50162e-05,2.84088e-09,-2.02521e-13,55427.1,-24.0726], Tmin=(1039.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]C=CC1[CH]C1(183)',
    structure = SMILES('CC=C[CH]C1[CH]C1'),
    E0 = (422.246,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,600,849.257,849.257,849.257,849.257,849.257,849.257,849.257,849.257,849.257,2373.29],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74747,0.0341555,4.52614e-05,-7.28226e-08,2.72356e-11,50879,24.0345], Tmin=(100,'K'), Tmax=(1020.26,'K')), NASAPolynomial(coeffs=[10.5777,0.0373332,-1.4981e-05,2.85311e-09,-2.0533e-13,47110,-28.384], Tmin=(1020.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(cyclopropane)"""),
)

species(
    label = 'CC=CC=C1CC1(184)',
    structure = SMILES('CC=CC=C1CC1'),
    E0 = (183.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39147,0.0431012,2.4198e-05,-5.5937e-08,2.26973e-11,22194.5,22.2187], Tmin=(100,'K'), Tmax=(1004.84,'K')), NASAPolynomial(coeffs=[12.0795,0.0343995,-1.33346e-05,2.4834e-09,-1.76809e-13,18337.9,-37.8988], Tmin=(1004.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 2"""),
)

species(
    label = '[CH2]C([CH2])C=CC=C(185)',
    structure = SMILES('[CH2]C([CH2])C=CC=C'),
    E0 = (413.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,358.23,358.525],'cm^-1')),
        HinderedRotor(inertia=(0.00130371,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143705,'amu*angstrom^2'), symmetry=1, barrier=(13.1622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144582,'amu*angstrom^2'), symmetry=1, barrier=(13.1697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.840124,'amu*angstrom^2'), symmetry=1, barrier=(76.6335,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736486,0.0603369,-1.494e-05,-2.88339e-08,1.77208e-11,49880.2,29.2977], Tmin=(100,'K'), Tmax=(918.927,'K')), NASAPolynomial(coeffs=[15.5206,0.0258831,-7.50641e-06,1.18159e-09,-7.82353e-14,45900.7,-47.646], Tmin=(918.927,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC[C]C1CC1(186)',
    structure = SMILES('C=CC[C]C1CC1'),
    E0 = (454.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34899,0.0413402,3.62459e-05,-7.27591e-08,2.96165e-11,54723.8,24.037], Tmin=(100,'K'), Tmax=(985.464,'K')), NASAPolynomial(coeffs=[13.8303,0.0318423,-1.19532e-05,2.23456e-09,-1.61466e-13,50265,-46.1365], Tmin=(985.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C=C[C]CC1CC1(187)',
    structure = SMILES('C=C[C]CC1CC1'),
    E0 = (454.022,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35614,0.0417958,3.34534e-05,-6.87377e-08,2.79034e-11,54715.9,24.741], Tmin=(100,'K'), Tmax=(990.9,'K')), NASAPolynomial(coeffs=[13.4052,0.032628,-1.24191e-05,2.32433e-09,-1.67404e-13,50390.2,-43.0566], Tmin=(990.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C[C]C=CC1CC1(188)',
    structure = SMILES('C[C]C=CC1CC1'),
    E0 = (426.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30705,0.0445434,2.28625e-05,-5.59038e-08,2.29196e-11,51463.4,23.9565], Tmin=(100,'K'), Tmax=(1004.54,'K')), NASAPolynomial(coeffs=[12.6896,0.0339622,-1.32166e-05,2.46975e-09,-1.76327e-13,47423.6,-39.7346], Tmin=(1004.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]CC=CC1CC1(189)',
    structure = SMILES('[CH]CC=CC1CC1'),
    E0 = (508.564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23923,0.0445235,2.78393e-05,-6.43545e-08,2.66984e-11,61279.9,23.6299], Tmin=(100,'K'), Tmax=(989.505,'K')), NASAPolynomial(coeffs=[13.9541,0.0319984,-1.21029e-05,2.25881e-09,-1.62549e-13,56860.5,-47.1941], Tmin=(989.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C1=CC(C1)C1CC1(190)',
    structure = SMILES('C1=CC(C1)C1CC1'),
    E0 = (226.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00349,0.0190918,0.00010255,-1.43633e-07,5.52654e-11,27292.7,20.7081], Tmin=(100,'K'), Tmax=(962.125,'K')), NASAPolynomial(coeffs=[14.8861,0.0282662,-9.55727e-06,1.81737e-09,-1.38097e-13,21910.2,-56.0347], Tmin=(962.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH2]C=C[CH]C1CC1(122)',
    structure = SMILES('[CH2]C=C[CH]C1CC1'),
    E0 = (347.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,600,889.969,889.969,889.969,889.969,889.969,889.969,889.969,889.969,889.969,889.969,2323.37],'cm^-1')),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3576.14,'J/mol'), sigma=(6.33006,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=558.58 K, Pc=31.99 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59545,0.0309575,7.2312e-05,-1.12596e-07,4.41233e-11,41940.9,22.5109], Tmin=(100,'K'), Tmax=(973.468,'K')), NASAPolynomial(coeffs=[15.1191,0.0301256,-1.07493e-05,2.0487e-09,-1.52953e-13,36714.3,-55.6904], Tmin=(973.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C[CH]C1C=CC1(123)',
    structure = SMILES('[CH2]C[CH]C1C=CC1'),
    E0 = (464.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7761,0.0318038,5.39804e-05,-8.39819e-08,3.16708e-11,55916,26.4958], Tmin=(100,'K'), Tmax=(1007.21,'K')), NASAPolynomial(coeffs=[11.6248,0.0350278,-1.38719e-05,2.66237e-09,-1.93823e-13,51784.6,-31.75], Tmin=(1007.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(464.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C1C=CCC1(124)',
    structure = SMILES('[CH2][CH]C1C=CCC1'),
    E0 = (347.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9669,0.0260017,7.01338e-05,-1.0027e-07,3.72998e-11,41926.5,25.3892], Tmin=(100,'K'), Tmax=(999.607,'K')), NASAPolynomial(coeffs=[11.4792,0.0349033,-1.36999e-05,2.64332e-09,-1.93939e-13,37678.4,-32.2351], Tmin=(999.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C=CC=CC=C(125)',
    structure = SMILES('[CH2]C=CC=CC=C'),
    E0 = (227.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.25172,'amu*angstrom^2'), symmetry=1, barrier=(28.7794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24928,'amu*angstrom^2'), symmetry=1, barrier=(28.7234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24986,'amu*angstrom^2'), symmetry=1, barrier=(28.7367,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00082,0.0483137,2.14331e-05,-6.83322e-08,3.15747e-11,27512.9,24.066], Tmin=(100,'K'), Tmax=(942.064,'K')), NASAPolynomial(coeffs=[17.7906,0.0213086,-6.07929e-06,1.0358e-09,-7.56134e-14,22384.4,-66.3624], Tmin=(942.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]CC=CC=C=C(126)',
    structure = SMILES('[CH2]CC=CC=C=C'),
    E0 = (381.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,540,610,2055,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.78762,'amu*angstrom^2'), symmetry=1, barrier=(18.1089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786465,'amu*angstrom^2'), symmetry=1, barrier=(18.0824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.785725,'amu*angstrom^2'), symmetry=1, barrier=(18.0654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.734247,0.0612921,-2.8976e-05,-6.54089e-09,6.98271e-12,45959.5,27.082], Tmin=(100,'K'), Tmax=(1020.39,'K')), NASAPolynomial(coeffs=[14.9233,0.0271344,-1.03164e-05,1.88298e-09,-1.31932e-13,41946.4,-47.134], Tmin=(1020.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(RCCJ)"""),
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
    label = '[CH]=CC=C[CH2](127)',
    structure = SMILES('[CH]=CC=C[CH2]'),
    E0 = (423.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(1.52957,'amu*angstrom^2'), symmetry=1, barrier=(35.1678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53261,'amu*angstrom^2'), symmetry=1, barrier=(35.2378,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22714,0.0272988,2.28188e-05,-5.2078e-08,2.29916e-11,50955.4,18.1255], Tmin=(100,'K'), Tmax=(937.472,'K')), NASAPolynomial(coeffs=[12.1493,0.0145305,-4.06026e-06,6.79535e-10,-4.91973e-14,47795.8,-36.0319], Tmin=(937.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=CC=C[CH]C(128)',
    structure = SMILES('[CH2]C=CC=C[CH]C'),
    E0 = (258.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16035,0.0465166,2.35889e-05,-6.26236e-08,2.70685e-11,31194,25.1629], Tmin=(100,'K'), Tmax=(970.334,'K')), NASAPolynomial(coeffs=[14.5416,0.0300014,-1.06227e-05,1.92702e-09,-1.3763e-13,26777.8,-48.366], Tmin=(970.334,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]CC=CC=[C]C(129)',
    structure = SMILES('[CH2]CC=CC=[C]C'),
    E0 = (442.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.520896,0.0677224,-4.89513e-05,1.84227e-08,-2.82292e-12,53330.4,29.2799], Tmin=(100,'K'), Tmax=(1528.13,'K')), NASAPolynomial(coeffs=[15.2041,0.0292873,-1.1223e-05,1.96293e-09,-1.30074e-13,48842.9,-47.784], Tmin=(1528.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC=[C]CC(130)',
    structure = SMILES('[CH2]C=CC=[C]CC'),
    E0 = (355.122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,952.275],'cm^-1')),
        HinderedRotor(inertia=(0.672776,'amu*angstrom^2'), symmetry=1, barrier=(15.4684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.674226,'amu*angstrom^2'), symmetry=1, barrier=(15.5018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400491,'amu*angstrom^2'), symmetry=1, barrier=(83.2973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.673519,'amu*angstrom^2'), symmetry=1, barrier=(15.4855,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.909141,0.0565739,-1.09454e-05,-2.38807e-08,1.27316e-11,42832.6,27.645], Tmin=(100,'K'), Tmax=(998.382,'K')), NASAPolynomial(coeffs=[13.4173,0.0318057,-1.18131e-05,2.12664e-09,-1.48019e-13,39071.9,-39.0058], Tmin=(998.382,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(Cds_S) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]CC=C[C]=CC(131)',
    structure = SMILES('[CH2]CC=C[C]=CC'),
    E0 = (403.466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,246.561,246.899],'cm^-1')),
        HinderedRotor(inertia=(0.282818,'amu*angstrom^2'), symmetry=1, barrier=(12.198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290513,'amu*angstrom^2'), symmetry=1, barrier=(12.2138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286326,'amu*angstrom^2'), symmetry=1, barrier=(12.2092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.7556,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.532197,0.0691234,-5.19104e-05,2.10378e-08,-3.51019e-12,48656.5,28.3542], Tmin=(100,'K'), Tmax=(1406.2,'K')), NASAPolynomial(coeffs=[13.7251,0.0315959,-1.188e-05,2.05992e-09,-1.36258e-13,44946.1,-39.7914], Tmin=(1406.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CJC=C) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C[C]=CCC(132)',
    structure = SMILES('[CH2]C=C[C]=CCC'),
    E0 = (316.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865961,0.0585054,-1.52556e-05,-2.01915e-08,1.1858e-11,38161.3,26.9224], Tmin=(100,'K'), Tmax=(978.011,'K')), NASAPolynomial(coeffs=[13.0041,0.0324953,-1.16116e-05,2.03327e-09,-1.39165e-13,34656.7,-37.1458], Tmin=(978.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]CC=[C]C=CC(133)',
    structure = SMILES('[CH2]CC=[C]C=CC'),
    E0 = (403.466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,245.966,245.99],'cm^-1')),
        HinderedRotor(inertia=(0.284261,'amu*angstrom^2'), symmetry=1, barrier=(12.2071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284338,'amu*angstrom^2'), symmetry=1, barrier=(12.207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284375,'amu*angstrom^2'), symmetry=1, barrier=(12.2069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.786,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.532197,0.0691234,-5.19104e-05,2.10378e-08,-3.51019e-12,48656.5,28.3542], Tmin=(100,'K'), Tmax=(1406.21,'K')), NASAPolynomial(coeffs=[13.7251,0.0315959,-1.188e-05,2.05992e-09,-1.36258e-13,44946.1,-39.7914], Tmin=(1406.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(RCCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=[C]C=CCC(134)',
    structure = SMILES('[CH2]C=[C]C=CCC'),
    E0 = (316.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865961,0.0585054,-1.52556e-05,-2.01915e-08,1.1858e-11,38161.3,26.9224], Tmin=(100,'K'), Tmax=(978.011,'K')), NASAPolynomial(coeffs=[13.0041,0.0324953,-1.16116e-05,2.03327e-09,-1.39165e-13,34656.7,-37.1458], Tmin=(978.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C[C]=CC=CC(135)',
    structure = SMILES('[CH2]C[C]=CC=CC'),
    E0 = (442.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,320.477,320.477],'cm^-1')),
        HinderedRotor(inertia=(0.158764,'amu*angstrom^2'), symmetry=1, barrier=(11.5711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158765,'amu*angstrom^2'), symmetry=1, barrier=(11.5711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158764,'amu*angstrom^2'), symmetry=1, barrier=(11.5711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158763,'amu*angstrom^2'), symmetry=1, barrier=(11.5711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.520896,0.0677224,-4.89513e-05,1.84227e-08,-2.82292e-12,53330.4,29.2799], Tmin=(100,'K'), Tmax=(1528.13,'K')), NASAPolynomial(coeffs=[15.2041,0.0292873,-1.1223e-05,1.96293e-09,-1.30074e-13,48842.9,-47.784], Tmin=(1528.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC=CCC(136)',
    structure = SMILES('C=[C]C=C[CH]CC'),
    E0 = (351.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,1685,370,180,446.19,446.272],'cm^-1')),
        HinderedRotor(inertia=(0.548575,'amu*angstrom^2'), symmetry=1, barrier=(77.478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.414351,'amu*angstrom^2'), symmetry=1, barrier=(13.6034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0962564,'amu*angstrom^2'), symmetry=1, barrier=(13.6034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167952,'amu*angstrom^2'), symmetry=1, barrier=(23.7282,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.708437,0.0622888,-2.60171e-05,-8.89255e-09,7.58583e-12,42412.4,25.2624], Tmin=(100,'K'), Tmax=(1012.75,'K')), NASAPolynomial(coeffs=[13.663,0.0319466,-1.19192e-05,2.12965e-09,-1.46745e-13,38720.5,-42.6722], Tmin=(1012.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Allyl_S)"""),
)

species(
    label = 'C2H4(T)(137)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,995.638,996.357,2246.15,2948.66],'cm^-1')),
        HinderedRotor(inertia=(0.00724744,'amu*angstrom^2'), symmetry=1, barrier=(5.084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40735,0.0100313,6.40897e-06,-1.41287e-08,5.92654e-12,38288.2,6.11706], Tmin=(100,'K'), Tmax=(954.266,'K')), NASAPolynomial(coeffs=[5.52251,0.00856169,-2.90741e-06,5.02348e-10,-3.44567e-14,37547.7,-5.75289], Tmin=(954.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C[CH2](138)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,227.3,229.293,232.862],'cm^-1')),
        HinderedRotor(inertia=(1.31489,'amu*angstrom^2'), symmetry=1, barrier=(50.4659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31911,0.00817972,3.34731e-05,-4.36187e-08,1.58211e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.758,'K')), NASAPolynomial(coeffs=[5.36759,0.0170742,-6.35103e-06,1.16619e-09,-8.27612e-14,44095,-3.4463], Tmin=(983.758,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC[CH2](139)',
    structure = SMILES('[CH]=CC[CH2]'),
    E0 = (435.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.207738,'amu*angstrom^2'), symmetry=1, barrier=(4.7763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204754,'amu*angstrom^2'), symmetry=1, barrier=(4.70771,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52058,0.0279427,-9.18462e-06,-4.19541e-09,2.67723e-12,52484.1,16.9976], Tmin=(100,'K'), Tmax=(1112.3,'K')), NASAPolynomial(coeffs=[7.71984,0.0177228,-6.83496e-06,1.24847e-09,-8.64351e-14,50803.1,-10.9966], Tmin=(1112.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=CC1[CH]CC1(141)',
    structure = SMILES('[CH2]C=CC1[CH]CC1'),
    E0 = (338.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77947,0.0293277,6.72206e-05,-1.01207e-07,3.86177e-11,40787.4,24.3257], Tmin=(100,'K'), Tmax=(988.511,'K')), NASAPolynomial(coeffs=[12.8161,0.0333065,-1.26227e-05,2.41589e-09,-1.77696e-13,36229.1,-40.812], Tmin=(988.511,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutane) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC1[CH]C=CC1(142)',
    structure = SMILES('[CH2]CC1[CH]C=CC1'),
    E0 = (279.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98149,0.0218491,9.08855e-05,-1.27232e-07,4.81861e-11,33706.5,21.4302], Tmin=(100,'K'), Tmax=(977.236,'K')), NASAPolynomial(coeffs=[13.4507,0.031766,-1.16166e-05,2.23722e-09,-1.67364e-13,28749.7,-47.5297], Tmin=(977.236,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-allyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1C=C[CH]CC1(66)',
    structure = SMILES('[CH2]C1[CH]C=CCC1'),
    E0 = (278.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98168,0.0222816,9.00546e-05,-1.27933e-07,4.92605e-11,33618,20.4273], Tmin=(100,'K'), Tmax=(962.164,'K')), NASAPolynomial(coeffs=[13.3362,0.0308834,-1.03562e-05,1.92089e-09,-1.42435e-13,28849.9,-47.3373], Tmin=(962.164,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclohexene-allyl) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC=CC=CC(143)',
    structure = SMILES('C=CC=CC=CC'),
    E0 = (109.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.707389,0.057383,-2.46618e-06,-4.14177e-08,2.1072e-11,13322.1,24.2541], Tmin=(100,'K'), Tmax=(964.432,'K')), NASAPolynomial(coeffs=[17.1632,0.0251951,-8.49311e-06,1.52032e-09,-1.08853e-13,8470.79,-63.2348], Tmin=(964.432,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=CC=CCC(144)',
    structure = SMILES('C=C=CC=CCC'),
    E0 = (175.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.701239,0.0605125,-1.83838e-05,-1.89409e-08,1.14648e-11,21276.7,25.4137], Tmin=(100,'K'), Tmax=(1007.05,'K')), NASAPolynomial(coeffs=[15.0975,0.0294814,-1.11145e-05,2.03249e-09,-1.43106e-13,17051.1,-50.725], Tmin=(1007.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = '[CH2]C=CC=C[CH2](145)',
    structure = SMILES('[CH2]C=CC=C[CH2]'),
    E0 = (257.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,467.253],'cm^-1')),
        HinderedRotor(inertia=(0.358058,'amu*angstrom^2'), symmetry=1, barrier=(55.8652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.35927,'amu*angstrom^2'), symmetry=1, barrier=(55.8751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.360823,'amu*angstrom^2'), symmetry=1, barrier=(55.8249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91532,0.0306751,4.13821e-05,-7.59867e-08,3.17688e-11,31116.9,21.5123], Tmin=(100,'K'), Tmax=(942.942,'K')), NASAPolynomial(coeffs=[12.739,0.0229472,-7.07087e-06,1.218e-09,-8.69952e-14,27378,-39.072], Tmin=(942.942,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC=CC[CH2](146)',
    structure = SMILES('[CH]=CC=CC[CH2]'),
    E0 = (487.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,266.637],'cm^-1')),
        HinderedRotor(inertia=(0.307685,'amu*angstrom^2'), symmetry=1, barrier=(15.2649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.306819,'amu*angstrom^2'), symmetry=1, barrier=(15.3068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308442,'amu*angstrom^2'), symmetry=1, barrier=(15.2691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30261,0.0492009,-1.47391e-05,-1.75105e-08,1.06853e-11,58750,24.5846], Tmin=(100,'K'), Tmax=(985.494,'K')), NASAPolynomial(coeffs=[13.7032,0.0217746,-7.859e-06,1.42057e-09,-1.00264e-13,55193.5,-40.7032], Tmin=(985.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC=C=CC=C(147)',
    structure = SMILES('[CH2]CC=C=CC=C'),
    E0 = (381.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,540,610,2055,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.78762,'amu*angstrom^2'), symmetry=1, barrier=(18.1089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786465,'amu*angstrom^2'), symmetry=1, barrier=(18.0824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.785725,'amu*angstrom^2'), symmetry=1, barrier=(18.0654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.734247,0.0612921,-2.8976e-05,-6.54089e-09,6.98271e-12,45959.5,27.082], Tmin=(100,'K'), Tmax=(1020.39,'K')), NASAPolynomial(coeffs=[14.9233,0.0271344,-1.03164e-05,1.88298e-09,-1.31932e-13,41946.4,-47.134], Tmin=(1020.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC=[C]CC=C(148)',
    structure = SMILES('[CH2]CC=[C]CC=C'),
    E0 = (473.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,197.473,198.421,4000],'cm^-1')),
        HinderedRotor(inertia=(0.781725,'amu*angstrom^2'), symmetry=1, barrier=(21.7963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.791016,'amu*angstrom^2'), symmetry=1, barrier=(21.7918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375545,'amu*angstrom^2'), symmetry=1, barrier=(10.3526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00818045,'amu*angstrom^2'), symmetry=1, barrier=(21.7868,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719683,0.0631224,-4.05446e-05,1.30994e-08,-1.71645e-12,57034.4,31.6017], Tmin=(100,'K'), Tmax=(1754.47,'K')), NASAPolynomial(coeffs=[15.7769,0.0287936,-1.11949e-05,1.94705e-09,-1.27305e-13,51750.9,-49.5052], Tmin=(1754.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC=CC[C]=C(149)',
    structure = SMILES('[CH2]CC=CC[C]=C'),
    E0 = (473.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,197.473,198.421,4000],'cm^-1')),
        HinderedRotor(inertia=(0.781725,'amu*angstrom^2'), symmetry=1, barrier=(21.7963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.791016,'amu*angstrom^2'), symmetry=1, barrier=(21.7918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375545,'amu*angstrom^2'), symmetry=1, barrier=(10.3526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00818045,'amu*angstrom^2'), symmetry=1, barrier=(21.7868,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719683,0.0631224,-4.05446e-05,1.30994e-08,-1.71645e-12,57034.4,31.6017], Tmin=(100,'K'), Tmax=(1754.47,'K')), NASAPolynomial(coeffs=[15.7769,0.0287936,-1.11949e-05,1.94705e-09,-1.27305e-13,51750.9,-49.5052], Tmin=(1754.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C[C]=CCC=C(150)',
    structure = SMILES('[CH2]C[C]=CCC=C'),
    E0 = (473.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,197.473,198.421,4000],'cm^-1')),
        HinderedRotor(inertia=(0.781725,'amu*angstrom^2'), symmetry=1, barrier=(21.7963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.791016,'amu*angstrom^2'), symmetry=1, barrier=(21.7918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375545,'amu*angstrom^2'), symmetry=1, barrier=(10.3526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00818045,'amu*angstrom^2'), symmetry=1, barrier=(21.7868,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719683,0.0631224,-4.05446e-05,1.30994e-08,-1.71645e-12,57034.4,31.6017], Tmin=(100,'K'), Tmax=(1754.47,'K')), NASAPolynomial(coeffs=[15.7769,0.0287936,-1.11949e-05,1.94705e-09,-1.27305e-13,51750.9,-49.5052], Tmin=(1754.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CCC=CC[CH2](151)',
    structure = SMILES('[CH]=CCC=CC[CH2]'),
    E0 = (482.421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,264.391,264.394],'cm^-1')),
        HinderedRotor(inertia=(0.00241164,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00241158,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27409,'amu*angstrom^2'), symmetry=1, barrier=(13.5965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274108,'amu*angstrom^2'), symmetry=1, barrier=(13.5965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.790519,0.0610148,-3.00894e-05,-9.34766e-11,3.21551e-12,58145.4,31.2573], Tmin=(100,'K'), Tmax=(1119.14,'K')), NASAPolynomial(coeffs=[13.2192,0.0324369,-1.30228e-05,2.39083e-09,-1.65477e-13,54371.3,-34.5362], Tmin=(1119.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]C=CCC=C(152)',
    structure = SMILES('[CH2]C=C[CH]CC=C'),
    E0 = (320.286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,381.612,382.11,382.446],'cm^-1')),
        HinderedRotor(inertia=(0.0011544,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00114613,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.299466,'amu*angstrom^2'), symmetry=1, barrier=(31.0678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.299606,'amu*angstrom^2'), symmetry=1, barrier=(31.0712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03015,0.052108,2.03912e-06,-3.57972e-08,1.61473e-11,38640,26.6827], Tmin=(100,'K'), Tmax=(1021.43,'K')), NASAPolynomial(coeffs=[13.6023,0.0325735,-1.28877e-05,2.41086e-09,-1.71396e-13,34522.4,-41.8212], Tmin=(1021.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C[CH]C=CCC(153)',
    structure = SMILES('[CH]C=CC=CCC'),
    E0 = (369.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.656756,0.060847,-9.09604e-06,-2.78403e-08,1.40427e-11,44621.2,27.6289], Tmin=(100,'K'), Tmax=(1011.62,'K')), NASAPolynomial(coeffs=[13.8763,0.0363662,-1.40024e-05,2.54771e-09,-1.77653e-13,40524.6,-43.3285], Tmin=(1011.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]CC1[CH]C1C=C(154)',
    structure = SMILES('[CH2]CC1[CH]C1C=C'),
    E0 = (497.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51908,0.0414961,2.29381e-05,-4.99253e-08,1.93249e-11,59910,27.2281], Tmin=(100,'K'), Tmax=(1040.16,'K')), NASAPolynomial(coeffs=[10.739,0.0368768,-1.48699e-05,2.80865e-09,-1.99996e-13,56323.9,-25.6346], Tmin=(1040.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(cyclopropane)"""),
)

species(
    label = 'C=CC=CCC=C(155)',
    structure = SMILES('C=CC=CCC=C'),
    E0 = (140.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.928238,0.0526928,5.52723e-06,-4.52814e-08,2.12059e-11,17024.7,26.4845], Tmin=(100,'K'), Tmax=(980.349,'K')), NASAPolynomial(coeffs=[15.6613,0.0277119,-1.00057e-05,1.83674e-09,-1.31767e-13,12447.7,-52.9123], Tmin=(980.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC=C=CCC(156)',
    structure = SMILES('C=CC=C=CCC'),
    E0 = (175.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.701239,0.0605125,-1.83838e-05,-1.89409e-08,1.14648e-11,21276.7,25.4137], Tmin=(100,'K'), Tmax=(1007.05,'K')), NASAPolynomial(coeffs=[15.0975,0.0294814,-1.11145e-05,2.03249e-09,-1.43106e-13,17051.1,-50.725], Tmin=(1007.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]CC1C=CC1[CH2](79)',
    structure = SMILES('[CH2]CC1C=CC1[CH2]'),
    E0 = (468.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,2231.58],'cm^-1')),
        HinderedRotor(inertia=(0.0951468,'amu*angstrom^2'), symmetry=1, barrier=(2.18761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0951468,'amu*angstrom^2'), symmetry=1, barrier=(2.18761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0951468,'amu*angstrom^2'), symmetry=1, barrier=(2.18761,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44248,0.0408177,3.32491e-05,-6.77816e-08,2.77761e-11,56447.6,26.1745], Tmin=(100,'K'), Tmax=(977.723,'K')), NASAPolynomial(coeffs=[12.6761,0.0324508,-1.15862e-05,2.11346e-09,-1.50784e-13,52454.1,-36.9565], Tmin=(977.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC[C]=CC=C(159)',
    structure = SMILES('[CH2]CC[C]=CC=C'),
    E0 = (454.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,208.485,449.907,455.182],'cm^-1')),
        HinderedRotor(inertia=(0.0851836,'amu*angstrom^2'), symmetry=1, barrier=(12.1932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0836752,'amu*angstrom^2'), symmetry=1, barrier=(12.1964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0833087,'amu*angstrom^2'), symmetry=1, barrier=(12.1951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0834059,'amu*angstrom^2'), symmetry=1, barrier=(12.1925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.538823,0.0668724,-4.27758e-05,9.85528e-09,5.459e-13,54803.1,29.6451], Tmin=(100,'K'), Tmax=(1112.5,'K')), NASAPolynomial(coeffs=[14.3354,0.0306229,-1.19084e-05,2.1493e-09,-1.47427e-13,50906.8,-42.101], Tmin=(1112.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCC=[C]C=C(160)',
    structure = SMILES('[CH2]CCC=[C]C=C'),
    E0 = (415.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,420.558,490.307,497.749],'cm^-1')),
        HinderedRotor(inertia=(0.0742585,'amu*angstrom^2'), symmetry=1, barrier=(12.962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.56337,'amu*angstrom^2'), symmetry=1, barrier=(12.953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.563327,'amu*angstrom^2'), symmetry=1, barrier=(12.952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.53499,'amu*angstrom^2'), symmetry=1, barrier=(104.268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.482421,0.068984,-4.78332e-05,1.46887e-08,-8.94686e-13,50132.3,28.9684], Tmin=(100,'K'), Tmax=(1096.54,'K')), NASAPolynomial(coeffs=[13.846,0.0314272,-1.1767e-05,2.06913e-09,-1.39605e-13,46528.7,-39.8025], Tmin=(1096.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCC=C[C]=C(161)',
    structure = SMILES('[CH2]CCC=C[C]=C'),
    E0 = (415.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,243.243,243.319,244.993],'cm^-1')),
        HinderedRotor(inertia=(0.378274,'amu*angstrom^2'), symmetry=1, barrier=(15.8794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00284945,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.379586,'amu*angstrom^2'), symmetry=1, barrier=(15.8798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.84744,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.482412,0.0689841,-4.78335e-05,1.46892e-08,-8.94868e-13,50132.3,28.9684], Tmin=(100,'K'), Tmax=(1096.54,'K')), NASAPolynomial(coeffs=[13.846,0.0314271,-1.1767e-05,2.06912e-09,-1.39604e-13,46528.7,-39.8029], Tmin=(1096.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CC=CCC[CH2](55)',
    structure = SMILES('[CH]=CC=CCC[CH2]'),
    E0 = (463.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,235.941,236.045],'cm^-1')),
        HinderedRotor(inertia=(0.00302811,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.390375,'amu*angstrom^2'), symmetry=1, barrier=(15.438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.390318,'amu*angstrom^2'), symmetry=1, barrier=(15.4384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.390266,'amu*angstrom^2'), symmetry=1, barrier=(15.4383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.549822,0.0652361,-3.29542e-05,-3.80032e-09,6.14858e-12,55917.2,29.532], Tmin=(100,'K'), Tmax=(1026.67,'K')), NASAPolynomial(coeffs=[15.217,0.0291836,-1.10964e-05,2.01692e-09,-1.40639e-13,51794,-47.0282], Tmin=(1026.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.253,0.021796,6.88126e-05,-9.4846e-08,3.49373e-11,31657.1,20.1973], Tmin=(100,'K'), Tmax=(994.535,'K')), NASAPolynomial(coeffs=[9.55786,0.0343077,-1.32408e-05,2.51009e-09,-1.81989e-13,28132.4,-25.4201], Tmin=(994.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28637,0.00959473,0.00012854,-1.73381e-07,6.69502e-11,20259.5,17.7359], Tmin=(100,'K'), Tmax=(947.918,'K')), NASAPolynomial(coeffs=[15.9797,0.0233646,-6.47538e-06,1.20573e-09,-9.64081e-14,14448.8,-64.5504], Tmin=(947.918,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CCJC=C)"""),
)

species(
    label = '[C]1=CC=CCCC1(44)',
    structure = SMILES('[C]1=CC=CCCC1'),
    E0 = (305.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94969,0.0261406,6.70518e-05,-1.00266e-07,3.85157e-11,36888.1,19.8877], Tmin=(100,'K'), Tmax=(980.34,'K')), NASAPolynomial(coeffs=[12.6586,0.029804,-1.10153e-05,2.09905e-09,-1.55011e-13,32512.7,-43.171], Tmin=(980.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89678,0.0281855,6.23559e-05,-9.61001e-08,3.74522e-11,32217.2,19.2001], Tmin=(100,'K'), Tmax=(970.584,'K')), NASAPolynomial(coeffs=[12.3028,0.0303972,-1.07588e-05,1.99277e-09,-1.45091e-13,28073,-41.6351], Tmin=(970.584,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CJC=C)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.40792,0.0142585,0.000100663,-1.29073e-07,4.69067e-11,34538.7,19.433], Tmin=(100,'K'), Tmax=(987.858,'K')), NASAPolynomial(coeffs=[9.96106,0.0367302,-1.40197e-05,2.69039e-09,-1.97948e-13,30457.7,-30.0171], Tmin=(987.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC) + radical(Allyl_S)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16762,0.024175,6.66823e-05,-9.1142e-08,3.29615e-11,46176.9,21.8767], Tmin=(100,'K'), Tmax=(1010.06,'K')), NASAPolynomial(coeffs=[8.84407,0.0385261,-1.5207e-05,2.88955e-09,-2.08316e-13,42747.4,-20.7004], Tmin=(1010.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC) + radical(Cds_S)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16762,0.024175,6.66823e-05,-9.1142e-08,3.29615e-11,46176.9,21.8767], Tmin=(100,'K'), Tmax=(1010.06,'K')), NASAPolynomial(coeffs=[8.84407,0.0385261,-1.5207e-05,2.88955e-09,-2.08316e-13,42747.4,-20.7004], Tmin=(1010.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC) + radical(Cds_S)"""),
)

species(
    label = '[CH]1C[CH]CC=CC1(49)',
    structure = SMILES('[CH]1C[CH]CC=CC1'),
    E0 = (339.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,301.663,801.854,801.854,801.854,801.854,801.854,801.854,801.854,801.854,801.854,801.854,801.854,801.854,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47137,0.0198095,6.85976e-05,-8.60544e-08,2.95889e-11,40945.9,22.1857], Tmin=(100,'K'), Tmax=(1032.63,'K')), NASAPolynomial(coeffs=[5.82997,0.0428902,-1.73551e-05,3.28286e-09,-2.33857e-13,38328.1,-3.44251], Tmin=(1032.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC) + radical(RCCJCC)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10288,0.018622,9.88448e-05,-1.34432e-07,5.04671e-11,39769.8,18.4366], Tmin=(100,'K'), Tmax=(977.308,'K')), NASAPolynomial(coeffs=[13.0785,0.0321989,-1.17787e-05,2.2757e-09,-1.7067e-13,34830.8,-48.5552], Tmin=(977.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(Allyl_S)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10288,0.018622,9.88448e-05,-1.34432e-07,5.04671e-11,39769.8,17.7434], Tmin=(100,'K'), Tmax=(977.308,'K')), NASAPolynomial(coeffs=[13.0785,0.0321989,-1.17787e-05,2.2757e-09,-1.7067e-13,34830.8,-49.2484], Tmin=(977.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = '[CH]1[CH]CCC=CC1(53)',
    structure = SMILES('[CH]1[CH]CCC=CC1'),
    E0 = (339.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,301.663,801.854,801.854,801.854,801.854,801.854,801.854,801.854,801.854,801.854,801.854,801.854,801.854,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39,1607.39],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47137,0.0198095,6.85976e-05,-8.60544e-08,2.95889e-11,40945.9,22.8788], Tmin=(100,'K'), Tmax=(1032.63,'K')), NASAPolynomial(coeffs=[5.82997,0.0428902,-1.73551e-05,3.28286e-09,-2.33857e-13,38328.1,-2.74936], Tmin=(1032.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC) + radical(RCCJCC)"""),
)

species(
    label = '[CH]=CCCCC=[CH](54)',
    structure = SMILES('[CH]=CCCCC=[CH]'),
    E0 = (534.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3115,3125,620,680,785,800,1600,1700,180,594.821],'cm^-1')),
        HinderedRotor(inertia=(0.533698,'amu*angstrom^2'), symmetry=1, barrier=(12.2708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0190427,'amu*angstrom^2'), symmetry=1, barrier=(12.2694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.533677,'amu*angstrom^2'), symmetry=1, barrier=(12.2703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.533606,'amu*angstrom^2'), symmetry=1, barrier=(12.2686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.475674,0.0683398,-4.89882e-05,1.80153e-08,-2.68644e-12,64373.3,30.0646], Tmin=(100,'K'), Tmax=(1568.66,'K')), NASAPolynomial(coeffs=[15.9415,0.0289025,-1.12769e-05,1.98826e-09,-1.32167e-13,59521.1,-51.5121], Tmin=(1568.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[C]1CC=CCCC1(57)',
    structure = SMILES('[C]1CC=CCCC1'),
    E0 = (347.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84214,0.0256056,8.14451e-05,-1.19213e-07,4.58824e-11,41870.3,19.8303], Tmin=(100,'K'), Tmax=(973.151,'K')), NASAPolynomial(coeffs=[13.886,0.0309084,-1.12076e-05,2.1332e-09,-1.58593e-13,36931,-51.2804], Tmin=(973.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[C]1C=CCCCC1(58)',
    structure = SMILES('[C]1C=CCCCC1'),
    E0 = (332.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8347,0.0259641,8.0182e-05,-1.17718e-07,4.53055e-11,40114.4,19.9202], Tmin=(100,'K'), Tmax=(973.885,'K')), NASAPolynomial(coeffs=[13.8091,0.0310733,-1.13075e-05,2.15189e-09,-1.5978e-13,35207.5,-50.7506], Tmin=(973.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]1CCCC2[CH]C12(64)',
    structure = SMILES('[CH]1CCCC2[CH]C12'),
    E0 = (400.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50382,0.0119786,0.000104284,-1.31069e-07,4.70483e-11,48196.3,23.3732], Tmin=(100,'K'), Tmax=(995.793,'K')), NASAPolynomial(coeffs=[9.68385,0.0370116,-1.4577e-05,2.83707e-09,-2.09723e-13,44095.2,-24.648], Tmin=(995.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[4.1.0]heptane-C3-7) + radical(bicyclo[4.1.0]heptane-C6-2)"""),
)

species(
    label = 'C1=CCCCCC=1(65)',
    structure = SMILES('C1=CCCCCC=1'),
    E0 = (181.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94212,0.0223748,9.06321e-05,-1.29301e-07,4.96932e-11,21928.9,18.5292], Tmin=(100,'K'), Tmax=(967.459,'K')), NASAPolynomial(coeffs=[14.0425,0.0299861,-1.05385e-05,1.99868e-09,-1.49545e-13,16890.1,-53.3892], Tmin=(967.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]1[CH]C2CCCC12(67)',
    structure = SMILES('[CH]1[CH]C2CCCC12'),
    E0 = (324.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66775,0.00592624,0.000123203,-1.51785e-07,5.46311e-11,39084.7,22.2493], Tmin=(100,'K'), Tmax=(986.716,'K')), NASAPolynomial(coeffs=[10.1915,0.0357275,-1.37706e-05,2.6963e-09,-2.01715e-13,34664.4,-28.8231], Tmin=(986.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-6) + radical(bicyclo[3.2.0]heptane-C5-6)"""),
)

species(
    label = 'C1=CCCC=CC1(68)',
    structure = SMILES('C1=CCCC=CC1'),
    E0 = (65.9993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96589,0.0219438,9.12463e-05,-1.29388e-07,4.95752e-11,8031.63,19.2218], Tmin=(100,'K'), Tmax=(968.516,'K')), NASAPolynomial(coeffs=[13.8278,0.0303786,-1.07544e-05,2.0418e-09,-1.52573e-13,3040.65,-51.5297], Tmin=(968.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.9993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C1=CC2CCC[C]12(69)',
    structure = SMILES('C1=CC2CCC[C]12'),
    E0 = (240.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.57983,0.00371919,0.000137207,-1.76667e-07,6.6575e-11,29002.2,16.7074], Tmin=(100,'K'), Tmax=(956.053,'K')), NASAPolynomial(coeffs=[13.9899,0.0262012,-8.23771e-06,1.57817e-09,-1.23486e-13,23611.3,-54.6103], Tmin=(956.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_T)"""),
)

species(
    label = '[CH]1CCC2C=CC12(70)',
    structure = SMILES('[CH]1CCC2C=CC12'),
    E0 = (303.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49761,0.00783687,0.000120574,-1.55992e-07,5.81223e-11,36527.3,20.6128], Tmin=(100,'K'), Tmax=(970.534,'K')), NASAPolynomial(coeffs=[13.2011,0.0280244,-1.00068e-05,1.96964e-09,-1.51578e-13,31421.3,-46.3067], Tmin=(970.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S)"""),
)

species(
    label = '[CH]1CC2C=CC2C1(71)',
    structure = SMILES('[CH]1CC2C=CC2C1'),
    E0 = (294.406,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.57722,0.00932315,0.000108826,-1.39235e-07,5.13164e-11,35479.7,19.642], Tmin=(100,'K'), Tmax=(973.889,'K')), NASAPolynomial(coeffs=[10.8741,0.031317,-1.14113e-05,2.19027e-09,-1.63631e-13,31204.6,-33.8179], Tmin=(973.889,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane)"""),
)

species(
    label = '[C]1=CC2CCCC12(72)',
    structure = SMILES('[C]1=CC2CCCC12'),
    E0 = (361.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27176,0.0136871,0.00010703,-1.44667e-07,5.49294e-11,43503.7,19.3405], Tmin=(100,'K'), Tmax=(965.454,'K')), NASAPolynomial(coeffs=[14.017,0.026744,-9.14678e-06,1.77017e-09,-1.3591e-13,38359.4,-51.8074], Tmin=(965.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH]1C[C]2CCCC12(73)',
    structure = SMILES('[CH]1C[C]2CCCC12'),
    E0 = (338.717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,303.444,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58889,0.012735,9.48241e-05,-1.16716e-07,4.09595e-11,40805.1,22.7077], Tmin=(100,'K'), Tmax=(1010.55,'K')), NASAPolynomial(coeffs=[7.78141,0.039824,-1.61028e-05,3.11628e-09,-2.27358e-13,37323,-14.4339], Tmin=(1010.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-tertiary) + radical(bicyclo[3.2.0]heptane-C5-6)"""),
)

species(
    label = '[CH]1CC2CCC[C]12(74)',
    structure = SMILES('[CH]1CC2CCC[C]12'),
    E0 = (338.717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,303.444,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58891,0.0127347,9.48252e-05,-1.16717e-07,4.09601e-11,40805.1,22.7076], Tmin=(100,'K'), Tmax=(1010.54,'K')), NASAPolynomial(coeffs=[7.7813,0.0398242,-1.61029e-05,3.11631e-09,-2.2736e-13,37323,-14.4332], Tmin=(1010.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-6) + radical(bicyclo[3.2.0]heptane-tertiary)"""),
)

species(
    label = '[CH]1CCC2[CH]CC12(75)',
    structure = SMILES('[CH]1CCC2[CH]CC12'),
    E0 = (319.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66775,0.00592624,0.000123203,-1.51785e-07,5.46311e-11,38531.1,22.9424], Tmin=(100,'K'), Tmax=(986.716,'K')), NASAPolynomial(coeffs=[10.1915,0.0357275,-1.37706e-05,2.6963e-09,-2.01715e-13,34110.9,-28.1299], Tmin=(986.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-6) + radical(bicyclo[3.2.0]heptane-C5-2)"""),
)

species(
    label = '[CH]1CCC2C[CH]C12(76)',
    structure = SMILES('[CH]1CCC2C[CH]C12'),
    E0 = (319.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66775,0.00592624,0.000123203,-1.51785e-07,5.46311e-11,38531.1,22.9424], Tmin=(100,'K'), Tmax=(986.716,'K')), NASAPolynomial(coeffs=[10.1915,0.0357275,-1.37706e-05,2.6963e-09,-2.01715e-13,34110.9,-28.1299], Tmin=(986.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-6) + radical(bicyclo[3.2.0]heptane-C5-2)"""),
)

species(
    label = '[CH]1CC2[CH]CC2C1(77)',
    structure = SMILES('[CH]1CC2[CH]CC2C1'),
    E0 = (326.489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66775,0.00592624,0.000123203,-1.51785e-07,5.46311e-11,39336.3,22.9424], Tmin=(100,'K'), Tmax=(986.716,'K')), NASAPolynomial(coeffs=[10.1915,0.0357275,-1.37706e-05,2.6963e-09,-2.01715e-13,34916,-28.1299], Tmin=(986.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-6) + radical(bicyclo[3.2.0]heptane-C5-3)"""),
)

species(
    label = '[CH]=CC1[CH]CCC1(78)',
    structure = SMILES('[CH]=CC1[CH]CCC1'),
    E0 = (403.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88819,0.0268691,7.11195e-05,-1.03575e-07,3.90066e-11,48579.9,25.2513], Tmin=(100,'K'), Tmax=(993.488,'K')), NASAPolynomial(coeffs=[12.4137,0.033516,-1.29354e-05,2.4985e-09,-1.84363e-13,44069.1,-37.6356], Tmin=(993.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CCC1[CH]C=C1(80)',
    structure = SMILES('[CH2]CCC1[CH]C=C1'),
    E0 = (434.123,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,363.369,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,2256.16],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79043,0.0277409,7.43956e-05,-1.1054e-07,4.23935e-11,52311,22.78], Tmin=(100,'K'), Tmax=(981.018,'K')), NASAPolynomial(coeffs=[13.5491,0.0319656,-1.18326e-05,2.26558e-09,-1.67938e-13,47493.6,-46.5185], Tmin=(981.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(434.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutene-allyl) + radical(RCCJ)"""),
)

species(
    label = '[C]1CC2CCCC12(81)',
    structure = SMILES('[C]1CC2CCCC12'),
    E0 = (337.886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2288,0.0110772,0.000127348,-1.69745e-07,6.44766e-11,40727.3,19.2909], Tmin=(100,'K'), Tmax=(960.68,'K')), NASAPolynomial(coeffs=[15.2612,0.027668,-9.18817e-06,1.77724e-09,-1.37976e-13,34953.7,-60.0772], Tmin=(960.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]1C=CCC1(211)',
    structure = SMILES('[CH]1C=CCC1'),
    E0 = (125.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22824,-0.00471828,0.000117614,-1.45678e-07,5.42093e-11,15155.9,11.3913], Tmin=(100,'K'), Tmax=(955.98,'K')), NASAPolynomial(coeffs=[11.0449,0.0181898,-5.59321e-06,1.09704e-09,-8.81308e-14,11120.1,-39.2592], Tmin=(955.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-allyl)"""),
)

species(
    label = 'C=C[C]1C=CCC1(212)',
    structure = SMILES('[CH2]C=C1C=CCC1'),
    E0 = (157.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03268,0.0218362,8.50901e-05,-1.22182e-07,4.72242e-11,19066.8,21.0911], Tmin=(100,'K'), Tmax=(962.414,'K')), NASAPolynomial(coeffs=[13.3923,0.0286435,-9.7147e-06,1.81172e-09,-1.34825e-13,14378.5,-46.2748], Tmin=(962.414,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC1[CH]CC=C1(213)',
    structure = SMILES('C=CC1[CH]CC=C1'),
    E0 = (263.676,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01242,0.0257893,6.44898e-05,-9.50585e-08,3.60183e-11,31800.3,22.8013], Tmin=(100,'K'), Tmax=(988.869,'K')), NASAPolynomial(coeffs=[11.682,0.0314455,-1.20007e-05,2.29264e-09,-1.68133e-13,27699,-34.8084], Tmin=(988.869,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-4)"""),
)

species(
    label = 'C=CC1C=C[CH]C1(214)',
    structure = SMILES('C=CC1[CH]C=CC1'),
    E0 = (201.842,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05624,0.019612,9.29599e-05,-1.31301e-07,5.04547e-11,24366.8,19.0084], Tmin=(100,'K'), Tmax=(966.64,'K')), NASAPolynomial(coeffs=[14.347,0.0271626,-9.39632e-06,1.80358e-09,-1.37029e-13,19261.7,-53.986], Tmin=(966.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-allyl)"""),
)

species(
    label = 'C=CC1[C]=CCC1(215)',
    structure = SMILES('C=CC1[C]=CCC1'),
    E0 = (334.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81267,0.0298364,5.77916e-05,-9.18499e-08,3.58604e-11,40287.7,21.6439], Tmin=(100,'K'), Tmax=(979.858,'K')), NASAPolynomial(coeffs=[13.0712,0.0292443,-1.07528e-05,2.03797e-09,-1.49965e-13,35903.4,-43.5566], Tmin=(979.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C=CC1C=[C]CC1(216)',
    structure = SMILES('C=CC1C=[C]CC1'),
    E0 = (334.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81267,0.0298364,5.77916e-05,-9.18499e-08,3.58604e-11,40287.7,21.6439], Tmin=(100,'K'), Tmax=(979.858,'K')), NASAPolynomial(coeffs=[13.0712,0.0292443,-1.07528e-05,2.03797e-09,-1.49965e-13,35903.4,-43.5566], Tmin=(979.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C=[C]C1C=CCC1(217)',
    structure = SMILES('C=[C]C1C=CCC1'),
    E0 = (314.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81267,0.0298364,5.77916e-05,-9.18499e-08,3.58604e-11,37922.6,21.6439], Tmin=(100,'K'), Tmax=(979.858,'K')), NASAPolynomial(coeffs=[13.0712,0.0292443,-1.07528e-05,2.03797e-09,-1.49965e-13,33538.3,-43.5566], Tmin=(979.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC1C=CCC1(218)',
    structure = SMILES('[CH]=CC1C=CCC1'),
    E0 = (323.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77166,0.0287655,6.58861e-05,-1.03671e-07,4.08748e-11,39039,21.7203], Tmin=(100,'K'), Tmax=(970.273,'K')), NASAPolynomial(coeffs=[14.4703,0.026962,-9.46966e-06,1.79684e-09,-1.34314e-13,34195.4,-51.4212], Tmin=(970.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_P)"""),
)

species(
    label = 'C=C[C]1C[CH]CC1(219)',
    structure = SMILES('[CH2]C=C1C[CH]CC1'),
    E0 = (260.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92025,0.0283869,6.30027e-05,-9.23789e-08,3.45631e-11,31370.7,25.2893], Tmin=(100,'K'), Tmax=(996.931,'K')), NASAPolynomial(coeffs=[10.7646,0.0362073,-1.39241e-05,2.63733e-09,-1.9107e-13,27455.2,-28.146], Tmin=(996.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(cyclopentane)"""),
)

species(
    label = '[CH2]C[C]1C=CCC1(220)',
    structure = SMILES('[CH2]CC1=C[CH]CC1'),
    E0 = (263.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92328,0.0247914,8.01805e-05,-1.14609e-07,4.33563e-11,31776,22.2282], Tmin=(100,'K'), Tmax=(983.374,'K')), NASAPolynomial(coeffs=[12.8077,0.0329823,-1.23415e-05,2.36903e-09,-1.75435e-13,27098.6,-42.9987], Tmin=(983.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-allyl) + radical(RCCJ)"""),
)

species(
    label = 'C=CC1[CH]C[CH]C1(221)',
    structure = SMILES('C=CC1[CH]C[CH]C1'),
    E0 = (341.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23266,0.0235714,6.49863e-05,-8.67847e-08,3.07003e-11,41201,26.1785], Tmin=(100,'K'), Tmax=(1026.87,'K')), NASAPolynomial(coeffs=[8.03455,0.0401068,-1.63354e-05,3.12573e-09,-2.24967e-13,37946,-12.013], Tmin=(1026.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(cyclopentane)"""),
)

species(
    label = '[CH2]CC1[CH]CC=C1(222)',
    structure = SMILES('[CH2]CC1[CH]CC=C1'),
    E0 = (340.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94011,0.0278148,6.32436e-05,-9.20792e-08,3.42244e-11,41019.5,25.2765], Tmin=(100,'K'), Tmax=(1003.67,'K')), NASAPolynomial(coeffs=[10.9005,0.0358327,-1.40926e-05,2.69961e-09,-1.96529e-13,37018.3,-28.9575], Tmin=(1003.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(340.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(cyclopentene-4)"""),
)

species(
    label = '[CH2]CC1[C]=CCC1(223)',
    structure = SMILES('[CH2]CC1[C]=CCC1'),
    E0 = (410.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73975,0.0318732,5.64848e-05,-8.87575e-08,3.40012e-11,49506.9,24.1211], Tmin=(100,'K'), Tmax=(994.041,'K')), NASAPolynomial(coeffs=[12.2715,0.033661,-1.28611e-05,2.44869e-09,-1.78664e-13,45231,-37.6015], Tmin=(994.041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C=[C]C1C[CH]CC1(224)',
    structure = SMILES('C=[C]C1C[CH]CC1'),
    E0 = (385.225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00247,0.0294826,5.11786e-05,-7.50358e-08,2.72869e-11,46416.3,24.921], Tmin=(100,'K'), Tmax=(1022.49,'K')), NASAPolynomial(coeffs=[8.81932,0.0388755,-1.55023e-05,2.93233e-09,-2.09786e-13,43137.2,-17.3358], Tmin=(1022.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(Cds_S)"""),
)

species(
    label = 'C=CC1[CH]CC[CH]1(225)',
    structure = SMILES('C=CC1[CH]CC[CH]1'),
    E0 = (350.586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15769,0.0220419,7.68239e-05,-1.03554e-07,3.74597e-11,42248.3,26.4388], Tmin=(100,'K'), Tmax=(1011.48,'K')), NASAPolynomial(coeffs=[10.2731,0.0369587,-1.50119e-05,2.92382e-09,-2.14442e-13,38201.8,-24.6933], Tmin=(1011.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]C1[CH]CC=C1(226)',
    structure = SMILES('C[CH]C1[CH]CC=C1'),
    E0 = (329.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09343,0.0258338,6.25648e-05,-8.67367e-08,3.11355e-11,39725.1,25.4979], Tmin=(100,'K'), Tmax=(1025.31,'K')), NASAPolynomial(coeffs=[9.27138,0.0386297,-1.58431e-05,3.05472e-09,-2.21198e-13,36108.7,-19.7689], Tmin=(1025.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-4) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]C1[C]=CCC1(227)',
    structure = SMILES('C[CH]C1[C]=CCC1'),
    E0 = (400.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89136,0.0299182,5.56858e-05,-8.32123e-08,3.08036e-11,48212.6,24.3483], Tmin=(100,'K'), Tmax=(1014.56,'K')), NASAPolynomial(coeffs=[10.6176,0.0364971,-1.46329e-05,2.80866e-09,-2.03725e-13,44332.6,-28.2722], Tmin=(1014.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C=[C]C1[CH]CCC1(228)',
    structure = SMILES('C=[C]C1[CH]CCC1'),
    E0 = (393.886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9288,0.027936,6.30855e-05,-9.19092e-08,3.40972e-11,47463.5,25.1769], Tmin=(100,'K'), Tmax=(1006.8,'K')), NASAPolynomial(coeffs=[11.0626,0.0357203,-1.4175e-05,2.7296e-09,-1.99197e-13,43390.7,-30.0433], Tmin=(1006.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(Cs_S)"""),
)

species(
    label = 'C=CC1C[CH][CH]C1(229)',
    structure = SMILES('C=CC1C[CH][CH]C1'),
    E0 = (333.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30517,0.0251212,5.31212e-05,-7.00495e-08,2.39883e-11,40153.8,25.9274], Tmin=(100,'K'), Tmax=(1051.65,'K')), NASAPolynomial(coeffs=[5.86148,0.0431494,-1.76005e-05,3.31424e-09,-2.34405e-13,37660.8,0.295044], Tmin=(1051.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(cyclopentane)"""),
)

species(
    label = 'C=C[C]1[CH]CCC1(230)',
    structure = SMILES('[CH2]C=C1[CH]CCC1'),
    E0 = (215.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85503,0.0228302,9.52286e-05,-1.35835e-07,5.21811e-11,25995.3,21.8514], Tmin=(100,'K'), Tmax=(970.078,'K')), NASAPolynomial(coeffs=[15.0515,0.0297948,-1.04482e-05,2.01252e-09,-1.5253e-13,20546.9,-56.2986], Tmin=(970.078,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC1C=[C]CC1(231)',
    structure = SMILES('[CH2]CC1C=[C]CC1'),
    E0 = (410.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73975,0.0318732,5.64848e-05,-8.87575e-08,3.40012e-11,49506.9,24.1211], Tmin=(100,'K'), Tmax=(994.041,'K')), NASAPolynomial(coeffs=[12.2715,0.033661,-1.28611e-05,2.44869e-09,-1.78664e-13,45231,-37.6015], Tmin=(994.041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[CH]=CC1C[CH]CC1(232)',
    structure = SMILES('[CH]=CC1C[CH]CC1'),
    E0 = (394.479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96395,0.028396,5.9255e-05,-8.67127e-08,3.21795e-11,47532.6,24.9877], Tmin=(100,'K'), Tmax=(1003.82,'K')), NASAPolynomial(coeffs=[10.1315,0.0367347,-1.42982e-05,2.70942e-09,-1.9562e-13,43832.9,-24.7072], Tmin=(1003.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(Cds_P)"""),
)

species(
    label = 'C[CH]C1C=C[CH]C1(233)',
    structure = SMILES('C[CH]C1[CH]C=CC1'),
    E0 = (268.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13513,0.0198802,9.00814e-05,-1.2159e-07,4.49033e-11,32412.1,21.6495], Tmin=(100,'K'), Tmax=(990.933,'K')), NASAPolynomial(coeffs=[11.7319,0.0347079,-1.34475e-05,2.61085e-09,-1.93537e-13,27880.2,-37.8319], Tmin=(990.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-allyl) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]C1C=[C]CC1(234)',
    structure = SMILES('C[CH]C1C=[C]CC1'),
    E0 = (400.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89136,0.0299182,5.56858e-05,-8.32123e-08,3.08036e-11,48212.6,24.3483], Tmin=(100,'K'), Tmax=(1014.56,'K')), NASAPolynomial(coeffs=[10.6176,0.0364971,-1.46329e-05,2.80866e-09,-2.03725e-13,44332.6,-28.2722], Tmin=(1014.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[CH]=CC(C=C)C[CH2](235)',
    structure = SMILES('[CH]=CC(C=C)C[CH2]'),
    E0 = (483.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,1188.85],'cm^-1')),
        HinderedRotor(inertia=(0.587691,'amu*angstrom^2'), symmetry=1, barrier=(13.5122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127981,'amu*angstrom^2'), symmetry=1, barrier=(2.94253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587665,'amu*angstrom^2'), symmetry=1, barrier=(13.5116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587664,'amu*angstrom^2'), symmetry=1, barrier=(13.5115,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.887803,0.0683422,-5.2066e-05,2.16267e-08,-3.81164e-12,58292.7,28.0129], Tmin=(100,'K'), Tmax=(1295.89,'K')), NASAPolynomial(coeffs=[11.0731,0.0369038,-1.56762e-05,2.90629e-09,-2.00181e-13,55652.8,-23.7654], Tmin=(1295.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CCC[CH]C=C(236)',
    structure = SMILES('[CH]=CCCC=C[CH2]'),
    E0 = (426.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,1002.14],'cm^-1')),
        HinderedRotor(inertia=(0.0982296,'amu*angstrom^2'), symmetry=1, barrier=(2.25849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.806,'amu*angstrom^2'), symmetry=1, barrier=(18.5315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.806081,'amu*angstrom^2'), symmetry=1, barrier=(18.5334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.80592,'amu*angstrom^2'), symmetry=1, barrier=(18.5297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.72659,0.0612087,-2.46998e-05,-8.65074e-09,6.80674e-12,51395.5,29.283], Tmin=(100,'K'), Tmax=(1059.96,'K')), NASAPolynomial(coeffs=[14.0378,0.0318345,-1.26495e-05,2.33586e-09,-1.63409e-13,47401.8,-41.2385], Tmin=(1059.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=CC([CH2])C=C(237)',
    structure = SMILES('[CH2]C=CC([CH2])C=C'),
    E0 = (375.736,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,3850.38],'cm^-1')),
        HinderedRotor(inertia=(0.53074,'amu*angstrom^2'), symmetry=1, barrier=(12.2028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205619,'amu*angstrom^2'), symmetry=1, barrier=(24.7637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.89346,'amu*angstrom^2'), symmetry=1, barrier=(89.5184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.90025,'amu*angstrom^2'), symmetry=1, barrier=(89.6744,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.547791,0.0685847,-5.08023e-05,2.03371e-08,-3.35664e-12,45320.9,28.6747], Tmin=(100,'K'), Tmax=(1419.12,'K')), NASAPolynomial(coeffs=[13.6133,0.0317576,-1.18766e-05,2.05085e-09,-1.3524e-13,41612.6,-38.9325], Tmin=(1419.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC1[C]CCC1(238)',
    structure = SMILES('C=CC1[C]CCC1'),
    E0 = (336.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64691,0.0310246,6.75956e-05,-1.06271e-07,4.16868e-11,40553.3,22.2452], Tmin=(100,'K'), Tmax=(973.123,'K')), NASAPolynomial(coeffs=[14.3651,0.0303726,-1.0977e-05,2.07423e-09,-1.53312e-13,35633.6,-51.3259], Tmin=(973.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C=CC1C[C]CC1(239)',
    structure = SMILES('C=CC1C[C]CC1'),
    E0 = (325.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62808,0.031314,6.7434e-05,-1.06589e-07,4.19243e-11,39303.8,22.3479], Tmin=(100,'K'), Tmax=(972.063,'K')), NASAPolynomial(coeffs=[14.5505,0.0300792,-1.081e-05,2.04105e-09,-1.51015e-13,34337.6,-52.2513], Tmin=(972.063,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C[C]C1C=CCC1(240)',
    structure = SMILES('C[C]C1C=CCC1'),
    E0 = (322.946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64362,0.0321932,6.20296e-05,-9.90666e-08,3.88046e-11,38943.8,21.4932], Tmin=(100,'K'), Tmax=(977.057,'K')), NASAPolynomial(coeffs=[13.7372,0.031404,-1.15567e-05,2.17894e-09,-1.59806e-13,34255,-48.4716], Tmin=(977.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]CC1C=CCC1(241)',
    structure = SMILES('[CH]CC1C=CCC1'),
    E0 = (397.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57457,0.0329999,6.27713e-05,-1.01665e-07,4.01308e-11,47857.2,21.9136], Tmin=(100,'K'), Tmax=(974.173,'K')), NASAPolynomial(coeffs=[14.4992,0.0304415,-1.10639e-05,2.08831e-09,-1.53958e-13,42942.2,-52.4035], Tmin=(974.173,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(397.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]1CC2[CH]C1CC2(242)',
    structure = SMILES('[CH]1CC2[CH]C1CC2'),
    E0 = (343.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.82119,-0.000636553,0.000145277,-1.77133e-07,6.42307e-11,41415.4,20.9102], Tmin=(100,'K'), Tmax=(975.933,'K')), NASAPolynomial(coeffs=[11.2195,0.033423,-1.23273e-05,2.42807e-09,-1.84964e-13,36514.9,-36.1105], Tmin=(975.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(7-norbornyl) + radical(2-norbornyl)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61263,-0.00100893,2.49899e-06,-1.43376e-09,2.58637e-13,-1051.1,2.6527], Tmin=(100,'K'), Tmax=(1817.03,'K')), NASAPolynomial(coeffs=[2.97589,0.00164143,-7.1973e-07,1.25379e-10,-7.91538e-15,-1025.84,5.53764], Tmin=(1817.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.69489,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: BurkeH2O2"""),
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
    E0 = (632.423,'kJ/mol'),
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
    E0 = (514.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (620.494,'kJ/mol'),
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
    E0 = (582.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (587.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (641.677,'kJ/mol'),
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
    E0 = (417.958,'kJ/mol'),
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
    E0 = (530.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (496.033,'kJ/mol'),
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
    E0 = (580.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (503.953,'kJ/mol'),
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
    E0 = (450.879,'kJ/mol'),
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
    E0 = (430.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (347.365,'kJ/mol'),
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
    E0 = (471.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (512.425,'kJ/mol'),
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
    E0 = (525.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (305.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (347.833,'kJ/mol'),
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
    E0 = (449.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (604.666,'kJ/mol'),
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
    E0 = (394.703,'kJ/mol'),
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
    E0 = (671.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (901.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (609.662,'kJ/mol'),
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
    E0 = (569.954,'kJ/mol'),
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
    E0 = (474.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (382.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (517.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (483.093,'kJ/mol'),
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
    E0 = (400.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS93',
    E0 = (381.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS94',
    E0 = (494.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS95',
    E0 = (411.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS96',
    E0 = (530.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS97',
    E0 = (471.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS98',
    E0 = (296.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS99',
    E0 = (438.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS100',
    E0 = (241.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS101',
    E0 = (546.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS102',
    E0 = (467.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS103',
    E0 = (533.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS104',
    E0 = (381.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS105',
    E0 = (324.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS106',
    E0 = (258.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS107',
    E0 = (452.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS108',
    E0 = (514.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS109',
    E0 = (506.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS110',
    E0 = (572.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS111',
    E0 = (361.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS112',
    E0 = (402.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS113',
    E0 = (383.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS114',
    E0 = (328.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS115',
    E0 = (334.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS116',
    E0 = (411.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS117',
    E0 = (475.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS118',
    E0 = (441.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS119',
    E0 = (429.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS120',
    E0 = (415.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS121',
    E0 = (369.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS122',
    E0 = (475.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS123',
    E0 = (416.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS124',
    E0 = (545.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS125',
    E0 = (545.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS126',
    E0 = (526.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS127',
    E0 = (535.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS128',
    E0 = (283.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS129',
    E0 = (285.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS130',
    E0 = (365.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS131',
    E0 = (403.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS132',
    E0 = (499.784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS133',
    E0 = (484.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS134',
    E0 = (413.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS135',
    E0 = (337.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS136',
    E0 = (439.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS137',
    E0 = (433.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS138',
    E0 = (358.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS139',
    E0 = (318.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS140',
    E0 = (254.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS141',
    E0 = (419.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS142',
    E0 = (428.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS143',
    E0 = (308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS144',
    E0 = (408.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS145',
    E0 = (428.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS146',
    E0 = (490.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS147',
    E0 = (433.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS148',
    E0 = (382.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS149',
    E0 = (427.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS150',
    E0 = (417.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS151',
    E0 = (414.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS152',
    E0 = (453.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS153',
    E0 = (438.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS154',
    E0 = (319.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction90',
    reactants = ['H(10)', 'C=CC=[C]C1CC1(167)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction1',
    reactants = ['CH2CHCHCH(163)', 'cC3H5(164)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.79546e+07,'m^3/(mol*s)'), n=-0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/NonDeC] for rate rule [Cd_pri_rad;C_rad/H/NonDeC]
Euclidian distance = 2.0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction88',
    reactants = ['H(10)', 'C=CC=C[C]1CC1(165)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction89',
    reactants = ['H(10)', 'C=CC=CC1[CH]C1(166)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction91',
    reactants = ['C2H3(32)', '[CH]=CC1CC1(168)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 89 used for Cd_pri_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction92',
    reactants = ['H(10)', 'C=C[C]=CC1CC1(169)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction93',
    reactants = ['H(10)', 'C=[C]C=CC1CC1(170)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction94',
    reactants = ['H(10)', '[CH]=CC=CC1CC1(171)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction95',
    reactants = ['CH2(S)(28)', 'C=CC=CC=C(158)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.13244e+08,'m^3/(mol*s)'), n=-0.441667, Ea=(7.08269,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;mb_db]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction96',
    reactants = ['C=C[CH]C[C]1CC1(172)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.64935e+10,'s^-1'), n=0.2551, Ea=(25.8154,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction97',
    reactants = ['[CH2]C[C]=CC1CC1(173)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction98',
    reactants = ['C=[C]C[CH]C1CC1(174)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction99',
    reactants = ['C=C[CH]CC1[CH]C1(175)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction100',
    reactants = ['[CH2]CC=[C]C1CC1(176)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction101',
    reactants = ['[CH]=CC[CH]C1CC1(177)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction102',
    reactants = ['C=CC[CH]C1[CH]C1(178)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction103',
    reactants = ['C[CH]C=[C]C1CC1(179)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction104',
    reactants = ['[CH]=C[CH]CC1CC1(180)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction105',
    reactants = ['[CH2]CC=C[C]1CC1(181)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction106',
    reactants = ['C[CH]C=C[C]1CC1(182)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction107',
    reactants = ['[CH2]CC=CC1[CH]C1(140)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction108',
    reactants = ['C[CH]C=CC1[CH]C1(183)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad_NDe;XH_Rrad] for rate rule [R6radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction109',
    reactants = ['CC=CC=C1CC1(184)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.02873e+09,'s^-1'), n=1.23767, Ea=(163.714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH3_1;Cd(C)C_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction110',
    reactants = ['[CH2]C([CH2])C=CC=C(185)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction86',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction111',
    reactants = ['C=CC[C]C1CC1(186)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;CsJ2C;CH2(C=C)] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C=C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction112',
    reactants = ['C=C[C]CC1CC1(187)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.22428e+12,'s^-1'), n=0.271335, Ea=(58.4032,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2(C=C);CH2(C)] + [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction113',
    reactants = ['C[C]C=CC1CC1(188)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.08343e+13,'s^-1'), n=0.129171, Ea=(99.5787,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction114',
    reactants = ['[CH]CC=CC1CC1(189)'],
    products = ['C=CC=CC1CC1(162)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;singletcarbene;CH2(C=C)] for rate rule [CsJ2-C;CsJ2H;CH2(C=C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction115',
    reactants = ['C=CC=CC1CC1(162)'],
    products = ['C1=CC(C1)C1CC1(190)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH2_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C=C[CH]C1CC1(122)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.34e+11,'s^-1'), n=0.21, Ea=(25.3061,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 124 used for R4_S_D;doublebond_intra_HCd_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_HCd_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 23.4 to 25.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C[CH]C1C=CC1(123)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra_HNd_pri;radadd_intra_cs] for rate rule [R5_SD_D;doublebond_intra_HNd_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2][CH]C1C=CCC1(124)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R6_SSM_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['H(10)', '[CH2]C=CC=CC=C(125)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.24e+08,'cm^3/(mol*s)'), n=1.64, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2580 used for Cds-CdH_Cds-HH;HJ
Exact match found for rate rule [Cds-CdH_Cds-HH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction48',
    reactants = ['H(10)', '[CH2]CC=CC=C=C(126)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C2H4(33)', '[CH]=CC=C[CH2](127)'],
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
    products = ['[CH2]C=CC=C[CH]C(128)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(9.08185e+06,'s^-1'), n=1.84946, Ea=(92.0377,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_single;Cs_H_out_H/(Cd-Cd-Cd)] + [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]CC=CC=[C]C(129)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C=CC=[C]CC(130)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]CC=C[C]=CC(131)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 288 used for R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C=C[C]=CCC(132)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.71035e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2]CC=[C]C=CC(133)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.33e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 284 used for R4H_SDS;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C=[C]C=CCC(134)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(484628,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;C_rad_out_2H;XH_out] for rate rule [R5H_SSMS;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2]C[C]=CC=CC(135)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSMS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2][C]=CC=CCC(136)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.27582e+06,'s^-1'), n=1.20683, Ea=(76.1033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_RSMSR;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C2H4(T)(137)', '[CH]=CC=C[CH2](127)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(7.76856e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_pri_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH]=C[CH2](138)', '[CH]=CC[CH2](139)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using an average for rate rule [Cd_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]CC=CC1[CH]C1(140)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 142 used for R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C=CC1[CH]CC1(141)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(6.76e+10,'s^-1'), n=0.19, Ea=(139.746,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 232 used for R4_Cs_HH_D;doublebond_intra_pri_HCd;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_HH_D;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]CC1[CH]C=CC1(142)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(9.95802e+09,'s^-1'), n=0.161323, Ea=(148.068,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 95 used for R5_SD_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H
Exact match found for rate rule [R5_SD_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]C1C=C[CH]CC1(66)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.15866e+09,'s^-1'), n=0.354662, Ea=(72.1768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSM_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R6_SSM_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['C=CC=CC=CC(143)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['C=C=CC=CCC(144)'],
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
    label = 'reaction67',
    reactants = ['CH2(T)(22)', '[CH2]C=CC=C[CH2](145)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(8.13728e+06,'m^3/(mol*s)'), n=0.369106, Ea=(32.7067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction68',
    reactants = ['CH2(T)(22)', '[CH]=CC=CC[CH2](146)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(4.06864e+06,'m^3/(mol*s)'), n=0.369106, Ea=(32.7067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction69',
    reactants = ['H(10)', '[CH2]CC=C=CC=C(147)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction70',
    reactants = ['[CH2]CC=[C]CC=C(148)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH2]CC=CC[C]=C(149)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH2]C[C]=CCC=C(150)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction73',
    reactants = ['[CH]=CCC=CC[CH2](151)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2][CH]C=CCC=C(152)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(6.92799e+06,'s^-1'), n=1.7075, Ea=(114.432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction75',
    reactants = ['[CH]=C[CH]C=CCC(153)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R7HJ_2;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction76',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]CC1[CH]C1C=C(154)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(5.21061e+14,'s^-1'), n=-0.479561, Ea=(247.427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;doublebond_intra_pri_HNd_Cs;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
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
    products = ['C=CC=CCC=C(155)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction78',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['C=CC=C=CCC(156)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction79',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R5_SSDS;C_rad_out_2H;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction80',
    reactants = ['[CH2]C=CC=CC[CH2](56)'],
    products = ['[CH2]CC1C=CC1[CH2](79)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra_csHNd] for rate rule [R5_SD_D;doublebond_intra_2H_pri;radadd_intra_csHNd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction81',
    reactants = ['CH2(T)(22)', 'C=CC=CC=C(158)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(0.18583,'m^3/(mol*s)'), n=2.36967, Ea=(33.8986,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-CdH;YJ] for rate rule [Cds-HH_Cds-CdH;CH2_triplet]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction82',
    reactants = ['[CH2]CC[C]=CC=C(159)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction83',
    reactants = ['[CH2]CCC=[C]C=C(160)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(6.1583e+09,'s^-1'), n=0.92705, Ea=(170.178,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] + [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction84',
    reactants = ['[CH2]CCC=C[C]=C(161)'],
    products = ['[CH2]C=CC=CC[CH2](56)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(2.22e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_Cd;Cs_H_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction85',
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
    kinetics = Arrhenius(A=(3.17358e+07,'m^3/(mol*s)'), n=0.285248, Ea=(3.35904,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
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
    label = 'reaction84',
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
    kinetics = Arrhenius(A=(2.22428e+12,'s^-1'), n=0.271335, Ea=(58.4032,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2(C=C);CH2(C)] + [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH2(C)]
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
    reactants = ['[CH]1C=C[CH]CCC1(51)'],
    products = ['[CH]1CCCC2[CH]C12(64)'],
    transitionState = 'TS92',
    kinetics = Arrhenius(A=(7.27808e+07,'s^-1'), n=0.867672, Ea=(166.949,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(10)', '[CH]1C=CC=CCC1(43)'],
    products = ['[CH]1C=C[CH]CCC1(51)'],
    transitionState = 'TS93',
    kinetics = Arrhenius(A=(2.7e+08,'cm^3/(mol*s)'), n=1.64, Ea=(1.58992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2557 used for Cds-CsH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CdH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(10)', '[C]1=CCCCC=C1(45)'],
    products = ['[CH]1C=C[CH]CCC1(51)'],
    transitionState = 'TS94',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction95',
    reactants = ['[CH]1C=CCC[CH]C1(46)'],
    products = ['[CH]1C=C[CH]CCC1(51)'],
    transitionState = 'TS95',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction96',
    reactants = ['[C]1=C[CH]CCCC1(50)'],
    products = ['[CH]1C=C[CH]CCC1(51)'],
    transitionState = 'TS96',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[C]1[CH]CCCCC=1(52)'],
    products = ['[CH]1C=C[CH]CCC1(51)'],
    transitionState = 'TS97',
    kinetics = Arrhenius(A=(6.64e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cd;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]1C=C[CH]CCC1(51)'],
    products = ['C1=CCCCCC=1(65)'],
    transitionState = 'TS98',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C1C=C[CH]CC1(66)'],
    products = ['[CH]1C=C[CH]CCC1(51)'],
    transitionState = 'TS99',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]1C=C[CH]CCC1(51)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS100',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Cpri_rad_out_single] + [R4;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SDS;C_rad_out_H/NonDeC;Cpri_rad_out_H/NonDeC]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[C]1=CCCC[CH]C1(47)'],
    products = ['[CH]1C=C[CH]CCC1(51)'],
    transitionState = 'TS101',
    kinetics = Arrhenius(A=(9.93038e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]1[CH]CCC=CC1(53)'],
    products = ['[CH]1C=C[CH]CCC1(51)'],
    transitionState = 'TS102',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] for rate rule [R3HJ;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[C]1=CC[CH]CCC1(48)'],
    products = ['[CH]1C=C[CH]CCC1(51)'],
    transitionState = 'TS103',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]1C[CH]CC=CC1(49)'],
    products = ['[CH]1C=C[CH]CCC1(51)'],
    transitionState = 'TS104',
    kinetics = Arrhenius(A=(1.004,'s^-1'), n=3.86, Ea=(41.6308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] for rate rule [R4HJ_2;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]1C=C[CH]CCC1(51)'],
    products = ['[CH]1[CH]C2CCCC12(67)'],
    transitionState = 'TS105',
    kinetics = Arrhenius(A=(1.53577e+10,'s^-1'), n=0.566412, Ea=(91.2011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra_pri;radadd_intra_csHCs] for rate rule [Rn0c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 90.4 to 91.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]1C=C[CH]CCC1(51)'],
    products = ['C1=CCCC=CC1(68)'],
    transitionState = 'TS106',
    kinetics = Arrhenius(A=(8.50442e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(10)', 'C1=CC2CCC[C]12(69)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS107',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(10)', '[CH]1CCC2C=CC12(70)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS108',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(10)', '[CH]1CC2C=CC2C1(71)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS109',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(10)', '[C]1=CC2CCCC12(72)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS110',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]1C[C]2CCCC12(73)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS111',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]1CC2CCC[C]12(74)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS112',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]1CCC2[CH]CC12(75)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS113',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]1CCC2C[CH]C12(76)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS114',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]1CC2[CH]CC2C1(77)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS115',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=CC1[CH]CCC1(78)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS116',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]CC1C=CC1[CH2](79)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS117',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for R5_SSSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R5_SSSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]CCC1[CH]C=C1(80)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS118',
    kinetics = Arrhenius(A=(2.49159e+11,'s^-1'), n=0.1555, Ea=(7.322,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R5_SSSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R5_SSSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[C]1CC2CCCC12(81)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS119',
    kinetics = Arrhenius(A=(1.83659e+12,'s^-1'), n=0.345437, Ea=(91.3338,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction149',
    reactants = ['[CH]1C=CCC1(211)', 'C2H3(32)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS120',
    kinetics = Arrhenius(A=(3.17358e+07,'m^3/(mol*s)'), n=0.285248, Ea=(3.35904,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;C_rad/H/CdCs] for rate rule [Cd_pri_rad;C_rad/H/CdCs]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction150',
    reactants = ['H(10)', 'C=C[C]1C=CCC1(212)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS121',
    kinetics = Arrhenius(A=(3.62e+13,'cm^3/(mol*s)'), n=0.228, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 13 used for C_rad/TwoDeCs;H_rad
Exact match found for rate rule [C_rad/TwoDeCs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction151',
    reactants = ['H(10)', 'C=CC1[CH]CC=C1(213)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS122',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction152',
    reactants = ['H(10)', 'C=CC1C=C[CH]C1(214)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS123',
    kinetics = Arrhenius(A=(1.58679e+07,'m^3/(mol*s)'), n=0.285248, Ea=(3.35904,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction153',
    reactants = ['H(10)', 'C=CC1[C]=CCC1(215)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS124',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction154',
    reactants = ['H(10)', 'C=CC1C=[C]CC1(216)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS125',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction155',
    reactants = ['H(10)', 'C=[C]C1C=CCC1(217)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS126',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction156',
    reactants = ['H(10)', '[CH]=CC1C=CCC1(218)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS127',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction157',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS128',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction158',
    reactants = ['[CH2]C[C]1C=CCC1(220)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS129',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction130',
    reactants = ['C=CC1[CH]C[CH]C1(221)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS130',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction160',
    reactants = ['[CH2]CC1[CH]CC=C1(222)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS131',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction161',
    reactants = ['[CH2]CC1[C]=CCC1(223)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS132',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction162',
    reactants = ['C=[C]C1C[CH]CC1(224)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS133',
    kinetics = Arrhenius(A=(6.94203e+09,'s^-1'), n=0.37, Ea=(99.6901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad_NDe] + [R3radExo;Y_rad;XH_Rrad_NDe] for rate rule [R3radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction163',
    reactants = ['C=CC1[CH]CC[CH]1(225)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS134',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction164',
    reactants = ['C[CH]C1[CH]CC=C1(226)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS135',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction165',
    reactants = ['C[CH]C1[C]=CCC1(227)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS136',
    kinetics = Arrhenius(A=(1.34494e+09,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction166',
    reactants = ['C=[C]C1[CH]CCC1(228)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS137',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction167',
    reactants = ['C=CC1C[CH][CH]C1(229)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS138',
    kinetics = Arrhenius(A=(2.05689e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction168',
    reactants = ['[CH2]CC1[CH]C=CC1(142)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS139',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction169',
    reactants = ['C=C[C]1[CH]CCC1(230)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS140',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction170',
    reactants = ['[CH2]CC1C=[C]CC1(231)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS141',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction171',
    reactants = ['[CH]=CC1C[CH]CC1(232)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS142',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction172',
    reactants = ['C[CH]C1C=C[CH]C1(233)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS143',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction173',
    reactants = ['C[CH]C1C=[C]CC1(234)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS144',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction174',
    reactants = ['[CH]=CC1[CH]CCC1(78)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS145',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction175',
    reactants = ['[CH]=CC(C=C)C[CH2](235)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS146',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Ypri_rad_out] for rate rule [R5_SSSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction176',
    reactants = ['[CH]=CCC[CH]C=C(236)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS147',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R5_SSSD;C_rad_out_H/OneDe;CdsinglepriH_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction177',
    reactants = ['[CH2]C=CC([CH2])C=C(237)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS148',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R5_SSDS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction178',
    reactants = ['C=CC1[C]CCC1(238)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS149',
    kinetics = Arrhenius(A=(1.83659e+12,'s^-1'), n=0.345437, Ea=(91.3338,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction179',
    reactants = ['C=CC1C[C]CC1(239)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS150',
    kinetics = Arrhenius(A=(1.83659e+12,'s^-1'), n=0.345437, Ea=(91.3338,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction180',
    reactants = ['C[C]C1C=CCC1(240)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS151',
    kinetics = Arrhenius(A=(2.75488e+12,'s^-1'), n=0.345437, Ea=(91.3338,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction181',
    reactants = ['[CH]CC1C=CCC1(241)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS152',
    kinetics = Arrhenius(A=(1.57948e+12,'s^-1'), n=0.325402, Ea=(56.3419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;singletcarbene;CH2(C)] + [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH2(C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction182',
    reactants = ['[CH]1CC2[CH]C1CC2(242)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS153',
    kinetics = Arrhenius(A=(1.28917e+12,'s^-1'), n=0.239384, Ea=(94.3194,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction183',
    reactants = ['[CH]1CCC2[CH]CC12(75)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS154',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R7JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

network(
    label = '18',
    isomers = [
        'C=CC=CC1CC1(162)',
        '[CH2]C=CC=CC[CH2](56)',
        'C7H10(1)',
        '[CH]1C=C[CH]CCC1(51)',
        'C1=CC2CCCC12(59)',
        'C=CC1C=CCC1(157)',
    ],
    reactants = [
        ('H(10)', 'C=CC=[C]C1CC1(167)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '18',
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

