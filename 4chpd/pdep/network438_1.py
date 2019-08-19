species(
    label = 'CC=C[CH]CC[O](2185)',
    structure = SMILES('CC=C[CH]CC[O]'),
    E0 = (130.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,302.918,302.935,1948.13,1948.14],'cm^-1')),
        HinderedRotor(inertia=(0.132823,'amu*angstrom^2'), symmetry=1, barrier=(8.64925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132814,'amu*angstrom^2'), symmetry=1, barrier=(8.64918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00183695,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594723,'amu*angstrom^2'), symmetry=1, barrier=(38.7274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84586,0.0550909,-2.77083e-05,6.11333e-09,-5.16057e-13,15784,24.2823], Tmin=(100,'K'), Tmax=(2710.28,'K')), NASAPolynomial(coeffs=[24.5141,0.0216361,-9.19304e-06,1.55904e-09,-9.59681e-14,3496.48,-107.68], Tmin=(2710.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(CCOJ)"""),
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
    label = 'C[CH]C1[CH]CCO1(3220)',
    structure = SMILES('C[CH]C1[CH]CCO1'),
    E0 = (141.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89095,0.0256874,8.30734e-05,-1.29417e-07,5.37453e-11,17058.3,24.2251], Tmin=(100,'K'), Tmax=(908.108,'K')), NASAPolynomial(coeffs=[14.2101,0.0252729,-5.18787e-06,6.75228e-10,-4.51772e-14,12600.5,-46.2452], Tmin=(908.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(141.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(CCJCO)"""),
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
    label = 'CC=CC=CC[O](3225)',
    structure = SMILES('CC=CC=CC[O]'),
    E0 = (94.4394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,267.439,267.617,267.847],'cm^-1')),
        HinderedRotor(inertia=(0.15287,'amu*angstrom^2'), symmetry=1, barrier=(7.76235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153387,'amu*angstrom^2'), symmetry=1, barrier=(7.76731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152973,'amu*angstrom^2'), symmetry=1, barrier=(7.76823,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53797,0.0560786,-3.37482e-05,1.02632e-08,-1.32006e-12,11444.5,25.0449], Tmin=(100,'K'), Tmax=(1652,'K')), NASAPolynomial(coeffs=[9.96629,0.0356711,-1.52186e-05,2.78559e-09,-1.88465e-13,8659.77,-19.8478], Tmin=(1652,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.4394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(CCOJ)"""),
)

species(
    label = 'CC=C[CH]CC=O(3226)',
    structure = SMILES('C[CH]C=CCC=O'),
    E0 = (-10.6972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,420.236,420.24],'cm^-1')),
        HinderedRotor(inertia=(0.000954561,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000954568,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181809,'amu*angstrom^2'), symmetry=1, barrier=(22.7843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181807,'amu*angstrom^2'), symmetry=1, barrier=(22.7843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43329,0.0437563,1.1295e-05,-3.70623e-08,1.41961e-11,-1183.58,24.5025], Tmin=(100,'K'), Tmax=(1111.99,'K')), NASAPolynomial(coeffs=[12.3431,0.0341521,-1.5732e-05,3.11137e-09,-2.24799e-13,-5442.44,-37.529], Tmin=(1111.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.6972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_S)"""),
)

species(
    label = 'CC=C=CCC[O](3227)',
    structure = SMILES('CC=C=CCC[O]'),
    E0 = (153.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.235383,'amu*angstrom^2'), symmetry=1, barrier=(5.41192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235301,'amu*angstrom^2'), symmetry=1, barrier=(5.41004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235278,'amu*angstrom^2'), symmetry=1, barrier=(5.40951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32897,0.0470352,3.86733e-05,-1.98074e-07,1.89277e-10,18560.2,22.6063], Tmin=(100,'K'), Tmax=(400.52,'K')), NASAPolynomial(coeffs=[3.96569,0.0469558,-2.19502e-05,4.23691e-09,-2.9809e-13,18298.6,14.5789], Tmin=(400.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(CCOJ)"""),
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
    label = 'CC=CC[CH]C[O](3228)',
    structure = SMILES('CC=CC[CH]C[O]'),
    E0 = (189.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,293.766,293.935,1732.7,4000],'cm^-1')),
        HinderedRotor(inertia=(0.121078,'amu*angstrom^2'), symmetry=1, barrier=(7.40925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00195614,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121114,'amu*angstrom^2'), symmetry=1, barrier=(7.41122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.520729,'amu*angstrom^2'), symmetry=1, barrier=(31.9228,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86022,0.0550888,-2.87413e-05,6.8352e-09,-6.44973e-13,22853.6,27.8179], Tmin=(100,'K'), Tmax=(2229.91,'K')), NASAPolynomial(coeffs=[13.6141,0.0340048,-1.45587e-05,2.59512e-09,-1.69611e-13,17611.5,-38.3143], Tmin=(2229.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCJCO) + radical(CCOJ)"""),
)

species(
    label = 'CC=C[CH]C[CH]O(3229)',
    structure = SMILES('CC=C[CH]C[CH]O'),
    E0 = (85.2554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.903358,0.0698374,-5.68938e-05,2.61067e-08,-5.14142e-12,10363.8,26.3919], Tmin=(100,'K'), Tmax=(1166.86,'K')), NASAPolynomial(coeffs=[10.0336,0.0385391,-1.66599e-05,3.11982e-09,-2.16495e-13,8233.08,-19.0652], Tmin=(1166.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.2554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(CCsJOH)"""),
)

species(
    label = 'CC=[C]CCC[O](3230)',
    structure = SMILES('CC=[C]CCC[O]'),
    E0 = (227.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3010,987.5,1337.5,450,1655,295.037,295.114,295.124,2568.46],'cm^-1')),
        HinderedRotor(inertia=(0.105522,'amu*angstrom^2'), symmetry=1, barrier=(6.516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105402,'amu*angstrom^2'), symmetry=1, barrier=(6.51458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105372,'amu*angstrom^2'), symmetry=1, barrier=(6.51487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105417,'amu*angstrom^2'), symmetry=1, barrier=(6.51614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06918,0.0710112,-8.22701e-05,7.19695e-08,-2.76333e-11,27448.1,28.6861], Tmin=(100,'K'), Tmax=(772.301,'K')), NASAPolynomial(coeffs=[2.96195,0.0510768,-2.38756e-05,4.57651e-09,-3.19319e-13,27457.9,21.9997], Tmin=(772.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'CC=CCC[CH][O](3231)',
    structure = SMILES('CC=CCC[CH][O]'),
    E0 = (169.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,330.858,330.868,330.868,2083.88],'cm^-1')),
        HinderedRotor(inertia=(0.100302,'amu*angstrom^2'), symmetry=1, barrier=(7.79203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100302,'amu*angstrom^2'), symmetry=1, barrier=(7.79195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100305,'amu*angstrom^2'), symmetry=1, barrier=(7.79206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100303,'amu*angstrom^2'), symmetry=1, barrier=(7.79201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.888619,0.0744432,-8.59409e-05,7.09628e-08,-2.58946e-11,20534.3,28.4126], Tmin=(100,'K'), Tmax=(757.317,'K')), NASAPolynomial(coeffs=[4.57537,0.0487228,-2.26225e-05,4.33013e-09,-3.02228e-13,20155,12.8339], Tmin=(757.317,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'CC=C[CH][CH]CO(3232)',
    structure = SMILES('C[CH]C=C[CH]CO'),
    E0 = (23.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.460684,0.0678298,-4.13819e-05,8.31869e-09,8.25374e-13,2903.12,26.3034], Tmin=(100,'K'), Tmax=(1152.6,'K')), NASAPolynomial(coeffs=[14.808,0.0320004,-1.29232e-05,2.36743e-09,-1.6326e-13,-1331.63,-48.9748], Tmin=(1152.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(C=CCJCO)"""),
)

species(
    label = 'C[C]=CCCC[O](3233)',
    structure = SMILES('C[C]=CCCC[O]'),
    E0 = (227.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06918,0.0710112,-8.22701e-05,7.19695e-08,-2.76333e-11,27448.1,28.6861], Tmin=(100,'K'), Tmax=(772.301,'K')), NASAPolynomial(coeffs=[2.96195,0.0510768,-2.38756e-05,4.57651e-09,-3.19319e-13,27457.9,21.9997], Tmin=(772.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CCCC[O](3234)',
    structure = SMILES('[CH2]C=CCCC[O]'),
    E0 = (141.049,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,335.541,335.688,2312.69,2313.14],'cm^-1')),
        HinderedRotor(inertia=(0.366019,'amu*angstrom^2'), symmetry=1, barrier=(29.2418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103051,'amu*angstrom^2'), symmetry=1, barrier=(8.24704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102774,'amu*angstrom^2'), symmetry=1, barrier=(8.25004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364814,'amu*angstrom^2'), symmetry=1, barrier=(29.2498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45667,0.0599447,-3.51885e-05,1.01817e-08,-1.23758e-12,17051.2,27.0784], Tmin=(100,'K'), Tmax=(1719.1,'K')), NASAPolynomial(coeffs=[10.4164,0.0390972,-1.69982e-05,3.12752e-09,-2.11734e-13,13970.6,-21.0018], Tmin=(1719.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(141.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCOJ)"""),
)

species(
    label = 'CC=[C][CH]CCO(3235)',
    structure = SMILES('CC=[C][CH]CCO'),
    E0 = (142.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18908,0.0651257,-4.85216e-05,2.06067e-08,-3.88164e-12,17273.1,26.2904], Tmin=(100,'K'), Tmax=(1179.83,'K')), NASAPolynomial(coeffs=[8.25154,0.0411818,-1.80803e-05,3.40584e-09,-2.3688e-13,15606.6,-8.94989], Tmin=(1179.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'C[C]=C[CH]CCO(3236)',
    structure = SMILES('C[C]=C[CH]CCO'),
    E0 = (142.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18908,0.0651257,-4.85216e-05,2.06067e-08,-3.88164e-12,17273.1,26.2904], Tmin=(100,'K'), Tmax=(1179.83,'K')), NASAPolynomial(coeffs=[8.25154,0.0411818,-1.80803e-05,3.40584e-09,-2.3688e-13,15606.6,-8.94989], Tmin=(1179.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C[CH]CCO(3237)',
    structure = SMILES('[CH2]C=C[CH]CCO'),
    E0 = (56.4567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615088,0.0650776,-3.92876e-05,8.44814e-09,4.41948e-13,6919.82,28.163], Tmin=(100,'K'), Tmax=(1171.88,'K')), NASAPolynomial(coeffs=[13.872,0.0321076,-1.28043e-05,2.32392e-09,-1.59111e-13,2969.52,-41.494], Tmin=(1171.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.4567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'CC1[CH]C1CC[O](3238)',
    structure = SMILES('CC1[CH]C1CC[O]'),
    E0 = (244.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32621,0.0508077,-1.2977e-05,-8.48881e-09,4.07928e-12,29502,27.9188], Tmin=(100,'K'), Tmax=(1241.4,'K')), NASAPolynomial(coeffs=[9.9593,0.0381347,-1.5963e-05,2.94191e-09,-2.0176e-13,26191.6,-20.2977], Tmin=(1241.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CCOJ) + radical(cyclopropane)"""),
)

species(
    label = 'CC1[CH][CH]CCO1(3157)',
    structure = SMILES('CC1[CH][CH]CCO1'),
    E0 = (113.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03294,0.0286987,6.12424e-05,-9.5356e-08,3.86768e-11,13768.3,22.8346], Tmin=(100,'K'), Tmax=(914.292,'K')), NASAPolynomial(coeffs=[9.21387,0.0347795,-1.0252e-05,1.63181e-09,-1.08638e-13,10887.9,-19.7366], Tmin=(914.292,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxane) + radical(RCCJCC) + radical(CCJCO)"""),
)

species(
    label = 'CC=CC=CCO(3239)',
    structure = SMILES('CC=CC=CCO'),
    E0 = (-131.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.663494,0.0646464,-4.1149e-05,1.03543e-08,2.38981e-15,-15660.2,27.3419], Tmin=(100,'K'), Tmax=(1137.77,'K')), NASAPolynomial(coeffs=[13.5032,0.0310097,-1.19691e-05,2.14282e-09,-1.45904e-13,-19326.5,-39.5317], Tmin=(1137.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H)"""),
)

species(
    label = 'CC=CCCC=O(3240)',
    structure = SMILES('CC=CCCC=O'),
    E0 = (-159.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4883,0.0578452,-3.2577e-05,8.7911e-09,-9.69802e-13,-19051.8,25.1553], Tmin=(100,'K'), Tmax=(1926.54,'K')), NASAPolynomial(coeffs=[12.3753,0.0352411,-1.49777e-05,2.70104e-09,-1.79526e-13,-23246.7,-34.5074], Tmin=(1926.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-159.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'CC=C=CCCO(3241)',
    structure = SMILES('CC=C=CCCO'),
    E0 = (-71.7961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01901,0.0667308,-5.03926e-05,2.10646e-08,-3.78755e-12,-8528.95,26.1452], Tmin=(100,'K'), Tmax=(1259.17,'K')), NASAPolynomial(coeffs=[10.0105,0.0381673,-1.63658e-05,3.04906e-09,-2.10663e-13,-10793.3,-19.3058], Tmin=(1259.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.7961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds)"""),
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
    label = 'C=C[CH]CC[O](3242)',
    structure = SMILES('[CH2]C=CCC[O]'),
    E0 = (164.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,245.422,245.638,2444.53],'cm^-1')),
        HinderedRotor(inertia=(0.00245214,'amu*angstrom^2'), symmetry=1, barrier=(10.4,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766442,'amu*angstrom^2'), symmetry=1, barrier=(32.7797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766241,'amu*angstrom^2'), symmetry=1, barrier=(32.778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98027,0.0464601,-2.5237e-05,6.38765e-09,-6.46696e-13,19894.2,22.9639], Tmin=(100,'K'), Tmax=(2131.54,'K')), NASAPolynomial(coeffs=[12.6767,0.0263877,-1.11119e-05,1.9699e-09,-1.28561e-13,15334.1,-36.7357], Tmin=(2131.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(164.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=CC)C[O](2186)',
    structure = SMILES('[CH2]C(C=CC)C[O]'),
    E0 = (186.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,1198.17,1198.42],'cm^-1')),
        HinderedRotor(inertia=(0.249587,'amu*angstrom^2'), symmetry=1, barrier=(5.73849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250246,'amu*angstrom^2'), symmetry=1, barrier=(5.75365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249709,'amu*angstrom^2'), symmetry=1, barrier=(5.7413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249535,'amu*angstrom^2'), symmetry=1, barrier=(5.7373,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32951,0.0627137,-4.78961e-05,2.42306e-08,-5.88253e-12,22542.8,28.8295], Tmin=(100,'K'), Tmax=(911.733,'K')), NASAPolynomial(coeffs=[5.2164,0.0456614,-1.98423e-05,3.718e-09,-2.58069e-13,21834,10.4366], Tmin=(911.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(CCOJ)"""),
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
    label = '[CH2]C[O](362)',
    structure = SMILES('[CH2]C[O]'),
    E0 = (188.892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1398.33],'cm^-1')),
        HinderedRotor(inertia=(0.00547724,'amu*angstrom^2'), symmetry=1, barrier=(7.58298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57172,0.0102136,5.90928e-06,-7.99883e-09,2.07083e-12,22733,11.7516], Tmin=(100,'K'), Tmax=(1490.83,'K')), NASAPolynomial(coeffs=[4.74087,0.0150202,-6.91924e-06,1.31181e-09,-8.98257e-14,21501.7,2.68369], Tmin=(1490.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CCOJ) + radical(CJCO)"""),
)

species(
    label = '[CH]C=CC(1185)',
    structure = SMILES('[CH]C=CC'),
    E0 = (342.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,213.913,213.918,213.995],'cm^-1')),
        HinderedRotor(inertia=(1.56238,'amu*angstrom^2'), symmetry=1, barrier=(51.0832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57548,'amu*angstrom^2'), symmetry=1, barrier=(51.0779,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51998,0.0276813,9.43783e-07,-1.25009e-08,4.67145e-12,41238.4,14.8113], Tmin=(100,'K'), Tmax=(1157.79,'K')), NASAPolynomial(coeffs=[6.04226,0.0260072,-1.04843e-05,1.90881e-09,-1.30601e-13,39719.4,-5.73539], Tmin=(1157.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
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
    label = '[CH2]C[CH]C=CC(1186)',
    structure = SMILES('[CH2]C[CH]C=CC'),
    E0 = (270.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,258.401,1598.77],'cm^-1')),
        HinderedRotor(inertia=(0.104608,'amu*angstrom^2'), symmetry=1, barrier=(4.96624,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.465896,'amu*angstrom^2'), symmetry=1, barrier=(21.9635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10458,'amu*angstrom^2'), symmetry=1, barrier=(4.96392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0121183,'amu*angstrom^2'), symmetry=1, barrier=(21.9632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43521,0.0478989,-1.06593e-05,-1.21243e-08,5.88414e-12,32576.7,24.0216], Tmin=(100,'K'), Tmax=(1133.02,'K')), NASAPolynomial(coeffs=[10.0362,0.0339424,-1.3905e-05,2.56703e-09,-1.77702e-13,29574.5,-23.1948], Tmin=(1133.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]C1CCO1(3243)',
    structure = SMILES('C[CH][CH]C1CCO1'),
    E0 = (214.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61033,0.0421349,2.25474e-05,-5.8026e-08,2.65233e-11,25839.5,27.6726], Tmin=(100,'K'), Tmax=(892.237,'K')), NASAPolynomial(coeffs=[9.963,0.0329561,-9.54341e-06,1.45945e-09,-9.32836e-14,23223.9,-17.9771], Tmin=(892.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(RCCJC) + radical(CCJCO)"""),
)

species(
    label = 'C=CC=CCC[O](3244)',
    structure = SMILES('C=CC=CCC[O]'),
    E0 = (101.128,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,288.464,288.469,288.469,288.472],'cm^-1')),
        HinderedRotor(inertia=(0.00202583,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211005,'amu*angstrom^2'), symmetry=1, barrier=(12.4597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210999,'amu*angstrom^2'), symmetry=1, barrier=(12.4597,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06875,0.0607627,-4.05668e-05,1.39196e-08,-1.97888e-12,12271,26.8012], Tmin=(100,'K'), Tmax=(1589.37,'K')), NASAPolynomial(coeffs=[12.6918,0.0315109,-1.296e-05,2.33994e-09,-1.57463e-13,8576.33,-34.6589], Tmin=(1589.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]CC=CCC[O](3245)',
    structure = SMILES('[CH2]CC=CCC[O]'),
    E0 = (195.931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,263.206,263.385,1390.36,3622.87],'cm^-1')),
        HinderedRotor(inertia=(0.00242938,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134905,'amu*angstrom^2'), symmetry=1, barrier=(6.63602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.6514,'amu*angstrom^2'), symmetry=1, barrier=(32.0083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134875,'amu*angstrom^2'), symmetry=1, barrier=(6.63659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80558,0.0567825,-3.14569e-05,8.39812e-09,-9.46501e-13,23635.1,27.9821], Tmin=(100,'K'), Tmax=(1743.78,'K')), NASAPolynomial(coeffs=[8.41213,0.0416279,-1.84209e-05,3.41432e-09,-2.3199e-13,21331.1,-7.56426], Tmin=(1743.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = 'CC[C]=CCC[O](3246)',
    structure = SMILES('CC[C]=CCC[O]'),
    E0 = (228.527,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3010,987.5,1337.5,450,1655,199.194,199.196,199.198,2472.12],'cm^-1')),
        HinderedRotor(inertia=(0.00424841,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319339,'amu*angstrom^2'), symmetry=1, barrier=(8.99234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319361,'amu*angstrom^2'), symmetry=1, barrier=(8.99234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00424846,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25565,0.0485394,3.77281e-05,-1.88653e-07,1.76218e-10,27537.4,25.3834], Tmin=(100,'K'), Tmax=(405.51,'K')), NASAPolynomial(coeffs=[3.70052,0.0496215,-2.2997e-05,4.43348e-09,-3.12236e-13,27294.1,18.1623], Tmin=(405.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'CCC=[C]CC[O](3247)',
    structure = SMILES('CCC=[C]CC[O]'),
    E0 = (228.527,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3010,987.5,1337.5,450,1655,199.194,199.196,199.198,2472.12],'cm^-1')),
        HinderedRotor(inertia=(0.00424841,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319339,'amu*angstrom^2'), symmetry=1, barrier=(8.99234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319361,'amu*angstrom^2'), symmetry=1, barrier=(8.99234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00424846,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25565,0.0485394,3.77281e-05,-1.88653e-07,1.76218e-10,27537.4,25.3834], Tmin=(100,'K'), Tmax=(405.51,'K')), NASAPolynomial(coeffs=[3.70052,0.0496215,-2.2997e-05,4.43348e-09,-3.12236e-13,27294.1,18.1623], Tmin=(405.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'C[CH]C=[C]CCO(3248)',
    structure = SMILES('C[CH]C=[C]CCO'),
    E0 = (143.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01195,0.0643142,-4.44698e-05,1.6304e-08,-2.52478e-12,17419.6,27.5879], Tmin=(100,'K'), Tmax=(1450.9,'K')), NASAPolynomial(coeffs=[11.3777,0.0357371,-1.49258e-05,2.72901e-09,-1.85736e-13,14411.6,-26.2787], Tmin=(1450.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = 'CCC=C[CH]C[O](3249)',
    structure = SMILES('CCC=C[CH]C[O]'),
    E0 = (107.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24416,0.0632969,-3.90518e-05,1.19708e-08,-1.53983e-12,13037.3,25.4339], Tmin=(100,'K'), Tmax=(1661.75,'K')), NASAPolynomial(coeffs=[11.3411,0.0389924,-1.71129e-05,3.16921e-09,-2.15688e-13,9681.58,-28.406], Tmin=(1661.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = 'CCC=CC[CH][O](3250)',
    structure = SMILES('CCC=CC[CH][O]'),
    E0 = (170.983,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,305.233,305.238,305.244,2073.07],'cm^-1')),
        HinderedRotor(inertia=(0.00180918,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143664,'amu*angstrom^2'), symmetry=1, barrier=(9.49855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143672,'amu*angstrom^2'), symmetry=1, barrier=(9.49869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143684,'amu*angstrom^2'), symmetry=1, barrier=(9.49875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14556,0.0681682,-6.11345e-05,3.71938e-08,-1.06666e-11,20662.5,28.1764], Tmin=(100,'K'), Tmax=(795.65,'K')), NASAPolynomial(coeffs=[5.30238,0.0472707,-2.17379e-05,4.1842e-09,-2.94807e-13,20001,9.07233], Tmin=(795.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'CCC=CCC=O(3251)',
    structure = SMILES('CCC=CCC=O'),
    E0 = (-151.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2195,0.0491283,2.00958e-06,-2.8429e-08,1.12906e-11,-18148.5,26.3787], Tmin=(100,'K'), Tmax=(1138.82,'K')), NASAPolynomial(coeffs=[12.6041,0.0361546,-1.64832e-05,3.22592e-09,-2.30998e-13,-22493.2,-37.7161], Tmin=(1138.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-151.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'C=CC=CCCO(3252)',
    structure = SMILES('C=CC=CCCO'),
    E0 = (-124.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576141,0.0651247,-3.44981e-05,-1.91078e-09,5.53915e-12,-14851.1,27.7061], Tmin=(100,'K'), Tmax=(1021.28,'K')), NASAPolynomial(coeffs=[14.9839,0.0288468,-1.08133e-05,1.94924e-09,-1.35287e-13,-18844.9,-47.2516], Tmin=(1021.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC1C=CCCO1(2191)',
    structure = SMILES('CC1C=CCCO1'),
    E0 = (-162.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63377,0.0305462,7.32927e-05,-1.17723e-07,4.79321e-11,-19446.9,19.6843], Tmin=(100,'K'), Tmax=(940.279,'K')), NASAPolynomial(coeffs=[15.2916,0.0268984,-7.75592e-06,1.33089e-09,-9.73101e-14,-24422.5,-58.1663], Tmin=(940.279,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro2hpyran)"""),
)

species(
    label = 'CHCH3(T)(359)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438699,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82363,-0.000909515,3.2138e-05,-3.7348e-08,1.3309e-11,41371.4,7.10948], Tmin=(100,'K'), Tmax=(960.812,'K')), NASAPolynomial(coeffs=[4.30487,0.00943069,-3.27559e-06,5.95121e-10,-4.27307e-14,40709.1,1.84202], Tmin=(960.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=CCC[O](3253)',
    structure = SMILES('[CH]=CCC[O]'),
    E0 = (296.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3120,650,792.5,1650,246.951,1525.25],'cm^-1')),
        HinderedRotor(inertia=(0.160242,'amu*angstrom^2'), symmetry=1, barrier=(6.9343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160245,'amu*angstrom^2'), symmetry=1, barrier=(6.93433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.02091,0.0321044,-1.57854e-05,3.34523e-09,-2.69928e-13,35678.3,17.4877], Tmin=(100,'K'), Tmax=(2745.35,'K')), NASAPolynomial(coeffs=[13.8715,0.0162945,-7.14691e-06,1.24744e-09,-7.8891e-14,29720.8,-45.8177], Tmin=(2745.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
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
    E0 = (130.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (148.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (307.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (218.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (381.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (255.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (179.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (314.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (214.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (427.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (294.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (205.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (388.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (559.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (294.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (320.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (237.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (432.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (366.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (164.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (208.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (219.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (155.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (583.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (346.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (138.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (531.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (513.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (252.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (312.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (309.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (386.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (377.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (245.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (549.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (238.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (139.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (155.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (138.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (640.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['CH2O(35)', 'C=CC=CC(381)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['C[CH]C1[CH]CCO1(3220)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.2, Ea=(17.9912,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R6;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', 'CC=CC=CC[O](3225)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.35e+08,'cm^3/(mol*s)'), n=1.64, Ea=(1.58992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2557 used for Cds-CsH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'CC=C[CH]CC=O(3226)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.6e+09,'cm^3/(mol*s)'), n=0.935, Ea=(17.4473,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2782 used for CO-CsH_O;HJ
Exact match found for rate rule [CO-CsH_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(19)', 'CC=C=CCC[O](3227)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=CC=CC(381)', '[CH2][O](354)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0267117,'m^3/(mol*s)'), n=2.459, Ea=(4.91982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CdH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2O(35)', '[CH2]C=C[CH]C(377)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(225.36,'m^3/(mol*s)'), n=0.996465, Ea=(58.8821,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-HH_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CC=CC[CH]C[O](3228)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC=C[CH]C[CH]O(3229)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4500,'s^-1'), n=2.62, Ea=(129.286,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 322 used for R2H_S;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC=[C]CCC[O](3230)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CC=CCC[CH][O](3231)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['CC=C[CH][CH]CO(3232)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['C[C]=CCCC[O](3233)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.96e+10,'s^-1'), n=0.83, Ea=(257.734,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 201 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C=CCCC[O](3234)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CC=[C][CH]CCO(3235)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[C]=C[CH]CCO(3236)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.44313e+09,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;XH_out] for rate rule [R6HJ_2;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['[CH2]C=C[CH]CCO(3237)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.85e+08,'s^-1'), n=0, Ea=(106.901,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7Hall;O_rad_out;Cs_H_out_2H] for rate rule [R7HJ_3;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C=C[CH]C(377)', '[CH2][O](354)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.35017e+06,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['CC1[CH]C1CC[O](3238)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.48125e+10,'s^-1'), n=0.454219, Ea=(235.628,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;doublebond_intra_pri_HNd_Cs;radadd_intra_cs] + [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['CC1[CH][CH]CCO1(3157)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.79146e+14,'s^-1'), n=-0.345723, Ea=(33.7118,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['CC=CC=CCO(3239)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['CC=CCCC=O(3240)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['CC=C=CCCO(3241)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(S)(40)', 'C=C[CH]CC[O](3242)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.97e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C=CC)C[O](2186)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['CC=CC1CCO1(2193)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C[O](362)', '[CH]C=CC(1185)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['O(20)', '[CH2]C[CH]C=CC(1186)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/Cs;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['C[CH][CH]C1CCO1(3243)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(19)', 'C=CC=CCC[O](3244)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]CC=CCC[O](3245)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CC[C]=CCC[O](3246)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CCC=[C]CC[O](3247)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 200 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C[CH]C=[C]CCO(3248)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['CCC=C[CH]C[O](3249)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CCC=CC[CH][O](3250)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(662081,'s^-1'), n=1.61604, Ea=(67.5001,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;Y_rad_out;Cs_H_out_H/NonDeC] + [R5H_SSMS;Y_rad_out;Cs_H_out] for rate rule [R5H_SSMS;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['CCC=CCC=O(3251)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['C=CC=CCCO(3252)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CC=C[CH]CC[O](2185)'],
    products = ['CC1C=CCCO1(2191)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CHCH3(T)(359)', '[CH]=CCC[O](3253)'],
    products = ['CC=C[CH]CC[O](2185)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

network(
    label = '438',
    isomers = [
        'CC=C[CH]CC[O](2185)',
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
    label = '438',
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

