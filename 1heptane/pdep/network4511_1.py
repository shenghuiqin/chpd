species(
    label = '[CH]=C(C=C)C([CH2])CC(26541)',
    structure = SMILES('[CH]=C(C=C)C([CH2])CC'),
    E0 = (430.331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.658357,'amu*angstrom^2'), symmetry=1, barrier=(15.1369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658289,'amu*angstrom^2'), symmetry=1, barrier=(15.1354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.413876,'amu*angstrom^2'), symmetry=1, barrier=(15.1482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0201855,'amu*angstrom^2'), symmetry=1, barrier=(15.14,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0201747,'amu*angstrom^2'), symmetry=1, barrier=(15.1319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.288468,0.082486,-4.70324e-05,-1.4862e-09,7.93038e-12,51921.9,33.3284], Tmin=(100,'K'), Tmax=(971.479,'K')), NASAPolynomial(coeffs=[18.0826,0.0331499,-1.14724e-05,1.984e-09,-1.35411e-13,47111.1,-61.1586], Tmin=(971.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = 'butene1(127)',
    structure = SMILES('C=CCC'),
    E0 = (-16.4325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,252.555],'cm^-1')),
        HinderedRotor(inertia=(0.178654,'amu*angstrom^2'), symmetry=1, barrier=(7.72883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185589,'amu*angstrom^2'), symmetry=1, barrier=(7.72103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58773,0.0232778,1.93412e-05,-3.55496e-08,1.36906e-11,-1918.73,14.5751], Tmin=(100,'K'), Tmax=(1007.28,'K')), NASAPolynomial(coeffs=[7.20517,0.0236362,-9.0315e-06,1.65393e-09,-1.16019e-13,-3797.34,-12.4426], Tmin=(1007.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.4325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene1""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2CHCCH(26391)',
    structure = SMILES('C#CC=C'),
    E0 = (274.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,3010,987.5,1337.5,450,1655,2175,525],'cm^-1')),
        HinderedRotor(inertia=(1.46338,'amu*angstrom^2'), symmetry=1, barrier=(33.6459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.87083,0.0182042,1.06711e-05,-2.72492e-08,1.19478e-11,33023.8,11.2934], Tmin=(100,'K'), Tmax=(955.249,'K')), NASAPolynomial(coeffs=[8.52653,0.0108962,-3.56564e-06,6.31243e-10,-4.51891e-14,31196.2,-19.6435], Tmin=(955.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""CH2CHCCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C1C=C1C([CH2])CC(26961)',
    structure = SMILES('[CH2]C1C=C1C([CH2])CC'),
    E0 = (534.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.366574,0.0871287,-7.25398e-05,3.35469e-08,-6.33013e-12,64391.6,32.7978], Tmin=(100,'K'), Tmax=(1269.18,'K')), NASAPolynomial(coeffs=[15.7091,0.0364649,-1.26632e-05,2.09595e-09,-1.35111e-13,60311,-48.5904], Tmin=(1269.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C1C([CH2])CC1CC(26962)',
    structure = SMILES('[CH]=C1C([CH2])CC1CC'),
    E0 = (465.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[9.80087,0.00423322,0.000119419,-1.21402e-07,2.82521e-11,55632.1,-10.6516], Tmin=(100,'K'), Tmax=(1719.04,'K')), NASAPolynomial(coeffs=[77.3893,0.0364629,-7.40585e-05,1.77572e-08,-1.31187e-12,4395.19,-454.783], Tmin=(1719.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Cds_P) + radical(Isobutyl)"""),
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
    label = '[CH]=C(C=C)C(=C)CC(26963)',
    structure = SMILES('[CH]=C(C=C)C(=C)CC'),
    E0 = (351.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,180,391.097],'cm^-1')),
        HinderedRotor(inertia=(0.788875,'amu*angstrom^2'), symmetry=1, barrier=(18.1378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788881,'amu*angstrom^2'), symmetry=1, barrier=(18.1379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0200846,'amu*angstrom^2'), symmetry=1, barrier=(18.1462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.78923,'amu*angstrom^2'), symmetry=1, barrier=(18.1459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.213862,0.0706338,-2.72841e-05,-1.35302e-08,9.86507e-12,42441.5,28.5391], Tmin=(100,'K'), Tmax=(1021.68,'K')), NASAPolynomial(coeffs=[15.9902,0.0350115,-1.33679e-05,2.43503e-09,-1.70049e-13,37853.3,-54.5886], Tmin=(1021.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C2H5(29)',
    structure = SMILES('C[CH2]'),
    E0 = (107.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1190.6,1642.82,1642.96,3622.23,3622.39],'cm^-1')),
        HinderedRotor(inertia=(0.866817,'amu*angstrom^2'), symmetry=1, barrier=(19.9298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24186,-0.00356905,4.82667e-05,-5.85401e-08,2.25805e-11,12969,4.44704], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.32196,0.0123931,-4.39681e-06,7.0352e-10,-4.18435e-14,12175.9,0.171104], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(107.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H5""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH]=C(C=C)C=C(26631)',
    structure = SMILES('[CH]=C(C=C)C=C'),
    E0 = (411.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.17526,'amu*angstrom^2'), symmetry=1, barrier=(27.0215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17459,'amu*angstrom^2'), symmetry=1, barrier=(27.0061,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75262,0.0382342,7.41856e-06,-3.70507e-08,1.69951e-11,49629.8,19.502], Tmin=(100,'K'), Tmax=(974.113,'K')), NASAPolynomial(coeffs=[12.41,0.0212121,-7.54597e-06,1.37118e-09,-9.78936e-14,46284.8,-38.1462], Tmin=(974.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C=C(4699)',
    structure = SMILES('[CH]=C=C[CH2]'),
    E0 = (451.584,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3000,3100,440,815,1455,1000,180,1024.85,1025.53,1026.61],'cm^-1')),
        HinderedRotor(inertia=(0.00938781,'amu*angstrom^2'), symmetry=1, barrier=(7.01846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76805,0.020302,8.75519e-06,-2.87666e-08,1.37354e-11,54363.7,13.5565], Tmin=(100,'K'), Tmax=(915.031,'K')), NASAPolynomial(coeffs=[9.46747,0.00887314,-1.78262e-06,2.38534e-10,-1.6263e-14,52390.1,-22.2544], Tmin=(915.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C=CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC(130)',
    structure = SMILES('[CH2][CH]CC'),
    E0 = (255.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2031.24],'cm^-1')),
        HinderedRotor(inertia=(0.244974,'amu*angstrom^2'), symmetry=1, barrier=(5.63244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00192352,'amu*angstrom^2'), symmetry=1, barrier=(5.63177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244928,'amu*angstrom^2'), symmetry=1, barrier=(5.63137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98997,0.0287412,-9.51469e-06,4.19232e-10,1.90526e-13,30780.1,16.8971], Tmin=(100,'K'), Tmax=(2154.56,'K')), NASAPolynomial(coeffs=[12.4231,0.0182241,-7.06316e-06,1.16769e-09,-7.11818e-14,25091.4,-39.6212], Tmin=(2154.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C2H3(30)',
    structure = SMILES('[CH]=C'),
    E0 = (286.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,677.08,1086.68,3788.01],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C#CC([CH2])CC(3268)',
    structure = SMILES('C#CC([CH2])CC'),
    E0 = (297.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,333.743],'cm^-1')),
        HinderedRotor(inertia=(0.00340305,'amu*angstrom^2'), symmetry=1, barrier=(11.8362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149684,'amu*angstrom^2'), symmetry=1, barrier=(11.8359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.927664,'amu*angstrom^2'), symmetry=1, barrier=(73.3293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.927823,'amu*angstrom^2'), symmetry=1, barrier=(73.3286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31625,0.0530696,-2.51599e-05,-5.62724e-09,7.0068e-12,35940.6,23.5452], Tmin=(100,'K'), Tmax=(915.067,'K')), NASAPolynomial(coeffs=[10.828,0.0273274,-8.9215e-06,1.45456e-09,-9.4854e-14,33536.8,-25.1218], Tmin=(915.067,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C(C=C)[C](C)CC(26964)',
    structure = SMILES('[CH]C(C=C)=C(C)CC'),
    E0 = (329.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.44218,0.0858936,-5.14325e-05,9.59426e-09,1.54289e-12,39774.2,31.1284], Tmin=(100,'K'), Tmax=(1106.5,'K')), NASAPolynomial(coeffs=[16.358,0.042711,-1.6684e-05,2.9922e-09,-2.04006e-13,34981.9,-56.4775], Tmin=(1106.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C=C)C(C)[CH]C(26965)',
    structure = SMILES('[CH]=C(C=C)C(C)[CH]C'),
    E0 = (419.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.189923,0.0803918,-4.84126e-05,6.88763e-09,2.7535e-12,50650.1,33.7531], Tmin=(100,'K'), Tmax=(1076.24,'K')), NASAPolynomial(coeffs=[17.0013,0.0361248,-1.40704e-05,2.55923e-09,-1.77086e-13,45813,-55.7286], Tmin=(1076.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cs_S)"""),
)

species(
    label = '[CH2][C](C=C)C(=C)CC(24623)',
    structure = SMILES('[CH2]C(C=C)=C([CH2])CC'),
    E0 = (228.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0497728,0.0733281,-1.6094e-05,-3.35123e-08,1.88363e-11,27601.1,30.4448], Tmin=(100,'K'), Tmax=(975.095,'K')), NASAPolynomial(coeffs=[18.3695,0.0342638,-1.21408e-05,2.16747e-09,-1.52112e-13,22274,-66.8493], Tmin=(975.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(CC)C(=C)[C]=C(26966)',
    structure = SMILES('[CH2]C(CC)C(=C)[C]=C'),
    E0 = (382.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.298592,0.0855799,-5.97426e-05,1.43922e-08,1.90134e-12,46134.5,32.5581], Tmin=(100,'K'), Tmax=(969.924,'K')), NASAPolynomial(coeffs=[16.3223,0.0360349,-1.25044e-05,2.12015e-09,-1.41251e-13,42016.5,-51.7276], Tmin=(969.924,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(382.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]C)C(=C)C=C(26967)',
    structure = SMILES('[CH2]C([CH]C)C(=C)C=C'),
    E0 = (377.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.057099,0.0781219,-4.3007e-05,7.91346e-10,5.38807e-12,45591.8,34.6518], Tmin=(100,'K'), Tmax=(1016.68,'K')), NASAPolynomial(coeffs=[16.0662,0.0363793,-1.3425e-05,2.3799e-09,-1.6304e-13,41192.3,-48.9139], Tmin=(1016.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cs_S)"""),
)

species(
    label = '[CH]=C(C=C)C(C)C[CH2](26968)',
    structure = SMILES('[CH]=C(C=C)C(C)C[CH2]'),
    E0 = (430.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,341.498,341.498],'cm^-1')),
        HinderedRotor(inertia=(0.187077,'amu*angstrom^2'), symmetry=1, barrier=(15.4819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00144552,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187078,'amu*angstrom^2'), symmetry=1, barrier=(15.4819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187078,'amu*angstrom^2'), symmetry=1, barrier=(15.4819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187076,'amu*angstrom^2'), symmetry=1, barrier=(15.4819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.299766,0.081879,-4.61112e-05,-3.82873e-10,6.57827e-12,51942.5,33.3747], Tmin=(100,'K'), Tmax=(1014.24,'K')), NASAPolynomial(coeffs=[18.315,0.0338468,-1.26123e-05,2.27196e-09,-1.57969e-13,46861.1,-63.1299], Tmin=(1014.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C([C]=C)C(C)CC(26969)',
    structure = SMILES('[CH]C(=C=C)C(C)CC'),
    E0 = (400.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.455209,0.0862452,-5.74997e-05,1.96918e-08,-2.76077e-12,48297.5,33.0617], Tmin=(100,'K'), Tmax=(1643.13,'K')), NASAPolynomial(coeffs=[18.3716,0.0404133,-1.56597e-05,2.71584e-09,-1.77876e-13,42110.6,-67.1158], Tmin=(1643.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC(=C)C([CH2])CC(26970)',
    structure = SMILES('[CH]=CC(=C)C([CH2])CC'),
    E0 = (430.331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.658357,'amu*angstrom^2'), symmetry=1, barrier=(15.1369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658289,'amu*angstrom^2'), symmetry=1, barrier=(15.1354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.413876,'amu*angstrom^2'), symmetry=1, barrier=(15.1482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0201855,'amu*angstrom^2'), symmetry=1, barrier=(15.14,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0201747,'amu*angstrom^2'), symmetry=1, barrier=(15.1319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.288468,0.082486,-4.70324e-05,-1.4862e-09,7.93038e-12,51921.9,33.3284], Tmin=(100,'K'), Tmax=(971.479,'K')), NASAPolynomial(coeffs=[18.0826,0.0331499,-1.14724e-05,1.984e-09,-1.35411e-13,47111.1,-61.1586], Tmin=(971.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC([CH2])C(=C)C=C(26971)',
    structure = SMILES('[CH2]CC([CH2])C(=C)C=C'),
    E0 = (388.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.175911,0.0796777,-4.07565e-05,-6.7114e-09,9.45239e-12,46884.7,34.3078], Tmin=(100,'K'), Tmax=(973.823,'K')), NASAPolynomial(coeffs=[17.6185,0.0337181,-1.1755e-05,2.04413e-09,-1.39997e-13,42132.5,-57.6734], Tmin=(973.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=CC(=[CH])C(C)CC(26972)',
    structure = SMILES('[CH]=CC(=[CH])C(C)CC'),
    E0 = (472.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,260.132],'cm^-1')),
        HinderedRotor(inertia=(0.307683,'amu*angstrom^2'), symmetry=1, barrier=(14.7752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307632,'amu*angstrom^2'), symmetry=1, barrier=(14.7756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307661,'amu*angstrom^2'), symmetry=1, barrier=(14.7753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307706,'amu*angstrom^2'), symmetry=1, barrier=(14.7745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307764,'amu*angstrom^2'), symmetry=1, barrier=(14.7746,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.415733,0.0847291,-5.25413e-05,5.05354e-09,4.96126e-12,56979.8,32.4074], Tmin=(100,'K'), Tmax=(1014.25,'K')), NASAPolynomial(coeffs=[18.7883,0.0332626,-1.23203e-05,2.20957e-09,-1.53195e-13,51835.9,-66.6668], Tmin=(1014.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(472.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(CC)C1[CH]CC=1(26973)',
    structure = SMILES('[CH2]C(CC)C1[CH]CC=1'),
    E0 = (396.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.712144,0.0531909,3.28174e-05,-7.88837e-08,3.37232e-11,47866.8,28.9959], Tmin=(100,'K'), Tmax=(966.88,'K')), NASAPolynomial(coeffs=[16.2153,0.0357154,-1.24613e-05,2.24914e-09,-1.60532e-13,42687.8,-56.555], Tmin=(966.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Isobutyl) + radical(cyclobutene-allyl)"""),
)

species(
    label = '[CH]=C1[CH]CCC1CC(26974)',
    structure = SMILES('[CH]C1=CCCC1CC'),
    E0 = (298.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.715496,0.0512987,4.88061e-05,-9.30531e-08,3.72597e-11,35989.3,27.8709], Tmin=(100,'K'), Tmax=(988.132,'K')), NASAPolynomial(coeffs=[15.1539,0.0436255,-1.66218e-05,3.09058e-09,-2.21353e-13,30657.1,-54.1567], Tmin=(988.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC(=C)C(=C)CC(24630)',
    structure = SMILES('C=CC(=C)C(=C)CC'),
    E0 = (104.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.296747,0.0670054,-1.0269e-05,-3.13536e-08,1.59568e-11,12721.3,27.8384], Tmin=(100,'K'), Tmax=(1010.3,'K')), NASAPolynomial(coeffs=[15.6889,0.0379462,-1.44599e-05,2.64736e-09,-1.86033e-13,7984.11,-54.6302], Tmin=(1010.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    label = '[CH]=C(C=C)C([CH2])C(26551)',
    structure = SMILES('[CH]=C(C=C)C([CH2])C'),
    E0 = (454.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,546.426],'cm^-1')),
        HinderedRotor(inertia=(0.184121,'amu*angstrom^2'), symmetry=1, barrier=(14.3488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.914455,'amu*angstrom^2'), symmetry=1, barrier=(21.0251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.623998,'amu*angstrom^2'), symmetry=1, barrier=(14.3469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0678068,'amu*angstrom^2'), symmetry=1, barrier=(14.3487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3398.3,'J/mol'), sigma=(6.02224,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=530.81 K, Pc=35.3 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.456642,0.0665153,-2.89033e-05,-1.53206e-08,1.26391e-11,54755,28.4102], Tmin=(100,'K'), Tmax=(943.627,'K')), NASAPolynomial(coeffs=[16.7291,0.0254789,-8.08814e-06,1.3537e-09,-9.22647e-14,50439.9,-55.7424], Tmin=(943.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C([CH]CCC)C=C(26975)',
    structure = SMILES('[CH]C(C=C)=CCCC'),
    E0 = (344.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.276246,0.0810984,-3.65024e-05,-6.62698e-09,7.38704e-12,41606.5,32.0183], Tmin=(100,'K'), Tmax=(1049.73,'K')), NASAPolynomial(coeffs=[16.2368,0.0426921,-1.66547e-05,3.01657e-09,-2.08216e-13,36788.9,-54.883], Tmin=(1049.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C=C)C[CH]CC(26537)',
    structure = SMILES('[CH]=C(C=C)C[CH]CC'),
    E0 = (427.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.132537,0.0817263,-6.15919e-05,2.51653e-08,-4.30676e-12,51579.7,32.8448], Tmin=(100,'K'), Tmax=(1348.33,'K')), NASAPolynomial(coeffs=[13.6805,0.0415346,-1.68793e-05,3.05769e-09,-2.07717e-13,47926.2,-36.5655], Tmin=(1348.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1=CCC1CC(26976)',
    structure = SMILES('C=CC1=CCC1CC'),
    E0 = (139.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.52436,0.0534286,4.26854e-05,-9.63812e-08,4.14791e-11,16880.9,26.7556], Tmin=(100,'K'), Tmax=(960.264,'K')), NASAPolynomial(coeffs=[19.7128,0.0303293,-1.0005e-05,1.8302e-09,-1.35107e-13,10575.5,-78.6823], Tmin=(960.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = 'CH2(19)',
    structure = SMILES('[CH2]'),
    E0 = (381.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1032.72,2936.3,3459],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8328,0.000224446,4.68033e-06,-6.04743e-09,2.59009e-12,45920.8,1.40666], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16229,0.00281798,-7.56235e-07,5.05446e-11,5.65236e-15,46099.1,4.77656], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(381.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH]=C([CH]CC)C=C(26786)',
    structure = SMILES('[CH]C(C=C)=CCC'),
    E0 = (368.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.486824,0.0649554,-1.79811e-05,-2.06189e-08,1.19882e-11,44438.8,27.0332], Tmin=(100,'K'), Tmax=(1009.03,'K')), NASAPolynomial(coeffs=[14.6011,0.0354818,-1.35284e-05,2.44586e-09,-1.6993e-13,40242.6,-47.8661], Tmin=(1009.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
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
    E0 = (430.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (534.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (489.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (574.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (542.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (454.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (556.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (605.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (532.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (577.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (622.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (652.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (474.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (513.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (541.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (833.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (463.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (505.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (707.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (568.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (544.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (493.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (873.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (590.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (590.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (438.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (749.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['butene1(127)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['[CH2]C1C=C1C([CH2])CC(26961)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(103.848,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 14 used for R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH
Exact match found for rate rule [R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['[CH]=C1C([CH2])CC1CC(26962)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.07e+06,'s^-1'), n=1.46, Ea=(58.9944,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 343 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH]=C(C=C)C(=C)CC(26963)'],
    products = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(72.3521,'m^3/(mol*s)'), n=1.66655, Ea=(10.8198,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeCs_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H5(29)', '[CH]=C(C=C)C=C(26631)'],
    products = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00238215,'m^3/(mol*s)'), n=2.41, Ea=(22.4155,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;CsJ-CsHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['butene1(127)', '[CH]=[C]C=C(4699)'],
    products = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH]CC(130)', 'CH2CHCCH(26391)'],
    products = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0345031,'m^3/(mol*s)'), n=2.33733, Ea=(26.7584,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-Cd_Ct-H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2H3(30)', 'C#CC([CH2])CC(3268)'],
    products = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(94600,'cm^3/(mol*s)'), n=2.41, Ea=(20.7945,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2340 used for Ct-Cs_Ct-H;CdsJ-H
Exact match found for rate rule [Ct-Cs_Ct-H;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(C=C)[C](C)CC(26964)'],
    products = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['[CH]=C(C=C)C(C)[CH]C(26965)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['[CH2][C](C=C)C(=C)CC(24623)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['[CH2]C(CC)C(=C)[C]=C(26966)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['[CH2]C([CH]C)C(=C)C=C(26967)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(C=C)C(C)C[CH2](26968)'],
    products = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['[CH]=C([C]=C)C(C)CC(26969)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(356.546,'s^-1'), n=2.79552, Ea=(110.814,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;Cd_H_out_doubleC] + [R4H_SSS;C_rad_out_2H;XH_out] for rate rule [R4H_SSS;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=CC(=C)C([CH2])CC(26970)'],
    products = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.27529e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_singleH] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['[CH2]CC([CH2])C(=C)C=C(26971)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC(=[CH])C(C)CC(26972)'],
    products = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]CC(130)', '[CH]=[C]C=C(4699)'],
    products = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['[CH2]C(CC)C1[CH]CC=1(26973)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R4_D_D;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['[CH]=C1[CH]CCC1CC(26974)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.12e+09,'s^-1'), n=0.63, Ea=(114.642,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['C=CC(=C)C(=C)CC(24630)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(S)(23)', '[CH]=C(C=C)C([CH2])C(26551)'],
    products = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['[CH]=C([CH]CCC)C=C(26975)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['[CH]=C(C=C)C[CH]CC(26537)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    products = ['C=CC1=CCC1CC(26976)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(19)', '[CH]=C([CH]CC)C=C(26786)'],
    products = ['[CH]=C(C=C)C([CH2])CC(26541)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4511',
    isomers = [
        '[CH]=C(C=C)C([CH2])CC(26541)',
    ],
    reactants = [
        ('butene1(127)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4511',
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

