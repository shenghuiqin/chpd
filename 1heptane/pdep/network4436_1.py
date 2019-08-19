species(
    label = 'C#C[CH]CC([O])=O(26470)',
    structure = SMILES('[CH]=C=CCC([O])=O'),
    E0 = (161.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,540,610,2055,412.243,412.514,415.53,415.869,417.685,418.142],'cm^-1')),
        HinderedRotor(inertia=(0.000976976,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000972391,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12031,0.041257,-2.25608e-05,4.69213e-09,-2.78647e-13,19544,21.4056], Tmin=(100,'K'), Tmax=(1935.31,'K')), NASAPolynomial(coeffs=[17.3725,0.0163551,-8.39266e-06,1.57959e-09,-1.04969e-13,12400.3,-65.4519], Tmin=(1935.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=C=CJ)"""),
)

species(
    label = 'CO2(13)',
    structure = SMILES('O=C=O'),
    E0 = (-403.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([459.166,1086.67,1086.68,2300.05],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2779,0.00275783,7.12787e-06,-1.07855e-08,4.14228e-12,-48475.6,5.97856], Tmin=(100,'K'), Tmax=(988.185,'K')), NASAPolynomial(coeffs=[4.55071,0.00290728,-1.14643e-06,2.25798e-10,-1.69526e-14,-48986,-1.45662], Tmin=(988.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.131,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: BurkeH2O2"""),
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
    label = '[CH]=[C]C1CC(=O)O1(27931)',
    structure = SMILES('[CH]=[C]C1CC(=O)O1'),
    E0 = (255.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23757,0.0248049,3.82149e-05,-7.2098e-08,3.12966e-11,30831.9,23.5421], Tmin=(100,'K'), Tmax=(919.293,'K')), NASAPolynomial(coeffs=[13.0293,0.013924,-2.89528e-06,4.03241e-10,-2.90003e-14,27323.3,-35.9048], Tmin=(919.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(Beta-Propiolactone) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[O]C1([O])C=C=CC1(27926)',
    structure = SMILES('[O]C1([O])C=C=CC1'),
    E0 = (516.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13919,0.0536308,-5.1175e-05,2.50556e-08,-4.74759e-12,62196.4,21.5157], Tmin=(100,'K'), Tmax=(1396.71,'K')), NASAPolynomial(coeffs=[14.1919,0.0127382,-3.48727e-06,4.93709e-10,-2.90379e-14,58892.7,-44.591], Tmin=(1396.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene) + radical(C=CC(C)(O)OJ) + radical(C=CC(C)(O)OJ)"""),
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
    label = '[CH]=C=C[CH]C(=O)O(27932)',
    structure = SMILES('[CH]=C=C[CH]C(=O)O'),
    E0 = (53.1469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20471,0.0521916,-3.40806e-05,6.01521e-09,1.0919e-12,6500.46,21.7033], Tmin=(100,'K'), Tmax=(1180.62,'K')), NASAPolynomial(coeffs=[15.3227,0.0183223,-8.78948e-06,1.75157e-09,-1.26518e-13,2193.73,-52.8732], Tmin=(1180.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.1469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(C=C=CJ)"""),
)

species(
    label = 'C=C=[C]CC([O])=O(27933)',
    structure = SMILES('[CH2]C#CCC([O])=O'),
    E0 = (140.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27752,0.0337833,-7.04604e-06,-7.0623e-09,2.75045e-12,17021.7,22.6995], Tmin=(100,'K'), Tmax=(1406.14,'K')), NASAPolynomial(coeffs=[10.1768,0.0250444,-1.23724e-05,2.40793e-09,-1.68029e-13,13442.6,-22.9296], Tmin=(1406.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CtHH) + group(Cs-CtHHH) + group(Cds-OdCsOs) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Propargyl) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C=[C]CC(=O)O(27934)',
    structure = SMILES('[CH]=C=[C]CC(=O)O'),
    E0 = (174.072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,3120,650,792.5,1650,3615,1277.5,1000,1685,370,180,490.529,1124.05,1125.93,1127.26],'cm^-1')),
        HinderedRotor(inertia=(1.34017,'amu*angstrom^2'), symmetry=1, barrier=(30.8131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34018,'amu*angstrom^2'), symmetry=1, barrier=(30.8134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34076,'amu*angstrom^2'), symmetry=1, barrier=(30.8266,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7799,0.0484264,-3.64069e-05,1.31212e-08,-1.91776e-12,21015.7,22.8992], Tmin=(100,'K'), Tmax=(1561.36,'K')), NASAPolynomial(coeffs=[12.3267,0.0214064,-1.04484e-05,2.03729e-09,-1.43019e-13,17722.3,-32.6821], Tmin=(1561.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C]C=CC([O])=O(27916)',
    structure = SMILES('C=C=C[CH]C([O])=O'),
    E0 = (124.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81609,0.0448875,-2.38578e-05,3.37656e-09,2.78035e-13,15040.1,19.7756], Tmin=(100,'K'), Tmax=(1645.55,'K')), NASAPolynomial(coeffs=[16.6486,0.0197535,-1.09017e-05,2.16064e-09,-1.49951e-13,8679.98,-63.6628], Tmin=(1645.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(CCOJ)"""),
)

species(
    label = 'C3H2(81)',
    structure = SMILES('[CH]C#C'),
    E0 = (530.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.621997,'amu*angstrom^2'), symmetry=1, barrier=(14.3009,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6874,0.0353211,-7.31245e-05,7.03801e-08,-2.45313e-11,63833.3,8.93418], Tmin=(100,'K'), Tmax=(908.454,'K')), NASAPolynomial(coeffs=[4.78002,0.0100611,-4.92178e-06,8.868e-10,-5.66859e-14,64115.2,2.68365], Tmin=(908.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCCH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C([O])=O(1440)',
    structure = SMILES('[CH2]C([O])=O'),
    E0 = (-8.20105,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,841.694,844.034,844.61,844.694,845.141],'cm^-1')),
        HinderedRotor(inertia=(0.00480388,'amu*angstrom^2'), symmetry=1, barrier=(2.45209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.60873,0.0103835,4.164e-06,-6.72683e-09,1.73621e-12,-973.928,12.8424], Tmin=(100,'K'), Tmax=(1564.29,'K')), NASAPolynomial(coeffs=[5.70762,0.0133434,-6.65897e-06,1.2886e-09,-8.86347e-14,-2649.38,-1.47903], Tmin=(1564.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.20105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJCO) + radical(CCOJ)"""),
)

species(
    label = '[O][C]=O(669)',
    structure = SMILES('[O][C]=O'),
    E0 = (31.9507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90478,-0.000175995,8.15126e-06,-1.13656e-08,4.4768e-12,3848.25,8.04855], Tmin=(100,'K'), Tmax=(975.388,'K')), NASAPolynomial(coeffs=[5.59398,-0.00122084,7.11747e-07,-9.7712e-11,3.97995e-15,3238.91,-1.49318], Tmin=(975.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.9507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(OJC=O) + radical((O)CJOH)"""),
)

species(
    label = '[CH]=C=CC[C]1OO1(27935)',
    structure = SMILES('C#C[CH]C[C]1OO1'),
    E0 = (511.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.327466,0.0606457,-6.13953e-05,3.0862e-08,-5.68531e-12,61658,27.2885], Tmin=(100,'K'), Tmax=(1619.06,'K')), NASAPolynomial(coeffs=[15.0636,0.00906015,4.59645e-07,-3.97974e-10,3.56606e-14,58875.8,-44.7613], Tmin=(1619.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(511.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(dioxirane) + radical(Cs_P) + radical(Sec_Propargyl)"""),
)

species(
    label = '[O]C(=O)CC1[C]=C1(27936)',
    structure = SMILES('[O]C(=O)CC1[C]=C1'),
    E0 = (322.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68785,0.0591263,-9.1612e-05,9.05215e-08,-3.54077e-11,38891.8,20.1689], Tmin=(100,'K'), Tmax=(810.18,'K')), NASAPolynomial(coeffs=[2.41837,0.0377299,-1.9061e-05,3.71958e-09,-2.59909e-13,39357.2,20.4016], Tmin=(810.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cyclopropene) + radical(cyclopropenyl-vinyl) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C1[CH]CC(=O)O1(27767)',
    structure = SMILES('[CH]C1=CCC(=O)O1'),
    E0 = (48.6702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49325,0.0109463,9.3203e-05,-1.2457e-07,4.62422e-11,5927.65,19.2864], Tmin=(100,'K'), Tmax=(991.276,'K')), NASAPolynomial(coeffs=[13.5866,0.0225692,-9.70978e-06,2.02655e-09,-1.57589e-13,957.958,-48.1094], Tmin=(991.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.6702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cyclopentane) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'O=C1C[CH][C]=CO1(27753)',
    structure = SMILES('O=C1CC=[C][CH]O1'),
    E0 = (38.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2422,0.00030625,9.60355e-05,-1.07419e-07,3.47913e-11,4645.62,15.4868], Tmin=(100,'K'), Tmax=(1086.12,'K')), NASAPolynomial(coeffs=[7.15805,0.0344329,-1.81432e-05,3.8186e-09,-2.85685e-13,931.484,-16.9108], Tmin=(1086.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cyclohexane) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=COC(=C)[O](27937)',
    structure = SMILES('[CH]=C=COC([CH2])=O'),
    E0 = (156.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,180,180,267.565,610.445,610.592,610.73],'cm^-1')),
        HinderedRotor(inertia=(0.0545268,'amu*angstrom^2'), symmetry=1, barrier=(14.4074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108481,'amu*angstrom^2'), symmetry=1, barrier=(28.6962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.46047,'amu*angstrom^2'), symmetry=1, barrier=(79.563,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.899354,0.0625266,-6.37084e-05,2.9838e-08,-4.58663e-12,18948.4,25.3793], Tmin=(100,'K'), Tmax=(977.681,'K')), NASAPolynomial(coeffs=[15.656,0.0130063,-4.38483e-06,7.41324e-10,-4.99926e-14,15544.2,-48.133], Tmin=(977.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CJCO)"""),
)

species(
    label = 'O(4)',
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
    label = 'C#C[CH]C[C]=O(27336)',
    structure = SMILES('[CH]=C=CC[C]=O'),
    E0 = (361.877,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,540,610,2055],'cm^-1')),
        HinderedRotor(inertia=(0.942019,'amu*angstrom^2'), symmetry=1, barrier=(21.6589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.942967,'amu*angstrom^2'), symmetry=1, barrier=(21.6807,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83042,0.0411588,-2.90855e-05,8.37138e-09,-5.06207e-13,43607.3,21.7665], Tmin=(100,'K'), Tmax=(1230.55,'K')), NASAPolynomial(coeffs=[12.1824,0.0151463,-6.68669e-06,1.2802e-09,-9.02369e-14,40481.4,-32.6727], Tmin=(1230.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C#CC1CC1([O])[O](27938)',
    structure = SMILES('C#CC1CC1([O])[O]'),
    E0 = (358.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78913,0.042521,-2.46961e-05,3.97885e-09,9.5633e-13,43216.2,21.4096], Tmin=(100,'K'), Tmax=(1122.09,'K')), NASAPolynomial(coeffs=[10.3022,0.0212593,-8.41904e-06,1.52409e-09,-1.04463e-13,40733.8,-23.1903], Tmin=(1122.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
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
    label = 'C#CC=CC([O])=O(27939)',
    structure = SMILES('C#CC=CC([O])=O'),
    E0 = (200.418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,750,770,3400,2100,180,180,180,180,342.569,1523.44],'cm^-1')),
        HinderedRotor(inertia=(0.0018168,'amu*angstrom^2'), symmetry=1, barrier=(2.97952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48377,'amu*angstrom^2'), symmetry=1, barrier=(34.1149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79352,0.0558365,-8.92923e-05,8.72628e-08,-3.39363e-11,24177.1,21.5848], Tmin=(100,'K'), Tmax=(792.612,'K')), NASAPolynomial(coeffs=[3.68808,0.031884,-1.67273e-05,3.32054e-09,-2.3441e-13,24328.8,15.7367], Tmin=(792.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cd-Cd(CO)H) + group(Cds-CdsCtH) + group(Cds-O2d(Cds-Cds)O2s) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CCOJ)"""),
)

species(
    label = 'C#CC[CH]C([O])=O(27940)',
    structure = SMILES('C#CCC=C([O])[O]'),
    E0 = (175.126,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,3010,987.5,1337.5,450,1655,350,440,435,1725,328.129,328.532,329.284],'cm^-1')),
        HinderedRotor(inertia=(0.30901,'amu*angstrom^2'), symmetry=1, barrier=(23.5647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310092,'amu*angstrom^2'), symmetry=1, barrier=(23.5682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22765,0.0578547,-5.86821e-05,3.03488e-08,-6.21543e-12,21165.2,22.4088], Tmin=(100,'K'), Tmax=(1185.92,'K')), NASAPolynomial(coeffs=[13.3449,0.016984,-6.98655e-06,1.28778e-09,-8.91006e-14,18291.2,-38.116], Tmin=(1185.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.126,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[C]#CCCC([O])=O(27941)',
    structure = SMILES('[C]#CCCC([O])=O'),
    E0 = (340.164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2175,525,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66353,0.0605995,-9.69397e-05,9.35472e-08,-3.45166e-11,40987.6,22.8011], Tmin=(100,'K'), Tmax=(869.601,'K')), NASAPolynomial(coeffs=[2.30221,0.0358339,-1.65691e-05,3.06723e-09,-2.06297e-13,41701.9,24.5545], Tmin=(869.601,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(340.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cds-OdCsOs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCOJ) + radical(Acetyl)"""),
)

species(
    label = '[C]#C[CH]CC(=O)O(27942)',
    structure = SMILES('[C]#C[CH]CC(=O)O'),
    E0 = (314.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2175,525,3615,1277.5,1000,3025,407.5,1350,352.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49729,0.0580984,-7.57984e-05,5.84835e-08,-1.82394e-11,37896.1,25.7279], Tmin=(100,'K'), Tmax=(875.226,'K')), NASAPolynomial(coeffs=[7.59044,0.025065,-1.02961e-05,1.81961e-09,-1.20045e-13,37028.2,-1.72116], Tmin=(875.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cds-OdCsOs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CCJCC=O)"""),
)

species(
    label = 'C#CC=CC(=O)O(26469)',
    structure = SMILES('C#CC=CC(=O)O'),
    E0 = (-25.2867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7061,0.0550211,-6.28539e-05,4.09845e-08,-1.12828e-11,-2962.62,21.7516], Tmin=(100,'K'), Tmax=(864.955,'K')), NASAPolynomial(coeffs=[7.91766,0.0262958,-1.30389e-05,2.58963e-09,-1.8547e-13,-4037.17,-7.31443], Tmin=(864.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.2867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cd-Cd(CO)H) + group(Cds-CdsCtH) + group(Cds-O2d(Cds-Cds)O2s) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC([CH2])C([O])=O(26471)',
    structure = SMILES('C#CC([CH2])C([O])=O'),
    E0 = (206.483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,180,1033.65,1036.93,1037.28,1037.57],'cm^-1')),
        HinderedRotor(inertia=(0.161355,'amu*angstrom^2'), symmetry=1, barrier=(3.70986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.46774,'amu*angstrom^2'), symmetry=1, barrier=(56.7382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.46282,'amu*angstrom^2'), symmetry=1, barrier=(56.625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27576,0.068391,-0.000110418,1.01886e-07,-3.64125e-11,24924.1,22.7699], Tmin=(100,'K'), Tmax=(856.44,'K')), NASAPolynomial(coeffs=[5.00114,0.0325915,-1.54906e-05,2.90735e-09,-1.97313e-13,24960.8,9.3141], Tmin=(856.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCOJ) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C#CC1CC(=O)O1(27741)',
    structure = SMILES('C#CC1CC(=O)O1'),
    E0 = (-63.2723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.323,0.0194823,6.11764e-05,-1.02923e-07,4.4702e-11,-7533.13,19.9777], Tmin=(100,'K'), Tmax=(895.745,'K')), NASAPolynomial(coeffs=[14.9013,0.0094636,6.70551e-07,-3.72268e-10,2.70386e-14,-11638,-49.6554], Tmin=(895.745,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.2723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsOs) + group(Ct-CtCs) + group(Ct-CtH) + ring(Beta-Propiolactone)"""),
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
    E0 = (161.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (283.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (516.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (161.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (321.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (294.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (218.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (206.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (523.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (483.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (512.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (324.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (244.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (179.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (470.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (604.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (358.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (417.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (309.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (281.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (487.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (388.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (240.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (329.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (169.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['CO2(13)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['[CH]=[C]C1CC(=O)O1(27931)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.448e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['[O]C1([O])C=C=CC1(27926)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.92551e+11,'s^-1'), n=0.201102, Ea=(354.987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R6;carbonylbond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CO2(13)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC([O])=O(26470)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.000178913,'m^3/(mol*s)'), n=2.74787, Ea=(113.483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CsJ-CdHH] for rate rule [CO2;CsJ-CdHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 108.4 to 113.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['[CH]=C=C[CH]C(=O)O(27932)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.64e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['C=C=[C]CC([O])=O(27933)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C=[C]CC(=O)O(27934)'],
    products = ['C#C[CH]CC([O])=O(26470)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SSS;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['C=[C]C=CC([O])=O(27916)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C3H2(81)', '[CH2]C([O])=O(1440)'],
    products = ['C#C[CH]CC([O])=O(26470)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.28206e+10,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_allenic;C_pri_rad] for rate rule [Cd_allenic;C_rad/H2/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O][C]=O(669)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC([O])=O(26470)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.56662e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H2/Cd;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['[CH]=C=CC[C]1OO1(27935)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(350.431,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['[O]C(=O)CC1[C]=C1(27936)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.2354e+12,'s^-1'), n=-0.1205, Ea=(162.597,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cdsingleH] + [R3;doublebond_intra_CdCdd;radadd_intra_cdsingle] for rate rule [R3;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['[CH]=C1[CH]CC(=O)O1(27767)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.33332e+10,'s^-1'), n=0.302034, Ea=(82.5645,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['O=C1C[CH][C]=CO1(27753)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.18439e+10,'s^-1'), n=0.253963, Ea=(17.5802,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R6_linear;carbonyl_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C=COC(=C)[O](27937)'],
    products = ['C#C[CH]CC([O])=O(26470)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O(4)', 'C#C[CH]C[C]=O(27336)'],
    products = ['C#C[CH]CC([O])=O(26470)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['C#CC1CC1([O])[O](27938)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.90568e+10,'s^-1'), n=0.237, Ea=(196.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHDe] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_csHCt]
Euclidian distance = 1.73205080757
family: Intra_R_Add_Exocyclic
Ea raised from 196.3 to 196.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', 'C#CC=CC([O])=O(27939)'],
    products = ['C#C[CH]CC([O])=O(26470)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.165e+08,'cm^3/(mol*s)'), n=1.533, Ea=(5.18398,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 2853 used for Cds-COH_Cds;HJ
Exact match found for rate rule [Cds-COH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O][C]=O(669)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CC([O])=O(26470)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0132398,'m^3/(mol*s)'), n=2.333, Ea=(2.89605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CtH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#CC[CH]C([O])=O(27940)'],
    products = ['C#C[CH]CC([O])=O(26470)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.12976e+09,'s^-1'), n=1.23766, Ea=(105.9,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/OneDe;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/OneDe;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[C]#CCCC([O])=O(27941)'],
    products = ['C#C[CH]CC([O])=O(26470)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[C]#C[CH]CC(=O)O(27942)'],
    products = ['C#C[CH]CC([O])=O(26470)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(19101.7,'s^-1'), n=1.79759, Ea=(74.4436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;XH_out] for rate rule [R6HJ_2;Ct_rad_out;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['C#CC=CC(=O)O(26469)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.01596e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#CC([CH2])C([O])=O(26471)'],
    products = ['C#C[CH]CC([O])=O(26470)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO] for rate rule [cCs(-HC)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C#C[CH]CC([O])=O(26470)'],
    products = ['C#CC1CC(=O)O1(27741)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.6e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

network(
    label = '4436',
    isomers = [
        'C#C[CH]CC([O])=O(26470)',
    ],
    reactants = [
        ('CO2(13)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4436',
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

