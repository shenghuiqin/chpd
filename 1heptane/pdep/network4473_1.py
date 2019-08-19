species(
    label = 'C#CC([CH2])O[C]=C(26505)',
    structure = SMILES('C#CC([CH2])O[C]=C'),
    E0 = (497.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2175,525,2950,3100,1380,975,1025,1650,750,770,3400,2100,3000,3100,440,815,1455,1000,1685,370,181.154,181.306,181.526],'cm^-1')),
        HinderedRotor(inertia=(0.998118,'amu*angstrom^2'), symmetry=1, barrier=(23.3231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.995985,'amu*angstrom^2'), symmetry=1, barrier=(23.3221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.001,'amu*angstrom^2'), symmetry=1, barrier=(23.3231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.999733,'amu*angstrom^2'), symmetry=1, barrier=(23.3219,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.193817,0.0790118,-8.82202e-05,4.75864e-08,-9.64007e-12,60050.7,29.7585], Tmin=(100,'K'), Tmax=(1356.94,'K')), NASAPolynomial(coeffs=[20.64,0.00884537,-9.812e-07,-2.76931e-11,8.00066e-15,55202.4,-74.1425], Tmin=(1356.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO) + radical(CJC(C)OC)"""),
)

species(
    label = 'CH2CO(28)',
    structure = SMILES('C=C=O'),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH]=C1CC1O[C]=C(27425)',
    structure = SMILES('[CH]=C1CC1O[C]=C'),
    E0 = (572.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.988261,0.0483532,7.22651e-06,-5.66946e-08,2.86108e-11,68925.6,25.9577], Tmin=(100,'K'), Tmax=(940.895,'K')), NASAPolynomial(coeffs=[20.8993,0.00871813,-1.34566e-06,2.23995e-10,-2.24833e-14,63186.3,-79.4763], Tmin=(940.895,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C(=C)OC1[CH2](27426)',
    structure = SMILES('[CH]=C1C(=C)OC1[CH2]'),
    E0 = (458.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.703617,0.0580206,-2.11626e-05,-2.59345e-08,1.72753e-11,55248.4,21.0219], Tmin=(100,'K'), Tmax=(949.636,'K')), NASAPolynomial(coeffs=[20.3208,0.0108753,-2.7446e-06,4.84311e-10,-3.86447e-14,49922.6,-81.0303], Tmin=(949.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(Cds_P) + radical(CJC(C)OC)"""),
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
    label = 'C#CC(=C)O[C]=C(27427)',
    structure = SMILES('C#CC(=C)O[C]=C'),
    E0 = (433.016,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2175,525,1685,370,750,770,3400,2100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,285.42,285.42,285.42],'cm^-1')),
        HinderedRotor(inertia=(0.426972,'amu*angstrom^2'), symmetry=1, barrier=(24.6829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.426973,'amu*angstrom^2'), symmetry=1, barrier=(24.6829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.426973,'amu*angstrom^2'), symmetry=1, barrier=(24.6829,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32435,0.0593907,-5.91578e-05,3.1416e-08,-6.80739e-12,52175.7,24.6897], Tmin=(100,'K'), Tmax=(1103.42,'K')), NASAPolynomial(coeffs=[11.2092,0.0235572,-1.04453e-05,1.98466e-09,-1.39167e-13,49994.3,-23.9717], Tmin=(1103.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJO)"""),
)

species(
    label = 'C#COC([CH2])C#C(27428)',
    structure = SMILES('C#COC([CH2])C#C'),
    E0 = (467.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.33385,'amu*angstrom^2'), symmetry=1, barrier=(30.6679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33673,'amu*angstrom^2'), symmetry=1, barrier=(30.734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34062,'amu*angstrom^2'), symmetry=1, barrier=(30.8234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33884,'amu*angstrom^2'), symmetry=1, barrier=(30.7825,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.438094,0.0833262,-9.90894e-05,5.4841e-08,-1.12547e-11,56402,24.9615], Tmin=(100,'K'), Tmax=(1354.7,'K')), NASAPolynomial(coeffs=[22.8673,0.00365289,1.15395e-06,-4.07658e-10,3.31158e-14,51084.2,-90.8707], Tmin=(1354.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtOs) + group(Ct-CtH) + group(Ct-CtH) + radical(CJC(C)OC)"""),
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
    label = 'C2H(33)',
    structure = SMILES('[C]#C'),
    E0 = (557.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89868,0.0132988,-2.80733e-05,2.89485e-08,-1.07502e-11,67061.6,6.18548], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.6627,0.00382492,-1.36633e-06,2.13455e-10,-1.23217e-14,67168.4,3.92206], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(557.301,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""C2H""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=[C]OC=C(2409)',
    structure = SMILES('C=[C]OC=C'),
    E0 = (209.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1685,370,343.748,343.774,344.322],'cm^-1')),
        HinderedRotor(inertia=(0.221848,'amu*angstrom^2'), symmetry=1, barrier=(18.6172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221784,'amu*angstrom^2'), symmetry=1, barrier=(18.6094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21804,0.0314889,-3.55286e-06,-1.9957e-08,1.05742e-11,25275.1,20.564], Tmin=(100,'K'), Tmax=(960.345,'K')), NASAPolynomial(coeffs=[10.9617,0.0135903,-4.52395e-06,7.9853e-10,-5.66264e-14,22741.7,-25.7118], Tmin=(960.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'C#C[C](C)O[C]=C(27429)',
    structure = SMILES('[CH]=C=C(C)O[C]=C'),
    E0 = (462.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17385,0.0642544,-7.22695e-05,4.68064e-08,-1.24742e-11,55764.5,27.8095], Tmin=(100,'K'), Tmax=(906.575,'K')), NASAPolynomial(coeffs=[9.50273,0.0275043,-1.14614e-05,2.08843e-09,-1.42188e-13,54254.4,-11.5553], Tmin=(906.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=COC([CH2])C#C(27430)',
    structure = SMILES('[CH]=COC([CH2])C#C'),
    E0 = (505.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.18862,'amu*angstrom^2'), symmetry=1, barrier=(27.3288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18826,'amu*angstrom^2'), symmetry=1, barrier=(27.3205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19023,'amu*angstrom^2'), symmetry=1, barrier=(27.3656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1888,'amu*angstrom^2'), symmetry=1, barrier=(27.3329,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.949101,0.0862643,-9.69121e-05,5.04515e-08,-9.61423e-12,60970.2,29.749], Tmin=(100,'K'), Tmax=(1513.44,'K')), NASAPolynomial(coeffs=[24.196,0.0029093,2.44956e-06,-6.93986e-10,5.29029e-14,55294.1,-95.5885], Tmin=(1513.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(505.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(CJC(C)OC)"""),
)

species(
    label = 'C#C[C]([CH2])OC=C(27431)',
    structure = SMILES('[CH]=C=C([CH2])OC=C'),
    E0 = (381.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.392465,0.0694589,-7.1504e-05,3.75166e-08,-7.58318e-12,46081.6,28.1668], Tmin=(100,'K'), Tmax=(1308.22,'K')), NASAPolynomial(coeffs=[16.8521,0.0148006,-3.86664e-06,5.17785e-10,-2.90762e-14,42145.7,-54.247], Tmin=(1308.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(O)CJ) + radical(C=C=CJ)"""),
)

species(
    label = '[C]#CC(C)O[C]=C(27432)',
    structure = SMILES('[C]#CC(C)O[C]=C'),
    E0 = (624.573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2175,525,1685,370,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,334.088,334.153,334.154,334.413],'cm^-1')),
        HinderedRotor(inertia=(0.239082,'amu*angstrom^2'), symmetry=1, barrier=(18.9341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238748,'amu*angstrom^2'), symmetry=1, barrier=(18.9308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23916,'amu*angstrom^2'), symmetry=1, barrier=(18.9274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05292,'amu*angstrom^2'), symmetry=1, barrier=(83.3987,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0271492,0.0715356,-7.31223e-05,3.69674e-08,-7.02437e-12,75275.4,29.9601], Tmin=(100,'K'), Tmax=(1478.95,'K')), NASAPolynomial(coeffs=[18.5956,0.0111341,-1.53535e-06,4.35318e-11,3.99985e-15,70896.5,-63.1239], Tmin=(1478.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(624.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO) + radical(Acetyl)"""),
)

species(
    label = '[CH]=[C]OC(C)C#C(27433)',
    structure = SMILES('[CH]=[C]OC(C)C#C'),
    E0 = (534.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06034,'amu*angstrom^2'), symmetry=1, barrier=(24.3794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06206,'amu*angstrom^2'), symmetry=1, barrier=(24.4189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06105,'amu*angstrom^2'), symmetry=1, barrier=(24.3956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06104,'amu*angstrom^2'), symmetry=1, barrier=(24.3955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.28091,0.0738561,-7.47253e-05,3.63391e-08,-6.57873e-12,64460.2,30.7372], Tmin=(100,'K'), Tmax=(1570.71,'K')), NASAPolynomial(coeffs=[20.5502,0.00853357,-6.22325e-07,-8.77828e-11,1.09183e-14,59430.3,-74.3469], Tmin=(1570.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[C]#CC([CH2])OC=C(27434)',
    structure = SMILES('[C]#CC([CH2])OC=C'),
    E0 = (595.339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2175,525,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.09868,'amu*angstrom^2'), symmetry=1, barrier=(25.2609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09809,'amu*angstrom^2'), symmetry=1, barrier=(25.2472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09786,'amu*angstrom^2'), symmetry=1, barrier=(25.242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09535,'amu*angstrom^2'), symmetry=1, barrier=(25.1843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.660415,0.0841386,-9.58533e-05,5.16244e-08,-1.02354e-11,71786.4,29.044], Tmin=(100,'K'), Tmax=(1445.79,'K')), NASAPolynomial(coeffs=[22.3199,0.00537872,1.61021e-06,-5.79629e-10,4.73584e-14,66728.1,-84.8071], Tmin=(1445.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(595.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=[C][O](173)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97566e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02476e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CsCJ=O) + radical(CJC=O)"""),
)

species(
    label = 'C=[C]OC1[C]=CC1(27423)',
    structure = SMILES('C=[C]OC1[C]=CC1'),
    E0 = (534.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37778,0.0385401,2.9688e-05,-7.55469e-08,3.40957e-11,64452.1,25.7541], Tmin=(100,'K'), Tmax=(947.707,'K')), NASAPolynomial(coeffs=[19.2898,0.0108953,-2.46131e-06,4.63802e-10,-4.08614e-14,58903.4,-71.0613], Tmin=(947.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1[C]=CC(=C)O1(27435)',
    structure = SMILES('[CH2]C1[C]=CC(=C)O1'),
    E0 = (361.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28002,0.0361555,4.92104e-05,-1.06025e-07,4.79303e-11,43638.9,20.5535], Tmin=(100,'K'), Tmax=(924.257,'K')), NASAPolynomial(coeffs=[22.8162,0.00407634,2.07141e-06,-4.69691e-10,2.45591e-14,37047.1,-95.7734], Tmin=(924.257,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CJC(C)OC) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(=C)OC=C(27436)',
    structure = SMILES('C#CC(=C)OC=C'),
    E0 = (193.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80412,0.0613455,-4.55098e-05,1.01724e-08,1.75293e-12,23368.1,23.4242], Tmin=(100,'K'), Tmax=(1017.68,'K')), NASAPolynomial(coeffs=[16.1308,0.0181054,-6.8354e-06,1.25308e-09,-8.85507e-14,19368.1,-55.1123], Tmin=(1017.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC[CH]O[C]=C(27437)',
    structure = SMILES('C#CC[CH]O[C]=C'),
    E0 = (494.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.046515,0.0792083,-9.13225e-05,5.15695e-08,-1.10296e-11,59627.5,28.6116], Tmin=(100,'K'), Tmax=(1256.44,'K')), NASAPolynomial(coeffs=[19.5014,0.010665,-1.95806e-06,1.55157e-10,-4.02288e-15,55213.4,-68.176], Tmin=(1256.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO) + radical(CCsJOC(O))"""),
)

species(
    label = 'C#CC1CC(=C)O1(27438)',
    structure = SMILES('C#CC1CC(=C)O1'),
    E0 = (177.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35964,0.0345112,5.61878e-05,-1.18712e-07,5.52381e-11,21520.4,17.5402], Tmin=(100,'K'), Tmax=(893.186,'K')), NASAPolynomial(coeffs=[23.0274,0.00102831,5.68911e-06,-1.35851e-09,9.40966e-14,15114.7,-98.7376], Tmin=(893.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(2methyleneoxetane)"""),
)

species(
    label = 'H2CC(41)',
    structure = SMILES('[C]=C'),
    E0 = (401.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2480.69,'J/mol'), sigma=(4.48499,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=387.48 K, Pc=62.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28155,0.00697643,-2.38528e-06,-1.21078e-09,9.82042e-13,48319.2,5.92036], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.27807,0.00475623,-1.63007e-06,2.54623e-10,-1.4886e-14,48014,0.639979], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(401.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C#CC([CH2])[O](26706)',
    structure = SMILES('C#CC([CH2])[O]'),
    E0 = (418.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,348.094],'cm^-1')),
        HinderedRotor(inertia=(0.289955,'amu*angstrom^2'), symmetry=1, barrier=(24.9403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.397163,'amu*angstrom^2'), symmetry=1, barrier=(34.191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46647,0.0444261,-4.38277e-05,2.13721e-08,-3.91288e-12,50443.4,20.5021], Tmin=(100,'K'), Tmax=(1545.21,'K')), NASAPolynomial(coeffs=[13.2214,0.00716331,-1.02144e-06,4.17503e-11,1.21996e-15,47626.5,-38.6839], Tmin=(1545.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CJCO)"""),
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
    label = '[CH]=C=CO[C]=C(27424)',
    structure = SMILES('[CH]=C=CO[C]=C'),
    E0 = (504.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21848,'amu*angstrom^2'), symmetry=1, barrier=(28.0152,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21747,'amu*angstrom^2'), symmetry=1, barrier=(27.992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3446,0.051008,-5.2679e-05,2.77941e-08,-5.63681e-12,60792.9,25.1839], Tmin=(100,'K'), Tmax=(1318.65,'K')), NASAPolynomial(coeffs=[13.3315,0.0108537,-2.68755e-06,3.38636e-10,-1.80053e-14,57961.4,-34.7112], Tmin=(1318.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO)"""),
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
    E0 = (497.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (572.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (578.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (654.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (694.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (497.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (787.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (605.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (610.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (639.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (720.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (567.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (639.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (611.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (628.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (556.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (586.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (656.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (506.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (819.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (886.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([CH2])O[C]=C(26505)'],
    products = ['CH2CO(28)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([CH2])O[C]=C(26505)'],
    products = ['[CH]=C1CC1O[C]=C(27425)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(74.1021,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for R4_S_T;triplebond_intra_H;radadd_intra_cs2H
Exact match found for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 70.8 to 74.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([CH2])O[C]=C(26505)'],
    products = ['[CH]=C1C(=C)OC1[CH2](27426)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.98674e+07,'s^-1'), n=1.31443, Ea=(80.4773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#CC(=C)O[C]=C(27427)'],
    products = ['C#CC([CH2])O[C]=C(26505)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(87.5179,'m^3/(mol*s)'), n=1.66467, Ea=(9.4546,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;HJ] for rate rule [Cds-CtOs_Cds;HJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#COC([CH2])C#C(27428)'],
    products = ['C#CC([CH2])O[C]=C(26505)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2CO(28)', '[CH]=[C]C=C(4699)'],
    products = ['C#CC([CH2])O[C]=C(26505)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(11.6997,'m^3/(mol*s)'), n=2.021, Ea=(107.174,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd;CJ]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond
Ea raised from 104.3 to 107.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H(33)', 'C=[C]OC=C(2409)'],
    products = ['C#CC([CH2])O[C]=C(26505)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0018005,'m^3/(mol*s)'), n=2.50408, Ea=(20.5143,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CJ] for rate rule [Cds-OsH_Cds;CtJ_Ct]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C#CC([CH2])O[C]=C(26505)'],
    products = ['C#C[C](C)O[C]=C(27429)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.89e+09,'s^-1'), n=0.81, Ea=(107.529,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out_OneDe] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Ct]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=COC([CH2])C#C(27430)'],
    products = ['C#CC([CH2])O[C]=C(26505)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#CC([CH2])O[C]=C(26505)'],
    products = ['C#C[C]([CH2])OC=C(27431)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[C]#CC(C)O[C]=C(27432)'],
    products = ['C#CC([CH2])O[C]=C(26505)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.39293e+07,'s^-1'), n=1.32074, Ea=(96.0416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_TSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]OC(C)C#C(27433)'],
    products = ['C#CC([CH2])O[C]=C(26505)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]#CC([CH2])OC=C(27434)'],
    products = ['C#CC([CH2])O[C]=C(26505)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_TSSS;Ct_rad_out;Cd_H_out_doubleC]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C][O](173)', '[CH]=[C]C=C(4699)'],
    products = ['C#CC([CH2])O[C]=C(26505)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#CC([CH2])O[C]=C(26505)'],
    products = ['C=[C]OC1[C]=CC1(27423)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(130.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#CC([CH2])O[C]=C(26505)'],
    products = ['[CH2]C1[C]=CC(=C)O1(27435)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 13 used for R5_SS_T;triplebond_intra_H;radadd_intra_cddouble
Exact match found for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#CC([CH2])O[C]=C(26505)'],
    products = ['C#CC(=C)OC=C(27436)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#CC([CH2])O[C]=C(26505)'],
    products = ['C#CC[CH]O[C]=C(27437)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#CC([CH2])O[C]=C(26505)'],
    products = ['C#CC1CC(=C)O1(27438)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H2CC(41)', 'C#CC([CH2])[O](26706)'],
    products = ['C#CC([CH2])O[C]=C(26505)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(19)', '[CH]=C=CO[C]=C(27424)'],
    products = ['C#CC([CH2])O[C]=C(26505)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4473',
    isomers = [
        'C#CC([CH2])O[C]=C(26505)',
    ],
    reactants = [
        ('CH2CO(28)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4473',
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

