species(
    label = 'C#C[CH]CC([CH2])C=C(26568)',
    structure = SMILES('C#C[CH]CC([CH2])C=C'),
    E0 = (525.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,750,770,3400,2100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.190021,0.0815974,-7.22381e-05,3.55589e-08,-7.00565e-12,63370.5,33.7853], Tmin=(100,'K'), Tmax=(1298.25,'K')), NASAPolynomial(coeffs=[15.6646,0.0298069,-9.00093e-06,1.34069e-09,-8.03084e-14,59501.7,-45.8878], Tmin=(1298.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'butadiene13(2459)',
    structure = SMILES('C=CC=C'),
    E0 = (96.4553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.30712,'amu*angstrom^2'), symmetry=1, barrier=(30.0532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80599,0.0102584,6.1726e-05,-9.01643e-08,3.59117e-11,11658.5,12.0621], Tmin=(100,'K'), Tmax=(946.047,'K')), NASAPolynomial(coeffs=[12.4694,0.0100554,-2.41207e-06,4.57077e-10,-3.93161e-14,8010.78,-43.6375], Tmin=(946.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""butadiene13""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2][CH]C=C(2458)',
    structure = SMILES('[CH2]C=C[CH2]'),
    E0 = (274.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.210311,'amu*angstrom^2'), symmetry=1, barrier=(25.2351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779031,'amu*angstrom^2'), symmetry=1, barrier=(93.4717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56318,0.0223429,1.87067e-05,-3.93099e-08,1.63982e-11,33100.5,13.4097], Tmin=(100,'K'), Tmax=(974.264,'K')), NASAPolynomial(coeffs=[9.82995,0.0151966,-5.22272e-06,9.67656e-10,-7.0786e-14,30607.7,-26.985], Tmin=(974.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C#C[CH]CC1CC1[CH2](27573)',
    structure = SMILES('C#C[CH]CC1CC1[CH2]'),
    E0 = (550.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.604087,0.0604483,5.77225e-07,-5.07869e-08,2.69588e-11,66317.3,29.0102], Tmin=(100,'K'), Tmax=(907.227,'K')), NASAPolynomial(coeffs=[16.8419,0.0269331,-6.96657e-06,1.02027e-09,-6.61605e-14,61804,-56.383], Tmin=(907.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(550.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(Sec_Propargyl) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC1CC([CH2])C1[CH2](27300)',
    structure = SMILES('C#CC1CC([CH2])C1[CH2]'),
    E0 = (595.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.915917,0.0498278,3.27845e-05,-8.42993e-08,3.8788e-11,71761.2,28.8754], Tmin=(100,'K'), Tmax=(912.609,'K')), NASAPolynomial(coeffs=[16.9832,0.0265389,-6.40905e-06,9.25698e-10,-6.18216e-14,66865.8,-57.9243], Tmin=(912.609,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(595.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C1[CH]CC(C=C)C1(27574)',
    structure = SMILES('[CH]C1=CCC(C=C)C1'),
    E0 = (425.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03756,0.0433507,5.89159e-05,-1.03232e-07,4.12797e-11,51347.3,27.1886], Tmin=(100,'K'), Tmax=(976.749,'K')), NASAPolynomial(coeffs=[15.6757,0.0370311,-1.37343e-05,2.5653e-09,-1.86325e-13,45929.6,-56.1828], Tmin=(976.749,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopentene) + radical(AllylJ2_triplet)"""),
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
    label = 'C#C[CH]CC(=C)C=C(27575)',
    structure = SMILES('C#C[CH]CC(=C)C=C'),
    E0 = (426.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,252.368,253.651,254.43],'cm^-1')),
        HinderedRotor(inertia=(0.379159,'amu*angstrom^2'), symmetry=1, barrier=(17.2126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.60486,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.63777,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.65163,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.560284,0.08356,-7.60547e-05,3.61026e-08,-6.63641e-12,51440.9,30.5272], Tmin=(100,'K'), Tmax=(1473.29,'K')), NASAPolynomial(coeffs=[19.2502,0.0221348,-5.73796e-06,7.64566e-10,-4.27438e-14,46432.7,-69.9099], Tmin=(1473.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC=CC([CH2])C=C(27576)',
    structure = SMILES('C#CC=CC([CH2])C=C'),
    E0 = (488.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,750,770,3400,2100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.807602,'amu*angstrom^2'), symmetry=1, barrier=(18.5684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.806442,'amu*angstrom^2'), symmetry=1, barrier=(18.5417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.169913,0.0781783,-6.75053e-05,3.11868e-08,-5.84533e-12,58934.4,29.4699], Tmin=(100,'K'), Tmax=(1274.87,'K')), NASAPolynomial(coeffs=[15.239,0.0308981,-1.18762e-05,2.09694e-09,-1.40879e-13,55092.2,-46.8894], Tmin=(1274.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Isobutyl)"""),
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
    label = 'C#C[CH]CC=C(26789)',
    structure = SMILES('C#C[CH]CC=C'),
    E0 = (376.014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,410.876,410.906],'cm^-1')),
        HinderedRotor(inertia=(0.0812548,'amu*angstrom^2'), symmetry=1, barrier=(9.73456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.63521,'amu*angstrom^2'), symmetry=1, barrier=(76.1042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.635141,'amu*angstrom^2'), symmetry=1, barrier=(76.1043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59773,0.0469098,-2.16011e-05,-7.36055e-09,7.51589e-12,45316,21.3454], Tmin=(100,'K'), Tmax=(911.508,'K')), NASAPolynomial(coeffs=[10.894,0.0220604,-6.94892e-06,1.11507e-09,-7.2543e-14,42958.8,-26.2758], Tmin=(911.508,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH]C[C](C)C=C(27577)',
    structure = SMILES('C#C[CH]C[C](C)C=C'),
    E0 = (452.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,360,370,350,2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.455224,0.069971,-3.84041e-05,-2.18379e-10,5.97447e-12,54553.4,27.9418], Tmin=(100,'K'), Tmax=(947.544,'K')), NASAPolynomial(coeffs=[13.3576,0.0343147,-1.17365e-05,1.96979e-09,-1.30516e-13,51263.9,-38.0656], Tmin=(947.544,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_T) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC[CH]C([CH2])C=C(27319)',
    structure = SMILES('C#CC[CH]C([CH2])C=C'),
    E0 = (573.895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,750,770,3400,2100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.320735,0.0757652,-6.17979e-05,2.83772e-08,-5.40345e-12,69160.4,34.6709], Tmin=(100,'K'), Tmax=(1243.11,'K')), NASAPolynomial(coeffs=[12.8711,0.0353817,-1.30694e-05,2.24475e-09,-1.48018e-13,66040.1,-28.6085], Tmin=(1243.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cs_S)"""),
)

species(
    label = 'C#C[CH][CH]C(C)C=C(27578)',
    structure = SMILES('[CH]=C=C[CH]C(C)C=C'),
    E0 = (468.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.390564,0.0679548,-2.85324e-05,-1.07414e-08,8.95187e-12,56441.3,30.375], Tmin=(100,'K'), Tmax=(1006.4,'K')), NASAPolynomial(coeffs=[15.2161,0.0330408,-1.22826e-05,2.20148e-09,-1.52442e-13,52241.2,-47.2857], Tmin=(1006.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Allyl_S)"""),
)

species(
    label = 'C#C[CH]CC(C)[C]=C(26580)',
    structure = SMILES('C#C[CH]CC(C)[C]=C'),
    E0 = (558.323,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0763857,0.0807808,-7.15244e-05,3.5528e-08,-7.21357e-12,67296.6,31.5209], Tmin=(100,'K'), Tmax=(1182.29,'K')), NASAPolynomial(coeffs=[13.9108,0.0339756,-1.21418e-05,2.04361e-09,-1.33186e-13,64025.4,-37.5386], Tmin=(1182.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CCC[C]([CH2])C=C(27579)',
    structure = SMILES('C#CCCC([CH2])=C[CH2]'),
    E0 = (457.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0021803,0.0779156,-5.14823e-05,9.20494e-09,2.89897e-12,55203.8,30.3564], Tmin=(100,'K'), Tmax=(1003.52,'K')), NASAPolynomial(coeffs=[16.5062,0.031646,-1.14907e-05,2.01525e-09,-1.37362e-13,50908.8,-54.2198], Tmin=(1003.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[C]#CCCC([CH2])C=C(27580)',
    structure = SMILES('[C]#CCCC([CH2])C=C'),
    E0 = (716.497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.112487,0.0810388,-7.49825e-05,4.00472e-08,-8.75619e-12,86318.4,33.6023], Tmin=(100,'K'), Tmax=(1103.45,'K')), NASAPolynomial(coeffs=[12.9111,0.0346448,-1.19169e-05,1.94582e-09,-1.24027e-13,83493.9,-29.4037], Tmin=(1103.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(716.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Acetyl)"""),
)

species(
    label = 'C#CCCC([CH2])[C]=C(27581)',
    structure = SMILES('C#CCCC([CH2])[C]=C'),
    E0 = (617.195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.122248,0.0813019,-7.43204e-05,3.85423e-08,-8.19287e-12,74374.4,33.3001], Tmin=(100,'K'), Tmax=(1130.2,'K')), NASAPolynomial(coeffs=[13.2613,0.0347998,-1.26022e-05,2.13636e-09,-1.39802e-13,71404.5,-31.6961], Tmin=(1130.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C)C[CH]C#C(27582)',
    structure = SMILES('[CH]=CC(C)C[CH]C#C'),
    E0 = (567.577,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,280.135,3724.29],'cm^-1')),
        HinderedRotor(inertia=(0.294469,'amu*angstrom^2'), symmetry=1, barrier=(16.3984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294469,'amu*angstrom^2'), symmetry=1, barrier=(16.3984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49856,'amu*angstrom^2'), symmetry=1, barrier=(83.4521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49857,'amu*angstrom^2'), symmetry=1, barrier=(83.4522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49857,'amu*angstrom^2'), symmetry=1, barrier=(83.4521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.17803,0.0821791,-7.18427e-05,3.4301e-08,-6.58623e-12,68422.4,32.366], Tmin=(100,'K'), Tmax=(1260.15,'K')), NASAPolynomial(coeffs=[16.1732,0.0302777,-1.00639e-05,1.61845e-09,-1.02503e-13,64301.3,-50.3006], Tmin=(1260.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])CCC#C(27583)',
    structure = SMILES('[CH]=CC([CH2])CCC#C'),
    E0 = (626.449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2175,525,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,750,770,3400,2100,318.944,318.945],'cm^-1')),
        HinderedRotor(inertia=(0.00128687,'amu*angstrom^2'), symmetry=1, barrier=(13.1149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181681,'amu*angstrom^2'), symmetry=1, barrier=(13.1149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17532,'amu*angstrom^2'), symmetry=1, barrier=(84.8423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17532,'amu*angstrom^2'), symmetry=1, barrier=(84.8422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17531,'amu*angstrom^2'), symmetry=1, barrier=(84.8423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.107654,0.0824421,-7.38667e-05,3.64754e-08,-7.26685e-12,75499,34.055], Tmin=(100,'K'), Tmax=(1213.9,'K')), NASAPolynomial(coeffs=[15.4308,0.0312394,-1.05953e-05,1.72657e-09,-1.10305e-13,71726.6,-43.921], Tmin=(1213.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[C]#C[CH]CC(C)C=C(27584)',
    structure = SMILES('[C]#C[CH]CC(C)C=C'),
    E0 = (657.625,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0710126,0.0804653,-7.20043e-05,3.68054e-08,-7.68601e-12,79240.5,31.8075], Tmin=(100,'K'), Tmax=(1154.28,'K')), NASAPolynomial(coeffs=[13.5043,0.0339134,-1.15089e-05,1.86525e-09,-1.1841e-13,76139.4,-34.9277], Tmin=(1154.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(657.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Acetyl)"""),
)

species(
    label = '[CH2]C(C=C)CC1[C]=C1(27585)',
    structure = SMILES('[CH2]C(C=C)CC1[C]=C1'),
    E0 = (699.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.213602,0.0785393,-6.52809e-05,2.97307e-08,-5.57663e-12,84219.1,30.7005], Tmin=(100,'K'), Tmax=(1261.88,'K')), NASAPolynomial(coeffs=[13.946,0.0350095,-1.35373e-05,2.39416e-09,-1.60827e-13,80753.4,-38.7449], Tmin=(1261.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(699.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(cyclopropenyl-vinyl) + radical(Isobutyl)"""),
)

species(
    label = 'C#C[CH]CC1[CH]CC1(27586)',
    structure = SMILES('C#C[CH]CC1[CH]CC1'),
    E0 = (537.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.955874,0.0537302,7.63044e-06,-4.53788e-08,2.08869e-11,64785.2,29.2975], Tmin=(100,'K'), Tmax=(966.75,'K')), NASAPolynomial(coeffs=[13.1312,0.0349545,-1.22687e-05,2.1554e-09,-1.49063e-13,60954.4,-36.6673], Tmin=(966.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane) + radical(Sec_Propargyl) + radical(cyclobutane)"""),
)

species(
    label = 'C#CC1C[CH]C([CH2])C1(27587)',
    structure = SMILES('C#CC1C[CH]C([CH2])C1'),
    E0 = (512.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38503,0.0387221,5.43728e-05,-9.83723e-08,4.12553e-11,61730.3,29.1296], Tmin=(100,'K'), Tmax=(934.046,'K')), NASAPolynomial(coeffs=[14.4524,0.030456,-8.94507e-06,1.48752e-09,-1.04318e-13,57208.7,-44.1583], Tmin=(934.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(512.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopentane) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC1C[CH][C]=CC1(27588)',
    structure = SMILES('C=CC1C[CH][C]=CC1'),
    E0 = (427.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46503,0.0333903,7.27858e-05,-1.13609e-07,4.42179e-11,51518.6,23.0074], Tmin=(100,'K'), Tmax=(979.574,'K')), NASAPolynomial(coeffs=[15.062,0.0334784,-1.25035e-05,2.38982e-09,-1.76867e-13,46186.7,-55.9282], Tmin=(979.574,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclohexene) + radical(cyclohexene-allyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CCCC(=C)C=C(27589)',
    structure = SMILES('C#CCCC(=C)C=C'),
    E0 = (280.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0606653,0.0755317,-4.29141e-05,-3.01587e-09,8.32733e-12,33829.3,28.4525], Tmin=(100,'K'), Tmax=(960.435,'K')), NASAPolynomial(coeffs=[17.4693,0.0286318,-9.653e-06,1.65262e-09,-1.12719e-13,29304.4,-60.9788], Tmin=(960.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC=CC(C)C=C(27590)',
    structure = SMILES('C#CC=CC(C)C=C'),
    E0 = (283.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0054291,0.0780359,-5.95267e-05,2.33134e-08,-3.68736e-12,34277.7,28.2905], Tmin=(100,'K'), Tmax=(1492.3,'K')), NASAPolynomial(coeffs=[17.5687,0.0309584,-1.22059e-05,2.1732e-09,-1.45787e-13,29035.8,-63.4731], Tmin=(1492.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#C[CH]CC[CH]C=C(26564)',
    structure = SMILES('C#C[CH]CCC=C[CH2]'),
    E0 = (467.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.158128,0.0765436,-5.38631e-05,1.44172e-08,7.8694e-13,56397.5,30.8415], Tmin=(100,'K'), Tmax=(986.639,'K')), NASAPolynomial(coeffs=[14.6957,0.0332245,-1.175e-05,2.00637e-09,-1.33839e-13,52768.6,-42.951], Tmin=(986.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Allyl_P)"""),
)

species(
    label = 'C#C[CH]C[CH]CC=C(27591)',
    structure = SMILES('C#C[CH]C[CH]CC=C'),
    E0 = (522.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615433,0.0762222,-7.04614e-05,4.10912e-08,-1.0382e-11,63012,31.9304], Tmin=(100,'K'), Tmax=(934.4,'K')), NASAPolynomial(coeffs=[8.37147,0.0430184,-1.71566e-05,3.05792e-09,-2.05641e-13,61562.6,-4.96139], Tmin=(934.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(RCCJCC)"""),
)

species(
    label = 'C#CC([CH2])C([CH2])C=C(26569)',
    structure = SMILES('C#CC([CH2])C([CH2])C=C'),
    E0 = (575.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,750,770,3400,2100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,335.718,335.875],'cm^-1')),
        HinderedRotor(inertia=(3.70374e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.904772,'amu*angstrom^2'), symmetry=1, barrier=(72.8523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.906402,'amu*angstrom^2'), symmetry=1, barrier=(72.8519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128953,'amu*angstrom^2'), symmetry=1, barrier=(10.327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.911175,'amu*angstrom^2'), symmetry=1, barrier=(72.8883,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3560.51,'J/mol'), sigma=(6.30174,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=556.14 K, Pc=32.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.297611,0.0815261,-7.15621e-05,3.45681e-08,-6.60905e-12,69359.6,35.8367], Tmin=(100,'K'), Tmax=(1387.22,'K')), NASAPolynomial(coeffs=[16.3164,0.0281223,-7.87185e-06,1.1031e-09,-6.32268e-14,65279.1,-47.8477], Tmin=(1387.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC1CC(C=C)C1(26574)',
    structure = SMILES('C#CC1CC(C=C)C1'),
    E0 = (323.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11767,0.0416371,5.59639e-05,-1.04462e-07,4.37956e-11,39010.3,26.0673], Tmin=(100,'K'), Tmax=(948.609,'K')), NASAPolynomial(coeffs=[17.313,0.0278985,-8.57331e-06,1.51652e-09,-1.11312e-13,33483.3,-64.1485], Tmin=(948.609,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane)"""),
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
    label = 'C#C[CH]C[CH]C=C(27592)',
    structure = SMILES('C#C[CH]CC=C[CH2]'),
    E0 = (491.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,433.817,438.145],'cm^-1')),
        HinderedRotor(inertia=(0.0632098,'amu*angstrom^2'), symmetry=1, barrier=(8.75477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.643662,'amu*angstrom^2'), symmetry=1, barrier=(84.8187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608377,'amu*angstrom^2'), symmetry=1, barrier=(84.7284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658565,'amu*angstrom^2'), symmetry=1, barrier=(84.9366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924384,0.0603159,-3.48015e-05,-6.7064e-10,6.04832e-12,59229.7,25.848], Tmin=(100,'K'), Tmax=(934.069,'K')), NASAPolynomial(coeffs=[13.2645,0.0256867,-8.44297e-06,1.39437e-09,-9.22156e-14,56129.7,-37.0982], Tmin=(934.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Allyl_P)"""),
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
    label = '[CH2]C([CH2])C=C(2950)',
    structure = SMILES('[CH2]C([CH2])C=C'),
    E0 = (361.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,445.142],'cm^-1')),
        HinderedRotor(inertia=(0.0443051,'amu*angstrom^2'), symmetry=1, barrier=(6.21855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0429448,'amu*angstrom^2'), symmetry=1, barrier=(6.01401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0394571,'amu*angstrom^2'), symmetry=1, barrier=(69.5556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03658,0.0357957,3.13672e-06,-3.00425e-08,1.51079e-11,43603,21.9964], Tmin=(100,'K'), Tmax=(906.69,'K')), NASAPolynomial(coeffs=[9.64728,0.0219245,-6.51397e-06,1.02238e-09,-6.65481e-14,41412.9,-18.4418], Tmin=(906.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=[C]C1CC(C=C)C1(27593)',
    structure = SMILES('[CH]=[C]C1CC(C=C)C1'),
    E0 = (643.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.98679,0.0465488,3.79495e-05,-8.17958e-08,3.42384e-11,77498.5,29.6063], Tmin=(100,'K'), Tmax=(976.344,'K')), NASAPolynomial(coeffs=[16.6898,0.0303046,-1.09761e-05,2.05981e-09,-1.51093e-13,72140.2,-57.5136], Tmin=(976.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(643.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1C=C=CCC1[CH2](27266)',
    structure = SMILES('[CH2]C1C=C=CCC1[CH2]'),
    E0 = (514.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19925,0.11752,-0.000142317,9.59849e-08,-2.61533e-11,62055.9,18.5313], Tmin=(100,'K'), Tmax=(894.296,'K')), NASAPolynomial(coeffs=[15.3905,0.0433164,-1.7853e-05,3.20015e-09,-2.14949e-13,59088.7,-59.651], Tmin=(894.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(514.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C=C)C[C]=C=C(27594)',
    structure = SMILES('[CH2]C#CCC([CH2])C=C'),
    E0 = (517.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.482028,0.0690453,-3.83492e-05,1.49867e-09,4.72443e-12,62358,33.768], Tmin=(100,'K'), Tmax=(978.37,'K')), NASAPolynomial(coeffs=[13.3106,0.0342632,-1.21088e-05,2.07494e-09,-1.39002e-13,59002.2,-32.1632], Tmin=(978.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(517.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Isobutyl) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]CC(C)C=C(27595)',
    structure = SMILES('[CH]=C=[C]CC(C)C=C'),
    E0 = (564.845,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,1685,370,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,278.703,279.28],'cm^-1')),
        HinderedRotor(inertia=(0.212971,'amu*angstrom^2'), symmetry=1, barrier=(11.7676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21297,'amu*angstrom^2'), symmetry=1, barrier=(11.7695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212816,'amu*angstrom^2'), symmetry=1, barrier=(11.7714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489768,'amu*angstrom^2'), symmetry=1, barrier=(27.0657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.129972,0.0780757,-6.30066e-05,2.74113e-08,-4.87569e-12,68080.3,32.8931], Tmin=(100,'K'), Tmax=(1332.5,'K')), NASAPolynomial(coeffs=[15.0379,0.0333244,-1.26307e-05,2.20792e-09,-1.47161e-13,64107.3,-43.3088], Tmin=(1332.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C=C)C=C[C]=C(27247)',
    structure = SMILES('[CH2]C(C=C)C=C[C]=C'),
    E0 = (511.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.549316,'amu*angstrom^2'), symmetry=1, barrier=(12.6298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.550066,'amu*angstrom^2'), symmetry=1, barrier=(12.6471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.14173,'amu*angstrom^2'), symmetry=1, barrier=(72.2345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05189,'amu*angstrom^2'), symmetry=1, barrier=(24.1849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00957926,0.0846348,-7.96539e-05,4.20261e-08,-9.04715e-12,61610.1,30.8151], Tmin=(100,'K'), Tmax=(1118.29,'K')), NASAPolynomial(coeffs=[13.9113,0.0348411,-1.28637e-05,2.20926e-09,-1.45848e-13,58496.6,-37.9017], Tmin=(1118.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(511.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C=C)CC=C=C(27596)',
    structure = SMILES('[CH2]C=C([CH2])CC=C=C'),
    E0 = (457.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0624414,0.0747112,-4.23516e-05,1.64994e-09,4.59371e-12,55146,30.943], Tmin=(100,'K'), Tmax=(1056.32,'K')), NASAPolynomial(coeffs=[17.0253,0.032143,-1.2669e-05,2.33302e-09,-1.63221e-13,50353.6,-57.5437], Tmin=(1056.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=C)CC=C=C(27597)',
    structure = SMILES('[CH2]C([C]=C)CC=C=C'),
    E0 = (615.451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1685,370,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,394.007,394.21,394.26],'cm^-1')),
        HinderedRotor(inertia=(0.0849566,'amu*angstrom^2'), symmetry=1, barrier=(9.34348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0845386,'amu*angstrom^2'), symmetry=1, barrier=(9.34271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0846963,'amu*angstrom^2'), symmetry=1, barrier=(9.34187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291402,'amu*angstrom^2'), symmetry=1, barrier=(32.0644,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.266665,0.0772376,-6.23557e-05,2.74322e-08,-4.98185e-12,74160,33.5559], Tmin=(100,'K'), Tmax=(1296.9,'K')), NASAPolynomial(coeffs=[13.8953,0.035203,-1.37381e-05,2.44045e-09,-1.64235e-13,70625,-35.7375], Tmin=(1296.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(615.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])CC=C=C(27598)',
    structure = SMILES('[CH]=CC([CH2])CC=C=C'),
    E0 = (624.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,298.396,608.066],'cm^-1')),
        HinderedRotor(inertia=(0.489598,'amu*angstrom^2'), symmetry=1, barrier=(11.2568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0430347,'amu*angstrom^2'), symmetry=1, barrier=(11.2639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48998,'amu*angstrom^2'), symmetry=1, barrier=(11.2656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54414,'amu*angstrom^2'), symmetry=1, barrier=(35.5028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0242847,0.0790152,-6.37821e-05,2.73747e-08,-4.75508e-12,75287.4,34.5358], Tmin=(100,'K'), Tmax=(1373.22,'K')), NASAPolynomial(coeffs=[16.4215,0.0311101,-1.14536e-05,1.97005e-09,-1.30016e-13,70770.7,-50.0216], Tmin=(1373.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(624.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1[CH]CC=C=CC1(27599)',
    structure = SMILES('[CH2]C1[CH]CC=C=CC1'),
    E0 = (563.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2584,0.048066,1.30348e-05,-4.07246e-08,1.6134e-11,67848.8,25.1817], Tmin=(100,'K'), Tmax=(1052.97,'K')), NASAPolynomial(coeffs=[10.4149,0.0410669,-1.65745e-05,3.08104e-09,-2.161e-13,64380.2,-26.7798], Tmin=(1052.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(563.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = 'C=C=CCC(=C)C=C(27600)',
    structure = SMILES('C=C=CCC(=C)C=C'),
    E0 = (279.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.119854,0.0723755,-3.41348e-05,-9.82397e-09,9.56411e-12,33771.4,29.0407], Tmin=(100,'K'), Tmax=(1004.98,'K')), NASAPolynomial(coeffs=[17.8141,0.0294074,-1.09848e-05,2.00542e-09,-1.41408e-13,28828.3,-63.3106], Tmin=(1004.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH]=C(C=C)C[CH]C=C(26567)',
    structure = SMILES('[CH]=C(C=C)CC=C[CH2]'),
    E0 = (501.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0152626,0.0736648,-3.28953e-05,-1.45596e-08,1.20924e-11,60475.6,31.0457], Tmin=(100,'K'), Tmax=(986.764,'K')), NASAPolynomial(coeffs=[19.0509,0.0275334,-9.94313e-06,1.80425e-09,-1.27941e-13,55208,-68.1918], Tmin=(986.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Allyl_P)"""),
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
    E0 = (525.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (562.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (620.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (594.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (638.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (703.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (569.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (676.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (551.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (655.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (698.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (643.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (710.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (726.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (863.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (661.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (611.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (665.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (732.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (726.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (751.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (649.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (595.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (585.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (614.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (603.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (685.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (685.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (735.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (533.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (873.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (892.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (643.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (651.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (657.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (609.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (569.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (558.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (860.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (815.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (566.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (550.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (672.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['butadiene13(2459)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C#C[CH]CC1CC1[CH2](27573)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C#CC1CC([CH2])C1[CH2](27300)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.43e+10,'s^-1'), n=0.21, Ea=(94.9768,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 161 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_csHCt
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_csHCt]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['[CH]=C1[CH]CC(C=C)C1(27574)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(69.2285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;triplebond_intra_H;radadd_intra_cs2H] for rate rule [R6;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#C[CH]CC(=C)C=C(27575)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.88727,'m^3/(mol*s)'), n=1.9687, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 18 used for Cds-CdCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -8.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C#CC=CC([CH2])C=C(27576)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(2.76144,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2560 used for Cds-CsH_Cds-CtH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CtH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['butadiene13(2459)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0136943,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2H3(30)', 'C#C[CH]CC=C(26789)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6870,'cm^3/(mol*s)'), n=2.41, Ea=(13.7235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 432 used for Cds-CsH_Cds-HH;CdsJ-H
Exact match found for rate rule [Cds-CsH_Cds-HH;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH]C=C(2458)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0264797,'m^3/(mol*s)'), n=2.333, Ea=(2.89605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CtH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#C[CH]C[C](C)C=C(27577)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#CC[CH]C([CH2])C=C(27319)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C#C[CH][CH]C(C)C=C(27578)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#C[CH]CC(C)[C]=C(26580)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C#CCC[C]([CH2])C=C(27579)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.07201e+08,'s^-1'), n=1.25033, Ea=(201.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/OneDe;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_H/Ct;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[C]#CCCC([CH2])C=C(27580)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#CCCC([CH2])[C]=C(27581)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/Ct]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=CC(C)C[CH]C#C(27582)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC([CH2])CCC#C(27583)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(13813.5,'s^-1'), n=1.88327, Ea=(38.7799,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/Ct]
Euclidian distance = 3.74165738677
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[C]#C[CH]CC(C)C=C(27584)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.31883e+06,'s^-1'), n=1.02765, Ea=(75.0925,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R6HJ_2;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]C=C(2458)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['[CH2]C(C=C)CC1[C]=C1(27585)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_T;triplebond_intra_H;radadd_intra_cs] for rate rule [R3_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C#C[CH]CC1[CH]CC1(27586)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C#CC1C[CH]C([CH2])C1(27587)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.45491e+07,'s^-1'), n=1.06599, Ea=(69.5416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_csHCt]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C=CC1C[CH][C]=CC1(27588)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.33254e+09,'s^-1'), n=0.487896, Ea=(59.5573,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6_linear;multiplebond_intra;radadd_intra_cs2H] + [R6_linear;triplebond_intra_H;radadd_intra] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C#CCCC(=C)C=C(27589)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C#CC=CC(C)C=C(27590)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C#C[CH]CC[CH]C=C(26564)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C#C[CH]C[CH]CC=C(27591)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#CC([CH2])C([CH2])C=C(26569)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C#CC1CC(C=C)C1(26574)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH2(19)', 'C#C[CH]C[CH]C=C(27592)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['C3H2(81)', '[CH2]C([CH2])C=C(2950)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.26928e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['[CH]=[C]C1CC(C=C)C1(27593)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.20551e+07,'s^-1'), n=1.225, Ea=(117.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 113.8 to 117.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['[CH2]C1C=C=CCC1[CH2](27266)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(125.52,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7;doublebond_intra_2H_pri;radadd_intra_cdsingleH] for rate rule [R7_MMSR;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['[CH2]C(C=C)C[C]=C=C(27594)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C=[C]CC(C)C=C(27595)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['[CH2]C(C=C)C=C[C]=C(27247)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['[CH2][C](C=C)CC=C=C(27596)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C([C]=C)CC=C=C(27597)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0.627, Ea=(245.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_singleH] for rate rule [R6H;Cd_rad_out_Cd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=CC([CH2])CC=C=C(27598)'],
    products = ['C#C[CH]CC([CH2])C=C(26568)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R7H;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['[CH2]C1[CH]CC=C=CC1(27599)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['C=C=CCC(=C)C=C(27600)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C#C[CH]CC([CH2])C=C(26568)'],
    products = ['[CH]=C(C=C)C[CH]C=C(26567)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

network(
    label = '4538',
    isomers = [
        'C#C[CH]CC([CH2])C=C(26568)',
    ],
    reactants = [
        ('butadiene13(2459)', 'CH2CHCCH(26391)'),
        ('[CH2][CH]C=C(2458)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4538',
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

