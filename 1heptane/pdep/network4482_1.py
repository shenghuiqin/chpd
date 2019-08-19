species(
    label = '[CH]=CC[CH]C#C(26514)',
    structure = SMILES('[CH]=CC[CH]C#C'),
    E0 = (623.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,348.902],'cm^-1')),
        HinderedRotor(inertia=(0.121235,'amu*angstrom^2'), symmetry=1, barrier=(10.5335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.821338,'amu*angstrom^2'), symmetry=1, barrier=(71.0376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.822725,'amu*angstrom^2'), symmetry=1, barrier=(71.0309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24553,0.0537539,-5.00573e-05,2.56039e-08,-5.18339e-12,75047.8,23.0097], Tmin=(100,'K'), Tmax=(1291.12,'K')), NASAPolynomial(coeffs=[11.912,0.0179055,-5.15333e-06,7.36621e-10,-4.27892e-14,72527,-30.2708], Tmin=(1291.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(623.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_P)"""),
)

species(
    label = 'C2H2(32)',
    structure = SMILES('C#C'),
    E0 = (217.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,584.389,584.389,2772.01],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80868,0.0233616,-3.55172e-05,2.80153e-08,-8.50075e-12,26429,13.9397], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.65878,0.00488397,-1.60829e-06,2.46975e-10,-1.38606e-14,25759.4,-3.99838], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(217.784,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH]=C1[CH]CC=C1(27498)',
    structure = SMILES('[CH]C1C=CCC=1'),
    E0 = (455.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17037,0.0218777,6.51392e-05,-9.93375e-08,3.92635e-11,54845.3,17.4239], Tmin=(100,'K'), Tmax=(958.344,'K')), NASAPolynomial(coeffs=[13.1224,0.0222441,-7.55663e-06,1.40445e-09,-1.04788e-13,50630.2,-45.9871], Tmin=(958.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(455.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclopentadiene) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]=CC=CC#C(27720)',
    structure = SMILES('[CH]=CC=CC#C'),
    E0 = (569.575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.4068,'amu*angstrom^2'), symmetry=1, barrier=(32.345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4036,'amu*angstrom^2'), symmetry=1, barrier=(32.2715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45494,0.0445831,-1.27783e-05,-2.36073e-08,1.45126e-11,68606.1,19.0625], Tmin=(100,'K'), Tmax=(948.617,'K')), NASAPolynomial(coeffs=[16.2341,0.0102641,-2.78635e-06,4.85884e-10,-3.71599e-14,64542.3,-58.0992], Tmin=(948.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Cds_P)"""),
)

species(
    label = 'C#C[CH]CC#C(27721)',
    structure = SMILES('C#C[CH]CC#C'),
    E0 = (542.117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,3025,407.5,1350,352.5,297.566],'cm^-1')),
        HinderedRotor(inertia=(1.05299,'amu*angstrom^2'), symmetry=1, barrier=(66.1802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05323,'amu*angstrom^2'), symmetry=1, barrier=(66.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05336,'amu*angstrom^2'), symmetry=1, barrier=(66.1797,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30954,0.0536634,-5.75938e-05,3.371e-08,-7.5542e-12,65303.2,20.2123], Tmin=(100,'K'), Tmax=(1263.61,'K')), NASAPolynomial(coeffs=[11.0488,0.0153147,-3.1457e-06,2.74927e-10,-7.55966e-15,63442.2,-26.6778], Tmin=(1263.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(542.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C2H2(T)(175)',
    structure = SMILES('[CH]=[CH]'),
    E0 = (590.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([597.295,1198.22,1198.22,1198.22,2685.93,3758.47],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88057,-0.000843296,2.04089e-05,-2.51937e-08,9.47543e-12,71059.3,4.72119], Tmin=(100,'K'), Tmax=(935.375,'K')), NASAPolynomial(coeffs=[4.92145,0.00358266,-9.24397e-07,1.57263e-10,-1.19532e-14,70476.2,-2.30676], Tmin=(935.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""C2H2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=C[CH]CC#C(27508)',
    structure = SMILES('[CH]C=CCC#C'),
    E0 = (592.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.527,0.0476137,-2.38156e-05,1.63159e-09,1.66515e-12,71390.7,21.7855], Tmin=(100,'K'), Tmax=(1137.09,'K')), NASAPolynomial(coeffs=[10.0409,0.0277971,-1.10414e-05,1.97921e-09,-1.34292e-13,68799.4,-23.2635], Tmin=(1137.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(592.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#C[CH]C[C]=C(27523)',
    structure = SMILES('C#C[CH]C[C]=C'),
    E0 = (613.855,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,1685,370,276.5,276.957],'cm^-1')),
        HinderedRotor(inertia=(0.00217879,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44516,'amu*angstrom^2'), symmetry=1, barrier=(78.3056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44431,'amu*angstrom^2'), symmetry=1, barrier=(78.2948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61257,0.0509025,-4.40529e-05,1.85449e-08,-1.84925e-12,73917.3,21.7687], Tmin=(100,'K'), Tmax=(854.459,'K')), NASAPolynomial(coeffs=[9.60387,0.0217389,-7.3325e-06,1.18963e-09,-7.60246e-14,72250.6,-17.2893], Tmin=(854.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=[C]CCC#C(27722)',
    structure = SMILES('[CH]=[C]CCC#C'),
    E0 = (714.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2175,525,3120,650,792.5,1650,750,770,3400,2100,1685,370,260.073],'cm^-1')),
        HinderedRotor(inertia=(1.46212,'amu*angstrom^2'), symmetry=1, barrier=(70.1777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272324,'amu*angstrom^2'), symmetry=1, barrier=(13.0708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272323,'amu*angstrom^2'), symmetry=1, barrier=(13.0708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51183,0.0539541,-5.36722e-05,3.03209e-08,-7.01343e-12,86053.8,22.6926], Tmin=(100,'K'), Tmax=(1042.05,'K')), NASAPolynomial(coeffs=[9.66317,0.0226643,-8.63131e-06,1.50521e-09,-1.00162e-13,84355,-16.9686], Tmin=(1042.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(714.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C=CC=C(27476)',
    structure = SMILES('[CH]=C=CC=C[CH2]'),
    E0 = (471.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69631,0.0401007,1.64624e-06,-3.65636e-08,1.91391e-11,56741.8,21.1146], Tmin=(100,'K'), Tmax=(915.452,'K')), NASAPolynomial(coeffs=[13.6109,0.0157687,-3.91799e-06,5.745e-10,-3.83951e-14,53398.5,-41.6597], Tmin=(915.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[C]#CCCC=[CH](27723)',
    structure = SMILES('[C]#CCCC=[CH]'),
    E0 = (814.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2175,525,325.912,326.299],'cm^-1')),
        HinderedRotor(inertia=(0.115462,'amu*angstrom^2'), symmetry=1, barrier=(8.69276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115376,'amu*angstrom^2'), symmetry=1, barrier=(8.68857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.876315,'amu*angstrom^2'), symmetry=1, barrier=(66.1171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49441,0.0537885,-5.47e-05,3.23177e-08,-7.78773e-12,97998.1,23.0217], Tmin=(100,'K'), Tmax=(1006.68,'K')), NASAPolynomial(coeffs=[9.40552,0.0223521,-7.85536e-06,1.29322e-09,-8.26069e-14,96405.4,-15.197], Tmin=(1006.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(814.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Cds_P)"""),
)

species(
    label = '[C]#C[CH]CC=C(27724)',
    structure = SMILES('[C]#C[CH]CC=C'),
    E0 = (713.157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2175,525,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,309.804,313.278,1939.81],'cm^-1')),
        HinderedRotor(inertia=(0.124211,'amu*angstrom^2'), symmetry=1, barrier=(8.87784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.948768,'amu*angstrom^2'), symmetry=1, barrier=(68.7186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.98275,'amu*angstrom^2'), symmetry=1, barrier=(68.6078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74185,0.0488309,-3.74503e-05,8.92627e-09,3.28566e-12,85855.3,21.5821], Tmin=(100,'K'), Tmax=(780.362,'K')), NASAPolynomial(coeffs=[9.22717,0.021645,-6.68908e-06,1.0101e-09,-6.1233e-14,84346.6,-14.8555], Tmin=(780.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(713.157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Acetyl)"""),
)

species(
    label = '[CH]=CCC1[C]=C1(27725)',
    structure = SMILES('[CH]=CCC1[C]=C1'),
    E0 = (796.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64412,0.0507429,-4.32187e-05,1.98824e-08,-3.78568e-12,95896.7,19.9439], Tmin=(100,'K'), Tmax=(1234.6,'K')), NASAPolynomial(coeffs=[10.151,0.0231813,-9.73238e-06,1.80027e-09,-1.24158e-13,93796.1,-22.8898], Tmin=(1234.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(796.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(cyclopropenyl-vinyl) + radical(Cds_P)"""),
)

species(
    label = '[C]1[CH]CC=CC=1(27481)',
    structure = SMILES('[C]1=CC[CH]C=C1'),
    E0 = (384.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30338,0.0306596,5.26427e-06,-1.95796e-08,7.17224e-12,46281.9,15.1439], Tmin=(100,'K'), Tmax=(1144.3,'K')), NASAPolynomial(coeffs=[7.24701,0.028685,-1.22115e-05,2.29106e-09,-1.59741e-13,44148.4,-13.7513], Tmin=(1144.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(384.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(C=CJC=C) + radical(Aromatic_pi_S_1_3)"""),
)

species(
    label = 'C#CC=CC=C(27726)',
    structure = SMILES('C#CC=CC=C'),
    E0 = (322.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52703,0.0410803,3.80947e-06,-4.09006e-08,2.03914e-11,38886.4,18.4006], Tmin=(100,'K'), Tmax=(951.737,'K')), NASAPolynomial(coeffs=[15.9979,0.0130904,-3.81659e-06,6.83769e-10,-5.19547e-14,34645.1,-58.508], Tmin=(951.737,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = '[CH]=CC([CH2])C#C(26515)',
    structure = SMILES('[CH]=CC([CH2])C#C'),
    E0 = (648.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435],'cm^-1')),
        HinderedRotor(inertia=(0.844145,'amu*angstrom^2'), symmetry=1, barrier=(19.4086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843649,'amu*angstrom^2'), symmetry=1, barrier=(19.3972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.471,'amu*angstrom^2'), symmetry=1, barrier=(56.8131,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3244.47,'J/mol'), sigma=(5.6374,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=506.78 K, Pc=41.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16694,0.060797,-6.72891e-05,4.05297e-08,-9.72226e-12,78149.1,22.0966], Tmin=(100,'K'), Tmax=(1019.62,'K')), NASAPolynomial(coeffs=[11.5905,0.0199042,-7.12909e-06,1.19394e-09,-7.73697e-14,76023.5,-28.3934], Tmin=(1019.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = 'C#CC1C=CC1(27509)',
    structure = SMILES('C#CC1C=CC1'),
    E0 = (357.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94316,0.0321813,2.08389e-05,-5.22416e-08,2.28935e-11,43109.6,15.6544], Tmin=(100,'K'), Tmax=(959.803,'K')), NASAPolynomial(coeffs=[13.3312,0.0168948,-5.55246e-06,1.0143e-09,-7.49195e-14,39441.6,-46.5392], Tmin=(959.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene)"""),
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
    label = '[CH]=C[CH2](2461)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,229.711,230.18,230.787],'cm^-1')),
        HinderedRotor(inertia=(1.33306,'amu*angstrom^2'), symmetry=1, barrier=(50.5153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,44095,-3.44606], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C1C=CC1(27727)',
    structure = SMILES('[CH]=[C]C1C=CC1'),
    E0 = (701.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93605,0.0368927,-3.47843e-06,-1.95547e-08,9.42062e-12,84405.4,19.2252], Tmin=(100,'K'), Tmax=(1023.17,'K')), NASAPolynomial(coeffs=[10.6445,0.021569,-8.45943e-06,1.57395e-09,-1.11491e-13,81643.4,-27.7763], Tmin=(1023.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(701.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC[C]=C=C(27728)',
    structure = SMILES('[CH]=CCC#C[CH2]'),
    E0 = (619.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75283,0.0431105,-2.43106e-05,3.42395e-09,1.10565e-12,74589,22.9402], Tmin=(100,'K'), Tmax=(1132.63,'K')), NASAPolynomial(coeffs=[10.4241,0.0219075,-8.70632e-06,1.58249e-09,-1.08712e-13,72020.5,-22.641], Tmin=(1132.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(619.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cds_P) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]CC=C(27729)',
    structure = SMILES('[CH]=C=[C]CC=C'),
    E0 = (621.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53045,0.0497898,-4.15108e-05,1.88964e-08,-3.50845e-12,74859.5,23.6481], Tmin=(100,'K'), Tmax=(1285.2,'K')), NASAPolynomial(coeffs=[10.8683,0.020727,-7.59074e-06,1.30122e-09,-8.57954e-14,72459.3,-23.7448], Tmin=(1285.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(621.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC=C[C]=C(27474)',
    structure = SMILES('[CH]C=CC=C=C'),
    E0 = (569.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38513,0.048153,-1.41305e-05,-1.58135e-08,9.49626e-12,68556.4,21.1552], Tmin=(100,'K'), Tmax=(995.315,'K')), NASAPolynomial(coeffs=[12.4634,0.0242385,-9.14606e-06,1.64942e-09,-1.14883e-13,65330.4,-37.3666], Tmin=(995.315,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC=C=C(27730)',
    structure = SMILES('[CH]=[C]CC=C=C'),
    E0 = (714.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1685,370,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.793719,'amu*angstrom^2'), symmetry=1, barrier=(18.2492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.794406,'amu*angstrom^2'), symmetry=1, barrier=(18.265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68641,0.0494647,-4.02886e-05,1.7498e-08,-3.14662e-12,85990.8,22.8642], Tmin=(100,'K'), Tmax=(1298.36,'K')), NASAPolynomial(coeffs=[10.2755,0.0230035,-9.71798e-06,1.80111e-09,-1.24189e-13,83760.4,-20.816], Tmin=(1298.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(714.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(Cds_S)"""),
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
    E0 = (623.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (638.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (789.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (774.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (867.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (696.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (749.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (728.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (839.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (815.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (961.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (896.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1042.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (849.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (645.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (701.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (808.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (630.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (907.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (701.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (755.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (667.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (667.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (826.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['C2H2(32)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['[CH]=C1[CH]CC=C1(27498)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.42e+11,'s^-1'), n=0.258, Ea=(15.8866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[CH]=CC=CC#C(27720)'],
    products = ['[CH]=CC[CH]C#C(26514)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(105.863,'m^3/(mol*s)'), n=1.66278, Ea=(8.08939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#C[CH]CC#C(27721)'],
    products = ['[CH]=CC[CH]C#C(26514)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.276e+10,'cm^3/(mol*s)'), n=1.103, Ea=(20.5476,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 139 used for Ct-Cs_Ct-H;HJ
Exact match found for rate rule [Ct-Cs_Ct-H;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H2(T)(175)', 'CH2CHCCH(26391)'],
    products = ['[CH]=CC[CH]C#C(26514)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0264797,'m^3/(mol*s)'), n=2.333, Ea=(2.89605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CtH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H2(32)', '[CH]=[C]C=C(4699)'],
    products = ['[CH]=CC[CH]C#C(26514)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(46.4627,'m^3/(mol*s)'), n=1.51997, Ea=(27.4714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-H;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['[CH]=C[CH]CC#C(27508)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.26084e+07,'s^-1'), n=1.66833, Ea=(126.789,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_H/OneDe;Cs_H_out_H/OneDe] + [R2H_S;C_rad_out_H/Ct;Cs_H_out_1H] for rate rule [R2H_S;C_rad_out_H/Ct;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['C#C[CH]C[C]=C(27523)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]CCC#C(27722)'],
    products = ['[CH]=CC[CH]C#C(26514)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['[CH]=[C]C=CC=C(27476)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[C]#CCCC=[CH](27723)'],
    products = ['[CH]=CC[CH]C#C(26514)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[C]#C[CH]CC=C(27724)'],
    products = ['[CH]=CC[CH]C#C(26514)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(640643,'s^-1'), n=2.07799, Ea=(182.911,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Y_rad_out;Cd_H_out_singleH] for rate rule [R6HJ_2;Ct_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C2H2(T)(175)', '[CH]=[C]C=C(4699)'],
    products = ['[CH]=CC[CH]C#C(26514)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['[CH]=CCC1[C]=C1(27725)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_T;triplebond_intra_H;radadd_intra_cs] for rate rule [R3_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['[C]1[CH]CC=CC=1(27481)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.926e+10,'s^-1'), n=0.198, Ea=(22.8237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['C#CC=CC=C(27726)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=CC([CH2])C#C(26515)'],
    products = ['[CH]=CC[CH]C#C(26514)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['C#CC1C=CC1(27509)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriH_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C3H2(81)', '[CH]=C[CH2](2461)'],
    products = ['[CH]=CC[CH]C#C(26514)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['[CH]=[C]C1C=CC1(27727)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(77.9976,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 76.0 to 78.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['[CH]=CC[C]=C=C(27728)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['[CH]=C=[C]CC=C(27729)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=CC[CH]C#C(26514)'],
    products = ['[CH]=CC=C[C]=C(27474)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]CC=C=C(27730)'],
    products = ['[CH]=CC[CH]C#C(26514)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(33613.8,'s^-1'), n=2.10442, Ea=(111.806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cd_H_out_singleH] for rate rule [R5H;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = '4482',
    isomers = [
        '[CH]=CC[CH]C#C(26514)',
    ],
    reactants = [
        ('C2H2(32)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4482',
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

