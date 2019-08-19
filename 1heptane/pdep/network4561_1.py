species(
    label = 'C#C[CH]CC([CH2])C#C(26589)',
    structure = SMILES('C#C[CH]CC([CH2])C#C'),
    E0 = (690.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,245.282,1846.91],'cm^-1')),
        HinderedRotor(inertia=(1.63004,'amu*angstrom^2'), symmetry=1, barrier=(69.5913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.63003,'amu*angstrom^2'), symmetry=1, barrier=(69.5912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.63002,'amu*angstrom^2'), symmetry=1, barrier=(69.5912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00280207,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.63005,'amu*angstrom^2'), symmetry=1, barrier=(69.5912,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.413,0.0775862,-6.82273e-05,2.00101e-08,4.75277e-12,83180.4,29.4531], Tmin=(100,'K'), Tmax=(770.958,'K')), NASAPolynomial(coeffs=[14.2582,0.0259108,-6.9066e-06,8.99538e-10,-4.79033e-14,80446.5,-37.6258], Tmin=(770.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(690.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Isobutyl) + radical(Sec_Propargyl)"""),
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
    label = '[CH]=C1CC1C[CH]C#C(28310)',
    structure = SMILES('[CH]=C1CC1C[CH]C#C'),
    E0 = (771.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.520488,0.0658238,-3.01389e-05,-1.49778e-08,1.32117e-11,92870.3,27.2635], Tmin=(100,'K'), Tmax=(921.222,'K')), NASAPolynomial(coeffs=[16.6621,0.023643,-6.89741e-06,1.08685e-09,-7.19489e-14,88712.1,-55.713], Tmin=(921.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(771.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Methylene_cyclopropane) + radical(Sec_Propargyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C([CH2])CC1C#C(28235)',
    structure = SMILES('[CH]=C1C([CH2])CC1C#C'),
    E0 = (734.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.7098,0.00121743,0.000103239,-1.06398e-07,2.44459e-11,87911.7,-18.5579], Tmin=(100,'K'), Tmax=(1768.54,'K')), NASAPolynomial(coeffs=[85.2364,0.0166553,-6.59147e-05,1.61936e-08,-1.19961e-12,32776.1,-501.95], Tmin=(1768.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(734.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(methylenecyclobutane) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1[CH]CC(C#C)C1(28303)',
    structure = SMILES('[CH]C1=CCC(C#C)C1'),
    E0 = (590.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06358,0.0465445,3.55102e-05,-7.92133e-08,3.39046e-11,71182.1,24.9026], Tmin=(100,'K'), Tmax=(954.738,'K')), NASAPolynomial(coeffs=[15.3395,0.0313083,-1.05839e-05,1.87424e-09,-1.33164e-13,66424.6,-53.9486], Tmin=(954.738,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopentene) + radical(AllylJ2_triplet)"""),
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
    label = 'C#C[CH]CC(=C)C#C(28311)',
    structure = SMILES('C#C[CH]CC(=C)C#C'),
    E0 = (600.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,750,756.667,763.333,770,3350,3450,2000,2200,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,362.966,362.979],'cm^-1')),
        HinderedRotor(inertia=(0.85887,'amu*angstrom^2'), symmetry=1, barrier=(80.294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858768,'amu*angstrom^2'), symmetry=1, barrier=(80.2934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162485,'amu*angstrom^2'), symmetry=1, barrier=(15.1923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.85907,'amu*angstrom^2'), symmetry=1, barrier=(80.294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.467001,0.0735506,-6.7945e-05,3.00024e-08,-3.77184e-12,72311.4,26.7014], Tmin=(100,'K'), Tmax=(904.235,'K')), NASAPolynomial(coeffs=[14.2698,0.0250314,-8.25919e-06,1.33368e-09,-8.56323e-14,69302.6,-41.3343], Tmin=(904.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCtCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC=CC([CH2])C#C(28312)',
    structure = SMILES('C#CC=CC([CH2])C#C'),
    E0 = (630.371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.36642,'amu*angstrom^2'), symmetry=1, barrier=(31.4168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36218,'amu*angstrom^2'), symmetry=1, barrier=(31.3192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36848,'amu*angstrom^2'), symmetry=1, barrier=(31.4641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36256,'amu*angstrom^2'), symmetry=1, barrier=(31.3279,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0459988,0.08189,-8.57925e-05,4.68765e-08,-1.01084e-11,75962.8,27.2457], Tmin=(100,'K'), Tmax=(1134.05,'K')), NASAPolynomial(coeffs=[16.6043,0.0234854,-8.54041e-06,1.46238e-09,-9.6824e-14,72207.3,-54.7212], Tmin=(1134.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(630.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + group(Ct-CtH) + radical(Isobutyl)"""),
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
    label = 'C#C[CH]C[C](C)C#C(28313)',
    structure = SMILES('C#C[CH]C[C](C)C#C'),
    E0 = (620.552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,2750,2850,1437.5,1250,1305,750,350,360,370,350,251.996,252.103],'cm^-1')),
        HinderedRotor(inertia=(0.0243252,'amu*angstrom^2'), symmetry=1, barrier=(69.3079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53851,'amu*angstrom^2'), symmetry=1, barrier=(69.307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198189,'amu*angstrom^2'), symmetry=1, barrier=(8.93972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53624,'amu*angstrom^2'), symmetry=1, barrier=(69.3101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53799,'amu*angstrom^2'), symmetry=1, barrier=(69.3078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00865188,0.086245,-0.000101648,6.67675e-08,-1.71282e-11,74781.1,28.5714], Tmin=(100,'K'), Tmax=(1070.31,'K')), NASAPolynomial(coeffs=[13.2721,0.0280954,-8.21885e-06,1.13897e-09,-6.2525e-14,72426,-34.1244], Tmin=(1070.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(620.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Tert_Propargyl)"""),
)

species(
    label = 'C#CC[CH]C([CH2])C#C(28291)',
    structure = SMILES('C#CC[CH]C([CH2])C#C'),
    E0 = (738.845,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,255.387,255.419],'cm^-1')),
        HinderedRotor(inertia=(0.010064,'amu*angstrom^2'), symmetry=1, barrier=(77.4601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67341,'amu*angstrom^2'), symmetry=1, barrier=(77.4598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6735,'amu*angstrom^2'), symmetry=1, barrier=(77.46,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0100641,'amu*angstrom^2'), symmetry=1, barrier=(77.4598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67326,'amu*angstrom^2'), symmetry=1, barrier=(77.4598,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.571046,0.0760411,-7.35918e-05,3.49874e-08,-4.09851e-12,88985.8,31.5974], Tmin=(100,'K'), Tmax=(773.851,'K')), NASAPolynomial(coeffs=[11.6824,0.0311181,-1.07655e-05,1.75463e-09,-1.1159e-13,86891.4,-21.5799], Tmin=(773.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(738.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Isobutyl) + radical(Cs_S)"""),
)

species(
    label = 'C#C[CH][CH]C(C)C#C(28314)',
    structure = SMILES('[CH]=C=C[CH]C(C)C#C'),
    E0 = (633.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.452216,0.070726,-5.04593e-05,1.13789e-08,2.36953e-12,76274.6,27.9614], Tmin=(100,'K'), Tmax=(945.801,'K')), NASAPolynomial(coeffs=[14.7019,0.0276189,-9.30489e-06,1.55108e-09,-1.02646e-13,72811.7,-44.0478], Tmin=(945.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtH) + radical(C=C=CJ) + radical(Allyl_S)"""),
)

species(
    label = 'C#CCC[C]([CH2])C#C(28315)',
    structure = SMILES('[CH]=C=C([CH2])CCC#C'),
    E0 = (637.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.103518,0.0804752,-7.50086e-05,3.42683e-08,-5.07279e-12,76794.9,29.3926], Tmin=(100,'K'), Tmax=(937.251,'K')), NASAPolynomial(coeffs=[15.57,0.0268725,-9.07526e-06,1.49184e-09,-9.69381e-14,73350.8,-47.1287], Tmin=(937.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(637.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(C=C=CJ)"""),
)

species(
    label = '[C]#CCCC([CH2])C#C(28316)',
    structure = SMILES('[C]#CCCC([CH2])C#C'),
    E0 = (881.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2100,2250,500,550,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.541114,0.0786974,-7.45226e-05,2.45641e-08,5.81388e-12,106136,29.9176], Tmin=(100,'K'), Tmax=(695.304,'K')), NASAPolynomial(coeffs=[12.4935,0.0291059,-8.89202e-06,1.28784e-09,-7.38281e-14,104011,-26.7343], Tmin=(695.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(881.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Isobutyl) + radical(Acetyl)"""),
)

species(
    label = '[C]#CC(C)C[CH]C#C(28317)',
    structure = SMILES('[C]#CC(C)C[CH]C#C'),
    E0 = (822.575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2250,500,550,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.090383,0.0837798,-9.60517e-05,6.20129e-08,-1.57417e-11,99075.5,29.5427], Tmin=(100,'K'), Tmax=(1069.34,'K')), NASAPolynomial(coeffs=[12.9485,0.0285423,-8.55252e-06,1.21851e-09,-6.88332e-14,96733.8,-31.4435], Tmin=(1069.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(822.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl)"""),
)

species(
    label = '[C]#CC([CH2])CCC#C(28318)',
    structure = SMILES('[C]#CC([CH2])CCC#C'),
    E0 = (881.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2100,2250,500,550,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.541149,0.0786968,-7.45203e-05,2.45602e-08,5.81617e-12,106136,29.9175], Tmin=(100,'K'), Tmax=(695.298,'K')), NASAPolynomial(coeffs=[12.4935,0.029106,-8.89206e-06,1.28785e-09,-7.3829e-14,104011,-26.7341], Tmin=(695.298,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(881.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Acetyl) + radical(Isobutyl)"""),
)

species(
    label = '[C]#C[CH]CC(C)C#C(28319)',
    structure = SMILES('[C]#C[CH]CC(C)C#C'),
    E0 = (822.575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2250,500,550,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.090383,0.0837798,-9.60517e-05,6.20129e-08,-1.57417e-11,99075.5,29.5427], Tmin=(100,'K'), Tmax=(1069.34,'K')), NASAPolynomial(coeffs=[12.9485,0.0285423,-8.55252e-06,1.21851e-09,-6.88332e-14,96733.8,-31.4435], Tmin=(1069.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(822.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Acetyl)"""),
)

species(
    label = 'C#CC([CH2])CC1[C]=C1(28320)',
    structure = SMILES('C#CC([CH2])CC1[C]=C1'),
    E0 = (864.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.227296,0.0818741,-8.92371e-05,5.46917e-08,-1.35077e-11,104055,28.4599], Tmin=(100,'K'), Tmax=(986.684,'K')), NASAPolynomial(coeffs=[12.9159,0.0304347,-1.10367e-05,1.85462e-09,-1.20116e-13,101551,-32.5852], Tmin=(986.684,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(864.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropene) + radical(Isobutyl) + radical(cyclopropenyl-vinyl)"""),
)

species(
    label = 'C#C[CH]CC1[C]=CC1(28321)',
    structure = SMILES('C#C[CH]CC1[C]=CC1'),
    E0 = (733.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.900158,0.056137,-8.17096e-06,-3.31069e-08,1.83468e-11,88397.2,27.0946], Tmin=(100,'K'), Tmax=(935.394,'K')), NASAPolynomial(coeffs=[15.0647,0.0257977,-7.99928e-06,1.32326e-09,-9.00366e-14,84424.7,-47.3649], Tmin=(935.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene) + radical(Sec_Propargyl) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'C#CC1C=[C]C([CH2])C1(28251)',
    structure = SMILES('C#CC1C=[C]C([CH2])C1'),
    E0 = (666.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.813743,0.0576012,-1.14906e-05,-2.87829e-08,1.61396e-11,80237.6,24.634], Tmin=(100,'K'), Tmax=(961.369,'K')), NASAPolynomial(coeffs=[15.4956,0.0262228,-8.88541e-06,1.55465e-09,-1.08497e-13,76041.8,-52.7591], Tmin=(961.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(666.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopentene) + radical(Isobutyl) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C#CC1C[CH][C]=CC1(28322)',
    structure = SMILES('C#CC1C[CH][C]=CC1'),
    E0 = (592.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48659,0.0366411,4.91591e-05,-8.92689e-08,3.66889e-11,71353.6,20.7371], Tmin=(100,'K'), Tmax=(959.453,'K')), NASAPolynomial(coeffs=[14.7287,0.0277497,-9.3492e-06,1.69777e-09,-1.23619e-13,66680.8,-53.7094], Tmin=(959.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(592.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclohexene) + radical(cyclohexene-allyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CCCC(=C)C#C(28323)',
    structure = SMILES('C#CCCC(=C)C#C'),
    E0 = (453.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.259764,0.0751668,-6.79627e-05,3.3074e-08,-6.4642e-12,54735.9,27.6059], Tmin=(100,'K'), Tmax=(1236.74,'K')), NASAPolynomial(coeffs=[15.1357,0.0270541,-9.60922e-06,1.61896e-09,-1.0583e-13,51056.3,-47.3229], Tmin=(1236.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCtCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC=CC(C)C#C(28324)',
    structure = SMILES('C#CC=CC(C)C#C'),
    E0 = (425.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0113703,0.0806111,-7.44027e-05,3.53183e-08,-6.66502e-12,51301.2,25.6733], Tmin=(100,'K'), Tmax=(1281.36,'K')), NASAPolynomial(coeffs=[17.6311,0.0255367,-9.93054e-06,1.77457e-09,-1.20418e-13,46779.9,-63.8151], Tmin=(1281.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + group(Ct-CtH)"""),
)

species(
    label = 'C#C[CH]CC[CH]C#C(26588)',
    structure = SMILES('C#C[CH]CC[CH]C#C'),
    E0 = (640.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00747775,0.0843511,-9.48732e-05,5.92416e-08,-1.44726e-11,77214,28.5753], Tmin=(100,'K'), Tmax=(1110.58,'K')), NASAPolynomial(coeffs=[13.9285,0.026959,-7.63376e-06,1.03634e-09,-5.62696e-14,74562.5,-38.1205], Tmin=(1110.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(640.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH]C[CH]CC#C(28325)',
    structure = SMILES('C#C[CH]C[CH]CC#C'),
    E0 = (689.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.426906,0.0817575,-0.000101943,7.58496e-08,-2.2516e-11,82995,30.4421], Tmin=(100,'K'), Tmax=(946.923,'K')), NASAPolynomial(coeffs=[9.04263,0.0354563,-1.2906e-05,2.11592e-09,-1.32437e-13,81807.4,-8.30869], Tmin=(946.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(RCCJCC) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC([CH2])C([CH2])C#C(26592)',
    structure = SMILES('C#CC([CH2])C([CH2])C#C'),
    E0 = (740.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,343.401],'cm^-1')),
        HinderedRotor(inertia=(0.0868022,'amu*angstrom^2'), symmetry=1, barrier=(7.26059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.831556,'amu*angstrom^2'), symmetry=1, barrier=(69.5361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.831076,'amu*angstrom^2'), symmetry=1, barrier=(69.5377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.831467,'amu*angstrom^2'), symmetry=1, barrier=(69.5362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0247503,'amu*angstrom^2'), symmetry=1, barrier=(69.536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3581.95,'J/mol'), sigma=(6.26163,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.49 K, Pc=33.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.422844,0.0761775,-6.30966e-05,1.35801e-08,7.2947e-12,89164.4,30.3868], Tmin=(100,'K'), Tmax=(783.028,'K')), NASAPolynomial(coeffs=[14.8814,0.0243001,-5.82893e-06,6.75401e-10,-3.20101e-14,86226.2,-40.1344], Tmin=(783.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(740.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC1CC(C#C)C1(26598)',
    structure = SMILES('C#CC1CC(C#C)C1'),
    E0 = (488.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14887,0.0447403,3.30347e-05,-8.13388e-08,3.69471e-11,58844.9,22.3781], Tmin=(100,'K'), Tmax=(922.648,'K')), NASAPolynomial(coeffs=[17.078,0.0220103,-5.33027e-06,8.04052e-10,-5.64045e-14,53933.6,-63.8744], Tmin=(922.648,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + ring(Cyclobutane)"""),
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
    label = 'C#C[CH]C[CH]C#C(28326)',
    structure = SMILES('C#C[CH]C[CH]C#C'),
    E0 = (664.547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,3000,3050,390,425,1340,1360,335,370,302.995,1819.84],'cm^-1')),
        HinderedRotor(inertia=(1.04658,'amu*angstrom^2'), symmetry=1, barrier=(69.1286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0415,'amu*angstrom^2'), symmetry=1, barrier=(69.0912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04252,'amu*angstrom^2'), symmetry=1, barrier=(69.179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05174,'amu*angstrom^2'), symmetry=1, barrier=(69.1917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.1225,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.659248,0.0693599,-8.04423e-05,5.06247e-08,-1.21856e-11,80050.4,23.9347], Tmin=(100,'K'), Tmax=(1175.95,'K')), NASAPolynomial(coeffs=[12.4762,0.0194204,-4.31177e-06,4.18405e-10,-1.40061e-14,77945,-32.1256], Tmin=(1175.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(664.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Sec_Propargyl)"""),
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
    label = 'C#CC([CH2])[CH2](26629)',
    structure = SMILES('C#CC([CH2])[CH2]'),
    E0 = (526.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00245333,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00268778,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19638,'amu*angstrom^2'), symmetry=1, barrier=(57.3712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59865,0.0445704,-4.03451e-05,2.09026e-08,-4.1805e-12,63457.8,21.3678], Tmin=(100,'K'), Tmax=(1443.27,'K')), NASAPolynomial(coeffs=[9.29133,0.0160372,-3.19385e-06,2.79115e-10,-8.33799e-15,61988.6,-15.9649], Tmin=(1443.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(526.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=[C]C1CC(C#C)C1(28327)',
    structure = SMILES('[CH]=[C]C1CC(C#C)C1'),
    E0 = (808.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01424,0.0497237,1.46207e-05,-5.78946e-08,2.69219e-11,97333.3,27.3153], Tmin=(100,'K'), Tmax=(951.089,'K')), NASAPolynomial(coeffs=[16.3554,0.0245793,-7.82448e-06,1.36849e-09,-9.79124e-14,92634.2,-55.29], Tmin=(951.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(808.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C=C=CCC1[CH2](28052)',
    structure = SMILES('[CH]=C1C=C=CCC1[CH2]'),
    E0 = (633.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20193,0.0380142,6.25694e-05,-1.10239e-07,4.50518e-11,76323.7,21.6618], Tmin=(100,'K'), Tmax=(965.943,'K')), NASAPolynomial(coeffs=[18.3152,0.0255846,-8.87497e-06,1.70037e-09,-1.29397e-13,70291.4,-74.4187], Tmin=(965.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclohexane) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = 'C#CC([CH2])C[C]=C=C(28328)',
    structure = SMILES('C#CC([CH2])CC#C[CH2]'),
    E0 = (682.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.619996,0.0707927,-5.60458e-05,1.7044e-08,1.50869e-12,82188.1,31.0883], Tmin=(100,'K'), Tmax=(852.701,'K')), NASAPolynomial(coeffs=[12.8399,0.028791,-9.11201e-06,1.42175e-09,-8.90814e-14,79547.1,-29.1845], Tmin=(852.701,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(682.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]CC(C)C#C(28329)',
    structure = SMILES('[CH]=C=[C]CC(C)C#C'),
    E0 = (729.795,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.637702,'amu*angstrom^2'), symmetry=1, barrier=(14.662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.637312,'amu*angstrom^2'), symmetry=1, barrier=(14.653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.33636,'amu*angstrom^2'), symmetry=1, barrier=(76.7095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.34313,'amu*angstrom^2'), symmetry=1, barrier=(76.8651,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.231316,0.0804322,-8.38018e-05,4.86381e-08,-1.13659e-11,87911.9,30.3344], Tmin=(100,'K'), Tmax=(1041.31,'K')), NASAPolynomial(coeffs=[13.4585,0.0296215,-1.06078e-05,1.77698e-09,-1.15194e-13,85157.2,-34.0141], Tmin=(1041.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.795,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtH) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = 'C#CC([CH2])C=C[C]=C(28257)',
    structure = SMILES('C#CC([CH2])C=C[C]=C'),
    E0 = (652.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.14611,'amu*angstrom^2'), symmetry=1, barrier=(26.3512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.99228,'amu*angstrom^2'), symmetry=1, barrier=(45.8063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13279,'amu*angstrom^2'), symmetry=1, barrier=(26.045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14769,'amu*angstrom^2'), symmetry=1, barrier=(26.3876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.106035,0.0879548,-9.62203e-05,5.49266e-08,-1.18564e-11,78637.3,28.4965], Tmin=(100,'K'), Tmax=(897.011,'K')), NASAPolynomial(coeffs=[15.6045,0.0269246,-9.25904e-06,1.51484e-09,-9.70545e-14,75455.7,-47.6144], Tmin=(897.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(652.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJC=C) + radical(Isobutyl)"""),
)

species(
    label = 'C#C[C]([CH2])CC=C=C(28330)',
    structure = SMILES('[CH]=C=C([CH2])CC=C=C'),
    E0 = (636.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0243356,0.0789117,-7.15758e-05,3.40036e-08,-6.43597e-12,76743,30.4793], Tmin=(100,'K'), Tmax=(1279.25,'K')), NASAPolynomial(coeffs=[16.9728,0.0259166,-9.43537e-06,1.61964e-09,-1.07243e-13,72406.8,-55.461], Tmin=(1279.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(C=C=CJ)"""),
)

species(
    label = '[C]#CC([CH2])CC=C=C(28331)',
    structure = SMILES('[C]#CC([CH2])CC=C=C'),
    E0 = (879.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2175,525,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,193.15,193.728,3107.67],'cm^-1')),
        HinderedRotor(inertia=(0.431029,'amu*angstrom^2'), symmetry=1, barrier=(11.4601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.88043,'amu*angstrom^2'), symmetry=1, barrier=(76.3994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.86825,'amu*angstrom^2'), symmetry=1, barrier=(76.4087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.87558,'amu*angstrom^2'), symmetry=1, barrier=(76.4096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.282991,0.0801761,-8.65715e-05,5.34669e-08,-1.33297e-11,105939,31.5722], Tmin=(100,'K'), Tmax=(978.204,'K')), NASAPolynomial(coeffs=[12.3933,0.0306611,-1.0653e-05,1.73286e-09,-1.09534e-13,103569,-26.5877], Tmin=(978.204,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(879.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtH) + radical(Acetyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[C]=CC=C=CC1(28019)',
    structure = SMILES('[CH2]C1[C]=CC=C=CC1'),
    E0 = (706.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21388,0.0464091,1.818e-05,-5.61421e-08,2.48128e-11,85136.7,21.4365], Tmin=(100,'K'), Tmax=(967.016,'K')), NASAPolynomial(coeffs=[14.255,0.0285038,-9.94751e-06,1.78807e-09,-1.27105e-13,80929.5,-49.7544], Tmin=(967.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(706.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC(=C)CC=C=C(28332)',
    structure = SMILES('C#CC(=C)CC=C=C'),
    E0 = (453.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.310443,0.0719905,-5.85343e-05,2.45797e-08,-4.14261e-12,54678.6,28.2324], Tmin=(100,'K'), Tmax=(1416.74,'K')), NASAPolynomial(coeffs=[16.5194,0.0262264,-1.00807e-05,1.77927e-09,-1.19205e-13,50085.8,-55.613], Tmin=(1416.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCtCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = '[CH]=C(C=C)C=[C]C=C(26596)',
    structure = SMILES('[CH]=C(C=C)C=[C]C=C'),
    E0 = (662.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.37975,'amu*angstrom^2'), symmetry=1, barrier=(31.7232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38561,'amu*angstrom^2'), symmetry=1, barrier=(31.8578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37601,'amu*angstrom^2'), symmetry=1, barrier=(31.6371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3608.42,'J/mol'), sigma=(5.96244,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.63 K, Pc=38.63 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371162,0.0693608,-3.99084e-05,-2.74722e-09,7.74147e-12,79839.3,26.692], Tmin=(100,'K'), Tmax=(961.871,'K')), NASAPolynomial(coeffs=[16.7573,0.0252327,-8.54203e-06,1.46909e-09,-1.00661e-13,75576.1,-57.4993], Tmin=(961.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(662.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJC=C)"""),
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
    E0 = (690.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (771.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (751.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (759.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (827.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (844.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (743.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (952.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (823.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (863.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (808.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (891.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1028.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (918.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (920.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (897.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (903.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (916.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (821.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (846.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (750.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (779.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (768.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (850.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (850.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (900.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (698.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (1046.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (1057.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (808.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (867.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (822.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (774.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (734.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (723.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (1012.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (821.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (715.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (739.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['CH2CHCCH(26391)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['[CH]=C1CC1C[CH]C#C(28310)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(80.5299,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for R4_S_T;triplebond_intra_H;radadd_intra_cs2H
Exact match found for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 78.7 to 80.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['[CH]=C1C([CH2])CC1C#C(28235)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(60.6686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_csHCt]
Euclidian distance = 3.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['[CH]=C1[CH]CC(C#C)C1(28303)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(69.2285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;triplebond_intra_H;radadd_intra_cs2H] for rate rule [R6;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#C[CH]CC(=C)C#C(28311)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(15.6063,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2632 used for Cds-CtCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CtCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C#CC=CC([CH2])C#C(28312)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(2.76144,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2560 used for Cds-CsH_Cds-CtH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CtH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]C=C(4699)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00294841,'m^3/(mol*s)'), n=2.48333, Ea=(17.6885,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CtH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2H(33)', 'C#C[CH]CC=C(26789)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CJ] for rate rule [Cds-CsH_Cds-HH;CtJ_Ct]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#C[CH]C[C](C)C#C(28313)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#CC[CH]C([CH2])C#C(28291)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#C[CH][CH]C(C)C#C(28314)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#CCC[C]([CH2])C#C(28315)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.07201e+08,'s^-1'), n=1.25033, Ea=(201.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/OneDe;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_H/Ct;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]#CCCC([CH2])C#C(28316)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[C]#CC(C)C[CH]C#C(28317)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.39293e+07,'s^-1'), n=1.32074, Ea=(96.0416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_TSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[C]#CC([CH2])CCC#C(28318)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(13813.5,'s^-1'), n=1.88327, Ea=(38.7799,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R5H_TSSS;Ct_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[C]#C[CH]CC(C)C#C(28319)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.31883e+06,'s^-1'), n=1.02765, Ea=(75.0925,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R6HJ_2;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]C=C(4699)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#CC([CH2])CC1[C]=C1(28320)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_T;triplebond_intra_H;radadd_intra_cs] for rate rule [R3_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#C[CH]CC1[C]=CC1(28321)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(130.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#CC1C=[C]C([CH2])C1(28251)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.8435e+15,'s^-1'), n=-1.17677, Ea=(156.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHDe] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_csHCt]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#CC1C[CH][C]=CC1(28322)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.33254e+09,'s^-1'), n=0.487896, Ea=(59.5573,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6_linear;multiplebond_intra;radadd_intra_cs2H] + [R6_linear;triplebond_intra_H;radadd_intra] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#CCCC(=C)C#C(28323)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#CC=CC(C)C#C(28324)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#C[CH]C[CH]CC#C(28325)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#CC([CH2])C([CH2])C#C(26592)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.31121e+11,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#CC1CC(C#C)C1(26598)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(19)', 'C#C[CH]C[CH]C#C(28326)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['C3H2(81)', 'C#CC([CH2])[CH2](26629)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.26928e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['[CH]=[C]C1CC(C#C)C1(28327)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.20551e+07,'s^-1'), n=1.225, Ea=(117.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 113.8 to 117.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['[CH]=C1C=C=CCC1[CH2](28052)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.05e+09,'s^-1'), n=0.155, Ea=(177.192,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;triplebond_intra_H;radadd_intra_cdsingleH] for rate rule [R7_MMSR;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#CC([CH2])C[C]=C=C(28328)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C=[C]CC(C)C#C(28329)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#CC([CH2])C=C[C]=C(28257)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#C[C]([CH2])CC=C=C(28330)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[C]#CC([CH2])CC=C=C(28331)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.24948e+06,'s^-1'), n=1.65108, Ea=(132.669,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;Cd_H_out_singleH] + [R7H;Y_rad_out;XH_out] for rate rule [R7H;Ct_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['[CH2]C1[C]=CC=C=CC1(28019)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.07e+10,'s^-1'), n=0.124, Ea=(130.708,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 830 used for R7_linear;triplebond_intra_H;radadd_intra_cdsingleH
Exact match found for rate rule [R7_linear;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#CC(=C)CC=C=C(28332)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.213e+11,'s^-1'), n=0.07, Ea=(76.6885,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for 1_4_5_hexatriene
Exact match found for rate rule [1_4_5_hexatriene]
Euclidian distance = 0
family: 6_membered_central_C-C_shift"""),
)

network(
    label = '4561',
    isomers = [
        'C#C[CH]CC([CH2])C#C(26589)',
    ],
    reactants = [
        ('CH2CHCCH(26391)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4561',
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

