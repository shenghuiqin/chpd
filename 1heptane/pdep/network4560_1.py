species(
    label = 'C#C[CH]CC[CH]C#C(26588)',
    structure = SMILES('C#C[CH]CC[CH]C#C'),
    E0 = (640.767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,3000,3050,390,425,1340,1360,335,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00747775,0.0843511,-9.48732e-05,5.92416e-08,-1.44726e-11,77214,28.5753], Tmin=(100,'K'), Tmax=(1110.58,'K')), NASAPolynomial(coeffs=[13.9285,0.026959,-7.63376e-06,1.03634e-09,-5.62696e-14,74562.5,-38.1205], Tmin=(1110.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(640.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Sec_Propargyl)"""),
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
    label = '[CH]=C1[CH]CCC1C#C(28246)',
    structure = SMILES('[CH]C1=CCCC1C#C'),
    E0 = (566.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.657797,0.0587742,-2.37209e-07,-4.02186e-08,1.91421e-11,68313.1,23.4987], Tmin=(100,'K'), Tmax=(993.964,'K')), NASAPolynomial(coeffs=[15.4165,0.0331881,-1.26438e-05,2.32173e-09,-1.6428e-13,63709.2,-56.0149], Tmin=(993.964,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(566.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopentene) + radical(AllylJ2_triplet)"""),
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
    label = 'C#C[CH]CC=CC#C(28420)',
    structure = SMILES('C#C[CH]CC=CC#C'),
    E0 = (604.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,756.667,763.333,770,3350,3450,2000,2200,2100,2250,500,550,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07912,'amu*angstrom^2'), symmetry=1, barrier=(47.803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06538,'amu*angstrom^2'), symmetry=1, barrier=(47.4873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34331,'amu*angstrom^2'), symmetry=1, barrier=(30.8853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34516,'amu*angstrom^2'), symmetry=1, barrier=(30.9279,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0418034,0.0757532,-7.14199e-05,3.52503e-08,-6.81494e-12,72865.3,28.4607], Tmin=(100,'K'), Tmax=(1334.35,'K')), NASAPolynomial(coeffs=[17.2956,0.0209051,-6.2486e-06,9.33726e-10,-5.65333e-14,68539,-58.7128], Tmin=(1334.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(604.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
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
    label = 'C#C[CH]C[CH]CC#C(28325)',
    structure = SMILES('C#C[CH]C[CH]CC#C'),
    E0 = (689.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,3000,3050,390,425,1340,1360,335,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.426906,0.0817575,-0.000101943,7.58496e-08,-2.2516e-11,82995,30.4421], Tmin=(100,'K'), Tmax=(946.923,'K')), NASAPolynomial(coeffs=[9.04263,0.0354563,-1.2906e-05,2.11592e-09,-1.32437e-13,81807.4,-8.30869], Tmin=(946.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(RCCJCC) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH][CH]CCC#C(28421)',
    structure = SMILES('[CH]=C=C[CH]CCC#C'),
    E0 = (642.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.248492,0.0755006,-6.83921e-05,3.38563e-08,-6.75008e-12,77378.7,28.9722], Tmin=(100,'K'), Tmax=(1213.27,'K')), NASAPolynomial(coeffs=[14.6526,0.0280124,-9.68146e-06,1.59629e-09,-1.02786e-13,73883.5,-43.3039], Tmin=(1213.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(642.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtH) + radical(Allyl_S) + radical(C=C=CJ)"""),
)

species(
    label = '[C]#CCCC[CH]C#C(28422)',
    structure = SMILES('[C]#CCCC[CH]C#C'),
    E0 = (831.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2100,2250,500,550,750,770,3400,2100,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.109183,0.0859056,-0.000104637,7.23027e-08,-1.96511e-11,100170,29.7573], Tmin=(100,'K'), Tmax=(1012.46,'K')), NASAPolynomial(coeffs=[11.9392,0.0305474,-9.85005e-06,1.47982e-09,-8.68067e-14,98216.5,-25.2801], Tmin=(1012.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(831.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Acetyl)"""),
)

species(
    label = '[C]#C[CH]CCCC#C(28423)',
    structure = SMILES('[C]#C[CH]CCCC#C'),
    E0 = (831.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2100,2250,500,550,750,770,3400,2100,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.109183,0.0859056,-0.000104637,7.23027e-08,-1.96511e-11,100170,29.7573], Tmin=(100,'K'), Tmax=(1012.46,'K')), NASAPolynomial(coeffs=[11.9392,0.0305474,-9.85005e-06,1.47982e-09,-8.68067e-14,98216.5,-25.2801], Tmin=(1012.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(831.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH]CCC1[C]=C1(28424)',
    structure = SMILES('C#C[CH]CCC1[C]=C1'),
    E0 = (814.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.255659,0.0828258,-9.2769e-05,5.91232e-08,-1.52776e-11,98069,26.6967], Tmin=(100,'K'), Tmax=(941.457,'K')), NASAPolynomial(coeffs=[12.1484,0.0322965,-1.22617e-05,2.11397e-09,-1.38955e-13,95829.7,-29.9614], Tmin=(941.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(814.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropene) + radical(cyclopropenyl-vinyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC1C=[C][CH]CC1(28266)',
    structure = SMILES('C#CC1C=[C][CH]CC1'),
    E0 = (568.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08751,0.0487853,1.37432e-05,-5.07551e-08,2.21559e-11,68484.4,20.0027], Tmin=(100,'K'), Tmax=(997.08,'K')), NASAPolynomial(coeffs=[14.8019,0.0296376,-1.14145e-05,2.14666e-09,-1.5486e-13,63966.5,-55.0624], Tmin=(997.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(568.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclohexene) + radical(Cds_S) + radical(cyclohexene-allyl)"""),
)

species(
    label = 'C#CC=CCCC#C(28425)',
    structure = SMILES('C#CC=CCCC#C'),
    E0 = (458.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.34182,0.0714313,-5.09124e-05,1.21057e-08,1.4912e-12,55267.7,27.543], Tmin=(100,'K'), Tmax=(1007.59,'K')), NASAPolynomial(coeffs=[15.9546,0.0265907,-9.67412e-06,1.70293e-09,-1.16473e-13,51251.4,-52.2151], Tmin=(1007.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + group(Ct-CtH)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3553.55,'J/mol'), sigma=(6.22586,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=555.06 K, Pc=33.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.413,0.0775862,-6.82273e-05,2.00101e-08,4.75277e-12,83180.4,29.4531], Tmin=(100,'K'), Tmax=(770.958,'K')), NASAPolynomial(coeffs=[14.2582,0.0259108,-6.9066e-06,8.99538e-10,-4.79033e-14,80446.5,-37.6258], Tmin=(770.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(690.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Isobutyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC1CCC1C#C(26599)',
    structure = SMILES('C#CC1CCC1C#C'),
    E0 = (488.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14887,0.0447403,3.30347e-05,-8.13388e-08,3.69471e-11,58844.9,23.0713], Tmin=(100,'K'), Tmax=(922.648,'K')), NASAPolynomial(coeffs=[17.078,0.0220103,-5.33027e-06,8.04052e-10,-5.64045e-14,53933.6,-63.1813], Tmin=(922.648,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + ring(Cyclobutane)"""),
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
    label = '[CH]=C=CC[CH2](26618)',
    structure = SMILES('C#C[CH]C[CH2]'),
    E0 = (477.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.0141456,'amu*angstrom^2'), symmetry=1, barrier=(6.90141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.94797,'amu*angstrom^2'), symmetry=1, barrier=(67.7796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0405762,'amu*angstrom^2'), symmetry=1, barrier=(67.6916,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76795,0.0430915,-3.64272e-05,1.77221e-08,-3.47093e-12,57486.3,19.7996], Tmin=(100,'K'), Tmax=(1346.56,'K')), NASAPolynomial(coeffs=[9.46877,0.0175599,-5.02768e-06,7.11764e-10,-4.08662e-14,55653.2,-18.7496], Tmin=(1346.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(RCCJ)"""),
)

species(
    label = '[CH]=[C]C1CCC1C#C(28426)',
    structure = SMILES('[CH]=[C]C1CCC1C#C'),
    E0 = (808.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01424,0.0497237,1.46207e-05,-5.78946e-08,2.69219e-11,97333.3,27.3153], Tmin=(100,'K'), Tmax=(951.089,'K')), NASAPolynomial(coeffs=[16.3554,0.0245793,-7.82448e-06,1.36849e-09,-9.79124e-14,92634.2,-55.29], Tmin=(951.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(808.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C1[CH]CCC=C=C1(28059)',
    structure = SMILES('[CH]C1C=C=CCCC=1'),
    E0 = (628.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06553,0.047489,2.97765e-05,-6.80472e-08,2.80264e-11,75753.2,20.1183], Tmin=(100,'K'), Tmax=(992.505,'K')), NASAPolynomial(coeffs=[14.154,0.0355073,-1.37281e-05,2.56046e-09,-1.83332e-13,71147.2,-53.0429], Tmin=(992.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(628.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#C[CH]CC[C]=C=C(28427)',
    structure = SMILES('C#C[CH]CCC#C[CH2]'),
    E0 = (632.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.590081,0.0724789,-6.2378e-05,2.54638e-08,-2.13479e-12,76205,29.5311], Tmin=(100,'K'), Tmax=(855.099,'K')), NASAPolynomial(coeffs=[12.2281,0.0303792,-1.01754e-05,1.64229e-09,-1.04662e-13,73763.5,-27.432], Tmin=(855.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]CCCC#C(28428)',
    structure = SMILES('[CH]=C=[C]CCCC#C'),
    E0 = (738.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,1685,370,2175,525,750,770,3400,2100,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,180,180],'cm^-1')),
        HinderedRotor(inertia=(3.73022,'amu*angstrom^2'), symmetry=1, barrier=(85.7651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.494465,'amu*angstrom^2'), symmetry=1, barrier=(11.3687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.494517,'amu*angstrom^2'), symmetry=1, barrier=(11.3699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.73743,'amu*angstrom^2'), symmetry=1, barrier=(85.9309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.270121,0.0823087,-9.14687e-05,5.77078e-08,-1.47568e-11,89005.6,30.4783], Tmin=(100,'K'), Tmax=(951.149,'K')), NASAPolynomial(coeffs=[12.2338,0.0319944,-1.21183e-05,2.08875e-09,-1.37366e-13,86729.8,-26.64], Tmin=(951.149,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(738.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtH) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH]CC=C[C]=C(28337)',
    structure = SMILES('C#C[CH]CC=C[C]=C'),
    E0 = (626.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(3.10668,'amu*angstrom^2'), symmetry=1, barrier=(71.4287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843074,'amu*angstrom^2'), symmetry=1, barrier=(19.3839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.09628,'amu*angstrom^2'), symmetry=1, barrier=(71.1895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.09656,'amu*angstrom^2'), symmetry=1, barrier=(71.196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.26904,0.0773468,-6.6166e-05,2.27803e-08,3.35569e-13,75523.2,28.3503], Tmin=(100,'K'), Tmax=(870.628,'K')), NASAPolynomial(coeffs=[15.2039,0.0261562,-7.9932e-06,1.22521e-09,-7.63792e-14,72262.2,-45.4257], Tmin=(870.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=CJC=C)"""),
)

species(
    label = 'C#C[CH][CH]CC=C=C(28429)',
    structure = SMILES('[CH]=C=C[CH]CC=C=C'),
    E0 = (640.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.330306,0.072092,-5.83686e-05,2.48269e-08,-4.26431e-12,77167.2,29.4586], Tmin=(100,'K'), Tmax=(1387.88,'K')), NASAPolynomial(coeffs=[15.7006,0.0277932,-1.04908e-05,1.82868e-09,-1.21606e-13,72900.8,-49.732], Tmin=(1387.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(640.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Allyl_S) + radical(C=C=CJ)"""),
)

species(
    label = '[C]#C[CH]CCC=C=C(28430)',
    structure = SMILES('[C]#C[CH]CCC=C=C'),
    E0 = (829.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2175,525,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,266.142,266.142,266.144,266.144],'cm^-1')),
        HinderedRotor(inertia=(0.00237995,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130048,'amu*angstrom^2'), symmetry=1, barrier=(6.53675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3889,'amu*angstrom^2'), symmetry=1, barrier=(69.8108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38886,'amu*angstrom^2'), symmetry=1, barrier=(69.8107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.473768,0.0790027,-8.15593e-05,4.49175e-08,-8.54346e-12,99946.5,29.2393], Tmin=(100,'K'), Tmax=(776.108,'K')), NASAPolynomial(coeffs=[11.3542,0.0330141,-1.21737e-05,2.06425e-09,-1.34488e-13,97953.8,-22.452], Tmin=(776.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(829.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Acetyl)"""),
)

species(
    label = '[C]1[CH]CCC=C=CC=1(28027)',
    structure = SMILES('[C]1=CCCC=[C]C=C1'),
    E0 = (559.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06344,0.0459548,3.40192e-05,-8.23559e-08,3.68929e-11,67472.2,17.1279], Tmin=(100,'K'), Tmax=(932.49,'K')), NASAPolynomial(coeffs=[17.2916,0.0237205,-6.42635e-06,1.04594e-09,-7.47385e-14,62385.9,-71.0739], Tmin=(932.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cyclooctatriene) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'C#CC=CCC=C=C(28431)',
    structure = SMILES('C#CC=CCC=C=C'),
    E0 = (457.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.404374,0.0681955,-4.16505e-05,4.349e-09,3.28597e-12,55209.8,28.1216], Tmin=(100,'K'), Tmax=(1065.65,'K')), NASAPolynomial(coeffs=[16.4876,0.0270667,-1.08414e-05,2.01827e-09,-1.42143e-13,50689.5,-55.6194], Tmin=(1065.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = '[CH]=C=[C]CCC=C=C(28432)',
    structure = SMILES('[CH]=C=[C]CCC=C=C'),
    E0 = (737.176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,563.333,586.667,610,1970,2140,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.824278,'amu*angstrom^2'), symmetry=1, barrier=(18.9518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.823687,'amu*angstrom^2'), symmetry=1, barrier=(18.9382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.82455,'amu*angstrom^2'), symmetry=1, barrier=(18.958,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.530648,0.0768795,-7.48119e-05,4.07263e-08,-9.1637e-12,88786.3,30.3181], Tmin=(100,'K'), Tmax=(1060.89,'K')), NASAPolynomial(coeffs=[11.9122,0.033966,-1.41361e-05,2.59731e-09,-1.78524e-13,86371.3,-25.2642], Tmin=(1060.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(737.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C(C=C)C(=[CH])C=C(26597)',
    structure = SMILES('[CH]=C(C=C)C(=[CH])C=C'),
    E0 = (729.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3115,3125,620,680,785,800,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,337.663],'cm^-1')),
        HinderedRotor(inertia=(0.355212,'amu*angstrom^2'), symmetry=1, barrier=(28.7257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279382,'amu*angstrom^2'), symmetry=1, barrier=(22.4526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.998459,'amu*angstrom^2'), symmetry=1, barrier=(80.0405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3547.11,'J/mol'), sigma=(5.9162,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=554.05 K, Pc=38.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.794749,0.0598441,-2.34502e-05,-1.05305e-08,7.82432e-12,87906.3,25.6889], Tmin=(100,'K'), Tmax=(1032.31,'K')), NASAPolynomial(coeffs=[14.0191,0.0301847,-1.17143e-05,2.14326e-09,-1.4969e-13,84026,-44.1016], Tmin=(1032.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
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
    E0 = (640.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (740.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (819.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (728.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (814.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (841.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (978.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (915.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (903.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (866.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (716.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (733.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (850.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (648.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1007.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (808.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (817.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (772.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (783.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (685.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (673.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1012.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (717.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (658.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (982.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (876.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['CH2CHCCH(26391)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['[CH]=C1[CH]CCC1C#C(28246)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.6886e+11,'s^-1'), n=0.474, Ea=(99.6629,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;triplebond_intra_H;radadd_intra_csHDe] for rate rule [R6;triplebond_intra_H;radadd_intra_csHCt]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C#C[CH]CC=CC#C(28420)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(2.76144,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2560 used for Cds-CsH_Cds-CtH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CtH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=[C]C=C(4699)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0132398,'m^3/(mol*s)'), n=2.333, Ea=(2.89605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CtH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C#C[CH]C[CH]CC#C(28325)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['C#C[CH][CH]CCC#C(28421)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.2288e+09,'s^-1'), n=1.25033, Ea=(201.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/OneDe;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_H/Ct;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[C]#CCCC[CH]C#C(28422)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[C]#C[CH]CCCC#C(28423)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(49528.1,'s^-1'), n=1.95205, Ea=(83.4098,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R6HJ_2;Ct_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]C=C(4699)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.73038e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['C#C[CH]CCC1[C]=C1(28424)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.1e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_T;triplebond_intra_H;radadd_intra_cs] for rate rule [R3_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['C#CC1C=[C][CH]CC1(28266)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.01802e+09,'s^-1'), n=0.463766, Ea=(75.7068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHDe] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_csHCt]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['C#CC=CCCC#C(28425)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.16e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#C[CH]CC([CH2])C#C(26589)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['C#CC1CCC1C#C(26599)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] + [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C3H2(81)', '[CH]=C=CC[CH2](26618)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['[CH]=[C]C1CCC1C#C(28426)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.06771e+06,'s^-1'), n=1.35044, Ea=(167.494,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHCt]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 163.3 to 167.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['[CH]=C1[CH]CCC=C=C1(28059)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.1e+09,'s^-1'), n=0.155, Ea=(177.192,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;triplebond_intra_H;radadd_intra_cdsingleH] for rate rule [R8;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['C#C[CH]CC[C]=C=C(28427)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(887826,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C=[C]CCCC#C(28428)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;Cd_rad_out_double;Cs_H_out_H/Ct]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['C#C[CH]CC=C[C]=C(28337)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(148400,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['C#C[CH][CH]CC=C=C(28429)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(544000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[C]#C[CH]CCC=C=C(28430)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(640643,'s^-1'), n=2.07799, Ea=(182.911,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Y_rad_out;Cd_H_out_singleH] for rate rule [R8Hall;Ct_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['[C]1[CH]CCC=C=CC=1(28027)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.53286e+10,'s^-1'), n=0.161, Ea=(76.7659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;triplebond_intra_H;radadd_intra_cdsingleH] for rate rule [R8_linear;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#C[CH]CC[CH]C#C(26588)'],
    products = ['C#CC=CCC=C=C(28431)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.28e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C=[C]CCC=C=C(28432)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0.627, Ea=(245.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_double;Cd_H_out_singleH] for rate rule [R6H;Cd_rad_out_double;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

network(
    label = '4560',
    isomers = [
        'C#C[CH]CC[CH]C#C(26588)',
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
    label = '4560',
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

