species(
    label = '[CH]=C(O)C([CH2])C#C(27563)',
    structure = SMILES('[CH]=C(O)C([CH2])C#C'),
    E0 = (434.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(0.952197,'amu*angstrom^2'), symmetry=1, barrier=(21.8929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.954405,'amu*angstrom^2'), symmetry=1, barrier=(21.9436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.951224,'amu*angstrom^2'), symmetry=1, barrier=(21.8705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.951369,'amu*angstrom^2'), symmetry=1, barrier=(21.8738,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.207352,0.0789663,-8.75161e-05,4.16385e-08,-5.02731e-12,52381.2,25.6027], Tmin=(100,'K'), Tmax=(867.225,'K')), NASAPolynomial(coeffs=[18.6789,0.011797,-2.52066e-06,2.71892e-10,-1.30063e-14,48499.4,-64.7894], Tmin=(867.225,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(434.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = 'HCCOH(50)',
    structure = SMILES('C#CO'),
    E0 = (80.0402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3615,1277.5,1000,287.694,287.7,1588.42],'cm^-1')),
        HinderedRotor(inertia=(0.002037,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05541,0.0252003,-3.80822e-05,3.09891e-08,-9.898e-12,9768.72,12.2272], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.3751,0.00549429,-1.88137e-06,2.93804e-10,-1.71772e-14,8932.78,-8.24498], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(80.0402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH]=C(O)C1CC1=[CH](30123)',
    structure = SMILES('[CH]=C(O)C1CC1=[CH]'),
    E0 = (538.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0883803,0.0727683,-7.35549e-05,3.58483e-08,-6.56726e-12,64902,25.3636], Tmin=(100,'K'), Tmax=(1512.04,'K')), NASAPolynomial(coeffs=[20.5103,0.00921099,-1.51108e-06,1.18854e-10,-4.21386e-15,59709,-79.1036], Tmin=(1512.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C=C(O)C1[CH2](30124)',
    structure = SMILES('[CH]=C1C=C(O)C1[CH2]'),
    E0 = (423.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688545,0.0572939,-1.41919e-05,-4.13942e-08,2.60498e-11,51102.1,21.537], Tmin=(100,'K'), Tmax=(899.769,'K')), NASAPolynomial(coeffs=[22.089,0.00488752,1.93714e-06,-5.62856e-10,3.94519e-14,45521.3,-89.0599], Tmin=(899.769,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(Isobutyl)"""),
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
    label = '[CH]=C(O)C(=C)C#C(30125)',
    structure = SMILES('[CH]=C(O)C(=C)C#C'),
    E0 = (356.611,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,3615,1277.5,1000,2175,525,750,770,3400,2100,325,375,415,465,420,450,1700,1750],'cm^-1')),
        HinderedRotor(inertia=(1.30009,'amu*angstrom^2'), symmetry=1, barrier=(29.8916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30123,'amu*angstrom^2'), symmetry=1, barrier=(29.9177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30112,'amu*angstrom^2'), symmetry=1, barrier=(29.9154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.361119,0.0768299,-8.39438e-05,4.26397e-08,-8.0043e-12,43063.8,24.8478], Tmin=(100,'K'), Tmax=(1506.02,'K')), NASAPolynomial(coeffs=[22.4893,0.00388216,9.20262e-07,-3.31018e-10,2.59341e-14,37571.1,-90.134], Tmin=(1506.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)Ct) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]O(172)',
    structure = SMILES('[CH]=[C]O'),
    E0 = (348.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,1058.91,1059.8],'cm^-1')),
        HinderedRotor(inertia=(0.315045,'amu*angstrom^2'), symmetry=1, barrier=(7.24351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22339,0.0198737,-3.18466e-05,3.10436e-08,-1.17427e-11,41917.3,13.7606], Tmin=(100,'K'), Tmax=(829.31,'K')), NASAPolynomial(coeffs=[3.78101,0.0111997,-5.33323e-06,1.02849e-09,-7.13893e-14,42030.6,12.4155], Tmin=(829.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(Cds_P) + radical(C=CJO)"""),
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
    label = '[CH]=C(O)C=C(6273)',
    structure = SMILES('[CH]=C(O)C=C'),
    E0 = (127.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.06241,'amu*angstrom^2'), symmetry=1, barrier=(24.427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06165,'amu*angstrom^2'), symmetry=1, barrier=(24.4094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55553,0.0432866,-1.7698e-05,-2.104e-08,1.52122e-11,15387.4,16.4164], Tmin=(100,'K'), Tmax=(907.962,'K')), NASAPolynomial(coeffs=[16.9959,0.00364253,9.14579e-07,-2.83671e-10,1.91369e-14,11413.8,-63.025], Tmin=(907.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
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
    label = 'OH(5)',
    structure = SMILES('[OH]'),
    E0 = (28.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3287.46],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0073,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4858,0.00133397,-4.70043e-06,5.64379e-09,-2.06318e-12,3411.96,1.99788], Tmin=(100,'K'), Tmax=(1005.25,'K')), NASAPolynomial(coeffs=[2.88225,0.00103869,-2.35652e-07,1.40229e-11,6.34581e-16,3669.56,5.59053], Tmin=(1005.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.372,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C#CC([CH2])C#C(27501)',
    structure = SMILES('C#CC([CH2])C#C'),
    E0 = (605.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(1.10375,'amu*angstrom^2'), symmetry=1, barrier=(60.786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.494061,'amu*angstrom^2'), symmetry=1, barrier=(27.2294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10575,'amu*angstrom^2'), symmetry=1, barrier=(60.6988,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34005,0.0508861,-4.74747e-05,2.33958e-08,-4.52454e-12,72935.7,20.5937], Tmin=(100,'K'), Tmax=(1340.28,'K')), NASAPolynomial(coeffs=[12.7012,0.0146548,-4.32414e-06,6.38341e-10,-3.82616e-14,70099,-36.7658], Tmin=(1340.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCtCsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C(O)[C](C)C#C(30126)',
    structure = SMILES('[CH]C(O)=C(C)C#C'),
    E0 = (350.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.490323,0.0696523,-5.71893e-05,1.86414e-08,-2.71193e-13,42270.7,24.0314], Tmin=(100,'K'), Tmax=(972.273,'K')), NASAPolynomial(coeffs=[15.9261,0.0218565,-7.68513e-06,1.3141e-09,-8.84829e-14,38526.7,-53.8212], Tmin=(972.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsOs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#CC([CH2])C(=C)[O](26501)',
    structure = SMILES('C#CC([CH2])C(=C)[O]'),
    E0 = (325.059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,310.359,311.072],'cm^-1')),
        HinderedRotor(inertia=(0.203268,'amu*angstrom^2'), symmetry=1, barrier=(13.9207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202694,'amu*angstrom^2'), symmetry=1, barrier=(13.923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29526,'amu*angstrom^2'), symmetry=1, barrier=(89.006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3830.04,'J/mol'), sigma=(6.36116,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=598.24 K, Pc=33.76 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.647829,0.0714799,-7.85135e-05,4.17373e-08,-7.43084e-12,39218.5,25.1189], Tmin=(100,'K'), Tmax=(863.516,'K')), NASAPolynomial(coeffs=[14.6554,0.0180653,-5.6555e-06,8.73045e-10,-5.41924e-14,36371.6,-42.8808], Tmin=(863.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = 'C#C[C]([CH2])C(=C)O(27562)',
    structure = SMILES('[CH]=C=C([CH2])C(=C)O'),
    E0 = (288.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.550383,0.0822665,-9.14748e-05,4.84112e-08,-9.48745e-12,34936.9,25.9359], Tmin=(100,'K'), Tmax=(1448.74,'K')), NASAPolynomial(coeffs=[21.9961,0.00662962,6.9751e-07,-3.813e-10,3.29097e-14,29808.8,-86.3477], Tmin=(1448.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C([O])C(C)C#C(27564)',
    structure = SMILES('[CH]=C([O])C(C)C#C'),
    E0 = (367.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,350,440,435,1725,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.819119,'amu*angstrom^2'), symmetry=1, barrier=(18.8331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.820063,'amu*angstrom^2'), symmetry=1, barrier=(18.8549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.817634,'amu*angstrom^2'), symmetry=1, barrier=(18.799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.506351,0.073971,-8.53075e-05,5.06325e-08,-1.17541e-11,44277,24.2439], Tmin=(100,'K'), Tmax=(1059.7,'K')), NASAPolynomial(coeffs=[15.2367,0.01837,-6.60585e-06,1.12137e-09,-7.37719e-14,41155,-47.6757], Tmin=(1059.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[C]#CC(C)C(=[CH])O(30127)',
    structure = SMILES('[C]#CC(C)C(=[CH])O'),
    E0 = (566.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.833875,'amu*angstrom^2'), symmetry=1, barrier=(19.1724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.832703,'amu*angstrom^2'), symmetry=1, barrier=(19.1455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834225,'amu*angstrom^2'), symmetry=1, barrier=(19.1805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.83334,'amu*angstrom^2'), symmetry=1, barrier=(19.1601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.228054,0.0806592,-9.71523e-05,5.57484e-08,-1.12624e-11,68261.7,24.4879], Tmin=(100,'K'), Tmax=(867.726,'K')), NASAPolynomial(coeffs=[17.132,0.0148847,-4.45167e-06,6.62037e-10,-4.00788e-14,64870.7,-57.3011], Tmin=(867.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(566.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(Acetyl)"""),
)

species(
    label = '[C]#CC([CH2])C(=C)O(27566)',
    structure = SMILES('[C]#CC([CH2])C(=C)O'),
    E0 = (524.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,250.47,250.472],'cm^-1')),
        HinderedRotor(inertia=(0.345811,'amu*angstrom^2'), symmetry=1, barrier=(15.3952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64352,'amu*angstrom^2'), symmetry=1, barrier=(73.1765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.345901,'amu*angstrom^2'), symmetry=1, barrier=(15.3951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.345879,'amu*angstrom^2'), symmetry=1, barrier=(15.3953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0771888,0.0818293,-0.000104219,6.64705e-08,-1.6102e-11,63215.7,26.3978], Tmin=(100,'K'), Tmax=(1119.09,'K')), NASAPolynomial(coeffs=[17.1404,0.013531,-2.87733e-06,2.63037e-10,-7.79849e-15,59854.2,-55.7975], Tmin=(1119.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C(O)C1[C]=CC1(30121)',
    structure = SMILES('[CH]=C(O)C1[C]=CC1'),
    E0 = (501.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.955597,0.0554075,-2.54453e-05,-1.52825e-08,1.22943e-11,60399.7,22.8002], Tmin=(100,'K'), Tmax=(950.901,'K')), NASAPolynomial(coeffs=[17.4939,0.0138456,-4.06298e-06,7.00343e-10,-5.10434e-14,55988.2,-62.8128], Tmin=(950.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH2]C1[C]=CC=C1O(30128)',
    structure = SMILES('[CH2]C1[C]=CC=C1O'),
    E0 = (350.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07219,0.0512362,-1.00603e-05,-3.48494e-08,2.06727e-11,42222.5,20.8627], Tmin=(100,'K'), Tmax=(918.367,'K')), NASAPolynomial(coeffs=[18.0459,0.0118138,-2.0331e-06,2.38547e-10,-1.70045e-14,37649.7,-67.503], Tmin=(918.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclopentadiene) + radical(1,3-cyclopentadiene-vinyl-1) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC(=C)C(=C)O(27568)',
    structure = SMILES('C#CC(=C)C(=C)O'),
    E0 = (109.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.464813,0.0645515,-3.7141e-05,-1.32674e-08,1.41409e-11,13311.1,21.4732], Tmin=(100,'K'), Tmax=(929.188,'K')), NASAPolynomial(coeffs=[21.3669,0.00842017,-1.16958e-06,1.28263e-10,-1.11712e-14,7965.45,-85.6949], Tmin=(929.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)Ct) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = '[CH]=C(O)C[CH]C#C(27515)',
    structure = SMILES('[CH]=C(O)C[CH]C#C'),
    E0 = (408.551,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3615,1277.5,1000,2175,525,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,280.857],'cm^-1')),
        HinderedRotor(inertia=(0.35025,'amu*angstrom^2'), symmetry=1, barrier=(19.2629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.355583,'amu*angstrom^2'), symmetry=1, barrier=(19.2677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.354969,'amu*angstrom^2'), symmetry=1, barrier=(19.2689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35889,'amu*angstrom^2'), symmetry=1, barrier=(73.6604,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0294285,0.0757923,-8.45921e-05,4.65376e-08,-9.57452e-12,49293.5,27.6378], Tmin=(100,'K'), Tmax=(1369.34,'K')), NASAPolynomial(coeffs=[18.6608,0.0102105,-7.19147e-07,-1.54675e-10,1.96469e-14,45204.7,-64.6462], Tmin=(1369.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=C(O)[CH]CC#C(30129)',
    structure = SMILES('[CH]C(O)=CCC#C'),
    E0 = (378.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.60563,0.0654853,-4.38111e-05,3.79783e-09,5.22203e-12,45621.1,25.1453], Tmin=(100,'K'), Tmax=(958.429,'K')), NASAPolynomial(coeffs=[16.2982,0.021037,-7.18297e-06,1.22969e-09,-8.39544e-14,41646.5,-54.9383], Tmin=(958.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#CC1CC=C1O(30130)',
    structure = SMILES('C#CC1CC=C1O'),
    E0 = (143.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.948863,0.0508256,-1.39399e-06,-4.79524e-08,2.59311e-11,17343.2,19.281], Tmin=(100,'K'), Tmax=(929.125,'K')), NASAPolynomial(coeffs=[20.412,0.00878939,-9.40384e-07,9.05591e-11,-1.03644e-14,11924.2,-82.8857], Tmin=(929.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene)"""),
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
    label = '[CH]=C=CC(=[CH])O(30122)',
    structure = SMILES('[CH]C(O)=CC#C'),
    E0 = (390.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,750,770,3400,2100,2175,525,3615,1277.5,1000,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.01858,'amu*angstrom^2'), symmetry=1, barrier=(46.4111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01877,'amu*angstrom^2'), symmetry=1, barrier=(46.4155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02012,'amu*angstrom^2'), symmetry=1, barrier=(46.4466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17377,0.0508357,-2.10838e-05,-1.91953e-08,1.39694e-11,47115.5,20.0058], Tmin=(100,'K'), Tmax=(934.364,'K')), NASAPolynomial(coeffs=[17.1724,0.0108897,-2.77925e-06,4.39366e-10,-3.1962e-14,42879.8,-62.7598], Tmin=(934.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
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
    E0 = (434.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (538.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (487.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (586.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (640.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (697.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (554.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (633.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (532.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (626.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (626.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (518.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (662.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (598.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (799.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (565.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (484.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (497.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (594.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (594.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (442.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (772.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['HCCOH(50)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['[CH]=C(O)C1CC1=[CH](30123)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(103.929,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for R4_S_T;triplebond_intra_H;radadd_intra_cs2H
Exact match found for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 102.0 to 103.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['[CH]=C1C=C(O)C1[CH2](30124)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH]=C(O)C(=C)C#C(30125)'],
    products = ['[CH]=C(O)C([CH2])C#C(27563)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(454.257,'m^3/(mol*s)'), n=1.51715, Ea=(17.7318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-TwoDe_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=[C]O(172)', 'CH2CHCCH(26391)'],
    products = ['[CH]=C(O)C([CH2])C#C(27563)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00294841,'m^3/(mol*s)'), n=2.48333, Ea=(17.6885,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CtH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H(33)', '[CH]=C(O)C=C(6273)'],
    products = ['[CH]=C(O)C([CH2])C#C(27563)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00505363,'m^3/(mol*s)'), n=2.41967, Ea=(13.2813,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeH_Cds;CJ] for rate rule [Cds-OneDeH_Cds;CtJ_Ct]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['HCCOH(50)', '[CH]=[C]C=C(4699)'],
    products = ['[CH]=C(O)C([CH2])C#C(27563)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;CJ] for rate rule [Ct-O_Ct;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OH(5)', 'C#CC([CH2])C#C(27501)'],
    products = ['[CH]=C(O)C([CH2])C#C(27563)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.296e+06,'cm^3/(mol*s)'), n=1.876, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 212 used for Ct-Cs_Ct-H;OJ_pri
Exact match found for rate rule [Ct-Cs_Ct-H;OJ_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -1.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['[CH]=C(O)[C](C)C#C(30126)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['C#C[C]([CH2])C(=C)O(27562)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['[CH]=C([O])C(C)C#C(27564)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]#CC(C)C(=[CH])O(30127)'],
    products = ['[CH]=C(O)C([CH2])C#C(27563)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.39293e+07,'s^-1'), n=1.32074, Ea=(96.0416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_TSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[C]#CC([CH2])C(=C)O(27566)'],
    products = ['[CH]=C(O)C([CH2])C#C(27563)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(544115,'s^-1'), n=1.7475, Ea=(73.8225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cd_H_out_singleH] for rate rule [R5H_TSSD;Ct_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]O(172)', '[CH]=[C]C=C(4699)'],
    products = ['[CH]=C(O)C([CH2])C#C(27563)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['[CH]=C(O)C1[C]=CC1(30121)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(130.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['[CH2]C1[C]=CC=C1O(30128)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.88e+10,'s^-1'), n=0.31, Ea=(50.6264,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleH
Exact match found for rate rule [R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['C#CC(=C)C(=C)O(27568)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['[CH]=C(O)C[CH]C#C(27515)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['[CH]=C(O)[CH]CC#C(30129)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['C#CC1CC=C1O(30130)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(19)', '[CH]=C=CC(=[CH])O(30122)'],
    products = ['[CH]=C(O)C([CH2])C#C(27563)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/TwoDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4697',
    isomers = [
        '[CH]=C(O)C([CH2])C#C(27563)',
    ],
    reactants = [
        ('HCCOH(50)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4697',
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

