species(
    label = 'C#C[CH]CC=[C]O(29677)',
    structure = SMILES('C#C[CH]CC=[C]O'),
    E0 = (406.965,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,290.846,290.862],'cm^-1')),
        HinderedRotor(inertia=(0.257311,'amu*angstrom^2'), symmetry=1, barrier=(15.4611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4935,'amu*angstrom^2'), symmetry=1, barrier=(89.6655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257411,'amu*angstrom^2'), symmetry=1, barrier=(15.462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49234,'amu*angstrom^2'), symmetry=1, barrier=(89.6678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.227396,0.0709476,-7.63671e-05,4.15004e-08,-8.48798e-12,49092.8,29.2436], Tmin=(100,'K'), Tmax=(1376.98,'K')), NASAPolynomial(coeffs=[16.9307,0.0123705,-1.60267e-06,-2.61708e-13,9.61192e-15,45446.1,-53.2217], Tmin=(1376.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=CJO)"""),
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
    label = '[CH]=C1[CH]CC=C1O(30106)',
    structure = SMILES('[CH]C1=CCC=C1O'),
    E0 = (241.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10895,0.0448339,2.61773e-05,-7.41806e-08,3.3969e-11,29158.7,20.4765], Tmin=(100,'K'), Tmax=(941.438,'K')), NASAPolynomial(coeffs=[18.7453,0.0170277,-4.60749e-06,7.9229e-10,-6.03094e-14,23749.5,-74.6364], Tmin=(941.438,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentadiene) + radical(AllylJ2_triplet)"""),
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
    label = 'C#C[CH]CC=C=O(27375)',
    structure = SMILES('C#C[CH]CC=C=O'),
    E0 = (290.852,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,328.462,328.468],'cm^-1')),
        HinderedRotor(inertia=(0.0734386,'amu*angstrom^2'), symmetry=1, barrier=(5.62232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.917907,'amu*angstrom^2'), symmetry=1, barrier=(70.2747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.917891,'amu*angstrom^2'), symmetry=1, barrier=(70.2746,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22884,0.0597315,-6.52525e-05,4.00127e-08,-9.8572e-12,35082.2,23.4039], Tmin=(100,'K'), Tmax=(990.974,'K')), NASAPolynomial(coeffs=[10.6571,0.0216744,-7.6463e-06,1.2584e-09,-8.02964e-14,33213.6,-21.9967], Tmin=(990.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CtCsHH) + group(Cds-(Cdd-O2d)CsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC=CC=[C]O(30178)',
    structure = SMILES('C#CC=CC=[C]O'),
    E0 = (353.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,750,770,3400,2100,2175,525,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.25832,'amu*angstrom^2'), symmetry=1, barrier=(28.9313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25659,'amu*angstrom^2'), symmetry=1, barrier=(28.8915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25705,'amu*angstrom^2'), symmetry=1, barrier=(28.9021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.702802,0.0585491,-2.73274e-05,-2.37329e-08,1.84169e-11,42639.6,24.3479], Tmin=(100,'K'), Tmax=(916.768,'K')), NASAPolynomial(coeffs=[21.6328,0.00423461,9.91569e-07,-2.95234e-10,1.83462e-14,37246.9,-83.2901], Tmin=(916.768,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJO)"""),
)

species(
    label = 'C#C[CH]CC#CO(30179)',
    structure = SMILES('C#C[CH]CC#CO'),
    E0 = (400.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2100,2175,2250,500,525,550,3615,1277.5,1000,3025,407.5,1350,352.5,301.986,301.99],'cm^-1')),
        HinderedRotor(inertia=(0.994076,'amu*angstrom^2'), symmetry=1, barrier=(64.3321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00184854,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.994102,'amu*angstrom^2'), symmetry=1, barrier=(64.3319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.994094,'amu*angstrom^2'), symmetry=1, barrier=(64.3321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33302,0.0578707,-5.51424e-05,2.23638e-08,-4.70314e-13,48284.2,23.487], Tmin=(100,'K'), Tmax=(788.535,'K')), NASAPolynomial(coeffs=[11.3929,0.0193919,-5.82236e-06,8.52802e-10,-5.03962e-14,46307.5,-25.1303], Tmin=(788.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CtH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtOs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
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
    label = 'C#C[CH]C[CH]C=O(27376)',
    structure = SMILES('C#C[CH]CC=C[O]'),
    E0 = (308.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,288.093,288.102,288.117],'cm^-1')),
        HinderedRotor(inertia=(0.583639,'amu*angstrom^2'), symmetry=1, barrier=(34.3834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.474064,'amu*angstrom^2'), symmetry=1, barrier=(27.9184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1292,'amu*angstrom^2'), symmetry=1, barrier=(66.4958,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0647921,0.0677002,-6.37812e-05,3.0016e-08,-5.2724e-12,37291.1,28.5977], Tmin=(100,'K'), Tmax=(1651.61,'K')), NASAPolynomial(coeffs=[17.2967,0.0123796,-1.48396e-06,4.11297e-12,6.9518e-15,33366.6,-58.3917], Tmin=(1651.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=COJ)"""),
)

species(
    label = 'C#CC[CH]C=[C]O(30177)',
    structure = SMILES('C#CC[CH]C=[C]O'),
    E0 = (401.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.244516,0.0668292,-6.4775e-05,3.15075e-08,-5.79981e-12,48482.2,28.2507], Tmin=(100,'K'), Tmax=(1521,'K')), NASAPolynomial(coeffs=[17.3844,0.0126407,-2.34724e-06,2.05666e-10,-7.37656e-15,44322.4,-58.1621], Tmin=(1521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO) + radical(Allyl_S)"""),
)

species(
    label = 'C#C[CH]C[C]=CO(30180)',
    structure = SMILES('C#C[CH]C[C]=CO'),
    E0 = (405.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.997599,'amu*angstrom^2'), symmetry=1, barrier=(22.9368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.997379,'amu*angstrom^2'), symmetry=1, barrier=(22.9317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.997833,'amu*angstrom^2'), symmetry=1, barrier=(22.9421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.55705,'amu*angstrom^2'), symmetry=1, barrier=(58.7917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.259851,0.0766703,-8.43788e-05,4.52067e-08,-8.94628e-12,48885.8,28.3377], Tmin=(100,'K'), Tmax=(1459.4,'K')), NASAPolynomial(coeffs=[19.3254,0.00848098,6.19949e-07,-4.33581e-10,3.8988e-14,44714.3,-68.2601], Tmin=(1459.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(405.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CCC[C]=[C]O(30181)',
    structure = SMILES('C#CCC[C]=[C]O'),
    E0 = (498.597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2175,525,3615,1277.5,1000,1670,1700,300,440,211.994,211.995],'cm^-1')),
        HinderedRotor(inertia=(0.48375,'amu*angstrom^2'), symmetry=1, barrier=(15.4274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48375,'amu*angstrom^2'), symmetry=1, barrier=(15.4274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.36009,'amu*angstrom^2'), symmetry=1, barrier=(107.158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48375,'amu*angstrom^2'), symmetry=1, barrier=(15.4274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.724781,0.0683187,-6.95465e-05,3.18e-08,-3.73453e-12,60088.8,28.104], Tmin=(100,'K'), Tmax=(876.952,'K')), NASAPolynomial(coeffs=[15.1455,0.016497,-4.77587e-06,7.06139e-10,-4.31902e-14,57023,-42.6333], Tmin=(876.952,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH][CH]C=CO(30182)',
    structure = SMILES('[CH]=C=C[CH]C=CO'),
    E0 = (276.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.915359,0.0477128,1.94657e-05,-8.03857e-08,4.11262e-11,33413.9,24.9138], Tmin=(100,'K'), Tmax=(897.161,'K')), NASAPolynomial(coeffs=[23.3913,0.00220392,4.09798e-06,-1.0069e-09,6.94772e-14,27179.6,-93.3495], Tmin=(897.161,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJC=C)"""),
)

species(
    label = '[C]#CCCC=[C]O(30183)',
    structure = SMILES('[C]#CCCC=[C]O'),
    E0 = (597.899,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2175,525,1685,370,3615,1277.5,1000,3010,987.5,1337.5,450,1655,306.177,314.217,321.416],'cm^-1')),
        HinderedRotor(inertia=(1.09706,'amu*angstrom^2'), symmetry=1, barrier=(76.045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182445,'amu*angstrom^2'), symmetry=1, barrier=(12.5879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179212,'amu*angstrom^2'), symmetry=1, barrier=(12.5806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174832,'amu*angstrom^2'), symmetry=1, barrier=(12.603,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.826086,0.0666472,-6.47444e-05,2.527e-08,-3.57714e-13,72028.1,28.0137], Tmin=(100,'K'), Tmax=(834.991,'K')), NASAPolynomial(coeffs=[14.749,0.0164325,-4.14764e-06,5.29873e-10,-2.86491e-14,69128.4,-40.0864], Tmin=(834.991,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(597.899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO) + radical(Acetyl)"""),
)

species(
    label = 'C#CCC[CH][C]=O(27377)',
    structure = SMILES('C#CCC[CH][C]=O'),
    E0 = (346.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2175,525,1855,455,950,3025,407.5,1350,352.5,366.456,366.46],'cm^-1')),
        HinderedRotor(inertia=(0.125201,'amu*angstrom^2'), symmetry=1, barrier=(11.9309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125203,'amu*angstrom^2'), symmetry=1, barrier=(11.9309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0179503,'amu*angstrom^2'), symmetry=1, barrier=(73.1117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.767218,'amu*angstrom^2'), symmetry=1, barrier=(73.1117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09492,0.062409,-6.55034e-05,3.86681e-08,-9.21788e-12,41805.7,25.4481], Tmin=(100,'K'), Tmax=(1020.24,'K')), NASAPolynomial(coeffs=[10.9832,0.0236419,-8.50803e-06,1.42609e-09,-9.23589e-14,39787.9,-22.4554], Tmin=(1020.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCCJ=O) + radical(CCJCHO)"""),
)

species(
    label = '[C]#C[CH]CC=CO(30184)',
    structure = SMILES('[C]#C[CH]CC=CO'),
    E0 = (504.365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2175,525,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,233.329,233.329,233.329],'cm^-1')),
        HinderedRotor(inertia=(2.03809,'amu*angstrom^2'), symmetry=1, barrier=(78.7388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.461589,'amu*angstrom^2'), symmetry=1, barrier=(17.8329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.46159,'amu*angstrom^2'), symmetry=1, barrier=(17.8329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03809,'amu*angstrom^2'), symmetry=1, barrier=(78.7388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.238248,0.0760668,-8.39879e-05,4.55392e-08,-9.08876e-12,60828.4,28.5254], Tmin=(100,'K'), Tmax=(1463.12,'K')), NASAPolynomial(coeffs=[18.4453,0.00914363,8.6589e-07,-5.25657e-10,4.6915e-14,57057.1,-62.927], Tmin=(1463.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'O[C]=CCC1[C]=C1(30185)',
    structure = SMILES('O[C]=CCC1[C]=C1'),
    E0 = (580.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.564752,0.0685911,-7.1507e-05,3.79199e-08,-7.83351e-12,69944.4,26.4021], Tmin=(100,'K'), Tmax=(1191.04,'K')), NASAPolynomial(coeffs=[16.2127,0.0160368,-5.31725e-06,8.69722e-10,-5.63294e-14,66217.1,-51.8252], Tmin=(1191.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopropene) + radical(cyclopropenyl-vinyl) + radical(C=CJO)"""),
)

species(
    label = 'OC1C=[C][CH]CC=1(30115)',
    structure = SMILES('OC1C=[C][CH]CC=1'),
    E0 = (209.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25233,0.0521136,-3.11029e-05,4.45427e-09,1.49225e-12,25268,19.0342], Tmin=(100,'K'), Tmax=(1134.35,'K')), NASAPolynomial(coeffs=[13.1754,0.0229422,-9.54994e-06,1.79124e-09,-1.25607e-13,21734.9,-43.641], Tmin=(1134.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Aromatic_pi_S_1_3) + radical(Cds_S)"""),
)

species(
    label = 'C#CC=CC=CO(30186)',
    structure = SMILES('C#CC=CC=CO'),
    E0 = (113.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.487413,0.0568503,-6.10814e-07,-6.23874e-08,3.46336e-11,13818.9,21.9927], Tmin=(100,'K'), Tmax=(913.374,'K')), NASAPolynomial(coeffs=[25.7977,6.67753e-05,3.86275e-06,-8.52696e-10,5.45476e-14,6940.35,-110.166], Tmin=(913.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CCCC=C=O(27384)',
    structure = SMILES('C#CCCC=C=O'),
    E0 = (144.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30526,0.0579761,-5.3373e-05,2.75257e-08,-5.85501e-12,17494.5,23.2928], Tmin=(100,'K'), Tmax=(1121.84,'K')), NASAPolynomial(coeffs=[10.3865,0.0255961,-1.00777e-05,1.79676e-09,-1.21338e-13,15457,-21.5629], Tmin=(1121.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CtCsHH) + group(Cds-(Cdd-O2d)CsH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC([CH2])C=[C]O(29678)',
    structure = SMILES('C#CC([CH2])C=[C]O'),
    E0 = (432.764,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.924649,'amu*angstrom^2'), symmetry=1, barrier=(21.2595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.921977,'amu*angstrom^2'), symmetry=1, barrier=(21.1981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.924042,'amu*angstrom^2'), symmetry=1, barrier=(21.2455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92162,'amu*angstrom^2'), symmetry=1, barrier=(21.1899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3889.19,'J/mol'), sigma=(6.40612,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=607.48 K, Pc=33.57 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.454518,0.0742274,-7.96125e-05,3.69258e-08,-4.02564e-12,52180.9,27.2436], Tmin=(100,'K'), Tmax=(857.207,'K')), NASAPolynomial(coeffs=[16.9634,0.0139291,-3.38705e-06,4.22091e-10,-2.26821e-14,48735.6,-53.4457], Tmin=(857.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(C=CJO)"""),
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
    label = '[CH2]C=[C]O(6303)',
    structure = SMILES('[CH2]C=[C]O'),
    E0 = (224.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.206965,'amu*angstrom^2'), symmetry=1, barrier=(4.75853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.69007,'amu*angstrom^2'), symmetry=1, barrier=(15.8661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.71481,0.0306825,-3.73315e-05,3.08792e-08,-1.09328e-11,27020.3,15.471], Tmin=(100,'K'), Tmax=(771.374,'K')), NASAPolynomial(coeffs=[4.68419,0.0180544,-8.07735e-06,1.53601e-09,-1.0692e-13,26788.3,6.94699], Tmin=(771.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=CJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=[C]C1CC=C1O(30187)',
    structure = SMILES('[CH]=[C]C1CC=C1O'),
    E0 = (486.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.955597,0.0554075,-2.54453e-05,-1.52825e-08,1.22943e-11,58638.4,22.8002], Tmin=(100,'K'), Tmax=(950.901,'K')), NASAPolynomial(coeffs=[17.4939,0.0138456,-4.06298e-06,7.00343e-10,-5.10434e-14,54226.9,-62.8128], Tmin=(950.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(486.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=C=[C]CC=[C]O(30188)',
    structure = SMILES('[CH2]C#CCC=[C]O'),
    E0 = (403.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05447,0.0565364,-3.74734e-05,2.34506e-09,4.99378e-12,48620.2,28.027], Tmin=(100,'K'), Tmax=(958.247,'K')), NASAPolynomial(coeffs=[15.0723,0.0171026,-5.61379e-06,9.59722e-10,-6.61209e-14,45057.6,-43.5743], Tmin=(958.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=CJO) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]CC=CO(30189)',
    structure = SMILES('[CH]=C=[C]CC=CO'),
    E0 = (412.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(1.0605,'amu*angstrom^2'), symmetry=1, barrier=(24.383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05781,'amu*angstrom^2'), symmetry=1, barrier=(24.321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06062,'amu*angstrom^2'), symmetry=1, barrier=(24.3857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.291434,0.0747409,-7.79872e-05,3.91973e-08,-7.30235e-12,49826.3,30.0512], Tmin=(100,'K'), Tmax=(1537.92,'K')), NASAPolynomial(coeffs=[20.2516,0.00802156,4.78583e-08,-2.48154e-10,2.30885e-14,45079.1,-72.7904], Tmin=(1537.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C]C=CC=[C]O(30135)',
    structure = SMILES('C=C=C[CH]C=[C]O'),
    E0 = (362.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12076,0.0484167,4.43351e-07,-4.68246e-08,2.50554e-11,43656.6,26.6762], Tmin=(100,'K'), Tmax=(922.347,'K')), NASAPolynomial(coeffs=[18.6637,0.011029,-1.67744e-06,1.89215e-10,-1.50551e-14,38774.7,-65.462], Tmin=(922.347,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJC=C) + radical(C=CJO)"""),
)

species(
    label = 'C=C=CC[C]=[C]O(30190)',
    structure = SMILES('C=C=CC[C]=[C]O'),
    E0 = (498.122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1670,1700,300,440,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.792482,'amu*angstrom^2'), symmetry=1, barrier=(18.2207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.793374,'amu*angstrom^2'), symmetry=1, barrier=(18.2412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792406,'amu*angstrom^2'), symmetry=1, barrier=(18.219,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.604151,0.0673297,-6.85762e-05,3.54844e-08,-7.16199e-12,60038.7,29.3342], Tmin=(100,'K'), Tmax=(1217.18,'K')), NASAPolynomial(coeffs=[16.1942,0.0160959,-5.43721e-06,9.01931e-10,-5.894e-14,56243.6,-48.9427], Tmin=(1217.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=C=CC[CH][C]=O(27388)',
    structure = SMILES('C=C=CC[CH][C]=O'),
    E0 = (344.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3025,407.5,1350,352.5,1855,455,950,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,368.703,371.035],'cm^-1')),
        HinderedRotor(inertia=(0.0922048,'amu*angstrom^2'), symmetry=1, barrier=(8.91028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0882183,'amu*angstrom^2'), symmetry=1, barrier=(8.89708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270103,'amu*angstrom^2'), symmetry=1, barrier=(26.3334,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29966,0.0576537,-5.12448e-05,2.4798e-08,-4.93281e-12,41588.6,25.4865], Tmin=(100,'K'), Tmax=(1193.79,'K')), NASAPolynomial(coeffs=[11.0631,0.0249405,-1.01416e-05,1.84455e-09,-1.26084e-13,39257.5,-23.3459], Tmin=(1193.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C=CCC=C=O(27389)',
    structure = SMILES('C=C=CCC=C=O'),
    E0 = (142.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4076,0.0601413,-5.99076e-05,3.42705e-08,-8.29957e-12,17271,22.2613], Tmin=(100,'K'), Tmax=(975.765,'K')), NASAPolynomial(coeffs=[8.76621,0.0299753,-1.35338e-05,2.58627e-09,-1.81627e-13,15835,-13.059], Tmin=(975.765,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    E0 = (406.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (422.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (533.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (573.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (627.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (625.394,'kJ/mol'),
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
    E0 = (569.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (533.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (639.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (623.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (557.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (745.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (558.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (692.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (799.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (632.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (429.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (470.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (446.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (592.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (414.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (754.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (486.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (539.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (807.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (451.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (609.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (558.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (431.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['HCCOH(50)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['[CH]=C1[CH]CC=C1O(30106)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.42e+11,'s^-1'), n=0.258, Ea=(15.8866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;triplebond_intra_H;radadd_intra_cdsingle] for rate rule [R6;triplebond_intra_H;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C#C[CH]CC=C=O(27375)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.185e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""Estimated using template [Od_R;HJ] for rate rule [Od_Cdd;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#CC=CC=[C]O(30178)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(105.863,'m^3/(mol*s)'), n=1.66278, Ea=(8.08939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#C[CH]CC#CO(30179)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]O(172)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0132398,'m^3/(mol*s)'), n=2.333, Ea=(2.89605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CtH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['HCCOH(50)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['C#C[CH]C[CH]C=O(27376)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['C#CC[CH]C=[C]O(30177)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.26084e+07,'s^-1'), n=1.66833, Ea=(126.789,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_H/OneDe;Cs_H_out_H/OneDe] + [R2H_S;C_rad_out_H/Ct;Cs_H_out_1H] for rate rule [R2H_S;C_rad_out_H/Ct;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#C[CH]C[C]=CO(30180)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#CCC[C]=[C]O(30181)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['C#C[CH][CH]C=CO(30182)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]#CCCC=[C]O(30183)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['C#CCC[CH][C]=O(27377)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(345676,'s^-1'), n=1.93175, Ea=(151.677,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_H/Ct;XH_out] for rate rule [R5HJ_3;C_rad_out_H/Ct;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[C]#C[CH]CC=CO(30184)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.79585e+08,'s^-1'), n=1.26031, Ea=(187.742,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Y_rad_out;Cd_H_out_singleNd] for rate rule [R6HJ_2;Ct_rad_out;Cd_H_out_singleNd]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C]O(172)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['O[C]=CCC1[C]=C1(30185)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_T;triplebond_intra_H;radadd_intra_cs] for rate rule [R3_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['OC1C=[C][CH]CC=1(30115)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(9.926e+10,'s^-1'), n=0.198, Ea=(22.8237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;triplebond_intra_H;radadd_intra_cdsingle] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['C#CC=CC=CO(30186)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['C#CCCC=C=O(27384)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#CC([CH2])C=[C]O(29678)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['C#CC1CC=C1O(30130)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriND_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C3H2(81)', '[CH2]C=[C]O(6303)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['[CH]=[C]C1CC=C1O(30187)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(79.5836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 78.3 to 79.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['C=C=[C]CC=[C]O(30188)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C=[C]CC=CO(30189)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(9.01194e+11,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_double;Cd_H_out_singleNd]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['C=[C]C=CC=[C]O(30135)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C=CC[C]=[C]O(30190)'],
    products = ['C#C[CH]CC=[C]O(29677)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(33613.8,'s^-1'), n=2.10442, Ea=(111.806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cd_H_out_singleH] for rate rule [R5H;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['C=C=CC[CH][C]=O(27388)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7HJ_5;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#C[CH]CC=[C]O(29677)'],
    products = ['C=C=CCC=C=O(27389)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

network(
    label = '4692',
    isomers = [
        'C#C[CH]CC=[C]O(29677)',
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
    label = '4692',
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

