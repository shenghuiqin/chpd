species(
    label = '[CH]=C(O)C=[C]C=C(27539)',
    structure = SMILES('[CH]=C(O)C=[C]C=C'),
    E0 = (377.892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.21986,'amu*angstrom^2'), symmetry=1, barrier=(28.0469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22111,'amu*angstrom^2'), symmetry=1, barrier=(28.0756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22084,'amu*angstrom^2'), symmetry=1, barrier=(28.0695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.48788,0.082274,-9.28211e-05,4.98511e-08,-9.92935e-12,45625.6,25.9778], Tmin=(100,'K'), Tmax=(1417.73,'K')), NASAPolynomial(coeffs=[21.8732,0.00657196,6.17923e-07,-3.62138e-10,3.16511e-14,40552.7,-85.2368], Tmin=(1417.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJC=C)"""),
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
    label = '[CH2]C1[C]=CC(O)=C1(30108)',
    structure = SMILES('[CH2]C1[C]=CC(O)=C1'),
    E0 = (350.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01123,0.0554724,-2.65111e-05,-1.43606e-08,1.25098e-11,42302,20.2669], Tmin=(100,'K'), Tmax=(922.57,'K')), NASAPolynomial(coeffs=[16.5681,0.0147375,-3.71634e-06,5.54969e-10,-3.74858e-14,38294.7,-59.6937], Tmin=(922.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclopentadiene) + radical(1,3-cyclopentadiene-vinyl-1) + radical(Isobutyl)"""),
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
    label = '[CH]=C(O)C#CC=C(30109)',
    structure = SMILES('[CH]=C(O)C#CC=C'),
    E0 = (360.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2100,2250,500,550,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.34028,'amu*angstrom^2'), symmetry=1, barrier=(30.8158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34789,'amu*angstrom^2'), symmetry=1, barrier=(30.9906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34292,'amu*angstrom^2'), symmetry=1, barrier=(30.8764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.820255,0.0565316,-2.58001e-05,-1.93973e-08,1.4781e-11,43468.4,23.4758], Tmin=(100,'K'), Tmax=(948.798,'K')), NASAPolynomial(coeffs=[19.954,0.00868282,-2.03469e-06,3.5803e-10,-2.96267e-14,38360.4,-75.612], Tmin=(948.798,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(360.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-Ct(Cds-Cds)) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)C=C=C=C(30110)',
    structure = SMILES('[CH]=C(O)C=C=C=C'),
    E0 = (408.282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,563.333,586.667,610,1970,2140,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(1.24577,'amu*angstrom^2'), symmetry=1, barrier=(28.6427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24474,'amu*angstrom^2'), symmetry=1, barrier=(28.6191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0219201,0.0775234,-9.03717e-05,5.03929e-08,-1.0624e-11,49256.8,23.6372], Tmin=(100,'K'), Tmax=(1252.1,'K')), NASAPolynomial(coeffs=[20.3836,0.00789299,-1.4657e-06,1.32974e-10,-5.26339e-15,44517,-77.7394], Tmin=(1252.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Cds_P)"""),
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
    label = 'C#CC=[C]C=C(27472)',
    structure = SMILES('C#CC=[C]C=C'),
    E0 = (521.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.62249,'amu*angstrom^2'), symmetry=1, barrier=(37.3042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61979,'amu*angstrom^2'), symmetry=1, barrier=(37.2422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45098,0.0475986,-2.51853e-05,-8.16809e-09,8.69319e-12,62818.4,18.2704], Tmin=(100,'K'), Tmax=(937.719,'K')), NASAPolynomial(coeffs=[14.4661,0.0131635,-3.8272e-06,6.24214e-10,-4.31862e-14,59450.6,-48.6255], Tmin=(937.719,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(521.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C(O)[C]=CC=C(30111)',
    structure = SMILES('[CH]C(O)=C=CC=C'),
    E0 = (354.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.411418,0.066564,-3.56082e-05,-1.23843e-08,1.28334e-11,42789.1,24.7075], Tmin=(100,'K'), Tmax=(936.931,'K')), NASAPolynomial(coeffs=[19.3918,0.0163888,-4.68006e-06,7.59953e-10,-5.31549e-14,37878.1,-72.8531], Tmin=(936.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(O)C=C[C]=C(30112)',
    structure = SMILES('[CH]C(O)=CC=C=C'),
    E0 = (354.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.411418,0.066564,-3.56082e-05,-1.23843e-08,1.28334e-11,42789.1,24.7075], Tmin=(100,'K'), Tmax=(936.931,'K')), NASAPolynomial(coeffs=[19.3918,0.0163888,-4.68006e-06,7.59953e-10,-5.31549e-14,37878.1,-72.8531], Tmin=(936.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C[C]=CC(=C)[O](26502)',
    structure = SMILES('C=C[C]=CC(=C)[O]'),
    E0 = (268.601,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21254,'amu*angstrom^2'), symmetry=1, barrier=(27.8787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21667,'amu*angstrom^2'), symmetry=1, barrier=(27.9737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3876.3,'J/mol'), sigma=(6.23234,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.47 K, Pc=36.33 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.118517,0.0727959,-7.66533e-05,4.02971e-08,-8.02712e-12,32455.7,24.9008], Tmin=(100,'K'), Tmax=(1380.29,'K')), NASAPolynomial(coeffs=[18.2599,0.0122707,-2.23671e-06,1.80618e-10,-5.15332e-15,28205.2,-65.7238], Tmin=(1380.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C[C]=[C]C(=C)O(27538)',
    structure = SMILES('C=C[C]=[C]C(=C)O'),
    E0 = (329.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,3615,1277.5,1000,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.41227,'amu*angstrom^2'), symmetry=1, barrier=(32.4708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41564,'amu*angstrom^2'), symmetry=1, barrier=(32.5484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4145,'amu*angstrom^2'), symmetry=1, barrier=(32.5222,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.290419,0.0830004,-9.76143e-05,5.59586e-08,-1.20009e-11,39829,24.4572], Tmin=(100,'K'), Tmax=(1300.02,'K')), NASAPolynomial(coeffs=[20.0388,0.0096718,-5.69708e-07,-1.84336e-10,2.20773e-14,35454.1,-75.4507], Tmin=(1300.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CC=CC(=[CH])O(30113)',
    structure = SMILES('[CH]=CC=CC(=[CH])O'),
    E0 = (425.993,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3115,3125,620,680,785,800,1600,1700,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680],'cm^-1')),
        HinderedRotor(inertia=(1.08243,'amu*angstrom^2'), symmetry=1, barrier=(24.8872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08189,'amu*angstrom^2'), symmetry=1, barrier=(24.8747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0833,'amu*angstrom^2'), symmetry=1, barrier=(24.9071,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.745607,0.0821774,-8.98903e-05,4.57356e-08,-8.54818e-12,51425,27.7205], Tmin=(100,'K'), Tmax=(1535.38,'K')), NASAPolynomial(coeffs=[23.4949,0.00372962,1.69359e-06,-5.19062e-10,3.9808e-14,45784.3,-93.7481], Tmin=(1535.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C=CC=C(27541)',
    structure = SMILES('[CH]=C([O])C=CC=C'),
    E0 = (316.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.00795,'amu*angstrom^2'), symmetry=1, barrier=(23.1748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00996,'amu*angstrom^2'), symmetry=1, barrier=(23.2209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.139478,0.0726991,-7.37084e-05,3.61498e-08,-6.62886e-12,38255.1,26.6447], Tmin=(100,'K'), Tmax=(1536.01,'K')), NASAPolynomial(coeffs=[20.1187,0.00908893,-9.88952e-07,-1.31715e-11,5.83859e-15,33312.3,-75.6153], Tmin=(1536.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C][C]=CC(=C)O(27542)',
    structure = SMILES('[CH2]C#CC=C([CH2])O'),
    E0 = (293.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2100,2250,500,550,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,289.39],'cm^-1')),
        HinderedRotor(inertia=(0.402479,'amu*angstrom^2'), symmetry=1, barrier=(23.9144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.402542,'amu*angstrom^2'), symmetry=1, barrier=(23.9159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.402469,'amu*angstrom^2'), symmetry=1, barrier=(23.914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.402374,'amu*angstrom^2'), symmetry=1, barrier=(23.9146,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.741107,0.0595161,-3.07576e-05,-1.44786e-08,1.32656e-11,35400.2,25.6163], Tmin=(100,'K'), Tmax=(936.059,'K')), NASAPolynomial(coeffs=[19.148,0.0114289,-2.68644e-06,4.17648e-10,-3.08622e-14,30615,-69.1236], Tmin=(936.059,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=C(O)CJ) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C[C]=CC(=C)O(27543)',
    structure = SMILES('[CH]C=C=CC(=C)O'),
    E0 = (355.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.98277,'amu*angstrom^2'), symmetry=1, barrier=(45.5879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98337,'amu*angstrom^2'), symmetry=1, barrier=(45.6015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98263,'amu*angstrom^2'), symmetry=1, barrier=(45.5845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.341226,0.0709099,-5.24418e-05,8.59105e-09,4.47144e-12,42869,24.1446], Tmin=(100,'K'), Tmax=(950.937,'K')), NASAPolynomial(coeffs=[17.9602,0.0192335,-6.31765e-06,1.06557e-09,-7.27387e-14,38503.7,-65.304], Tmin=(950.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(O)C=C1[CH]C1(30114)',
    structure = SMILES('[CH]=C(O)C=C1[CH]C1'),
    E0 = (413.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.821036,0.0531112,-3.85183e-06,-4.8058e-08,2.66304e-11,49806.3,20.7618], Tmin=(100,'K'), Tmax=(925.87,'K')), NASAPolynomial(coeffs=[21.4289,0.0074744,-2.1881e-07,-5.27016e-11,-3.92104e-16,44130.3,-87.1161], Tmin=(925.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_P) + radical(Allyl_S)"""),
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
    label = 'C=C=C=CC(=C)O(27547)',
    structure = SMILES('C=C=C=CC(=C)O'),
    E0 = (161.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0825154,0.07602,-8.04058e-05,4.11841e-08,-7.98507e-12,19544.9,23.6143], Tmin=(100,'K'), Tmax=(1390.17,'K')), NASAPolynomial(coeffs=[20.6151,0.00993236,-2.04726e-06,2.25925e-10,-1.14283e-14,14421.6,-80.7874], Tmin=(1390.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC1C=C(O)C=1(30116)',
    structure = SMILES('C=CC1C=C(O)C=1'),
    E0 = (280.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25915,0.0458984,2.80428e-06,-4.42562e-08,2.24876e-11,33848.6,18.1882], Tmin=(100,'K'), Tmax=(944.332,'K')), NASAPolynomial(coeffs=[17.2,0.0142359,-3.86284e-06,6.62923e-10,-4.99032e-14,29239,-66.2702], Tmin=(944.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(cyclobutadiene_13)"""),
)

species(
    label = '[CH2]C=[C]C1C=C1O(30117)',
    structure = SMILES('[CH2]C=[C]C1C=C1O'),
    E0 = (489.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.255008,0.0697028,-6.97843e-05,3.45638e-08,-6.54912e-12,59032.6,24.9472], Tmin=(100,'K'), Tmax=(1403.18,'K')), NASAPolynomial(coeffs=[18.6451,0.0123648,-3.23718e-06,4.50751e-10,-2.66661e-14,54355.4,-68.2803], Tmin=(1403.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(O)C=C=[C]C(30118)',
    structure = SMILES('[CH]C(O)=CC#CC'),
    E0 = (348.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.777683,0.0587545,-1.95155e-05,-2.49005e-08,1.65215e-11,42057.3,24.8483], Tmin=(100,'K'), Tmax=(928.695,'K')), NASAPolynomial(coeffs=[17.1269,0.0192228,-5.55158e-06,8.86635e-10,-6.0489e-14,37688.6,-59.9885], Tmin=(928.695,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.626,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(O)[C]=C=CC(30119)',
    structure = SMILES('[CH]=C(O)C#C[CH]C'),
    E0 = (385.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2250,500,550,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.08117,'amu*angstrom^2'), symmetry=1, barrier=(24.8582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07754,'amu*angstrom^2'), symmetry=1, barrier=(24.7747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.02246,'amu*angstrom^2'), symmetry=1, barrier=(69.4924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07441,'amu*angstrom^2'), symmetry=1, barrier=(24.7029,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0284261,0.0715453,-7.34925e-05,3.75465e-08,-7.19784e-12,46515.2,27.1878], Tmin=(100,'K'), Tmax=(1475.84,'K')), NASAPolynomial(coeffs=[18.1846,0.0114738,-1.39749e-06,-7.40173e-12,8.38829e-15,42339,-63.4646], Tmin=(1475.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Cds_P) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=C([O])C=C=CC(27550)',
    structure = SMILES('[CH]=C([O])C=C=CC'),
    E0 = (369.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.880157,'amu*angstrom^2'), symmetry=1, barrier=(20.2365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.879483,'amu*angstrom^2'), symmetry=1, barrier=(20.221,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.64381,0.0704199,-7.66356e-05,4.3236e-08,-9.59867e-12,44562.3,23.8541], Tmin=(100,'K'), Tmax=(1103.45,'K')), NASAPolynomial(coeffs=[14.7526,0.0192758,-7.1119e-06,1.23233e-09,-8.22471e-14,41448.6,-45.6016], Tmin=(1103.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C1[CH]C(O)=C1(30120)',
    structure = SMILES('[CH2]C=C1[CH]C(O)=C1'),
    E0 = (204.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04798,0.0389025,5.372e-05,-1.17026e-07,5.33244e-11,24679.5,19.0371], Tmin=(100,'K'), Tmax=(920.092,'K')), NASAPolynomial(coeffs=[25.0252,0.00222584,3.36924e-06,-7.37728e-10,4.30279e-14,17407.5,-110.183], Tmin=(920.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclobutene) + radical(C=CCJCO) + radical(C=CC=CCJ)"""),
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
    E0 = (377.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (391.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (585.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (636.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (554.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (632.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (549.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (590.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (592.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (570.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (570.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (648.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (547.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (686.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (569.101,'kJ/mol'),
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
    E0 = (515.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (390.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (402.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (386.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (489.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (585.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (619.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (457.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (515.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (501.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (772.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['HCCOH(50)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['[CH2]C1[C]=CC(O)=C1(30108)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.4226e+09,'s^-1'), n=0.3735, Ea=(13.1942,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[CH]=C(O)C#CC=C(30109)'],
    products = ['[CH]=C(O)C=[C]C=C(27539)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2371.94,'m^3/(mol*s)'), n=1.49517, Ea=(13.5032,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-De_Ct-De;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH]=C(O)C=C=C=C(30110)'],
    products = ['[CH]=C(O)C=[C]C=C(27539)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['HCCOH(50)', '[CH]=[C]C=C(4699)'],
    products = ['[CH]=C(O)C=[C]C=C(27539)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;CJ] for rate rule [Ct-O_Ct;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]O(172)', 'CH2CHCCH(26391)'],
    products = ['[CH]=C(O)C=[C]C=C(27539)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0117197,'m^3/(mol*s)'), n=2.33233, Ea=(9.74423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(5)', 'C#CC=[C]C=C(27472)'],
    products = ['[CH]=C(O)C=[C]C=C(27539)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(498.593,'m^3/(mol*s)'), n=1.168, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;OJ_pri] for rate rule [Ct-De_Ct-H;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -0.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['[CH]=C(O)[C]=CC=C(30111)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.0945e+07,'s^-1'), n=1.951, Ea=(212.263,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_single] for rate rule [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['[CH]=C(O)C=C[C]=C(30112)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['C=C[C]=CC(=C)[O](26502)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['C=C[C]=[C]C(=C)O(27538)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC=CC(=[CH])O(30113)'],
    products = ['[CH]=C(O)C=[C]C=C(27539)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['[CH]=C([O])C=CC=C(27541)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.66043e+07,'s^-1'), n=1.66125, Ea=(169.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_RSR;Cd_rad_out_singleDe_Cd;XH_out] + [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['C=[C][C]=CC(=C)O(27542)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['[CH]=C[C]=CC(=C)O(27543)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R6HJ_3;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C]O(172)', '[CH]=[C]C=C(4699)'],
    products = ['[CH]=C(O)C=[C]C=C(27539)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['[CH]=C(O)C=C1[CH]C1(30114)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['OC1C=[C][CH]CC=1(30115)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.16959e+10,'s^-1'), n=0.31, Ea=(12.393,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['C=C=C=CC(=C)O(27547)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['C=CC1C=C(O)C=1(30116)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleH_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['[CH2]C=[C]C1C=C1O(30117)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(111.722,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 110.1 to 111.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['[CH]=C(O)C=C=[C]C(30118)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(O)[C]=C=CC(30119)'],
    products = ['[CH]=C(O)C=[C]C=C(27539)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.992e-05,'s^-1'), n=4.805, Ea=(234.476,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R4H_MMS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['[CH]=C([O])C=C=CC(27550)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1763.2,'s^-1'), n=2.17098, Ea=(79.9643,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;C_rad_out_2H;XH_out] for rate rule [R6H;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['[CH2]C=C1[CH]C(O)=C1(30120)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(O)C=[C]C=C(27539)'],
    products = ['[CH]=C(O)C1[C]=CC1(30121)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(123.301,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 120.5 to 123.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(19)', '[CH]=C=CC(=[CH])O(30122)'],
    products = ['[CH]=C(O)C=[C]C=C(27539)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4698',
    isomers = [
        '[CH]=C(O)C=[C]C=C(27539)',
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
    label = '4698',
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

