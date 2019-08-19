species(
    label = '[CH]=C(C=C)O[CH]CCC=O(29722)',
    structure = SMILES('[CH]=C(C=C)O[CH]CCC=O'),
    E0 = (171.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62126,0.121661,-0.000138283,8.01189e-08,-1.82936e-11,20870.4,37.2825], Tmin=(100,'K'), Tmax=(1070.41,'K')), NASAPolynomial(coeffs=[22.1614,0.0327869,-1.37401e-05,2.55144e-09,-1.77182e-13,15779,-79.0733], Tmin=(1070.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCsJOC(O))"""),
)

species(
    label = 'O=CCCC=O(5767)',
    structure = SMILES('O=CCCC=O'),
    E0 = (-309.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.18601,'amu*angstrom^2'), symmetry=1, barrier=(4.27673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186498,'amu*angstrom^2'), symmetry=1, barrier=(4.28796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186203,'amu*angstrom^2'), symmetry=1, barrier=(4.28117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3653.08,'J/mol'), sigma=(5.8998,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.60 K, Pc=40.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98042,0.0494518,-6.06096e-05,5.52155e-08,-2.1639e-11,-37204.9,20.2227], Tmin=(100,'K'), Tmax=(774.657,'K')), NASAPolynomial(coeffs=[2.99129,0.0355848,-1.7014e-05,3.28724e-09,-2.30123e-13,-37102,17.2786], Tmin=(774.657,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-309.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH)"""),
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
    label = '[CH2]C1C=C1O[CH]CCC=O(30339)',
    structure = SMILES('[CH2]C1C=C1O[CH]CCC=O'),
    E0 = (276.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.44883,0.119854,-0.00013978,8.52755e-08,-2.06264e-11,33428.6,36.5112], Tmin=(100,'K'), Tmax=(1010.24,'K')), NASAPolynomial(coeffs=[19.8479,0.0355321,-1.45824e-05,2.65848e-09,-1.82007e-13,29125.6,-66.4509], Tmin=(1010.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(Isobutyl) + radical(CCsJOC(O))"""),
)

species(
    label = '[CH]=C(C=C)OC1CCC1[O](30340)',
    structure = SMILES('[CH]=C(C=C)OC1CCC1[O]'),
    E0 = (264.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.941351,0.085396,-1.73408e-05,-5.25554e-08,2.9832e-11,31965.9,34.0264], Tmin=(100,'K'), Tmax=(959.667,'K')), NASAPolynomial(coeffs=[27.6341,0.0225191,-6.94955e-06,1.28046e-09,-9.77199e-14,23892.1,-116.148], Tmin=(959.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C1OC(CCC=O)C1[CH2](30341)',
    structure = SMILES('[CH]=C1OC(CCC=O)C1[CH2]'),
    E0 = (198.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.660615,0.0873596,-4.03786e-05,-2.26652e-08,1.95143e-11,24044.6,34.2495], Tmin=(100,'K'), Tmax=(915.457,'K')), NASAPolynomial(coeffs=[22.4089,0.026784,-7.03248e-06,1.04795e-09,-6.88368e-14,18135.2,-84.2157], Tmin=(915.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1=CC([O])CC[CH]O1(30342)',
    structure = SMILES('C=CC1=CC([O])CC[CH]O1'),
    E0 = (134.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.817829,0.0419438,0.000186287,-3.30881e-07,1.47689e-10,16367.1,33.8313], Tmin=(100,'K'), Tmax=(904,'K')), NASAPolynomial(coeffs=[53.7184,-0.0301793,2.52291e-05,-5.07673e-09,3.34121e-13,-406.221,-262.007], Tmin=(904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(CC(C)OJ) + radical(CCsJOC(O))"""),
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
    label = '[CH]=C(C=C)OC=CCC=O(30343)',
    structure = SMILES('[CH]=C(C=C)OC=CCC=O'),
    E0 = (128.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.156,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.986905,0.0956366,-8.02161e-05,3.26143e-08,-5.21705e-12,15621.3,37.6914], Tmin=(100,'K'), Tmax=(1501.1,'K')), NASAPolynomial(coeffs=[24.7261,0.0271188,-1.17485e-05,2.20658e-09,-1.52806e-13,7901.75,-96.8036], Tmin=(1501.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
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
    label = 'CH2CHO(40)',
    structure = SMILES('[CH2]C=O'),
    E0 = (1.22925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,526.75,532.597,975.94,1639.13,1641.45],'cm^-1')),
        HinderedRotor(inertia=(0.00114821,'amu*angstrom^2'), symmetry=1, barrier=(2.1986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66874,0.0096233,1.60617e-05,-2.87682e-08,1.2503e-11,219.438,12.5694], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91637,0.0088465,-3.14955e-06,5.05413e-10,-3.01305e-14,-1047.8,-6.1065], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(1.22925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2CHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH]=C(C=C)OC=C(26672)',
    structure = SMILES('[CH]=C(C=C)OC=C'),
    E0 = (263.582,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.94833,'amu*angstrom^2'), symmetry=1, barrier=(21.804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.948084,'amu*angstrom^2'), symmetry=1, barrier=(21.7983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.94931,'amu*angstrom^2'), symmetry=1, barrier=(21.8265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605443,0.0650354,-4.59656e-05,6.2025e-09,4.2548e-12,31832.4,25.6455], Tmin=(100,'K'), Tmax=(976.472,'K')), NASAPolynomial(coeffs=[17.1697,0.0181788,-6.24141e-06,1.10251e-09,-7.73555e-14,27596.4,-58.9987], Tmin=(976.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
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
    label = 'C#CO[CH]CCC=O(21985)',
    structure = SMILES('C#CO[CH]CCC=O'),
    E0 = (73.9157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2175,525,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.325155,0.0917069,-0.000104941,6.00961e-08,-1.34362e-11,9048.75,27.6332], Tmin=(100,'K'), Tmax=(1097.94,'K')), NASAPolynomial(coeffs=[19.043,0.0211448,-8.53866e-06,1.56097e-09,-1.07734e-13,4795.74,-67.6166], Tmin=(1097.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.9157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCsJOCs)"""),
)

species(
    label = '[CH]=C(C=C)OC[CH]CC=O(30344)',
    structure = SMILES('[CH]=C(C=C)OC[CH]CC=O'),
    E0 = (177.803,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20303,0.108783,-0.000108957,5.60698e-08,-1.14614e-11,21576.9,38.8008], Tmin=(100,'K'), Tmax=(1186.23,'K')), NASAPolynomial(coeffs=[21.147,0.0334168,-1.36551e-05,2.50914e-09,-1.7325e-13,16274.5,-72.8417], Tmin=(1186.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)OCC[CH]C=O(30345)',
    structure = SMILES('[CH]=C(C=C)OCCC=C[O]'),
    E0 = (121.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60696,0.105234,-7.19285e-05,-1.19572e-10,1.26257e-11,14789.6,37.4631], Tmin=(100,'K'), Tmax=(953.425,'K')), NASAPolynomial(coeffs=[28.8912,0.0202775,-5.91051e-06,1.01641e-09,-7.42681e-14,7019.9,-118.466], Tmin=(953.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C(=C)O[CH]CCC=O(30346)',
    structure = SMILES('C=[C]C(=C)O[CH]CCC=O'),
    E0 = (123.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60558,0.124398,-0.000149531,9.38695e-08,-2.33429e-11,15082,36.4235], Tmin=(100,'K'), Tmax=(983.663,'K')), NASAPolynomial(coeffs=[20.0725,0.0362466,-1.51106e-05,2.76882e-09,-1.8984e-13,10817.1,-67.8041], Tmin=(983.663,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C(C=C)OCCC[C]=O(30347)',
    structure = SMILES('[CH]=C(C=C)OCCC[C]=O'),
    E0 = (137.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33366,0.115319,-0.000126803,7.22573e-08,-1.63339e-11,16774.8,37.5947], Tmin=(100,'K'), Tmax=(1077.82,'K')), NASAPolynomial(coeffs=[20.4664,0.0344154,-1.42095e-05,2.61533e-09,-1.80564e-13,12075.4,-69.2121], Tmin=(1077.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([C]=C)OCCCC=O(30348)',
    structure = SMILES('[CH]C(=C=C)OCCCC=O'),
    E0 = (153.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.856238,0.108115,-0.000100938,5.11556e-08,-1.07525e-11,18648.2,36.6149], Tmin=(100,'K'), Tmax=(1124.68,'K')), NASAPolynomial(coeffs=[15.894,0.0485413,-2.14819e-05,4.05674e-09,-2.83017e-13,14880.5,-46.1633], Tmin=(1124.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC(=C)O[CH]CCC=O(30349)',
    structure = SMILES('[CH]=CC(=C)O[CH]CCC=O'),
    E0 = (171.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62126,0.121661,-0.000138283,8.01189e-08,-1.82936e-11,20870.4,37.2825], Tmin=(100,'K'), Tmax=(1070.41,'K')), NASAPolynomial(coeffs=[22.1614,0.0327869,-1.37401e-05,2.55144e-09,-1.77182e-13,15779,-79.0733], Tmin=(1070.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCsJOC(O))"""),
)

species(
    label = 'C=CC(=C)O[CH][CH]CC=O(30350)',
    structure = SMILES('[CH2]C=C([CH2])OC=CCC=O'),
    E0 = (67.1244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.873596,0.0925905,-6.56913e-05,1.83279e-08,-6.98134e-13,8260.82,39.6339], Tmin=(100,'K'), Tmax=(1186.41,'K')), NASAPolynomial(coeffs=[22.0939,0.0336776,-1.46245e-05,2.79139e-09,-1.97109e-13,1507.47,-80.5906], Tmin=(1186.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.1244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=C(O)CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=CC(=[CH])OCCCC=O(30351)',
    structure = SMILES('[CH]=CC(=[CH])OCCCC=O'),
    E0 = (224.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3010,987.5,1337.5,450,1655,3115,3125,620,680,785,800,1600,1700,2782.5,750,1395,475,1775,1000,350,440,435,1725,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25614,0.114926,-0.000124977,7.0758e-08,-1.59886e-11,27250.9,36.9029], Tmin=(100,'K'), Tmax=(1074.17,'K')), NASAPolynomial(coeffs=[19.7565,0.0366786,-1.57083e-05,2.94158e-09,-2.04924e-13,22736.7,-65.9739], Tmin=(1074.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(=C)O[CH]C[CH]C=O(30352)',
    structure = SMILES('C=CC(=C)O[CH]CC=C[O]'),
    E0 = (67.9838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.07093,0.12455,-0.000127565,6.19101e-08,-1.12385e-11,8457.66,41.8112], Tmin=(100,'K'), Tmax=(1531.51,'K')), NASAPolynomial(coeffs=[33.6147,0.0123149,-1.5569e-06,5.8052e-11,8.23794e-16,-853.583,-144.527], Tmin=(1531.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.9838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=COJ) + radical(CCsJOC(O))"""),
)

species(
    label = 'C=CC(=C)O[CH]CC[C]=O(30353)',
    structure = SMILES('C=CC(=C)O[CH]CC[C]=O'),
    E0 = (84.6928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69865,0.122052,-0.000140104,8.1612e-08,-1.86364e-11,10394.3,37.9739], Tmin=(100,'K'), Tmax=(1073.75,'K')), NASAPolynomial(coeffs=[22.8706,0.0305251,-1.2242e-05,2.22538e-09,-1.52837e-13,5118.08,-82.3069], Tmin=(1073.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.6928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(CCsJOC(O))"""),
)

species(
    label = '[O][CH]CCC=O(5764)',
    structure = SMILES('[O][CH]CCC=O'),
    E0 = (19.0721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,180,180,2092.49],'cm^-1')),
        HinderedRotor(inertia=(0.166661,'amu*angstrom^2'), symmetry=1, barrier=(3.83187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164989,'amu*angstrom^2'), symmetry=1, barrier=(3.79342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166144,'amu*angstrom^2'), symmetry=1, barrier=(3.81999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56305,0.0641173,-0.000107371,1.07453e-07,-4.10375e-11,2371.43,23.494], Tmin=(100,'K'), Tmax=(845.992,'K')), NASAPolynomial(coeffs=[1.74033,0.0387928,-1.90537e-05,3.64333e-09,-2.50242e-13,3217.68,27.8472], Tmin=(845.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.0721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'O=CCC[CH]OC1[CH]CC=1(30354)',
    structure = SMILES('O=CCC[CH]OC1[CH]CC=1'),
    E0 = (174.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.918236,0.0923951,-5.75084e-05,1.39484e-09,7.72728e-12,21180.2,38.2744], Tmin=(100,'K'), Tmax=(1012.88,'K')), NASAPolynomial(coeffs=[23.2525,0.028989,-1.10678e-05,2.06513e-09,-1.48069e-13,14639.9,-86.7596], Tmin=(1012.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(CCsJOC(O)) + radical(CCJCO)"""),
)

species(
    label = '[CH]=C(C=C)OC1CC[CH]O1(30355)',
    structure = SMILES('[CH]=C(C=C)OC1CC[CH]O1'),
    E0 = (146.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.866732,0.0844471,-1.03048e-05,-6.91872e-08,4.01923e-11,17815.7,30.8811], Tmin=(100,'K'), Tmax=(901.379,'K')), NASAPolynomial(coeffs=[28.121,0.0173968,-1.21274e-06,-1.11079e-10,1.03045e-14,10088,-119.836], Tmin=(901.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Tetrahydrofuran) + radical(CCsJOCs) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1[CH]CC(CCC=O)O1(30356)',
    structure = SMILES('[CH]C1=CCC(CCC=O)O1'),
    E0 = (33.5862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.208552,0.0772182,-1.59575e-05,-3.40024e-08,1.88545e-11,4204.71,32.2842], Tmin=(100,'K'), Tmax=(972.822,'K')), NASAPolynomial(coeffs=[17.2509,0.0408989,-1.46473e-05,2.57886e-09,-1.77891e-13,-870.669,-60.0929], Tmin=(972.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.5862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(2,3-Dihydrofuran) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC1=CO[CH]CC[CH]O1(30357)',
    structure = SMILES('C=CC1=CO[CH]CC[CH]O1'),
    E0 = (117.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14205,0.0776998,3.50671e-05,-1.28874e-07,6.26342e-11,14308.9,28.2268], Tmin=(100,'K'), Tmax=(930.703,'K')), NASAPolynomial(coeffs=[36.2601,0.00762021,1.88505e-06,-4.33276e-10,1.67807e-14,3419.94,-170.627], Tmin=(930.703,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclooctane) + radical(CCsJOC(O)) + radical(CCsJOC(O))"""),
)

species(
    label = 'C=CC(=C)OC=CCC=O(30358)',
    structure = SMILES('C=CC(=C)OC=CCC=O'),
    E0 = (-118.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.833684,0.0913366,-6.15821e-05,1.36982e-08,9.37223e-13,-14102.1,36.7278], Tmin=(100,'K'), Tmax=(1157.4,'K')), NASAPolynomial(coeffs=[22.1391,0.0335569,-1.4712e-05,2.83638e-09,-2.01931e-13,-20867.6,-83.7146], Tmin=(1157.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CO(12)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(762.44,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35308e-07,1.51269e-10,-9.88872e-15,-14292.7,6.51157], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH]=C(C=C)O[CH]CC(26561)',
    structure = SMILES('[CH]=C(C=C)O[CH]CC'),
    E0 = (278.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3521.87,'J/mol'), sigma=(6.18816,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=550.11 K, Pc=33.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33347,0.102542,-0.000105763,5.3938e-08,-1.05971e-11,33689.3,31.9723], Tmin=(100,'K'), Tmax=(1288.78,'K')), NASAPolynomial(coeffs=[24.9585,0.0187802,-5.76036e-06,9.08133e-10,-5.80399e-14,27091.7,-100.846], Tmin=(1288.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)OC([CH2])CC=O(30359)',
    structure = SMILES('[CH]=C(C=C)OC([CH2])CC=O'),
    E0 = (175.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57292,0.122403,-0.00014228,8.52254e-08,-2.01905e-11,21341,37.045], Tmin=(100,'K'), Tmax=(1031.15,'K')), NASAPolynomial(coeffs=[21.0394,0.0346862,-1.46791e-05,2.72787e-09,-1.89125e-13,16677.7,-72.7398], Tmin=(1031.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1=CC(CCC=O)O1(30360)',
    structure = SMILES('C=CC1=CC(CCC=O)O1'),
    E0 = (-121.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.709431,0.0834987,-2.51991e-05,-3.55362e-08,2.15278e-11,-14387.7,33.2208], Tmin=(100,'K'), Tmax=(982.97,'K')), NASAPolynomial(coeffs=[24.4472,0.0272984,-9.8926e-06,1.86617e-09,-1.37677e-13,-21563.9,-99.0593], Tmin=(982.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-121.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH]=C(C=C)O[CH]COC=C(30361)',
    structure = SMILES('[CH]=C(C=C)O[CH]COC=C'),
    E0 = (213.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.54252,0.13591,-0.000149386,7.6622e-08,-1.46109e-11,25966.3,41.2305], Tmin=(100,'K'), Tmax=(1458.82,'K')), NASAPolynomial(coeffs=[36.5109,0.00808986,5.45559e-07,-3.5158e-10,2.92987e-14,16195.1,-160.566], Tmin=(1458.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)O[CH]CC=CO(30362)',
    structure = SMILES('[CH]=C(C=C)O[CH]CC=CO'),
    E0 = (173.617,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.56537,0.135367,-0.000149738,7.71713e-08,-1.4729e-11,21180.3,42.5641], Tmin=(100,'K'), Tmax=(1471.49,'K')), NASAPolynomial(coeffs=[36.2652,0.00709576,1.4047e-06,-5.40397e-10,4.30013e-14,11623.4,-157.624], Tmin=(1471.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCsJOC(O))"""),
)

species(
    label = '[CH]=C([O])C=C(5816)',
    structure = SMILES('[CH]=C([O])C=C'),
    E0 = (264.929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.979636,'amu*angstrom^2'), symmetry=1, barrier=(22.5238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41648,0.0452197,-4.57103e-05,2.25525e-08,-4.15485e-12,31966.6,18.4199], Tmin=(100,'K'), Tmax=(1540.92,'K')), NASAPolynomial(coeffs=[13.6856,0.00610864,-5.68369e-07,-3.63413e-11,6.19928e-15,29047.6,-43.2788], Tmin=(1540.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]CCC=O(18484)',
    structure = SMILES('[CH]CCC=O'),
    E0 = (221.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,338.309,338.317,1520,1520.03],'cm^-1')),
        HinderedRotor(inertia=(0.00147232,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09193,'amu*angstrom^2'), symmetry=1, barrier=(7.47087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0919401,'amu*angstrom^2'), symmetry=1, barrier=(7.47082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.39867,0.0383489,-2.66828e-05,1.05264e-08,-1.89984e-12,26648.5,18.3896], Tmin=(100,'K'), Tmax=(1182.31,'K')), NASAPolynomial(coeffs=[5.79639,0.0268536,-1.20987e-05,2.30283e-09,-1.60963e-13,25845.1,1.42861], Tmin=(1182.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ2_triplet)"""),
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
    E0 = (171.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (276.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (264.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (226.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (320.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (346.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (177.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (300.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (376.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (286.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (299.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (394.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (268.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (230.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (575.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (317.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (268.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (216.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (323.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (470.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (309.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (229.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (229.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (226.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (196.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (516.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (334.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (179.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (527.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (317.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (486.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['O=CCCC=O(5767)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['[CH2]C1C=C1O[CH]CCC=O(30339)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(104.481,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 14 used for R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH
Exact match found for rate rule [R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 104.1 to 104.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['[CH]=C(C=C)OC1CCC1[O](30340)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13771e+06,'s^-1'), n=1.58803, Ea=(92.2965,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 87.4 to 92.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['[CH]=C1OC(CCC=O)C1[CH2](30341)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.73e+06,'s^-1'), n=1.31, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 344 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_csHNd
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['C=CC1=CC([O])CC[CH]O1(30342)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.77692e+10,'s^-1'), n=0.185556, Ea=(148.486,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R8;carbonylbond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]=C(C=C)OC=CCC=O(30343)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]C=C(4699)', 'O=CCCC=O(5767)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2CHO(40)', '[CH]=C(C=C)OC=C(26672)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0114756,'m^3/(mol*s)'), n=2.44484, Ea=(36.1431,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds;CsJ-OneDeHH] for rate rule [Cds-HH_Cds-OsH;CsJ-COHH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C2H3(30)', 'C#CO[CH]CCC=O(21985)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.324566,'m^3/(mol*s)'), n=2.2487, Ea=(16.299,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;CdsJ-H] for rate rule [Ct-O_Ct;CdsJ-H]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(C=C)OC[CH]CC=O(30344)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['[CH]=C(C=C)OCC[CH]C=O(30345)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_1H;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['C=[C]C(=C)O[CH]CCC=O(30346)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['[CH]=C(C=C)OCCC[C]=O(30347)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(27900,'s^-1'), n=1.97, Ea=(96.6504,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeO;XH_out] for rate rule [R4H_SSS;C_rad_out_H/NonDeO;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['[CH]=C([C]=C)OCCCC=O(30348)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.164999,'s^-1'), n=3.35583, Ea=(59.1095,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;XH_out] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC(=C)O[CH]CCC=O(30349)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.27529e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_singleH] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['C=CC(=C)O[CH][CH]CC=O(30350)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=CC(=[CH])OCCCC=O(30351)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(366176,'s^-1'), n=1.54456, Ea=(43.4053,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.74165738677
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['C=CC(=C)O[CH]C[CH]C=O(30352)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R6HJ_3;Cd_rad_out_singleH;Cs_H_out_H/CO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['C=CC(=C)O[CH]CC[C]=O(30353)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7HJ_3;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]C=C(4699)', '[O][CH]CCC=O(5764)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['O=CCC[CH]OC1[CH]CC=1(30354)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R4_D_D;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['[CH]=C(C=C)OC1CC[CH]O1(30355)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_csHO]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['[CH]=C1[CH]CC(CCC=O)O1(30356)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['C=CC1=CO[CH]CC[CH]O1(30357)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.25459e+11,'s^-1'), n=0.208106, Ea=(54.2604,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;multiplebond_intra;radadd_intra_cdsingleH] + [R6plus;carbonyl_intra_H;radadd_intra] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['C=CC(=C)OC=CCC=O(30358)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CO(12)', '[CH]=C(C=C)O[CH]CC(26561)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for CO;C_pri/NonDeC
Exact match found for rate rule [CO;C_pri/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C(C=C)OC([CH2])CC=O(30359)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    products = ['C=CC1=CC(CCC=O)O1(30360)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C(C=C)O[CH]COC=C(30361)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C(C=C)O[CH]CC=CO(30362)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C([O])C=C(5816)', '[CH]CCC=O(18484)'],
    products = ['[CH]=C(C=C)O[CH]CCC=O(29722)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '4741',
    isomers = [
        '[CH]=C(C=C)O[CH]CCC=O(29722)',
    ],
    reactants = [
        ('O=CCCC=O(5767)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4741',
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

