species(
    label = 'C=C([CH]C)C=[C]O(26141)',
    structure = SMILES('[CH2]C(C=[C]O)=CC'),
    E0 = (202.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.888381,'amu*angstrom^2'), symmetry=1, barrier=(20.4256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888947,'amu*angstrom^2'), symmetry=1, barrier=(20.4386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888946,'amu*angstrom^2'), symmetry=1, barrier=(20.4386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888298,'amu*angstrom^2'), symmetry=1, barrier=(20.4237,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.358423,0.0676042,-3.68312e-05,-1.2555e-08,1.33823e-11,24528.5,26.741], Tmin=(100,'K'), Tmax=(930.222,'K')), NASAPolynomial(coeffs=[19.9379,0.0152641,-3.79478e-06,5.79085e-10,-4.04404e-14,19507.7,-73.7099], Tmin=(930.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CJO)"""),
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
    label = 'CH3CHCCH2(18175)',
    structure = SMILES('C=C=CC'),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.759584,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74635,0.0218189,8.22353e-06,-2.14768e-08,8.55624e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.6,'K')), NASAPolynomial(coeffs=[6.82078,0.0192338,-7.45622e-06,1.36536e-09,-9.53195e-14,16028,-10.4333], Tmin=(1025.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C1([CH]C)C=C1O(26242)',
    structure = SMILES('[CH2]C1([CH]C)C=C1O'),
    E0 = (368.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.505468,0.0648724,-3.69532e-05,-4.4877e-09,7.78475e-12,44515.7,26.661], Tmin=(100,'K'), Tmax=(999.125,'K')), NASAPolynomial(coeffs=[18.2421,0.0194783,-7.25787e-06,1.35752e-09,-9.83105e-14,39693,-65.2905], Tmin=(999.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsCs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(Neopentyl) + radical(Cs_S)"""),
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
    label = '[CH2]C(C=C=O)=CC(24884)',
    structure = SMILES('C=C([CH]C)C=C=O'),
    E0 = (104.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,219.591,921.681],'cm^-1')),
        HinderedRotor(inertia=(0.230739,'amu*angstrom^2'), symmetry=1, barrier=(7.89542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.528809,'amu*angstrom^2'), symmetry=1, barrier=(18.0945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.437466,'amu*angstrom^2'), symmetry=1, barrier=(14.969,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32624,0.0602416,-5.14001e-05,2.42382e-08,-4.83978e-12,12702.7,22.2525], Tmin=(100,'K'), Tmax=(1161.93,'K')), NASAPolynomial(coeffs=[9.8195,0.0310033,-1.36548e-05,2.58159e-09,-1.80182e-13,10729,-19.9973], Tmin=(1161.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C(C#CO)=CC(26243)',
    structure = SMILES('[CH2]C(C#CO)=CC'),
    E0 = (204.259,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2100,2250,500,550,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1615.65],'cm^-1')),
        HinderedRotor(inertia=(0.750107,'amu*angstrom^2'), symmetry=1, barrier=(20.9025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.98292,'amu*angstrom^2'), symmetry=1, barrier=(27.4429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.452565,'amu*angstrom^2'), symmetry=1, barrier=(12.5923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.985502,'amu*angstrom^2'), symmetry=1, barrier=(27.4442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.075,0.0580731,-4.3661e-05,1.67495e-08,-2.60459e-12,24677.1,24.3937], Tmin=(100,'K'), Tmax=(1506.76,'K')), NASAPolynomial(coeffs=[13.8818,0.0240749,-9.81536e-06,1.77445e-09,-1.19955e-13,20817.7,-42.642], Tmin=(1506.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CtH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtOs) + radical(Allyl_P)"""),
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
    label = 'C=[C][CH]C(18176)',
    structure = SMILES('[CH2][C]=CC'),
    E0 = (361.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.352622,'amu*angstrom^2'), symmetry=1, barrier=(8.10748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828631,'amu*angstrom^2'), symmetry=1, barrier=(19.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12013e-13,43485.7,14.8304], Tmin=(100,'K'), Tmax=(2065.83,'K')), NASAPolynomial(coeffs=[10.7464,0.014324,-5.20136e-06,8.69079e-10,-5.48385e-14,40045.6,-31.3799], Tmin=(2065.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH]C=O)=CC(24885)',
    structure = SMILES('[CH2]C(C=C[O])=CC'),
    E0 = (104.474,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11053,'amu*angstrom^2'), symmetry=1, barrier=(25.5333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11085,'amu*angstrom^2'), symmetry=1, barrier=(25.5406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11105,'amu*angstrom^2'), symmetry=1, barrier=(25.5452,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.560222,0.0588774,-6.56668e-06,-4.50636e-08,2.48168e-11,12704.6,24.298], Tmin=(100,'K'), Tmax=(939.246,'K')), NASAPolynomial(coeffs=[20.7761,0.0149073,-3.61815e-06,5.93117e-10,-4.52202e-14,7048.99,-81.8559], Tmin=(939.246,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([C]=CO)=CC(26244)',
    structure = SMILES('[CH2]C([C]=CO)=CC'),
    E0 = (162.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.12407,0.0862253,-8.64173e-05,4.18271e-08,-7.49422e-12,19692.2,28.547], Tmin=(100,'K'), Tmax=(1613.96,'K')), NASAPolynomial(coeffs=[22.6969,0.0103442,-2.39602e-07,-2.3593e-10,2.28963e-14,14196.8,-90.9818], Tmin=(1613.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=C(C)C=[C]O(26245)',
    structure = SMILES('C[C]=C(C)C=[C]O'),
    E0 = (289.098,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1670,1700,300,440,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.71008,'amu*angstrom^2'), symmetry=1, barrier=(16.3261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709768,'amu*angstrom^2'), symmetry=1, barrier=(16.319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.711185,'amu*angstrom^2'), symmetry=1, barrier=(16.3515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.711337,'amu*angstrom^2'), symmetry=1, barrier=(16.355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.292245,0.0751187,-7.17929e-05,3.24463e-08,-4.74726e-12,34909.7,27.1703], Tmin=(100,'K'), Tmax=(975.102,'K')), NASAPolynomial(coeffs=[16.6338,0.0205124,-6.9112e-06,1.15901e-09,-7.70805e-14,31131.9,-54.286], Tmin=(975.102,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'CC=C(C)[C]=[C]O(26246)',
    structure = SMILES('CC=C(C)[C]=[C]O'),
    E0 = (250.252,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1670,1700,300,440,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.844101,'amu*angstrom^2'), symmetry=1, barrier=(19.4075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.844623,'amu*angstrom^2'), symmetry=1, barrier=(19.4195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843955,'amu*angstrom^2'), symmetry=1, barrier=(19.4042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84427,'amu*angstrom^2'), symmetry=1, barrier=(19.4114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.287236,0.0765641,-7.42238e-05,3.34187e-08,-4.32734e-12,30236.8,26.3131], Tmin=(100,'K'), Tmax=(916.498,'K')), NASAPolynomial(coeffs=[16.1701,0.0212962,-6.76758e-06,1.07992e-09,-6.94493e-14,26735.3,-52.1475], Tmin=(916.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C(=[C]C)C=CO(26247)',
    structure = SMILES('[CH2]C(=[C]C)C=CO'),
    E0 = (200.853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.01107,'amu*angstrom^2'), symmetry=1, barrier=(23.2466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00701,'amu*angstrom^2'), symmetry=1, barrier=(23.1531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01146,'amu*angstrom^2'), symmetry=1, barrier=(23.2554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00583,'amu*angstrom^2'), symmetry=1, barrier=(23.1259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.122488,0.0703583,-3.43982e-05,-2.25609e-08,1.88748e-11,24310.6,24.9334], Tmin=(100,'K'), Tmax=(918.261,'K')), NASAPolynomial(coeffs=[22.8639,0.0106806,-1.24998e-06,8.21586e-11,-6.52429e-15,18473.6,-91.8829], Tmin=(918.261,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C[C](C)C=[C]O(23557)',
    structure = SMILES('[CH2]C=C(C)C=[C]O'),
    E0 = (169.312,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.94182,'amu*angstrom^2'), symmetry=1, barrier=(21.6543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.942638,'amu*angstrom^2'), symmetry=1, barrier=(21.6731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.943764,'amu*angstrom^2'), symmetry=1, barrier=(21.699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.941309,'amu*angstrom^2'), symmetry=1, barrier=(21.6425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.639459,0.0611824,-2.20546e-05,-2.52954e-08,1.74848e-11,20496.3,27.01], Tmin=(100,'K'), Tmax=(923.008,'K')), NASAPolynomial(coeffs=[18.4153,0.0171899,-4.25839e-06,6.34914e-10,-4.33134e-14,15807.4,-64.9486], Tmin=(923.008,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][C](C=C)C=CO(23561)',
    structure = SMILES('[CH2]C=C([CH2])C=CO'),
    E0 = (81.0673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.431061,0.0568846,1.37032e-05,-7.81804e-08,4.02148e-11,9898.83,24.9112], Tmin=(100,'K'), Tmax=(911.594,'K')), NASAPolynomial(coeffs=[24.8221,0.00705615,1.57764e-06,-4.8334e-10,3.06844e-14,3075.32,-103.539], Tmin=(911.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.0673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(Allyl_P)"""),
)

species(
    label = 'CC=C(C)[CH][C]=O(24886)',
    structure = SMILES('C[CH]C(C)=C[C]=O'),
    E0 = (122.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26793,0.0662242,-5.92344e-05,3.02763e-08,-6.82967e-12,14847,22.1539], Tmin=(100,'K'), Tmax=(1010.02,'K')), NASAPolynomial(coeffs=[8.09387,0.039191,-1.90864e-05,3.77622e-09,-2.70288e-13,13468.1,-10.8452], Tmin=(1010.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_S) + radical(C=CCJ=O)"""),
)

species(
    label = 'CC1C[C]1C=[C]O(26248)',
    structure = SMILES('CC1C[C]1C=[C]O'),
    E0 = (244.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38647,0.0363236,5.01287e-05,-1.00764e-07,4.43568e-11,29512.2,24.2178], Tmin=(100,'K'), Tmax=(926.272,'K')), NASAPolynomial(coeffs=[19.0935,0.014245,-2.19185e-06,2.82826e-10,-2.45222e-14,23898.7,-72.4466], Tmin=(926.272,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.433,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(Allyl_T) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]1C=C(O)C1C(26249)',
    structure = SMILES('C=C1[CH][C](O)C1C'),
    E0 = (189.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.8692,-0.00682173,0.000116105,-1.15371e-07,2.65104e-11,22439,-21.2754], Tmin=(100,'K'), Tmax=(1746.23,'K')), NASAPolynomial(coeffs=[81.4901,0.0196438,-6.83191e-05,1.67666e-08,-1.24454e-12,-30924.3,-483.523], Tmin=(1746.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(C2CsJOH) + radical(C=CCJCO)"""),
)

species(
    label = 'CC=C(C)C=C=O(24891)',
    structure = SMILES('CC=C(C)C=C=O'),
    E0 = (-49.6667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35874,0.0653297,-5.81465e-05,1.60105e-08,9.53542e-12,-5885.28,22.2017], Tmin=(100,'K'), Tmax=(566.768,'K')), NASAPolynomial(coeffs=[6.67586,0.0389649,-1.79089e-05,3.42614e-09,-2.39795e-13,-6667.26,-2.01271], Tmin=(566.768,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.6667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = 'CH2(S)(23)',
    structure = SMILES('[CH2]'),
    E0 = (419.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.36,2789.41,2993.36],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19195,-0.00230793,8.0509e-06,-6.60123e-09,1.95638e-12,50484.3,-0.754589], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.28556,0.00460255,-1.97412e-06,4.09548e-10,-3.34695e-14,50922.4,8.67684], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(419.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]C(=C)C=[C]O(26250)',
    structure = SMILES('[CH2]C(=C)C=[C]O'),
    E0 = (238.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.07431,'amu*angstrom^2'), symmetry=1, barrier=(24.7004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07356,'amu*angstrom^2'), symmetry=1, barrier=(24.6833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07643,'amu*angstrom^2'), symmetry=1, barrier=(24.7493,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02812,0.0517311,-1.50839e-05,-3.27879e-08,2.10487e-11,28838.5,21.8394], Tmin=(100,'K'), Tmax=(910.478,'K')), NASAPolynomial(coeffs=[19.841,0.00561651,6.96708e-07,-2.68812e-10,1.76819e-14,23898.4,-75.474], Tmin=(910.478,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = 'CC=[C]CC=[C]O(23615)',
    structure = SMILES('CC=[C]CC=[C]O'),
    E0 = (321.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,191.912,193.426],'cm^-1')),
        HinderedRotor(inertia=(0.465503,'amu*angstrom^2'), symmetry=1, barrier=(12.5844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.461016,'amu*angstrom^2'), symmetry=1, barrier=(12.582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.462048,'amu*angstrom^2'), symmetry=1, barrier=(12.5768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.472012,'amu*angstrom^2'), symmetry=1, barrier=(12.5972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3932.86,'J/mol'), sigma=(6.48948,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=614.30 K, Pc=32.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.526236,0.0680975,-6.07589e-05,2.83443e-08,-5.25795e-12,38801.7,30.5631], Tmin=(100,'K'), Tmax=(1306.12,'K')), NASAPolynomial(coeffs=[15.5725,0.0220181,-7.83935e-06,1.33311e-09,-8.78251e-14,34871.2,-46.0445], Tmin=(1306.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'CC=C1C=C(O)C1(26226)',
    structure = SMILES('CC=C1C=C(O)C1'),
    E0 = (-30.8599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.983063,0.0439093,3.78662e-05,-9.41853e-08,4.32257e-11,-3582.25,21.2053], Tmin=(100,'K'), Tmax=(929.73,'K')), NASAPolynomial(coeffs=[22.1666,0.0103026,-7.33302e-07,4.94991e-11,-1.0475e-14,-10007.7,-92.8216], Tmin=(929.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.8599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclobutene)"""),
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
    label = 'CC=[C]C=[C]O(26251)',
    structure = SMILES('CC=[C]C=[C]O'),
    E0 = (287.843,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.00051,'amu*angstrom^2'), symmetry=1, barrier=(23.0037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00297,'amu*angstrom^2'), symmetry=1, barrier=(23.0602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00307,'amu*angstrom^2'), symmetry=1, barrier=(23.0626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.633672,0.0624463,-6.43657e-05,3.33842e-08,-6.55829e-12,34750.6,24.4668], Tmin=(100,'K'), Tmax=(1411.62,'K')), NASAPolynomial(coeffs=[16.0373,0.0112722,-1.99037e-06,1.493e-10,-3.4447e-15,31151.6,-52.5011], Tmin=(1411.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC(=C)C=[C]O(23555)',
    structure = SMILES('C=CC(=C)C=[C]O'),
    E0 = (195.745,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.00796,'amu*angstrom^2'), symmetry=1, barrier=(23.1751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00435,'amu*angstrom^2'), symmetry=1, barrier=(23.0921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00503,'amu*angstrom^2'), symmetry=1, barrier=(23.1076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00392,0.0521759,-7.13778e-06,-3.70004e-08,2.07346e-11,23663.2,24.7741], Tmin=(100,'K'), Tmax=(933.265,'K')), NASAPolynomial(coeffs=[17.7212,0.0153279,-3.85051e-06,6.093e-10,-4.39679e-14,19027.2,-62.8426], Tmin=(933.265,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]CC(=C)C=[C]O(26252)',
    structure = SMILES('[CH2]CC(=C)C=[C]O'),
    E0 = (269.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,269.035,269.04],'cm^-1')),
        HinderedRotor(inertia=(0.319508,'amu*angstrom^2'), symmetry=1, barrier=(16.411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319507,'amu*angstrom^2'), symmetry=1, barrier=(16.4109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319506,'amu*angstrom^2'), symmetry=1, barrier=(16.411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319513,'amu*angstrom^2'), symmetry=1, barrier=(16.411,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.39755,0.0670841,-3.75642e-05,-1.10843e-08,1.27513e-11,32600.3,29.2129], Tmin=(100,'K'), Tmax=(929.662,'K')), NASAPolynomial(coeffs=[19.6989,0.015076,-3.73048e-06,5.66614e-10,-3.94321e-14,27670.3,-69.7113], Tmin=(929.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = 'C=C([C]=[C]O)CC(26253)',
    structure = SMILES('C=C([C]=[C]O)CC'),
    E0 = (263.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.888265,'amu*angstrom^2'), symmetry=1, barrier=(20.423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.886148,'amu*angstrom^2'), symmetry=1, barrier=(20.3743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.887264,'amu*angstrom^2'), symmetry=1, barrier=(20.3999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.887683,'amu*angstrom^2'), symmetry=1, barrier=(20.4096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.308171,0.0798289,-8.03493e-05,4.07437e-08,-7.86838e-12,31875.5,29.5571], Tmin=(100,'K'), Tmax=(1418.08,'K')), NASAPolynomial(coeffs=[19.6322,0.0152512,-3.22834e-06,3.44633e-10,-1.58455e-14,27057.8,-70.6553], Tmin=(1418.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C(C=[C]O)CC(26254)',
    structure = SMILES('[CH]=C(C=[C]O)CC'),
    E0 = (311.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.79293,'amu*angstrom^2'), symmetry=1, barrier=(18.231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.791094,'amu*angstrom^2'), symmetry=1, barrier=(18.1888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792524,'amu*angstrom^2'), symmetry=1, barrier=(18.2217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.793242,'amu*angstrom^2'), symmetry=1, barrier=(18.2382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.585644,0.0799298,-7.79655e-05,3.71689e-08,-6.65893e-12,37675.8,31.3733], Tmin=(100,'K'), Tmax=(1567.01,'K')), NASAPolynomial(coeffs=[21.5036,0.0120347,-1.95662e-06,1.44677e-10,-4.31783e-15,32166.1,-80.6077], Tmin=(1567.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]C)C=CO(26255)',
    structure = SMILES('[CH]C(C=CO)=CC'),
    E0 = (182.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.154739,0.0681908,-1.78357e-05,-3.89529e-08,2.40392e-11,22066.8,25.172], Tmin=(100,'K'), Tmax=(924.429,'K')), NASAPolynomial(coeffs=[21.6766,0.0173702,-4.01727e-06,5.856e-10,-4.12019e-14,16280.1,-86.7444], Tmin=(924.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C([CH][C]=O)CC(24895)',
    structure = SMILES('[CH2]C(=C[C]=O)CC'),
    E0 = (133.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.99178,0.069998,-6.37891e-05,3.13318e-08,-6.48901e-12,16108,24.5225], Tmin=(100,'K'), Tmax=(1122.43,'K')), NASAPolynomial(coeffs=[10.9256,0.0345966,-1.64789e-05,3.23175e-09,-2.30217e-13,13878,-24.5498], Tmin=(1122.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CCJ=O)"""),
)

species(
    label = 'C[CH][C]1C=C(O)C1(26256)',
    structure = SMILES('CC=C1[CH][C](O)C1'),
    E0 = (185.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[11.3995,-0.00935621,0.000118132,-1.15477e-07,2.62856e-11,21897.9,-22.6991], Tmin=(100,'K'), Tmax=(1760.66,'K')), NASAPolynomial(coeffs=[83.1581,0.0186221,-6.84312e-05,1.67804e-08,-1.24302e-12,-32975.5,-493.559], Tmin=(1760.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(methylenecyclobutane) + radical(C2CsJOH) + radical(C=CCJCO)"""),
)

species(
    label = 'C=CC(=C)C=CO(23567)',
    structure = SMILES('C=CC(=C)C=CO'),
    E0 = (-43.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.79274,0.0504278,1.97475e-05,-7.58616e-08,3.70311e-11,-5157.75,22.4038], Tmin=(100,'K'), Tmax=(923.765,'K')), NASAPolynomial(coeffs=[21.8621,0.0112011,-1.00305e-06,5.74549e-11,-8.23324e-15,-11269.3,-89.583], Tmin=(923.765,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(C=C=O)CC(24899)',
    structure = SMILES('C=C(C=C=O)CC'),
    E0 = (-36.2866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16173,0.0650063,-5.84974e-05,3.0053e-08,-6.59037e-12,-4264.22,23.9544], Tmin=(100,'K'), Tmax=(1063.82,'K')), NASAPolynomial(coeffs=[9.34634,0.0342318,-1.51046e-05,2.85982e-09,-1.99876e-13,-6005.6,-16.0379], Tmin=(1063.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-36.2866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=[C]C(C)C=[C]O(23559)',
    structure = SMILES('C=[C]C(C)C=[C]O'),
    E0 = (323.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.603266,'amu*angstrom^2'), symmetry=1, barrier=(13.8703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602863,'amu*angstrom^2'), symmetry=1, barrier=(13.861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60326,'amu*angstrom^2'), symmetry=1, barrier=(13.8701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602705,'amu*angstrom^2'), symmetry=1, barrier=(13.8574,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3876.19,'J/mol'), sigma=(6.45323,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.45 K, Pc=32.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.528279,0.0735384,-7.46492e-05,4.05178e-08,-8.80462e-12,39092.4,28.3331], Tmin=(100,'K'), Tmax=(1116.7,'K')), NASAPolynomial(coeffs=[14.0405,0.0251369,-9.63296e-06,1.70258e-09,-1.14738e-13,36074.7,-38.3468], Tmin=(1116.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = 'C=C1C=C(O)C1C(26231)',
    structure = SMILES('C=C1C=C(O)C1C'),
    E0 = (-28.4041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.727621,0.0523261,1.38765e-05,-6.93916e-08,3.45105e-11,-3280.25,18.5163], Tmin=(100,'K'), Tmax=(927.958,'K')), NASAPolynomial(coeffs=[21.9711,0.0112355,-1.30157e-06,1.35511e-10,-1.42474e-14,-9396.31,-94.0938], Tmin=(927.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.4041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = 'CHCH3(T)(95)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438699,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82363,-0.000909515,3.2138e-05,-3.7348e-08,1.3309e-11,41371.4,7.10948], Tmin=(100,'K'), Tmax=(960.812,'K')), NASAPolynomial(coeffs=[4.30487,0.00943069,-3.27559e-06,5.95121e-10,-4.27307e-14,40709.1,1.84202], Tmin=(960.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=[C]C=[C]O(8474)',
    structure = SMILES('C=[C]C=[C]O'),
    E0 = (323.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1670,1700,300,440,3615,1277.5,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.2815,'amu*angstrom^2'), symmetry=1, barrier=(29.4642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2808,'amu*angstrom^2'), symmetry=1, barrier=(29.4482,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03271,0.0497411,-5.36705e-05,2.76028e-08,-5.16708e-12,39072.6,21.2319], Tmin=(100,'K'), Tmax=(1581.42,'K')), NASAPolynomial(coeffs=[14.4057,0.00378978,1.41685e-06,-4.68662e-10,3.71413e-14,36359.2,-44.6198], Tmin=(1581.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJO)"""),
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
    E0 = (225.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (368.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (347.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (431.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (525.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (463.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (365.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (399.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (451.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (435.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (604.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (352.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (235.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (276.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (709.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (296.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (340.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (227.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (658.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (491.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (211.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (669.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (407.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (383.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (436.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (456.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (606.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (229.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (340.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (211.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (211.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (418.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (210.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (667.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['HCCOH(50)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(22.8996,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 22.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['[CH2]C1([CH]C)C=C1O(26242)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.69117e+12,'s^-1'), n=0.1675, Ea=(166.235,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R4_D_D;doublebond_intra_HNd;radadd_intra_cdsingleNd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic
Ea raised from 165.6 to 166.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[CH2]C(C=C=O)=CC(24884)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.185e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""Estimated using template [Od_R;HJ] for rate rule [Od_Cdd;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(C#CO)=CC(26243)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=[C]O(172)', 'CH3CHCCH2(18175)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['HCCOH(50)', 'C=[C][CH]C(18176)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['[CH2]C([CH]C=O)=CC(24885)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['[CH2]C([C]=CO)=CC(26244)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[C]=C(C)C=[C]O(26245)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC=C(C)[C]=[C]O(26246)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.18351e+08,'s^-1'), n=1.55833, Ea=(185.212,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Cd_rad_out;Cs_H_out_2H] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out_2H] for rate rule [R3H_SS_2Cd;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(=[C]C)C=CO(26247)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_single] for rate rule [R4H_DSD;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['C=C[C](C)C=[C]O(23557)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['[CH2][C](C=C)C=CO(23561)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSMS;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['CC=C(C)[CH][C]=O(24886)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_3;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]O(172)', 'C=[C][CH]C(18176)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['CC1C[C]1C=[C]O(26248)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.68393e+12,'s^-1'), n=-0.105173, Ea=(93.5715,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe_HNd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['[CH2][C]1C=C(O)C1C(26249)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['CC=C(C)C=C=O(24891)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH2(S)(23)', '[CH2]C(=C)C=[C]O(26250)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC=[C]CC=[C]O(23615)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['CC=C1C=C(O)C1(26226)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriND_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(19)', 'CC=[C]C=[C]O(26251)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', 'C=CC(=C)C=[C]O(23555)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]CC(=C)C=[C]O(26252)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C([C]=[C]O)CC(26253)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.44583e+10,'s^-1'), n=0.916667, Ea=(172.729,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Cd_rad_out;Cs_H_out_H/NonDeC] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_2Cd;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(C=[C]O)CC(26254)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['[CH]=C([CH]C)C=CO(26255)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.27529e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_singleH] for rate rule [R4H_DSD;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['C=C([CH][C]=O)CC(24895)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R5HJ_3;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['C[CH][C]1C=C(O)C1(26256)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['C=CC(=C)C=CO(23567)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['C=C(C=C=O)CC(24899)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]C(C)C=[C]O(23559)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C([CH]C)C=[C]O(26141)'],
    products = ['C=C1C=C(O)C1C(26231)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriND_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CHCH3(T)(95)', 'C=[C]C=[C]O(8474)'],
    products = ['C=C([CH]C)C=[C]O(26141)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4311',
    isomers = [
        'C=C([CH]C)C=[C]O(26141)',
    ],
    reactants = [
        ('HCCOH(50)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4311',
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

