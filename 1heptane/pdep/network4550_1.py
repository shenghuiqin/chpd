species(
    label = 'C=C[C]=CC(=C)[CH]C(26578)',
    structure = SMILES('[CH2]C(C=[C]C=C)=CC'),
    E0 = (422.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11314,'amu*angstrom^2'), symmetry=1, barrier=(25.5932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11308,'amu*angstrom^2'), symmetry=1, barrier=(25.592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11281,'amu*angstrom^2'), symmetry=1, barrier=(25.5856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1124,'amu*angstrom^2'), symmetry=1, barrier=(25.5762,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.20964,0.0813771,-5.33489e-05,4.86344e-09,6.0616e-12,50985.3,28.0233], Tmin=(100,'K'), Tmax=(965.094,'K')), NASAPolynomial(coeffs=[18.7617,0.0279536,-9.49238e-06,1.631e-09,-1.11448e-13,46149.6,-68.9099], Tmin=(965.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJC=C)"""),
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
    label = '[CH2]C1([CH]C)C=C1C=C(28131)',
    structure = SMILES('[CH2]C1([CH]C)C=C1C=C'),
    E0 = (633.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00303269,0.0751271,-4.21966e-05,-2.91308e-10,5.59422e-12,76378.1,29.8398], Tmin=(100,'K'), Tmax=(1052.29,'K')), NASAPolynomial(coeffs=[18.2145,0.0299114,-1.1969e-05,2.24144e-09,-1.58888e-13,71216,-65.2647], Tmin=(1052.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsCs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Cs_S) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C1[C]=CC(=CC)C1(27993)',
    structure = SMILES('[CH2]C1[C]=CC(=CC)C1'),
    E0 = (491.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.961058,0.0518598,1.47115e-05,-5.37401e-08,2.38393e-11,59268.8,27.5104], Tmin=(100,'K'), Tmax=(974.944,'K')), NASAPolynomial(coeffs=[14.1969,0.0333638,-1.19238e-05,2.14542e-09,-1.51251e-13,54986.2,-44.7369], Tmin=(974.944,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(3-Methylenecyclopentene) + radical(Isobutyl) + radical(cyclopentene-vinyl)"""),
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
    label = '[CH2]C(C#CC=C)=CC(28207)',
    structure = SMILES('[CH2]C(C#CC=C)=CC'),
    E0 = (402.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,342.844,343.23],'cm^-1')),
        HinderedRotor(inertia=(0.289921,'amu*angstrom^2'), symmetry=1, barrier=(24.2104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289345,'amu*angstrom^2'), symmetry=1, barrier=(24.2136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290229,'amu*angstrom^2'), symmetry=1, barrier=(24.2134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289719,'amu*angstrom^2'), symmetry=1, barrier=(24.2127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.594822,0.0622587,-1.90943e-05,-1.84882e-08,1.10347e-11,48504.3,28.3008], Tmin=(100,'K'), Tmax=(1027.39,'K')), NASAPolynomial(coeffs=[15.5118,0.0308409,-1.21468e-05,2.26042e-09,-1.60065e-13,44032.3,-50.9153], Tmin=(1027.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-Ct(Cds-Cds)) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C=C=C)=CC(28208)',
    structure = SMILES('[CH2]C(C=C=C=C)=CC'),
    E0 = (452.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.16475,'amu*angstrom^2'), symmetry=1, barrier=(26.7799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16556,'amu*angstrom^2'), symmetry=1, barrier=(26.7986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16181,'amu*angstrom^2'), symmetry=1, barrier=(26.7124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.248219,0.0829739,-7.25099e-05,3.25411e-08,-5.80457e-12,54640.4,27.6575], Tmin=(100,'K'), Tmax=(1353.23,'K')), NASAPolynomial(coeffs=[18.8223,0.026604,-1.00269e-05,1.75924e-09,-1.17902e-13,49479,-70.1157], Tmin=(1353.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Allyl_P)"""),
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
    label = '[CH2]C([C]=CC=C)=CC(28209)',
    structure = SMILES('[CH2]C([C]=CC=C)=CC'),
    E0 = (422.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11314,'amu*angstrom^2'), symmetry=1, barrier=(25.5932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11308,'amu*angstrom^2'), symmetry=1, barrier=(25.592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11281,'amu*angstrom^2'), symmetry=1, barrier=(25.5856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1124,'amu*angstrom^2'), symmetry=1, barrier=(25.5762,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.20964,0.0813771,-5.33489e-05,4.86344e-09,6.0616e-12,50985.3,28.0233], Tmin=(100,'K'), Tmax=(965.094,'K')), NASAPolynomial(coeffs=[18.7617,0.0279536,-9.49238e-06,1.631e-09,-1.11448e-13,46149.6,-68.9099], Tmin=(965.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C([CH]C=C=C)=CC(28184)',
    structure = SMILES('[CH2]C([CH]C)=CC=C=C'),
    E0 = (395.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.35884,0.0661228,-1.58305e-05,-2.7435e-08,1.55174e-11,47764.3,27.8502], Tmin=(100,'K'), Tmax=(987.888,'K')), NASAPolynomial(coeffs=[16.6374,0.0317915,-1.16554e-05,2.10819e-09,-1.48297e-13,43007,-58.286], Tmin=(987.888,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC=CCJ) + radical(Allyl_S)"""),
)

species(
    label = 'C=C[C]=CC(C)=[C]C(28210)',
    structure = SMILES('C=C[C]=CC(C)=[C]C'),
    E0 = (508.914,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1670,1700,300,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.872917,'amu*angstrom^2'), symmetry=1, barrier=(20.0701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.87146,'amu*angstrom^2'), symmetry=1, barrier=(20.0366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.871532,'amu*angstrom^2'), symmetry=1, barrier=(20.0382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.871566,'amu*angstrom^2'), symmetry=1, barrier=(20.039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.168926,0.0875971,-8.35837e-05,4.33645e-08,-9.09611e-12,61361.8,28.071], Tmin=(100,'K'), Tmax=(1149.79,'K')), NASAPolynomial(coeffs=[15.4706,0.0331894,-1.26052e-05,2.21063e-09,-1.48097e-13,57765.3,-49.564], Tmin=(1149.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=[C]C(C)=CC(28211)',
    structure = SMILES('C=C[C]=[C]C(C)=CC'),
    E0 = (470.068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1670,1700,300,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06445,'amu*angstrom^2'), symmetry=1, barrier=(24.4739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06636,'amu*angstrom^2'), symmetry=1, barrier=(24.5176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06572,'amu*angstrom^2'), symmetry=1, barrier=(24.5029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06526,'amu*angstrom^2'), symmetry=1, barrier=(24.4923,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.253269,0.0900129,-8.95961e-05,4.92993e-08,-1.09542e-11,56692.3,27.4963], Tmin=(100,'K'), Tmax=(1090.65,'K')), NASAPolynomial(coeffs=[15.0491,0.0338909,-1.24098e-05,2.1186e-09,-1.39349e-13,53354.4,-47.6564], Tmin=(1090.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CC=CC([CH2])=CC(28212)',
    structure = SMILES('[CH]=CC=CC([CH2])=CC'),
    E0 = (470.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.969803,'amu*angstrom^2'), symmetry=1, barrier=(22.2977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.969787,'amu*angstrom^2'), symmetry=1, barrier=(22.2973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.969862,'amu*angstrom^2'), symmetry=1, barrier=(22.299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.969357,'amu*angstrom^2'), symmetry=1, barrier=(22.2874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.200559,0.0782961,-4.06861e-05,-1.09497e-08,1.20612e-11,56772.7,28.7973], Tmin=(100,'K'), Tmax=(967.639,'K')), NASAPolynomial(coeffs=[20.525,0.0250634,-8.45729e-06,1.49412e-09,-1.05547e-13,51242.9,-78.3577], Tmin=(967.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=[C]C)C=CC=C(28213)',
    structure = SMILES('[CH2]C(=[C]C)C=CC=C'),
    E0 = (461.418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.916024,'amu*angstrom^2'), symmetry=1, barrier=(21.0612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.915971,'amu*angstrom^2'), symmetry=1, barrier=(21.06,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.918178,'amu*angstrom^2'), symmetry=1, barrier=(21.1107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.917008,'amu*angstrom^2'), symmetry=1, barrier=(21.0838,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.171783,0.0795145,-4.93103e-05,1.57438e-09,6.74168e-12,55656.8,28.7645], Tmin=(100,'K'), Tmax=(986.408,'K')), NASAPolynomial(coeffs=[19.1738,0.027264,-9.69315e-06,1.72406e-09,-1.20269e-13,50565.8,-70.7624], Tmin=(986.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[C]=C[C](C)C=C(27245)',
    structure = SMILES('[CH2]C=C(C)C=[C]C=C'),
    E0 = (389.128,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.18545,'amu*angstrom^2'), symmetry=1, barrier=(27.2558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18514,'amu*angstrom^2'), symmetry=1, barrier=(27.2487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66756,'amu*angstrom^2'), symmetry=1, barrier=(38.3405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18422,'amu*angstrom^2'), symmetry=1, barrier=(27.2275,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0776901,0.074885,-3.83502e-05,-8.11734e-09,1.02395e-11,46952.8,28.2695], Tmin=(100,'K'), Tmax=(952.817,'K')), NASAPolynomial(coeffs=[17.1906,0.0299609,-1.00027e-05,1.69779e-09,-1.15227e-13,42469.9,-59.8753], Tmin=(952.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C](C=C)C=CC=C(27249)',
    structure = SMILES('[CH2]C=CC=C([CH2])C=C'),
    E0 = (308.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.450634,0.0592502,1.47434e-05,-6.80801e-08,3.26041e-11,37211.5,28.8927], Tmin=(100,'K'), Tmax=(939.517,'K')), NASAPolynomial(coeffs=[19.3591,0.0259869,-7.57054e-06,1.27112e-09,-9.062e-14,31573.7,-72.2455], Tmin=(939.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C][C]=CC(C)=CC(28214)',
    structure = SMILES('[CH2]C#CC=C(C)[CH]C'),
    E0 = (428.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2100,2250,500,550,350,440,435,1725,3000,3100,440,815,1455,1000,180,1074.92],'cm^-1')),
        HinderedRotor(inertia=(1.01012,'amu*angstrom^2'), symmetry=1, barrier=(23.2247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219021,'amu*angstrom^2'), symmetry=1, barrier=(5.03573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0412056,'amu*angstrom^2'), symmetry=1, barrier=(33.8705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01006,'amu*angstrom^2'), symmetry=1, barrier=(23.2234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0898588,'amu*angstrom^2'), symmetry=1, barrier=(23.245,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.402073,0.06774,-3.13604e-05,-3.39656e-09,4.84397e-12,51654.1,28.6556], Tmin=(100,'K'), Tmax=(1108.93,'K')), NASAPolynomial(coeffs=[14.9207,0.0353255,-1.45072e-05,2.69881e-09,-1.88456e-13,47207.1,-48.4218], Tmin=(1108.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Propargyl) + radical(Allyl_S)"""),
)

species(
    label = '[CH]=C[C]=CC(C)=CC(28215)',
    structure = SMILES('[CH]C=C=CC(C)=CC'),
    E0 = (495.543,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.04561,'amu*angstrom^2'), symmetry=1, barrier=(47.0325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05568,'amu*angstrom^2'), symmetry=1, barrier=(47.2641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05326,'amu*angstrom^2'), symmetry=1, barrier=(47.2084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04621,'amu*angstrom^2'), symmetry=1, barrier=(47.0464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0377573,0.0828825,-6.18866e-05,2.45295e-08,-4.03035e-12,59750,28.6712], Tmin=(100,'K'), Tmax=(1410.72,'K')), NASAPolynomial(coeffs=[15.1567,0.0397999,-1.60777e-05,2.88172e-09,-1.94072e-13,55463,-49.8619], Tmin=(1410.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C[C]=C[C]1CC1C(28216)',
    structure = SMILES('C=C[C]=C[C]1CC1C'),
    E0 = (464.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.841198,0.0498247,3.45677e-05,-8.45774e-08,3.75499e-11,55967.9,25.4186], Tmin=(100,'K'), Tmax=(941.119,'K')), NASAPolynomial(coeffs=[17.8081,0.0271207,-7.99707e-06,1.3602e-09,-9.76441e-14,50586.2,-67.0323], Tmin=(941.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(464.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C=CJC=C) + radical(Allyl_T)"""),
)

species(
    label = '[CH2]C(=CC)C=C1[CH]C1(28217)',
    structure = SMILES('[CH2]C(=CC)C=C1[CH]C1'),
    E0 = (457.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.445251,0.0599684,8.26802e-06,-5.71543e-08,2.70939e-11,55194.4,25.1511], Tmin=(100,'K'), Tmax=(969.139,'K')), NASAPolynomial(coeffs=[18.8586,0.027748,-9.62082e-06,1.76211e-09,-1.28057e-13,49569.5,-73.7121], Tmin=(969.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(Methylene_cyclopropane) + radical(Allyl_P) + radical(Allyl_S)"""),
)

species(
    label = '[CH2][C]1C=C(C=C)C1C(28218)',
    structure = SMILES('[CH2]C=C1C=C([CH2])C1C'),
    E0 = (347.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.787685,0.0515025,3.13582e-05,-7.99618e-08,3.53475e-11,41889.8,24.6145], Tmin=(100,'K'), Tmax=(948.511,'K')), NASAPolynomial(coeffs=[17.2705,0.0296335,-9.39864e-06,1.63826e-09,-1.17067e-13,36619.9,-65.3313], Tmin=(948.511,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclobutene) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'CC=C1C=[C][CH]CC1(28106)',
    structure = SMILES('CC=C1[CH][C]=CCC1'),
    E0 = (369.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73933,0.0239625,0.000104374,-1.44144e-07,5.3867e-11,44542.1,22.3173], Tmin=(100,'K'), Tmax=(987.575,'K')), NASAPolynomial(coeffs=[14.4487,0.037447,-1.4776e-05,2.89494e-09,-2.16226e-13,38864,-54.8777], Tmin=(987.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexane) + radical(Cds_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C=C=CC(C)=CC(28219)',
    structure = SMILES('C=C=C=CC(C)=CC'),
    E0 = (301.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0673712,0.081858,-6.89193e-05,3.07885e-08,-5.61863e-12,36402.9,26.2374], Tmin=(100,'K'), Tmax=(1296.74,'K')), NASAPolynomial(coeffs=[15.461,0.0343734,-1.39912e-05,2.54917e-09,-1.74274e-13,32410.6,-52.0279], Tmin=(1296.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    label = '[CH2]C(=C)C=[C]C=C(28220)',
    structure = SMILES('[CH2]C(=C)C=[C]C=C'),
    E0 = (458.597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.38097,'amu*angstrom^2'), symmetry=1, barrier=(31.7511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3789,'amu*angstrom^2'), symmetry=1, barrier=(31.7035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38774,'amu*angstrom^2'), symmetry=1, barrier=(31.907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.481321,0.0652591,-3.07872e-05,-1.63299e-08,1.40788e-11,55294.3,23.0451], Tmin=(100,'K'), Tmax=(931.274,'K')), NASAPolynomial(coeffs=[18.5292,0.0185358,-5.13312e-06,8.14283e-10,-5.59085e-14,50597.4,-69.9103], Tmin=(931.274,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[C]=CC[C]=CC(26586)',
    structure = SMILES('C=C[C]=CC[C]=CC'),
    E0 = (541.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.713501,'amu*angstrom^2'), symmetry=1, barrier=(16.4048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714239,'amu*angstrom^2'), symmetry=1, barrier=(16.4218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714825,'amu*angstrom^2'), symmetry=1, barrier=(16.4352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714532,'amu*angstrom^2'), symmetry=1, barrier=(16.4285,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3664.75,'J/mol'), sigma=(6.21721,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=572.42 K, Pc=34.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283233,0.0780811,-6.41802e-05,2.88981e-08,-5.38961e-12,65244.1,30.6758], Tmin=(100,'K'), Tmax=(1261.83,'K')), NASAPolynomial(coeffs=[13.5051,0.036168,-1.43563e-05,2.5746e-09,-1.74308e-13,61907.4,-36.1871], Tmin=(1261.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = 'C=CC1=CC(=CC)C1(27983)',
    structure = SMILES('C=CC1=CC(=CC)C1'),
    E0 = (233.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.515583,0.0538086,3.35545e-05,-9.06582e-08,4.10441e-11,28278.6,24.2553], Tmin=(100,'K'), Tmax=(945.902,'K')), NASAPolynomial(coeffs=[21.6779,0.0214945,-5.87183e-06,1.03264e-09,-7.91736e-14,21717.3,-90.1843], Tmin=(945.902,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
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
    label = 'C=C[C]=C[C]=CC(28221)',
    structure = SMILES('[CH2]C=[C]C=C=CC'),
    E0 = (479.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.61077,'amu*angstrom^2'), symmetry=1, barrier=(37.0347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61317,'amu*angstrom^2'), symmetry=1, barrier=(37.09,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61523,'amu*angstrom^2'), symmetry=1, barrier=(37.1372,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.706156,0.0641221,-4.98443e-05,2.10319e-08,-3.6216e-12,57796,25.4495], Tmin=(100,'K'), Tmax=(1380.35,'K')), NASAPolynomial(coeffs=[13.4754,0.0271192,-9.63391e-06,1.61155e-09,-1.04311e-13,54270.8,-40.2707], Tmin=(1380.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(479.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C1[C]=CC(=C)C1C(28157)',
    structure = SMILES('[CH2]C1[C]=CC(=C)C1C'),
    E0 = (496.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96216,0.0484096,3.24666e-05,-7.793e-08,3.39382e-11,59786.1,26.3133], Tmin=(100,'K'), Tmax=(953.031,'K')), NASAPolynomial(coeffs=[16.3375,0.0294919,-9.55254e-06,1.685e-09,-1.20706e-13,54783.9,-57.9922], Tmin=(953.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(3-Methylenecyclopentene) + radical(Isobutyl) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C=C[C]=CC(=C)C=C(27242)',
    structure = SMILES('C=C[C]=CC(=C)C=C'),
    E0 = (415.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,2980,3010,3040,3070,3100,1330,1380,1430,900,975,1050,1000,1025,1050,1600,1650,1700,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.37039,'amu*angstrom^2'), symmetry=1, barrier=(31.5079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36362,'amu*angstrom^2'), symmetry=1, barrier=(31.3523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36465,'amu*angstrom^2'), symmetry=1, barrier=(31.3761,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.443814,0.0658506,-2.32901e-05,-2.00871e-08,1.36437e-11,50119.6,26.0281], Tmin=(100,'K'), Tmax=(963.071,'K')), NASAPolynomial(coeffs=[16.5216,0.028058,-9.57185e-06,1.66689e-09,-1.1545e-13,45678.6,-57.9117], Tmin=(963.071,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]CC(=C)C=[C]C=C(28222)',
    structure = SMILES('[CH2]CC(=C)C=[C]C=C'),
    E0 = (489.699,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.915199,'amu*angstrom^2'), symmetry=1, barrier=(21.0422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.914235,'amu*angstrom^2'), symmetry=1, barrier=(21.0201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.91411,'amu*angstrom^2'), symmetry=1, barrier=(21.0172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.915715,'amu*angstrom^2'), symmetry=1, barrier=(21.0541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.171361,0.0808675,-5.41203e-05,6.38684e-09,5.40683e-12,59057.1,30.4981], Tmin=(100,'K'), Tmax=(965.796,'K')), NASAPolynomial(coeffs=[18.5254,0.0277608,-9.42534e-06,1.61788e-09,-1.10385e-13,54311,-64.9266], Tmin=(965.796,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(RCCJ)"""),
)

species(
    label = 'C=C[C]=[C]C(=C)CC(28223)',
    structure = SMILES('C=C[C]=[C]C(=C)CC'),
    E0 = (483.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.14274,'amu*angstrom^2'), symmetry=1, barrier=(26.2738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14299,'amu*angstrom^2'), symmetry=1, barrier=(26.2795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14131,'amu*angstrom^2'), symmetry=1, barrier=(26.241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14076,'amu*angstrom^2'), symmetry=1, barrier=(26.2284,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.593324,0.0902659,-8.51475e-05,4.27855e-08,-8.51279e-12,58319.8,29.8238], Tmin=(100,'K'), Tmax=(1236.21,'K')), NASAPolynomial(coeffs=[18.5575,0.0278746,-8.92721e-06,1.40312e-09,-8.77658e-14,53617.4,-66.4976], Tmin=(1236.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C(C=[C]C=C)CC(28224)',
    structure = SMILES('[CH]=C(C=[C]C=C)CC'),
    E0 = (531.549,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.966536,'amu*angstrom^2'), symmetry=1, barrier=(22.2226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.964927,'amu*angstrom^2'), symmetry=1, barrier=(22.1856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.967055,'amu*angstrom^2'), symmetry=1, barrier=(22.2345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966455,'amu*angstrom^2'), symmetry=1, barrier=(22.2207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.283111,0.0836655,-6.03568e-05,1.15554e-08,3.9116e-12,64094.2,29.5159], Tmin=(100,'K'), Tmax=(962.692,'K')), NASAPolynomial(coeffs=[18.9884,0.0271947,-9.144e-06,1.55804e-09,-1.05824e-13,59290,-68.406], Tmin=(962.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(531.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C([CH]C)C=CC=C(28225)',
    structure = SMILES('[CH]C(C=CC=C)=CC'),
    E0 = (442.761,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.13907,0.077334,-3.266e-05,-1.50056e-08,1.20255e-11,53413.1,29.0019], Tmin=(100,'K'), Tmax=(987.472,'K')), NASAPolynomial(coeffs=[18.0164,0.0339045,-1.24328e-05,2.22109e-09,-1.54422e-13,48359.3,-65.7928], Tmin=(987.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C][C]=CC(=C)CC(28226)',
    structure = SMILES('[CH2]C#CC=C([CH2])CC'),
    E0 = (408.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.284135,0.0717819,-4.48084e-05,1.05984e-08,2.2619e-13,49288.1,29.8352], Tmin=(100,'K'), Tmax=(1147.72,'K')), NASAPolynomial(coeffs=[14.6385,0.034857,-1.36738e-05,2.46026e-09,-1.67787e-13,45130.2,-45.1532], Tmin=(1147.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(CTCC=CCJ) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C[C]=CC(=C)CC(28227)',
    structure = SMILES('[CH]C=C=CC(=C)CC'),
    E0 = (508.923,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.195375,0.081013,-5.02743e-05,9.29527e-09,1.82515e-12,61369.9,30.344], Tmin=(100,'K'), Tmax=(1080.64,'K')), NASAPolynomial(coeffs=[16.4166,0.0374271,-1.46253e-05,2.63394e-09,-1.80559e-13,56734.2,-55.9244], Tmin=(1080.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC1=C[C]([CH]C)C1(28228)',
    structure = SMILES('[CH2]C=C1[CH]C(=CC)C1'),
    E0 = (377.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3376,0.0309278,9.44757e-05,-1.47593e-07,5.94761e-11,45500.8,26.8349], Tmin=(100,'K'), Tmax=(949.608,'K')), NASAPolynomial(coeffs=[19.1155,0.0257866,-7.57045e-06,1.38996e-09,-1.07729e-13,38979.8,-74.5711], Tmin=(949.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(C=CCJC=C) + radical(Allyl_P)"""),
)

species(
    label = 'C=C1C=[C][CH]CC1C(28123)',
    structure = SMILES('[CH2]C1=C[C]=CCC1C'),
    E0 = (334.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16713,0.0431696,4.74017e-05,-9.05025e-08,3.75003e-11,40364.4,23.3485], Tmin=(100,'K'), Tmax=(957.789,'K')), NASAPolynomial(coeffs=[14.9268,0.0331464,-1.11983e-05,2.0004e-09,-1.43036e-13,35552.6,-53.8], Tmin=(957.789,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC=CC(=C)C=C(27260)',
    structure = SMILES('C=CC=CC(=C)C=C'),
    E0 = (216.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.526365,0.0592493,6.02658e-06,-5.32876e-08,2.55662e-11,26187.3,26.1353], Tmin=(100,'K'), Tmax=(966.002,'K')), NASAPolynomial(coeffs=[18.0461,0.0279986,-9.56972e-06,1.72854e-09,-1.24399e-13,20875.8,-67.7542], Tmin=(966.002,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C=CC(=C)CC(28229)',
    structure = SMILES('C=C=C=CC(=C)CC'),
    E0 = (314.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.153163,0.0807073,-5.96724e-05,1.83521e-08,-8.26442e-13,38025.4,28.1369], Tmin=(100,'K'), Tmax=(1082.51,'K')), NASAPolynomial(coeffs=[17.6361,0.0305443,-1.17382e-05,2.11883e-09,-1.46018e-13,33261.7,-63.3101], Tmin=(1082.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=[C]C(C)C=[C]C=C(26582)',
    structure = SMILES('C=[C]C(C)C=[C]C=C'),
    E0 = (543.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.7818,'amu*angstrom^2'), symmetry=1, barrier=(17.9751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781839,'amu*angstrom^2'), symmetry=1, barrier=(17.976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781771,'amu*angstrom^2'), symmetry=1, barrier=(17.9745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781781,'amu*angstrom^2'), symmetry=1, barrier=(17.9747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3606.27,'J/mol'), sigma=(6.17911,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.29 K, Pc=34.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.194054,0.0844721,-8.08726e-05,4.40912e-08,-1.00072e-11,65539.2,28.7822], Tmin=(100,'K'), Tmax=(1049.12,'K')), NASAPolynomial(coeffs=[12.0739,0.039178,-1.61131e-05,2.93984e-09,-2.0115e-13,63046.5,-29.1009], Tmin=(1049.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC1=CC(=C)C1C(28073)',
    structure = SMILES('C=CC1=CC(=C)C1C'),
    E0 = (236.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.257097,0.0622633,9.42245e-06,-6.56652e-08,3.2237e-11,28580.8,21.5772], Tmin=(100,'K'), Tmax=(946.791,'K')), NASAPolynomial(coeffs=[21.4906,0.0224133,-6.43173e-06,1.11665e-09,-8.27781e-14,22325.4,-91.5021], Tmin=(946.791,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
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
    label = 'C=[C]C=[C]C=C(27475)',
    structure = SMILES('[CH2]C=[C]C=C=C'),
    E0 = (515.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.89636,'amu*angstrom^2'), symmetry=1, barrier=(43.601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89613,'amu*angstrom^2'), symmetry=1, barrier=(43.5957,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60107,0.0457135,-1.98528e-05,-8.88694e-09,7.69604e-12,62096.2,20.4258], Tmin=(100,'K'), Tmax=(942.99,'K')), NASAPolynomial(coeffs=[11.6413,0.0202885,-6.7122e-06,1.12524e-09,-7.5614e-14,59439.5,-31.4698], Tmin=(942.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(515.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=[C]C1CC1=CC(27965)',
    structure = SMILES('[CH2]C=[C]C1CC1=CC'),
    E0 = (575.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.263379,0.0739729,-5.27636e-05,1.94392e-08,-2.92523e-12,69410.2,27.6047], Tmin=(100,'K'), Tmax=(1543.9,'K')), NASAPolynomial(coeffs=[15.9958,0.033213,-1.31631e-05,2.33965e-09,-1.56359e-13,64552.3,-55.1282], Tmin=(1543.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Methylene_cyclopropane) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1([CH]C)C=C=CC1(28181)',
    structure = SMILES('[CH2]C1([CH]C)C=C=CC1'),
    E0 = (747.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.791425,0.0568324,-2.46446e-06,-3.08311e-08,1.36261e-11,89993.3,28.1857], Tmin=(100,'K'), Tmax=(1075.01,'K')), NASAPolynomial(coeffs=[14.2247,0.0365383,-1.55743e-05,2.98977e-09,-2.13653e-13,85389.6,-45.5727], Tmin=(1075.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(747.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene) + radical(Cs_S) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C(C=C=[C]C)=CC(28230)',
    structure = SMILES('[CH2]C([CH]C)=CC#CC'),
    E0 = (393.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.549269,0.0657551,-2.65495e-05,-9.636e-09,7.99114e-12,47439.7,26.7139], Tmin=(100,'K'), Tmax=(1005.61,'K')), NASAPolynomial(coeffs=[13.3842,0.0357533,-1.3199e-05,2.33064e-09,-1.59138e-13,43793.9,-40.5714], Tmin=(1005.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(CTCC=CCJ) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C([C]=C=CC)=CC(28231)',
    structure = SMILES('[CH2]C(C#C[CH]C)=CC'),
    E0 = (427.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.486715,0.0692701,-3.90669e-05,2.85776e-09,4.09536e-12,51521.3,29.5558], Tmin=(100,'K'), Tmax=(980.736,'K')), NASAPolynomial(coeffs=[12.9294,0.0352066,-1.2487e-05,2.13667e-09,-1.42716e-13,48278.3,-34.322], Tmin=(980.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Allyl_P) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH2]C(=[C]C)C=C=CC(28232)',
    structure = SMILES('[CH2]C(=[C]C)C=C=CC'),
    E0 = (514.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.817417,'amu*angstrom^2'), symmetry=1, barrier=(18.794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.821067,'amu*angstrom^2'), symmetry=1, barrier=(18.8779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.818313,'amu*angstrom^2'), symmetry=1, barrier=(18.8146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.818047,'amu*angstrom^2'), symmetry=1, barrier=(18.8085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.078409,0.083347,-7.27894e-05,3.41193e-08,-6.55927e-12,61987.4,27.898], Tmin=(100,'K'), Tmax=(1232.44,'K')), NASAPolynomial(coeffs=[14.7831,0.0356211,-1.47017e-05,2.6974e-09,-1.85283e-13,58362.9,-46.1165], Tmin=(1232.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(514.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C](C=C)C=C=CC(27269)',
    structure = SMILES('[CH2]C(C=C)=C[C]=CC'),
    E0 = (389.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0777078,0.0748848,-3.83494e-05,-8.11842e-09,1.02399e-11,46952.8,28.2694], Tmin=(100,'K'), Tmax=(952.813,'K')), NASAPolynomial(coeffs=[17.1906,0.029961,-1.00027e-05,1.69781e-09,-1.15228e-13,42469.9,-59.875], Tmin=(952.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C(=CC)C1[C]=CC1(28197)',
    structure = SMILES('[CH2]C(=CC)C1[C]=CC1'),
    E0 = (545.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.448191,0.0672177,-3.16541e-05,-2.52403e-09,4.55203e-12,65694.3,27.01], Tmin=(100,'K'), Tmax=(1100.31,'K')), NASAPolynomial(coeffs=[14.4049,0.0354718,-1.42665e-05,2.62744e-09,-1.82514e-13,61473.3,-46.8815], Tmin=(1100.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Allyl_P) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'C#C[CH]C([CH2])=CC(28200)',
    structure = SMILES('C#CC=C([CH2])[CH]C'),
    E0 = (435.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.417901,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292485,'amu*angstrom^2'), symmetry=1, barrier=(20.3659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.418421,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.421287,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.951004,0.0577526,-2.77367e-05,-4.57001e-09,5.78378e-12,52497.7,21.8521], Tmin=(100,'K'), Tmax=(1022.91,'K')), NASAPolynomial(coeffs=[13.4773,0.0273447,-1.03855e-05,1.87403e-09,-1.29861e-13,48963.2,-43.6142], Tmin=(1022.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CTCC=CCJ) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C=[C]C1C(=C)C1C(28089)',
    structure = SMILES('[CH2]C=[C]C1C(=C)C1C'),
    E0 = (580.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281391,0.0704847,-3.57132e-05,-2.51723e-09,5.61998e-12,69926.7,26.3372], Tmin=(100,'K'), Tmax=(1048.18,'K')), NASAPolynomial(coeffs=[15.6153,0.0332581,-1.29068e-05,2.35498e-09,-1.63834e-13,65542.6,-53.9404], Tmin=(1048.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C([CH]C)C=C=CC(26418)',
    structure = SMILES('[CH]C(C=C=CC)=CC'),
    E0 = (495.543,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.04561,'amu*angstrom^2'), symmetry=1, barrier=(47.0325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05568,'amu*angstrom^2'), symmetry=1, barrier=(47.2641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05326,'amu*angstrom^2'), symmetry=1, barrier=(47.2084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04621,'amu*angstrom^2'), symmetry=1, barrier=(47.0464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0377573,0.0828825,-6.18866e-05,2.45295e-08,-4.03035e-12,59750,28.6712], Tmin=(100,'K'), Tmax=(1410.72,'K')), NASAPolynomial(coeffs=[15.1567,0.0397999,-1.60777e-05,2.88172e-09,-1.94072e-13,55463,-49.8619], Tmin=(1410.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC(=C)C=C=CC(27275)',
    structure = SMILES('C=CC(=C)C=C=CC'),
    E0 = (269.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.473627,0.0666108,-2.96241e-05,-5.13706e-09,5.66211e-12,32531.2,26.3588], Tmin=(100,'K'), Tmax=(1074.63,'K')), NASAPolynomial(coeffs=[14.2838,0.0353047,-1.39804e-05,2.5622e-09,-1.77929e-13,28402.5,-46.661], Tmin=(1074.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    E0 = (422.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (633.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (542.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (626.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (681.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (629.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (644.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (634.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (637.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (670.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (655.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (692.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (864.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (572.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (455.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (570.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (627.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (812.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (516.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (560.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (560.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (510.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (447.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (878.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (711.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (430.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (861.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (496.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (627.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (603.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (656.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (676.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (846.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (600.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (654.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (560.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (448.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (461.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (430.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (638.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (430.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (859.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (575.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (747.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (630.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (660.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (676.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (485.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (545.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (817.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (580.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (627.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (447.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['CH3CHCCH2(18175)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2]C1([CH]C)C=C1C=C(28131)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.126e+14,'s^-1'), n=-0.355, Ea=(211.184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingleDe] for rate rule [R4_D_D;doublebond_intra_HNd;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 210.1 to 211.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2]C1[C]=CC(=CC)C1(27993)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.40382e+09,'s^-1'), n=0.352, Ea=(120.148,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(C#CC=C)=CC(28207)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.9e+08,'cm^3/(mol*s)'), n=1.64, Ea=(12.7612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2707 used for Ct-Cd_Ct-Cd;HJ
Exact match found for rate rule [Ct-Cd_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C(C=C=C=C)=CC(28208)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]C=C(4699)', 'CH3CHCCH2(18175)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C][CH]C(18176)', 'CH2CHCCH(26391)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0117197,'m^3/(mol*s)'), n=2.33233, Ea=(9.74423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C([C]=CC=C)=CC(28209)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.0945e+07,'s^-1'), n=1.951, Ea=(212.263,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_single] for rate rule [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2]C([CH]C=C=C)=CC(28184)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C[C]=CC(C)=[C]C(28210)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C[C]=[C]C(C)=CC(28211)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.18351e+08,'s^-1'), n=1.55833, Ea=(185.212,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Cd_rad_out;Cs_H_out_2H] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out_2H] for rate rule [R3H_SS_2Cd;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC=CC([CH2])=CC(28212)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(=[C]C)C=CC=C(28213)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_single] for rate rule [R4H_DSD;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['C=C[C]=C[C](C)C=C(27245)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2][C](C=C)C=CC=C(27249)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSMS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C][C]=CC(C)=CC(28214)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C[C]=CC(C)=CC(28215)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]C=C(4699)', 'C=[C][CH]C(18176)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['C=C[C]=C[C]1CC1C(28216)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.68393e+12,'s^-1'), n=-0.105173, Ea=(93.5715,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe_HNd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2]C(=CC)C=C1[CH]C1(28217)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2][C]1C=C(C=C)C1C(28218)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['CC=C1C=[C][CH]CC1(28106)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.63396e+08,'s^-1'), n=0.483333, Ea=(87.4777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['C=C=C=CC(C)=CC(28219)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(S)(23)', '[CH2]C(=C)C=[C]C=C(28220)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C[C]=CC[C]=CC(26586)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['C=CC1=CC(=CC)C1(27983)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(19)', 'C=C[C]=C[C]=CC(28221)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2]C1[C]=CC(=C)C1C(28157)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(187000,'s^-1'), n=1.48, Ea=(73.4697,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 69.0 to 73.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', 'C=C[C]=CC(=C)C=C(27242)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]CC(=C)C=[C]C=C(28222)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C[C]=[C]C(=C)CC(28223)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.44583e+10,'s^-1'), n=0.916667, Ea=(172.729,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Cd_rad_out;Cs_H_out_H/NonDeC] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_2Cd;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C(C=[C]C=C)CC(28224)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C([CH]C)C=CC=C(28225)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_single] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['C=[C][C]=CC(=C)CC(28226)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.37e+07,'s^-1'), n=1.41, Ea=(177.82,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;C_rad_out_H/NonDeC;Cd_H_out_doubleC] for rate rule [R5HJ_3;C_rad_out_H/NonDeC;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C[C]=CC(=C)CC(28227)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['C=CC1=C[C]([CH]C)C1(28228)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['C=C1C=[C][CH]CC1C(28123)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(487000,'s^-1'), n=1.17, Ea=(26.3592,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['C=CC=CC(=C)C=C(27260)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['C=C=C=CC(=C)CC(28229)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['C=CC1=CC(=C)C1C(28073)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriDe_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CHCH3(T)(95)', 'C=[C]C=[C]C=C(27475)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2]C=[C]C1CC1=CC(27965)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.50867e+09,'s^-1'), n=0.788889, Ea=(153.364,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R4_S_(Cd)_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 151.9 to 153.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2]C1([CH]C)C=C=CC1(28181)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(4.14729e+08,'s^-1'), n=0.705, Ea=(324.618,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;doublebond_intra_HNd;radadd_intra_cs2H] + [R6_SMM;doublebond_intra;radadd_intra_cs2H] for rate rule [R6_SMM;doublebond_intra_HNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 320.9 to 324.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2]C(C=C=[C]C)=CC(28230)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2]C([C]=C=CC)=CC(28231)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(0.0029655,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C(=[C]C)C=C=CC(28232)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R6H;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2][C](C=C)C=C=CC(27269)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(64.2,'s^-1'), n=2.1, Ea=(63.1784,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 115 used for R7H;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R7H;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2]C(=CC)C1[C]=CC1(28197)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(122.502,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 120.2 to 122.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction50',
    reactants = ['CH2(19)', 'C#C[CH]C([CH2])=CC(28200)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['[CH2]C=[C]C1C(=C)C1C(28089)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(157.638,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_csHNd] for rate rule [R4_S_(Cd)_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 155.9 to 157.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=C([CH]C)C=C=CC(26418)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6H;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C=C[C]=CC(=C)[CH]C(26578)'],
    products = ['C=CC(=C)C=C=CC(27275)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

network(
    label = '4550',
    isomers = [
        'C=C[C]=CC(=C)[CH]C(26578)',
    ],
    reactants = [
        ('CH3CHCCH2(18175)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4550',
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

