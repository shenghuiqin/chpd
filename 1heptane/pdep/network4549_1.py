species(
    label = 'C#CC([CH2])C(=C)[CH]C(26577)',
    structure = SMILES('C#CC([CH2])C([CH2])=CC'),
    E0 = (478.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,750,770,3400,2100,273.559],'cm^-1')),
        HinderedRotor(inertia=(0.20589,'amu*angstrom^2'), symmetry=1, barrier=(10.9289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.2534,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25263,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205847,'amu*angstrom^2'), symmetry=1, barrier=(10.929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25256,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.119174,0.0887879,-8.75458e-05,4.76961e-08,-1.05787e-11,57667.9,29.155], Tmin=(100,'K'), Tmax=(1086.04,'K')), NASAPolynomial(coeffs=[14.3672,0.0354333,-1.38548e-05,2.46114e-09,-1.65913e-13,54521.3,-41.9294], Tmin=(1086.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Allyl_P)"""),
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
    label = 'C#CC1CC1([CH2])[CH]C(28161)',
    structure = SMILES('C#CC1CC1([CH2])[CH]C'),
    E0 = (586.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579915,0.0597992,-9.59294e-07,-4.24954e-08,2.08555e-11,70661,29.7444], Tmin=(100,'K'), Tmax=(979.017,'K')), NASAPolynomial(coeffs=[16.7225,0.0300351,-1.08047e-05,1.96693e-09,-1.40153e-13,65765.9,-56.6494], Tmin=(979.017,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(Neopentyl) + radical(Cs_S)"""),
)

species(
    label = '[CH]=C1CC1C([CH2])=CC(28188)',
    structure = SMILES('[CH]=C1CC1C([CH2])=CC'),
    E0 = (582.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.145889,0.0792655,-6.11168e-05,2.4292e-08,-3.88012e-12,70177,27.9595], Tmin=(100,'K'), Tmax=(1488.97,'K')), NASAPolynomial(coeffs=[18.2896,0.0297397,-1.12238e-05,1.95292e-09,-1.29321e-13,64687,-68.3197], Tmin=(1488.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C1CC(=CC)C1[CH2](27973)',
    structure = SMILES('[CH]=C1CC(=CC)C1[CH2]'),
    E0 = (578.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.675307,0.0598105,-6.09802e-06,-3.33792e-08,1.68688e-11,69735.4,28.5336], Tmin=(100,'K'), Tmax=(985.123,'K')), NASAPolynomial(coeffs=[14.787,0.0331075,-1.20261e-05,2.15999e-09,-1.51096e-13,65470.4,-46.8712], Tmin=(985.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(578.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_P) + radical(Isobutyl)"""),
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
    label = 'C#CC(=C)C([CH2])=CC(28189)',
    structure = SMILES('C#CC(=C)C([CH2])=CC'),
    E0 = (401.291,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,750,770,3400,2100,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.2066,'amu*angstrom^2'), symmetry=1, barrier=(27.7422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20834,'amu*angstrom^2'), symmetry=1, barrier=(27.7822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20664,'amu*angstrom^2'), symmetry=1, barrier=(27.7429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20778,'amu*angstrom^2'), symmetry=1, barrier=(27.7692,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0056467,0.0751701,-4.24627e-05,-4.06682e-09,8.35668e-12,48419.8,26.6058], Tmin=(100,'K'), Tmax=(994.465,'K')), NASAPolynomial(coeffs=[19.1877,0.0256074,-9.39261e-06,1.71028e-09,-1.21159e-13,43235.7,-72.7558], Tmin=(994.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Ct) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Allyl_P)"""),
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
    label = '[CH2]C(C=C)=CC(24210)',
    structure = SMILES('[CH2]C(C=C)=CC'),
    E0 = (171.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.902009,'amu*angstrom^2'), symmetry=1, barrier=(20.739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.902046,'amu*angstrom^2'), symmetry=1, barrier=(20.7398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.901955,'amu*angstrom^2'), symmetry=1, barrier=(20.7377,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17187,0.0502492,-6.01379e-06,-2.94563e-08,1.53252e-11,20775.8,20.8331], Tmin=(100,'K'), Tmax=(976.501,'K')), NASAPolynomial(coeffs=[14.417,0.0239289,-8.49401e-06,1.53256e-09,-1.08637e-13,16857.1,-49.5715], Tmin=(976.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P)"""),
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
    label = 'C#C[C](C)C([CH2])=CC(28190)',
    structure = SMILES('C#CC(C)=C([CH2])[CH]C'),
    E0 = (395.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.433095,0.0746052,-5.68201e-05,2.37605e-08,-4.16177e-12,47645.5,25.2835], Tmin=(100,'K'), Tmax=(1324.72,'K')), NASAPolynomial(coeffs=[12.6988,0.0375688,-1.4883e-05,2.65558e-09,-1.78847e-13,44395.8,-37.3403], Tmin=(1324.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.054,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtCs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CTCC=CCJ) + radical(Allyl_S)"""),
)

species(
    label = 'C#C[C]([CH2])C(C)=CC(28191)',
    structure = SMILES('C#CC([CH2])=C(C)[CH]C'),
    E0 = (425.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0652111,0.0779655,-5.89652e-05,2.30559e-08,-3.6375e-12,51288.1,27.7197], Tmin=(100,'K'), Tmax=(1501.02,'K')), NASAPolynomial(coeffs=[17.7225,0.0305637,-1.15955e-05,2.01693e-09,-1.33366e-13,45948.1,-65.3199], Tmin=(1501.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtCs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'C#CC([CH2])C(C)=[C]C(28192)',
    structure = SMILES('C#CC([CH2])C(C)=[C]C'),
    E0 = (564.574,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,750,770,3400,2100,217.043],'cm^-1')),
        HinderedRotor(inertia=(0.324075,'amu*angstrom^2'), symmetry=1, barrier=(10.8268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.323976,'amu*angstrom^2'), symmetry=1, barrier=(10.827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.323962,'amu*angstrom^2'), symmetry=1, barrier=(10.8268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.323935,'amu*angstrom^2'), symmetry=1, barrier=(10.8266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97991,'amu*angstrom^2'), symmetry=1, barrier=(66.1633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.144112,0.0957746,-0.000120541,9.01514e-08,-2.77085e-11,68047.4,29.4398], Tmin=(100,'K'), Tmax=(823.938,'K')), NASAPolynomial(coeffs=[10.489,0.04168,-1.75569e-05,3.18119e-09,-2.14291e-13,66379.2,-19.29], Tmin=(823.938,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(C)C([CH2])=[C]C(28193)',
    structure = SMILES('C#CC(C)C([CH2])=[C]C'),
    E0 = (510.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,750,770,3400,2100,214.155],'cm^-1')),
        HinderedRotor(inertia=(0.465745,'amu*angstrom^2'), symmetry=1, barrier=(15.1577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.465749,'amu*angstrom^2'), symmetry=1, barrier=(15.1577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.533594,'amu*angstrom^2'), symmetry=1, barrier=(17.3658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.21014,'amu*angstrom^2'), symmetry=1, barrier=(71.9291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.465747,'amu*angstrom^2'), symmetry=1, barrier=(15.1577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0602892,0.0888894,-8.95969e-05,5.0727e-08,-1.19088e-11,61598,27.2103], Tmin=(100,'K'), Tmax=(1016.84,'K')), NASAPolynomial(coeffs=[12.5549,0.039738,-1.70901e-05,3.18921e-09,-2.21058e-13,59057.1,-33.2775], Tmin=(1016.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C#CC([CH2])[C](C)C=C(27304)',
    structure = SMILES('C#CC([CH2])C(C)=C[CH2]'),
    E0 = (478.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,750,770,3400,2100,273.549],'cm^-1')),
        HinderedRotor(inertia=(0.205842,'amu*angstrom^2'), symmetry=1, barrier=(10.9286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25261,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25479,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205694,'amu*angstrom^2'), symmetry=1, barrier=(10.9292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25325,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.11917,0.0887879,-8.75457e-05,4.7696e-08,-1.05786e-11,57667.9,29.155], Tmin=(100,'K'), Tmax=(1086.06,'K')), NASAPolynomial(coeffs=[14.3673,0.0354333,-1.38548e-05,2.46114e-09,-1.65912e-13,54521.3,-41.9295], Tmin=(1086.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(Isobutyl)"""),
)

species(
    label = '[C]#CC(C)C([CH2])=CC(28194)',
    structure = SMILES('[C]#CC(C)C([CH2])=CC'),
    E0 = (610.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,272.945,272.959],'cm^-1')),
        HinderedRotor(inertia=(0.199624,'amu*angstrom^2'), symmetry=1, barrier=(10.5537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199622,'amu*angstrom^2'), symmetry=1, barrier=(10.5537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199618,'amu*angstrom^2'), symmetry=1, barrier=(10.5537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551799,'amu*angstrom^2'), symmetry=1, barrier=(29.1704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49121,'amu*angstrom^2'), symmetry=1, barrier=(78.8337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0356133,0.0888033,-9.08766e-05,5.30193e-08,-1.27982e-11,73542.7,27.5658], Tmin=(100,'K'), Tmax=(993.307,'K')), NASAPolynomial(coeffs=[12.296,0.039431,-1.63184e-05,2.97846e-09,-2.0362e-13,71107.1,-31.501], Tmin=(993.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Allyl_P)"""),
)

species(
    label = '[C]#CC([CH2])C(C)=CC(28195)',
    structure = SMILES('[C]#CC([CH2])C(C)=CC'),
    E0 = (663.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,256.105,256.106],'cm^-1')),
        HinderedRotor(inertia=(0.184878,'amu*angstrom^2'), symmetry=1, barrier=(8.60541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184881,'amu*angstrom^2'), symmetry=1, barrier=(8.60543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184882,'amu*angstrom^2'), symmetry=1, barrier=(8.60541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184888,'amu*angstrom^2'), symmetry=1, barrier=(8.60541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41244,'amu*angstrom^2'), symmetry=1, barrier=(65.7405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.135906,0.0953045,-0.000120467,9.06027e-08,-2.77377e-11,79990.6,29.6769], Tmin=(100,'K'), Tmax=(867.038,'K')), NASAPolynomial(coeffs=[10.2997,0.0412444,-1.67065e-05,2.95102e-09,-1.95189e-13,78403.4,-17.8977], Tmin=(867.038,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(663.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC(C)[C]([CH2])C=C(27307)',
    structure = SMILES('C#CC(C)C([CH2])=C[CH2]'),
    E0 = (424.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.420985,0.0878918,-7.7617e-05,3.56373e-08,-6.53844e-12,51240.4,28.7396], Tmin=(100,'K'), Tmax=(1313.17,'K')), NASAPolynomial(coeffs=[18.6246,0.0298784,-1.13508e-05,1.99585e-09,-1.33912e-13,46238.3,-68.3337], Tmin=(1313.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C#CC([CH2])[C]1CC1C(28196)',
    structure = SMILES('C#CC([CH2])[C]1CC1C'),
    E0 = (580.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.669857,0.064943,-2.93956e-05,-6.52415e-09,7.46031e-12,69927.8,30.4812], Tmin=(100,'K'), Tmax=(959.768,'K')), NASAPolynomial(coeffs=[12.4246,0.0349174,-1.21084e-05,2.05553e-09,-1.37177e-13,66797.9,-30.2962], Tmin=(959.768,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(Isobutyl) + radical(Tertalkyl)"""),
)

species(
    label = '[CH]=[C]C1CC(C)C1=C(28127)',
    structure = SMILES('[CH]=[C]C1CC(C)C1=C'),
    E0 = (625.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.9244,-0.001505,0.000117043,-1.16756e-07,2.67166e-11,74817.7,-18.6413], Tmin=(100,'K'), Tmax=(1755.11,'K')), NASAPolynomial(coeffs=[83.1059,0.0264574,-7.13474e-05,1.72839e-08,-1.2761e-12,19836.3,-491.931], Tmin=(1755.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(625.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Cds_P) + radical(Cds_S)"""),
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
    label = '[CH2]C1[C]=CCC1=CC(27951)',
    structure = SMILES('[CH2]C1[C]=CCC1=CC'),
    E0 = (508.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.873141,0.0601265,-2.56443e-05,-1.26391e-10,2.11445e-12,61261.7,28.1928], Tmin=(100,'K'), Tmax=(1226.28,'K')), NASAPolynomial(coeffs=[11.0884,0.0393858,-1.56625e-05,2.81296e-09,-1.90323e-13,57810.5,-27.0305], Tmin=(1226.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(4-Methylenecyclopentene) + radical(cyclopentene-vinyl) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC(=C)C(C)=CC(28198)',
    structure = SMILES('C#CC(=C)C(C)=CC'),
    E0 = (249.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0629506,0.0783113,-5.31144e-05,1.17514e-08,1.42569e-12,30198.7,26.5329], Tmin=(100,'K'), Tmax=(1062.22,'K')), NASAPolynomial(coeffs=[17.4482,0.0307437,-1.18891e-05,2.16216e-09,-1.50029e-13,25442,-63.8846], Tmin=(1062.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Ct) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
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
    label = 'C#CC([CH2])C([CH2])=C(28199)',
    structure = SMILES('C#CC([CH2])C([CH2])=C'),
    E0 = (514.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,414.588],'cm^-1')),
        HinderedRotor(inertia=(3.33367,'amu*angstrom^2'), symmetry=1, barrier=(76.6477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361083,'amu*angstrom^2'), symmetry=1, barrier=(14.3649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.629114,'amu*angstrom^2'), symmetry=1, barrier=(76.6461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194461,'amu*angstrom^2'), symmetry=1, barrier=(23.7057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.333188,0.0755357,-7.5244e-05,4.0144e-08,-8.53195e-12,61987.2,25.0291], Tmin=(100,'K'), Tmax=(1146.09,'K')), NASAPolynomial(coeffs=[14.8868,0.0247408,-8.76231e-06,1.47157e-09,-9.60599e-14,58651.3,-47.1677], Tmin=(1146.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(514.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = 'C#C[CH]CC(=C)[CH]C(26576)',
    structure = SMILES('C#C[CH]CC([CH2])=CC'),
    E0 = (452.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,280.958,4000],'cm^-1')),
        HinderedRotor(inertia=(0.31491,'amu*angstrom^2'), symmetry=1, barrier=(17.4021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.69638,'amu*angstrom^2'), symmetry=1, barrier=(93.7549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.69715,'amu*angstrom^2'), symmetry=1, barrier=(93.852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.69066,'amu*angstrom^2'), symmetry=1, barrier=(93.8213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.68589,'amu*angstrom^2'), symmetry=1, barrier=(93.8855,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3617.13,'J/mol'), sigma=(6.33257,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=564.99 K, Pc=32.32 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0773253,0.0821768,-7.17924e-05,3.45969e-08,-6.76576e-12,54568.1,30.1996], Tmin=(100,'K'), Tmax=(1230.35,'K')), NASAPolynomial(coeffs=[15.1429,0.032694,-1.14643e-05,1.90791e-09,-1.23502e-13,50822.9,-46.3841], Tmin=(1230.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC[CH]C([CH2])=CC(28163)',
    structure = SMILES('C#CC[CH]C([CH2])=CC'),
    E0 = (447.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0596728,0.0780646,-6.02488e-05,2.46605e-08,-4.08945e-12,53957.4,29.2035], Tmin=(100,'K'), Tmax=(1434.13,'K')), NASAPolynomial(coeffs=[16.7108,0.0312897,-1.13261e-05,1.91866e-09,-1.25094e-13,49147.2,-57.7514], Tmin=(1434.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'C#CC([CH2])C[C]=CC(26585)',
    structure = SMILES('C#CC([CH2])C[C]=CC'),
    E0 = (603.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,750,770,3400,2100,218.855,219.126],'cm^-1')),
        HinderedRotor(inertia=(0.00350253,'amu*angstrom^2'), symmetry=1, barrier=(0.119786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161609,'amu*angstrom^2'), symmetry=1, barrier=(5.50317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.9251,'amu*angstrom^2'), symmetry=1, barrier=(101.848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170685,'amu*angstrom^2'), symmetry=1, barrier=(5.50478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.01599,'amu*angstrom^2'), symmetry=1, barrier=(102.07,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3619.53,'J/mol'), sigma=(6.34117,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=565.36 K, Pc=32.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.301188,0.0802179,-7.47905e-05,4.08001e-08,-9.26938e-12,72754.1,32.1875], Tmin=(100,'K'), Tmax=(1051.89,'K')), NASAPolynomial(coeffs=[11.4468,0.0378349,-1.43523e-05,2.4957e-09,-1.65688e-13,70409.3,-22.1473], Tmin=(1051.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CC1CCC1=CC(28109)',
    structure = SMILES('C#CC1CCC1=CC'),
    E0 = (277.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.7388,-0.000635209,0.000116702,-1.17013e-07,2.68619e-11,33013.3,-20.3185], Tmin=(100,'K'), Tmax=(1749.67,'K')), NASAPolynomial(coeffs=[82.4395,0.0271591,-7.14814e-05,1.73124e-08,-1.27893e-12,-21422.2,-490.203], Tmin=(1749.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(methylenecyclobutane)"""),
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
    label = 'C#CC([CH2])[C]=CC(28201)',
    structure = SMILES('C#CC([CH2])[C]=CC'),
    E0 = (603.629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.508653,'amu*angstrom^2'), symmetry=1, barrier=(11.6949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.571988,'amu*angstrom^2'), symmetry=1, barrier=(13.1511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.84644,'amu*angstrom^2'), symmetry=1, barrier=(65.4453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.508352,'amu*angstrom^2'), symmetry=1, barrier=(11.688,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.673053,0.0760649,-9.08781e-05,6.39507e-08,-1.85024e-11,72716.9,25.751], Tmin=(100,'K'), Tmax=(838.771,'K')), NASAPolynomial(coeffs=[9.75295,0.0327664,-1.34507e-05,2.41384e-09,-1.62102e-13,71193.7,-16.4584], Tmin=(838.771,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C1C([CH2])C(=C)C1C(28062)',
    structure = SMILES('[CH]=C1C([CH2])C(=C)C1C'),
    E0 = (581.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.396928,0.0684962,-3.10124e-05,-7.426e-09,7.68343e-12,70038.3,25.9269], Tmin=(100,'K'), Tmax=(1006.96,'K')), NASAPolynomial(coeffs=[14.7201,0.0338243,-1.24707e-05,2.21695e-09,-1.52469e-13,66027,-48.8681], Tmin=(1006.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = 'C#CC([CH2])C(=C)C=C(27302)',
    structure = SMILES('C#CC([CH2])C(=C)C=C'),
    E0 = (452.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03125,'amu*angstrom^2'), symmetry=1, barrier=(23.7105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03144,'amu*angstrom^2'), symmetry=1, barrier=(23.7148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03087,'amu*angstrom^2'), symmetry=1, barrier=(23.7018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03116,'amu*angstrom^2'), symmetry=1, barrier=(23.7085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.254927,0.0862255,-7.86394e-05,3.29108e-08,-3.8019e-12,54525.2,28.2262], Tmin=(100,'K'), Tmax=(957.701,'K')), NASAPolynomial(coeffs=[18.0421,0.0256535,-8.59123e-06,1.42886e-09,-9.44497e-14,50293.8,-63.05], Tmin=(957.701,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC([CH2])C(=C)C[CH2](28202)',
    structure = SMILES('C#CC([CH2])C(=C)C[CH2]'),
    E0 = (545.359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,750,770,3400,2100,180,3134.32],'cm^-1')),
        HinderedRotor(inertia=(0.357705,'amu*angstrom^2'), symmetry=1, barrier=(16.3837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0278674,'amu*angstrom^2'), symmetry=1, barrier=(16.3818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.47977,'amu*angstrom^2'), symmetry=1, barrier=(80.0067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712543,'amu*angstrom^2'), symmetry=1, barrier=(16.3828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.47974,'amu*angstrom^2'), symmetry=1, barrier=(80.006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0685429,0.0881258,-8.77578e-05,4.84677e-08,-1.09054e-11,65739.2,31.5861], Tmin=(100,'K'), Tmax=(1070.9,'K')), NASAPolynomial(coeffs=[14.0556,0.0353696,-1.38626e-05,2.46576e-09,-1.66326e-13,62714.1,-37.5223], Tmin=(1070.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = 'C#C[C]([CH2])C(=C)CC(28203)',
    structure = SMILES('C#CC([CH2])=C([CH2])CC'),
    E0 = (405.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0383366,0.0795971,-6.48939e-05,2.85317e-08,-5.11843e-12,48912.1,28.0908], Tmin=(100,'K'), Tmax=(1324.78,'K')), NASAPolynomial(coeffs=[15.3594,0.0333369,-1.2515e-05,2.17303e-09,-1.44251e-13,44852.8,-50.1334], Tmin=(1324.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(405.441,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtCs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Allyl_P) + radical(CTCC=CCJ)"""),
)

species(
    label = '[CH]=C(CC)C([CH2])C#C(28204)',
    structure = SMILES('[CH]=C(CC)C([CH2])C#C'),
    E0 = (587.209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,750,770,3400,2100,296.774],'cm^-1')),
        HinderedRotor(inertia=(0.176524,'amu*angstrom^2'), symmetry=1, barrier=(11.1472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1761,'amu*angstrom^2'), symmetry=1, barrier=(11.1657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14745,'amu*angstrom^2'), symmetry=1, barrier=(73.406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165531,'amu*angstrom^2'), symmetry=1, barrier=(11.1453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.35132,'amu*angstrom^2'), symmetry=1, barrier=(22.1888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.182989,0.0909487,-9.40559e-05,5.36951e-08,-1.24226e-11,70776.4,30.614], Tmin=(100,'K'), Tmax=(1045.61,'K')), NASAPolynomial(coeffs=[14.4735,0.0348822,-1.36277e-05,2.41706e-09,-1.627e-13,67711.3,-40.7491], Tmin=(1045.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(587.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C([CH]C)C(C)C#C(26415)',
    structure = SMILES('[CH]C(=CC)C(C)C#C'),
    E0 = (492.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,750,770,3400,2100,350,440,435,1725,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0939173,0.0867805,-7.35219e-05,3.52309e-08,-7.18591e-12,59354,27.4382], Tmin=(100,'K'), Tmax=(1140.2,'K')), NASAPolynomial(coeffs=[11.8233,0.0456323,-1.93896e-05,3.58051e-09,-2.46327e-13,56679.2,-30.6885], Tmin=(1140.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(492.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]#CC([CH2])C(=C)CC(28205)',
    structure = SMILES('[C]#CC([CH2])C(=C)CC'),
    E0 = (677.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0955444,0.0910583,-0.000100214,6.34548e-08,-1.64035e-11,81601.7,30.6407], Tmin=(100,'K'), Tmax=(938.16,'K')), NASAPolynomial(coeffs=[12.4815,0.0374338,-1.44751e-05,2.52746e-09,-1.67591e-13,79241.8,-29.2332], Tmin=(938.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(677.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=[C]C1CCC1=CC(28110)',
    structure = SMILES('[CH]=[C]C1CCC1=CC'),
    E0 = (620.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[11.4585,-0.00407591,0.000119164,-1.16947e-07,2.65163e-11,74276.4,-19.3863], Tmin=(100,'K'), Tmax=(1769.07,'K')), NASAPolynomial(coeffs=[84.8065,0.0254037,-7.14488e-05,1.72963e-08,-1.27453e-12,17760.4,-501.475], Tmin=(1769.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(620.968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1[C]=CC(C)C1=C(28074)',
    structure = SMILES('[CH2]C1[C]=CC(C)C1=C'),
    E0 = (510.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.890564,0.0655118,-3.97256e-05,1.25488e-08,-1.67154e-12,61551.2,24.5099], Tmin=(100,'K'), Tmax=(1645.37,'K')), NASAPolynomial(coeffs=[11.8064,0.0389747,-1.55333e-05,2.74659e-09,-1.82185e-13,57959,-33.5884], Tmin=(1645.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(4-Methylenecyclopentene) + radical(Isobutyl) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C#CC(=C)C(=C)CC(28206)',
    structure = SMILES('C#CC(=C)C(=C)CC'),
    E0 = (263.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00405388,0.0739208,-3.27625e-05,-1.48196e-08,1.21474e-11,31809,27.4258], Tmin=(100,'K'), Tmax=(989.432,'K')), NASAPolynomial(coeffs=[19.1617,0.0277024,-1.00903e-05,1.83893e-09,-1.30637e-13,26486,-72.5674], Tmin=(989.432,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Ct) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC(C)C(=C)C=C(27316)',
    structure = SMILES('C#CC(C)C(=C)C=C'),
    E0 = (246.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.331557,0.0852241,-6.84851e-05,2.33815e-08,-1.42364e-12,29864.4,26.7198], Tmin=(100,'K'), Tmax=(1025.41,'K')), NASAPolynomial(coeffs=[18.5959,0.028449,-1.03867e-05,1.83279e-09,-1.25409e-13,25085.9,-69.442], Tmin=(1025.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC([CH2])C(C)[C]=C(26581)',
    structure = SMILES('C#CC([CH2])C(C)[C]=C'),
    E0 = (608.07,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,1685,370,346.654,3567.82],'cm^-1')),
        HinderedRotor(inertia=(0.156181,'amu*angstrom^2'), symmetry=1, barrier=(13.1159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976502,'amu*angstrom^2'), symmetry=1, barrier=(82.4596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.979884,'amu*angstrom^2'), symmetry=1, barrier=(82.5209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156716,'amu*angstrom^2'), symmetry=1, barrier=(13.0809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.993785,'amu*angstrom^2'), symmetry=1, barrier=(82.5432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3560.51,'J/mol'), sigma=(6.30174,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=556.14 K, Pc=32.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0301223,0.0807078,-7.08754e-05,3.45879e-08,-6.83665e-12,73285.7,33.5674], Tmin=(100,'K'), Tmax=(1221.04,'K')), NASAPolynomial(coeffs=[14.9207,0.0317301,-1.07075e-05,1.73688e-09,-1.10546e-13,69634.6,-41.5471], Tmin=(1221.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(608.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC1CC(C)C1=C(28126)',
    structure = SMILES('C#CC1CC(C)C1=C'),
    E0 = (281.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.2139,0.00184847,0.000114806,-1.17026e-07,2.71212e-11,33554.2,-19.6082], Tmin=(100,'K'), Tmax=(1735.08,'K')), NASAPolynomial(coeffs=[80.7972,0.0281619,-7.13661e-05,1.72989e-08,-1.28054e-12,-19393.6,-481.024], Tmin=(1735.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(methylenecyclobutane)"""),
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
    label = 'C#CC([CH2])[C]=C(27503)',
    structure = SMILES('C#CC([CH2])[C]=C'),
    E0 = (639.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2175,525,2950,3100,1380,975,1025,1650,750,770,3400,2100,3000,3100,440,815,1455,1000,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.564352,'amu*angstrom^2'), symmetry=1, barrier=(12.9756,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2595,0.0612351,-7.30673e-05,4.92288e-08,-1.33647e-11,77030.5,21.837], Tmin=(100,'K'), Tmax=(899.573,'K')), NASAPolynomial(coeffs=[9.94523,0.0226112,-8.6597e-06,1.494e-09,-9.79377e-14,75467.9,-19.1468], Tmin=(899.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_S)"""),
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
    E0 = (478.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (586.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (582.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (578.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (631.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (652.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (751.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (629.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (576.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (680.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (726.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (555.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (628.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (706.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (703.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (533.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (812.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (709.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (627.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (608.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (536.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (541.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (934.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (638.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (638.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (773.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (486.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (817.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (985.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (581.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (663.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (659.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (687.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (732.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (536.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (720.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (623.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (536.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (541.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (503.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (702.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (486.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (983.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['CH3CHCCH2(18175)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#CC1CC1([CH2])[CH]C(28161)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(108.139,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 104.1 to 108.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['[CH]=C1CC1C([CH2])=CC(28188)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(103.929,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for R4_S_T;triplebond_intra_H;radadd_intra_cs2H
Exact match found for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 102.0 to 103.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['[CH]=C1CC(=CC)C1[CH2](27973)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.46159e+06,'s^-1'), n=1.55572, Ea=(100.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 96.4 to 100.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#CC(=C)C([CH2])=CC(28189)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.36e+08,'cm^3/(mol*s)'), n=1.64, Ea=(17.9912,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2645 used for Cds-CdCt_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCt_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C][CH]C(18176)', 'CH2CHCCH(26391)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00294841,'m^3/(mol*s)'), n=2.48333, Ea=(17.6885,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CtH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H(33)', '[CH2]C(C=C)=CC(24210)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00684716,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-HH;CJ] for rate rule [Cds-CdH_Cds-HH;CtJ_Ct]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]C=C(4699)', 'CH3CHCCH2(18175)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#C[C](C)C([CH2])=CC(28190)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#C[C]([CH2])C(C)=CC(28191)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.4947e+07,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#CC([CH2])C(C)=[C]C(28192)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#CC(C)C([CH2])=[C]C(28193)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#CC([CH2])[C](C)C=C(27304)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[C]#CC(C)C([CH2])=CC(28194)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.39293e+07,'s^-1'), n=1.32074, Ea=(96.0416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_TSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[C]#CC([CH2])C(C)=CC(28195)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(263079,'s^-1'), n=1.73643, Ea=(39.8993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_TSSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#CC(C)[C]([CH2])C=C(27307)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(121000,'s^-1'), n=1.9, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 92 used for R5H_SSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]C=C(4699)', 'C=[C][CH]C(18176)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#CC([CH2])[C]1CC1C(28196)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['[CH]=[C]C1CC(C)C1=C(28127)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(149.076,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['[CH2]C(=CC)C1[C]=CC1(28197)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(130.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['[CH2]C1[C]=CCC1=CC(27951)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra_H;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#CC(=C)C(C)=CC(28198)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(S)(23)', 'C#CC([CH2])C([CH2])=C(28199)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#CC[CH]C([CH2])=CC(28163)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#CC1CCC1=CC(28109)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(19)', 'C#C[CH]C([CH2])=CC(28200)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/TwoDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(19)', 'C#CC([CH2])[C]=CC(28201)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['[CH]=C1C([CH2])C(=C)C1C(28062)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.13771e+06,'s^-1'), n=1.58803, Ea=(102.944,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 100.0 to 102.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', 'C#CC([CH2])C(=C)C=C(27302)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#CC([CH2])C(=C)C[CH2](28202)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#C[C]([CH2])C(=C)CC(28203)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.18083e+09,'s^-1'), n=1.04667, Ea=(209.2,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C(CC)C([CH2])C#C(28204)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C([CH]C)C(C)C#C(26415)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[C]#CC([CH2])C(=C)CC(28205)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(366176,'s^-1'), n=1.54456, Ea=(43.4053,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_TSSS;Ct_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['[CH]=[C]C1CCC1=CC(28110)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(145.101,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['[CH2]C1[C]=CC(C)C1=C(28074)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHCs] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#CC(=C)C(=C)CC(28206)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#CC(C)C(=C)C=C(27316)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#CC([CH2])C(C)[C]=C(26581)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#CC1CC(C)C1=C(28126)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['CHCH3(T)(95)', 'C#CC([CH2])[C]=C(27503)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4549',
    isomers = [
        'C#CC([CH2])C(=C)[CH]C(26577)',
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
    label = '4549',
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

