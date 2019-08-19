species(
    label = 'C#C[CH]CC([CH2])(O)C(=O)CC(29693)',
    structure = SMILES('C#C[CH]CC([CH2])(O)C(=O)CC'),
    E0 = (92.9663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,180,981.09,1600,1816.5,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151096,'amu*angstrom^2'), symmetry=1, barrier=(3.47399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151096,'amu*angstrom^2'), symmetry=1, barrier=(3.47399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151096,'amu*angstrom^2'), symmetry=1, barrier=(3.47399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151096,'amu*angstrom^2'), symmetry=1, barrier=(3.47399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151096,'amu*angstrom^2'), symmetry=1, barrier=(3.47399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151096,'amu*angstrom^2'), symmetry=1, barrier=(3.47399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151096,'amu*angstrom^2'), symmetry=1, barrier=(3.47399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151096,'amu*angstrom^2'), symmetry=1, barrier=(3.47399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.31723,0.146596,-0.000187289,1.18181e-07,-2.50021e-11,11402,41.5806], Tmin=(100,'K'), Tmax=(707.552,'K')), NASAPolynomial(coeffs=[19.8269,0.0439193,-1.73374e-05,3.01477e-09,-1.97612e-13,7704.89,-61.5731], Tmin=(707.552,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.9663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)(C=O)O) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C=C(O)C(=O)CC(4626)',
    structure = SMILES('C=C(O)C(=O)CC'),
    E0 = (-358.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,332.629,332.633],'cm^-1')),
        HinderedRotor(inertia=(0.152178,'amu*angstrom^2'), symmetry=1, barrier=(11.9524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152152,'amu*angstrom^2'), symmetry=1, barrier=(11.9478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152157,'amu*angstrom^2'), symmetry=1, barrier=(11.9534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152209,'amu*angstrom^2'), symmetry=1, barrier=(11.954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4071.48,'J/mol'), sigma=(6.49965,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.96 K, Pc=33.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.785843,0.0702816,-6.54766e-05,3.23543e-08,-6.53762e-12,-42997.7,22.9463], Tmin=(100,'K'), Tmax=(1175.98,'K')), NASAPolynomial(coeffs=[12.9689,0.0288414,-1.26178e-05,2.38822e-09,-1.6711e-13,-45863,-37.8049], Tmin=(1175.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(O)=C([O])CC(4557)',
    structure = SMILES('[CH2]C(O)=C([O])CC'),
    E0 = (-184.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,402.631,402.631],'cm^-1')),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100034,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4303.36,'J/mol'), sigma=(7.05412,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=672.18 K, Pc=27.82 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.516209,0.0797578,-7.96062e-05,3.90546e-08,-7.19184e-12,-22064.1,30.4076], Tmin=(100,'K'), Tmax=(1528.49,'K')), NASAPolynomial(coeffs=[20.9225,0.0119894,-1.65417e-06,6.23996e-11,2.30461e-15,-27255.3,-77.6604], Tmin=(1528.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C#C[CH]CC1(O)CC1([O])CC(29917)',
    structure = SMILES('C#C[CH]CC1(O)CC1([O])CC'),
    E0 = (187.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43992,0.126511,-0.000127244,6.63368e-08,-1.35109e-11,22795.1,41.2678], Tmin=(100,'K'), Tmax=(1233.89,'K')), NASAPolynomial(coeffs=[26.0408,0.0320587,-9.83906e-06,1.5082e-09,-9.3176e-14,15928.4,-101.466], Tmin=(1233.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC1CC([CH2])(O)C1([O])CC(29918)',
    structure = SMILES('C#CC1CC([CH2])(O)C1([O])CC'),
    E0 = (239.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.04865,0.119994,-0.000102635,3.75099e-08,-2.68987e-12,29065.3,38.9046], Tmin=(100,'K'), Tmax=(991.798,'K')), NASAPolynomial(coeffs=[25.0627,0.0355061,-1.24442e-05,2.15211e-09,-1.46258e-13,22465.1,-97.8313], Tmin=(991.798,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH]=C1[CH]CC(O)(C1)C(=O)CC(29911)',
    structure = SMILES('[CH]C1=CCC(O)(C1)C(=O)CC'),
    E0 = (-13.2498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.691,0.112043,-8.72385e-05,3.52105e-08,-5.75949e-12,-1377.94,39.2745], Tmin=(100,'K'), Tmax=(1444.52,'K')), NASAPolynomial(coeffs=[22.4668,0.0451474,-1.77731e-05,3.15093e-09,-2.10977e-13,-8357.17,-86.1575], Tmin=(1444.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.2498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(AllylJ2_triplet)"""),
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
    label = 'C#CC=CC([CH2])(O)C(=O)CC(29919)',
    structure = SMILES('C#CC=CC([CH2])(O)C(=O)CC'),
    E0 = (57.2484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,205.038,1544.36,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149498,'amu*angstrom^2'), symmetry=1, barrier=(3.43724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149498,'amu*angstrom^2'), symmetry=1, barrier=(3.43724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149498,'amu*angstrom^2'), symmetry=1, barrier=(3.43724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149498,'amu*angstrom^2'), symmetry=1, barrier=(3.43724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149498,'amu*angstrom^2'), symmetry=1, barrier=(3.43724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149498,'amu*angstrom^2'), symmetry=1, barrier=(3.43724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149498,'amu*angstrom^2'), symmetry=1, barrier=(3.43724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.182,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.03423,0.140039,-0.000184961,1.32681e-07,-3.82745e-11,7096.11,39.6664], Tmin=(100,'K'), Tmax=(845.941,'K')), NASAPolynomial(coeffs=[17.4277,0.0480233,-2.18165e-05,4.12379e-09,-2.85759e-13,3803.05,-50.9721], Tmin=(845.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.2484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CC(O)(C=O)CJ)"""),
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
    label = 'C2H5CO(71)',
    structure = SMILES('CC[C]=O'),
    E0 = (-44.1874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.192008,'amu*angstrom^2'), symmetry=1, barrier=(4.41465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191563,'amu*angstrom^2'), symmetry=1, barrier=(4.40441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91313,0.0225012,-8.59328e-06,4.80319e-10,2.48743e-13,-5274.24,14.655], Tmin=(100,'K'), Tmax=(1673.18,'K')), NASAPolynomial(coeffs=[7.46306,0.0158512,-6.42141e-06,1.12499e-09,-7.32055e-14,-7388.54,-11.406], Tmin=(1673.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.1874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2CO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C#C[CH]CC(=C)O(27848)',
    structure = SMILES('C#C[CH]CC(=C)O'),
    E0 = (161.455,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,3615,1277.5,1000,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,347.607,349.16],'cm^-1')),
        HinderedRotor(inertia=(0.188638,'amu*angstrom^2'), symmetry=1, barrier=(16.6837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.85356,'amu*angstrom^2'), symmetry=1, barrier=(75.7529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194696,'amu*angstrom^2'), symmetry=1, barrier=(16.7365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267148,'amu*angstrom^2'), symmetry=1, barrier=(23.9723,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.162327,0.0745926,-7.55548e-05,3.83664e-08,-7.31337e-12,19582.9,27.7193], Tmin=(100,'K'), Tmax=(1489.94,'K')), NASAPolynomial(coeffs=[18.482,0.0128342,-1.59595e-06,1.40068e-12,8.63281e-15,15326.3,-65.3034], Tmin=(1489.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
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
    label = 'C#C[CH]CC(=C)C(=O)CC(29920)',
    structure = SMILES('C#C[CH]CC(=C)C(=O)CC'),
    E0 = (174.636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.472931,0.104688,-0.000118836,8.47842e-08,-2.58926e-11,21159.1,35.7195], Tmin=(100,'K'), Tmax=(783.856,'K')), NASAPolynomial(coeffs=[9.17872,0.0554319,-2.45697e-05,4.60376e-09,-3.17903e-13,19646.2,-8.49284], Tmin=(783.856,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC[CH]C([CH2])(O)C(=O)CC(29750)',
    structure = SMILES('C#CC[CH]C([CH2])(O)C(=O)CC'),
    E0 = (146.658,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,303.556,840.416,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153563,'amu*angstrom^2'), symmetry=1, barrier=(3.53071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153563,'amu*angstrom^2'), symmetry=1, barrier=(3.53071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153563,'amu*angstrom^2'), symmetry=1, barrier=(3.53071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153563,'amu*angstrom^2'), symmetry=1, barrier=(3.53071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153563,'amu*angstrom^2'), symmetry=1, barrier=(3.53071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153563,'amu*angstrom^2'), symmetry=1, barrier=(3.53071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153563,'amu*angstrom^2'), symmetry=1, barrier=(3.53071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153563,'amu*angstrom^2'), symmetry=1, barrier=(3.53071,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39932,0.144432,-0.000187287,1.29833e-07,-3.56911e-11,17866,44.4976], Tmin=(100,'K'), Tmax=(892.596,'K')), NASAPolynomial(coeffs=[20.1266,0.0434908,-1.76631e-05,3.14948e-09,-2.10716e-13,13844.5,-61.6184], Tmin=(892.596,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)(C=O)O) + radical(CCJCO)"""),
)

species(
    label = 'C#C[CH]CC(C)([O])C(=O)CC(29921)',
    structure = SMILES('C#C[CH]CC(C)([O])C(=O)CC'),
    E0 = (123.338,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,750,770,3400,2100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,180,811.56,1427.46,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156179,'amu*angstrom^2'), symmetry=1, barrier=(3.59086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156179,'amu*angstrom^2'), symmetry=1, barrier=(3.59086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156179,'amu*angstrom^2'), symmetry=1, barrier=(3.59086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156179,'amu*angstrom^2'), symmetry=1, barrier=(3.59086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156179,'amu*angstrom^2'), symmetry=1, barrier=(3.59086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156179,'amu*angstrom^2'), symmetry=1, barrier=(3.59086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156179,'amu*angstrom^2'), symmetry=1, barrier=(3.59086,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.87803,0.137794,-0.000192028,1.53127e-07,-4.86561e-11,15037.6,42.7602], Tmin=(100,'K'), Tmax=(871.289,'K')), NASAPolynomial(coeffs=[13.0547,0.054374,-2.28207e-05,4.07639e-09,-2.70146e-13,12999.7,-23.9863], Tmin=(871.289,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'C#C[CH][CH]C(C)(O)C(=O)CC(29922)',
    structure = SMILES('[CH]=C=C[CH]C(C)(O)C(=O)CC'),
    E0 = (4.79372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.11486,0.138317,-0.00016547,1.07262e-07,-2.79539e-11,793.48,39.967], Tmin=(100,'K'), Tmax=(934.185,'K')), NASAPolynomial(coeffs=[18.978,0.0480012,-2.04525e-05,3.77298e-09,-2.58996e-13,-3147.46,-60.358], Tmin=(934.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.79372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJCO)"""),
)

species(
    label = '[C]#CCCC([CH2])(O)C(=O)CC(29923)',
    structure = SMILES('[C]#CCCC([CH2])(O)C(=O)CC'),
    E0 = (283.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,375,552.5,462.5,1710,180,180,180,180,878.496,1349.81,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157407,'amu*angstrom^2'), symmetry=1, barrier=(3.6191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157407,'amu*angstrom^2'), symmetry=1, barrier=(3.6191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157407,'amu*angstrom^2'), symmetry=1, barrier=(3.6191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157407,'amu*angstrom^2'), symmetry=1, barrier=(3.6191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157407,'amu*angstrom^2'), symmetry=1, barrier=(3.6191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157407,'amu*angstrom^2'), symmetry=1, barrier=(3.6191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157407,'amu*angstrom^2'), symmetry=1, barrier=(3.6191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157407,'amu*angstrom^2'), symmetry=1, barrier=(3.6191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.68747,0.154662,-0.000224119,1.74224e-07,-5.30781e-11,34378.9,43.7708], Tmin=(100,'K'), Tmax=(889.146,'K')), NASAPolynomial(coeffs=[18.3294,0.046592,-1.89919e-05,3.31961e-09,-2.16276e-13,31175.9,-52.1484], Tmin=(889.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)(C=O)O) + radical(Acetyl)"""),
)

species(
    label = 'C#CCCC([CH2])([O])C(=O)CC(29924)',
    structure = SMILES('C#CCCC([CH2])([O])C(=O)CC'),
    E0 = (188.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,2175,525,750,770,3400,2100,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,180,180,180,180,991.019,1239.97,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160402,'amu*angstrom^2'), symmetry=1, barrier=(3.68796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160402,'amu*angstrom^2'), symmetry=1, barrier=(3.68796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160402,'amu*angstrom^2'), symmetry=1, barrier=(3.68796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160402,'amu*angstrom^2'), symmetry=1, barrier=(3.68796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160402,'amu*angstrom^2'), symmetry=1, barrier=(3.68796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160402,'amu*angstrom^2'), symmetry=1, barrier=(3.68796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160402,'amu*angstrom^2'), symmetry=1, barrier=(3.68796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.32275,0.148011,-0.000213437,1.70333e-07,-5.38887e-11,22919.2,44.1354], Tmin=(100,'K'), Tmax=(859.719,'K')), NASAPolynomial(coeffs=[15.4725,0.0518424,-2.23157e-05,4.03597e-09,-2.69393e-13,20353.6,-36.1524], Tmin=(859.719,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)(C=O)O) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'C#C[CH]CC(C)(O)C(=O)[CH]C(29925)',
    structure = SMILES('C#C[CH]CC(C)(O)C([O])=CC'),
    E0 = (17.8271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96825,0.125961,-0.000132204,7.05947e-08,-1.37422e-11,2363.88,43.5487], Tmin=(100,'K'), Tmax=(912.375,'K')), NASAPolynomial(coeffs=[21.8073,0.0374459,-1.25257e-05,2.03237e-09,-1.30321e-13,-2628.89,-72.5602], Tmin=(912.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.8271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C#CCCC([CH2])(O)C(=O)[CH]C(29926)',
    structure = SMILES('C#CCCC([CH2])(O)C([O])=CC'),
    E0 = (85.0603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.0707,0.130851,-0.000148362,8.97004e-08,-2.14895e-11,10451.3,44.9539], Tmin=(100,'K'), Tmax=(1022.12,'K')), NASAPolynomial(coeffs=[21.2018,0.039774,-1.47025e-05,2.52123e-09,-1.6618e-13,5693.92,-67.8316], Tmin=(1022.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.0603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'C#C[CH]CC(C)(O)C(=O)C[CH2](29927)',
    structure = SMILES('C#C[CH]CC(C)(O)C(=O)C[CH2]'),
    E0 = (92.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39323,0.145784,-0.000198742,1.47745e-07,-4.35366e-11,11403.9,42.8285], Tmin=(100,'K'), Tmax=(877.317,'K')), NASAPolynomial(coeffs=[18.1989,0.046113,-1.84378e-05,3.21788e-09,-2.10702e-13,8013.33,-52.5527], Tmin=(877.317,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJCC=O) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CCCC([CH2])(O)C(=O)C[CH2](29928)',
    structure = SMILES('C#CCCC([CH2])(O)C(=O)C[CH2]'),
    E0 = (158.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,552.5,462.5,1710,2175,525,750,770,3400,2100,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,289.477,857.364,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153369,'amu*angstrom^2'), symmetry=1, barrier=(3.52625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153369,'amu*angstrom^2'), symmetry=1, barrier=(3.52625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153369,'amu*angstrom^2'), symmetry=1, barrier=(3.52625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153369,'amu*angstrom^2'), symmetry=1, barrier=(3.52625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153369,'amu*angstrom^2'), symmetry=1, barrier=(3.52625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153369,'amu*angstrom^2'), symmetry=1, barrier=(3.52625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153369,'amu*angstrom^2'), symmetry=1, barrier=(3.52625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153369,'amu*angstrom^2'), symmetry=1, barrier=(3.52625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52125,0.151687,-0.000201836,1.35286e-07,-3.26977e-11,19272,43.1017], Tmin=(100,'K'), Tmax=(714.578,'K')), NASAPolynomial(coeffs=[20.239,0.0442821,-1.83622e-05,3.28336e-09,-2.19019e-13,15508.5,-62.6272], Tmin=(714.578,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJCC=O) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[C]#C[CH]CC(C)(O)C(=O)CC(29929)',
    structure = SMILES('[C]#C[CH]CC(C)(O)C(=O)CC'),
    E0 = (218.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,180,698.407,1535.23,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153648,'amu*angstrom^2'), symmetry=1, barrier=(3.53267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153648,'amu*angstrom^2'), symmetry=1, barrier=(3.53267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153648,'amu*angstrom^2'), symmetry=1, barrier=(3.53267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153648,'amu*angstrom^2'), symmetry=1, barrier=(3.53267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153648,'amu*angstrom^2'), symmetry=1, barrier=(3.53267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153648,'amu*angstrom^2'), symmetry=1, barrier=(3.53267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153648,'amu*angstrom^2'), symmetry=1, barrier=(3.53267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153648,'amu*angstrom^2'), symmetry=1, barrier=(3.53267,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.24284,0.14445,-0.000202747,1.57113e-07,-4.7915e-11,26497.3,42.3957], Tmin=(100,'K'), Tmax=(901.446,'K')), NASAPolynomial(coeffs=[15.8926,0.0491566,-1.95162e-05,3.36465e-09,-2.17415e-13,23829.9,-39.8754], Tmin=(901.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Acetyl)"""),
)

species(
    label = '[CH2]C(O)(CC1[C]=C1)C(=O)CC(29930)',
    structure = SMILES('[CH2]C(O)(CC1[C]=C1)C(=O)CC'),
    E0 = (266.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35268,0.148996,-0.000201204,1.43151e-07,-3.90396e-11,32269.8,40.0564], Tmin=(100,'K'), Tmax=(713.043,'K')), NASAPolynomial(coeffs=[18.0018,0.0493062,-2.19832e-05,4.09482e-09,-2.80391e-13,28998.6,-53.8419], Tmin=(713.043,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(cyclopropenyl-vinyl) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'C#C[CH]CC1(O)CO[C]1CC(29931)',
    structure = SMILES('C#C[CH]CC1(O)CO[C]1CC'),
    E0 = (178.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.9908,0.126872,-0.000144639,9.15487e-08,-2.26647e-11,21712.8,40.3125], Tmin=(100,'K'), Tmax=(1096.36,'K')), NASAPolynomial(coeffs=[18.7013,0.0404974,-1.15785e-05,1.58623e-09,-8.67424e-14,17829.5,-58.4363], Tmin=(1096.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane) + radical(Sec_Propargyl) + radical(C2CsJOCs)"""),
)

species(
    label = 'C#CC1CC([CH2])(O)[C](CC)O1(29932)',
    structure = SMILES('C#CC1CC([CH2])(O)[C](CC)O1'),
    E0 = (152.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.2697,0.123065,-0.000124137,6.70051e-08,-1.40781e-11,18629.6,38.7606], Tmin=(100,'K'), Tmax=(1272.34,'K')), NASAPolynomial(coeffs=[23.1386,0.0342893,-8.98824e-06,1.17489e-09,-6.33515e-14,12884.1,-87.1098], Tmin=(1272.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Tetrahydrofuran) + radical(CJC(C)2O) + radical(C2CsJOCs)"""),
)

species(
    label = 'CC[C]([O])C1(O)CC=C=CC1(29853)',
    structure = SMILES('CC[C]([O])C1(O)CC=C=CC1'),
    E0 = (141.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.68475,0.179909,-0.000264156,2.09669e-07,-6.64308e-11,17292.3,28.3715], Tmin=(100,'K'), Tmax=(804.032,'K')), NASAPolynomial(coeffs=[19.8059,0.0579119,-2.69837e-05,5.07697e-09,-3.47837e-13,13680.8,-78.8017], Tmin=(804.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(141.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = 'CH3CH2OH(54)',
    structure = SMILES('CCO'),
    E0 = (-248.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00248522,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168417,'amu*angstrom^2'), symmetry=1, barrier=(8.09857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0684,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3014.83,'J/mol'), sigma=(4.53,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.8587,-0.00374017,6.95554e-05,-8.86548e-08,3.51688e-11,-29996.1,4.80185], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.56244,0.0152042,-5.38968e-06,8.6225e-10,-5.12898e-14,-31525.6,-9.47302], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-248.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2OH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C#C[CH]CC([CH2])=C=O(29933)',
    structure = SMILES('C#C[CH]CC(=C)[C]=O'),
    E0 = (408.811,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,1855,455,950,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,628.246,629.819],'cm^-1')),
        HinderedRotor(inertia=(0.131521,'amu*angstrom^2'), symmetry=1, barrier=(13.0511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.56831,'amu*angstrom^2'), symmetry=1, barrier=(13.0666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.36201,'amu*angstrom^2'), symmetry=1, barrier=(77.2992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.03756,'amu*angstrom^2'), symmetry=1, barrier=(92.8315,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728768,0.0752909,-8.70382e-05,5.66776e-08,-1.5156e-11,49283.6,27.6613], Tmin=(100,'K'), Tmax=(902.929,'K')), NASAPolynomial(coeffs=[10.6835,0.0311908,-1.37765e-05,2.58558e-09,-1.79189e-13,47485.9,-19.3482], Tmin=(902.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=C(C)CJ=O)"""),
)

species(
    label = 'C#CC=CC(C)(O)C(=O)CC(29934)',
    structure = SMILES('C#CC=CC(C)(O)C(=O)CC'),
    E0 = (-155.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6317,0.128295,-0.000143693,8.79406e-08,-2.19275e-11,-18560.3,38.1677], Tmin=(100,'K'), Tmax=(967.871,'K')), NASAPolynomial(coeffs=[17.3537,0.0498313,-2.20888e-05,4.17933e-09,-2.91803e-13,-22235.4,-52.8056], Tmin=(967.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-155.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
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
    label = 'C#C[CH]CC([CH2])(O)C(C)=O(29935)',
    structure = SMILES('C#C[CH]CC([CH2])(O)C(C)=O'),
    E0 = (118.383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69259,0.12675,-0.000159232,1.00892e-07,-2.36605e-11,14442,37.977], Tmin=(100,'K'), Tmax=(807.094,'K')), NASAPolynomial(coeffs=[20.3805,0.0323281,-1.15764e-05,1.91503e-09,-1.22362e-13,10391.3,-66.8038], Tmin=(807.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)(C=O)O) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH]CC[C](O)C(=O)CC(29689)',
    structure = SMILES('C#C[CH]CCC(O)=C([O])CC'),
    E0 = (2.24988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.5816,0.129806,-0.000136622,7.40012e-08,-1.55322e-11,520.336,44.9195], Tmin=(100,'K'), Tmax=(1231.54,'K')), NASAPolynomial(coeffs=[26.7278,0.0295261,-8.28994e-06,1.17937e-09,-6.90462e-14,-6313.25,-101.02], Tmin=(1231.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.24988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH]C[C](O)CC(=O)CC(29936)',
    structure = SMILES('C#C[CH]C[C](O)CC(=O)CC'),
    E0 = (65.5214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07638,0.146593,-0.000222607,1.89698e-07,-6.31902e-11,8086.85,43.0398], Tmin=(100,'K'), Tmax=(874.65,'K')), NASAPolynomial(coeffs=[11.3372,0.0585108,-2.56929e-05,4.65652e-09,-3.09683e-13,6763.19,-14.0297], Tmin=(874.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.5214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C2CsJOH)"""),
)

species(
    label = 'C#CC([CH2])C([CH2])(O)C(=O)CC(29694)',
    structure = SMILES('C#CC([CH2])C([CH2])(O)C(=O)CC'),
    E0 = (142.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,750,770,3400,2100,180,180,180,407.214,1321.03,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.144703,'amu*angstrom^2'), symmetry=1, barrier=(3.32701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144703,'amu*angstrom^2'), symmetry=1, barrier=(3.32701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144703,'amu*angstrom^2'), symmetry=1, barrier=(3.32701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144703,'amu*angstrom^2'), symmetry=1, barrier=(3.32701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144703,'amu*angstrom^2'), symmetry=1, barrier=(3.32701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144703,'amu*angstrom^2'), symmetry=1, barrier=(3.32701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144703,'amu*angstrom^2'), symmetry=1, barrier=(3.32701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144703,'amu*angstrom^2'), symmetry=1, barrier=(3.32701,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4669.55,'J/mol'), sigma=(7.77953,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=729.37 K, Pc=22.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35725,0.145917,-0.000185555,1.17857e-07,-2.61417e-11,17388.1,43.3784], Tmin=(100,'K'), Tmax=(732.643,'K')), NASAPolynomial(coeffs=[20.4558,0.0422935,-1.62485e-05,2.78754e-09,-1.81433e-13,13483.6,-63.4178], Tmin=(732.643,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CtCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)(C=O)O) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC1CC(O)(C1)C(=O)CC(29697)',
    structure = SMILES('C#CC1CC(O)(C1)C(=O)CC'),
    E0 = (-115.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77856,0.112379,-9.7768e-05,4.45052e-08,-8.07364e-12,-13707.6,38.7506], Tmin=(100,'K'), Tmax=(1332.67,'K')), NASAPolynomial(coeffs=[23.2204,0.0373453,-1.33129e-05,2.25674e-09,-1.48127e-13,-20370.6,-89.0343], Tmin=(1332.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-115.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane)"""),
)

species(
    label = 'C#C[CH]CC([CH2])=C(CC)OO(29937)',
    structure = SMILES('C#C[CH]CC([CH2])=C(CC)OO'),
    E0 = (339.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,750,770,3400,2100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62636,0.123397,-0.00012916,7.40647e-08,-1.7222e-11,41062.2,42.7791], Tmin=(100,'K'), Tmax=(1039.29,'K')), NASAPolynomial(coeffs=[18.1916,0.0471212,-1.90702e-05,3.44572e-09,-2.34499e-13,36942.9,-53.5945], Tmin=(1039.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Allyl_P)"""),
)

species(
    label = 'C#C[CH]CO[C](CC)C(=C)O(29685)',
    structure = SMILES('[CH]=C=CCOC(CC)=C([CH2])O'),
    E0 = (63.0028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.02063,0.134104,-0.000138916,7.17967e-08,-1.42419e-11,7847.57,45.943], Tmin=(100,'K'), Tmax=(1322.6,'K')), NASAPolynomial(coeffs=[30.7091,0.0247715,-6.61546e-06,9.24044e-10,-5.43087e-14,-434.237,-123.793], Tmin=(1322.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.0028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C#C[CH]CC([CH2])(O)C(=C)OC(29938)',
    structure = SMILES('C#C[CH]CC([CH2])(O)C(=C)OC'),
    E0 = (148.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,657.238,1541.23,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152088,'amu*angstrom^2'), symmetry=1, barrier=(3.4968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152088,'amu*angstrom^2'), symmetry=1, barrier=(3.4968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152088,'amu*angstrom^2'), symmetry=1, barrier=(3.4968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152088,'amu*angstrom^2'), symmetry=1, barrier=(3.4968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152088,'amu*angstrom^2'), symmetry=1, barrier=(3.4968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152088,'amu*angstrom^2'), symmetry=1, barrier=(3.4968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152088,'amu*angstrom^2'), symmetry=1, barrier=(3.4968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152088,'amu*angstrom^2'), symmetry=1, barrier=(3.4968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.34956,0.132324,-0.000140303,7.29747e-08,-1.34511e-11,18056.7,44.1196], Tmin=(100,'K'), Tmax=(928.059,'K')), NASAPolynomial(coeffs=[24.885,0.0336261,-1.09793e-05,1.76983e-09,-1.13878e-13,12197,-89.5732], Tmin=(928.059,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'C#C[CH]CC([CH2])(O)C(O)=CC(29939)',
    structure = SMILES('C#C[CH]CC([CH2])(O)C(O)=CC'),
    E0 = (93.4657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42026,0.135425,-0.000148366,7.88772e-08,-1.45089e-11,11478,44.5477], Tmin=(100,'K'), Tmax=(892.973,'K')), NASAPolynomial(coeffs=[25.1115,0.0326713,-1.03207e-05,1.61662e-09,-1.01855e-13,5740.69,-89.7534], Tmin=(892.973,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(93.4657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CC(C)(O)CJ) + radical(Sec_Propargyl)"""),
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
    label = '[CH]=C=CCC(O)=C([O])CC(29862)',
    structure = SMILES('C#C[CH]CC(O)=C([O])CC'),
    E0 = (26.0301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95454,0.11526,-0.000123645,6.71349e-08,-1.39355e-11,3358.55,40.4228], Tmin=(100,'K'), Tmax=(1295.87,'K')), NASAPolynomial(coeffs=[25.2479,0.0220177,-4.97946e-06,5.63241e-10,-2.68802e-14,-2912.84,-94.859], Tmin=(1295.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.0301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(Sec_Propargyl)"""),
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
    label = '[CH2]C([CH2])(O)C(=O)CC(5988)',
    structure = SMILES('[CH2]C([CH2])(O)C(=O)CC'),
    E0 = (-64.1774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24882,0.119041,-0.00016428,1.17739e-07,-3.31244e-11,-7533.04,31.7431], Tmin=(100,'K'), Tmax=(875.483,'K')), NASAPolynomial(coeffs=[18.1268,0.0305211,-1.26247e-05,2.26293e-09,-1.51404e-13,-10925.9,-59.1577], Tmin=(875.483,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.1774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH]=[C]C1CC(O)(C1)C(=O)CC(29940)',
    structure = SMILES('[CH]=[C]C1CC(O)(C1)C(=O)CC'),
    E0 = (204.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47073,0.112212,-9.84e-05,4.50425e-08,-8.31496e-12,24761.3,40.7087], Tmin=(100,'K'), Tmax=(1294.57,'K')), NASAPolynomial(coeffs=[21.2925,0.0418776,-1.69039e-05,3.07412e-09,-2.10233e-13,18867.6,-74.9876], Tmin=(1294.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1(O)CC=C=CC1([O])CC(29941)',
    structure = SMILES('[CH2]C1(O)CC=C=CC1([O])CC'),
    E0 = (150.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.64751,0.178062,-0.000253016,1.92995e-07,-5.87755e-11,18345,29.2576], Tmin=(100,'K'), Tmax=(804.835,'K')), NASAPolynomial(coeffs=[20.8621,0.0562362,-2.59377e-05,4.87728e-09,-3.34997e-13,14400.2,-83.6625], Tmin=(804.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(C=CC(C)2OJ) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C(O)(C[C]=C=C)C(=O)CC(29942)',
    structure = SMILES('[CH2]C#CCC([CH2])(O)C(=O)CC'),
    E0 = (84.7623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.19968,0.141166,-0.00018171,1.27424e-07,-3.5732e-11,10413.4,43.5194], Tmin=(100,'K'), Tmax=(872.554,'K')), NASAPolynomial(coeffs=[18.378,0.0468369,-1.95578e-05,3.53931e-09,-2.38904e-13,6822.23,-52.9519], Tmin=(872.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.7623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CJC(C)(C=O)O) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]CC(C)(O)C(=O)CC(29943)',
    structure = SMILES('[CH]=C=[C]CC(C)(O)C(=O)CC'),
    E0 = (125.719,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,180,1600,1657.14,2857.31,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153175,'amu*angstrom^2'), symmetry=1, barrier=(3.52178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153175,'amu*angstrom^2'), symmetry=1, barrier=(3.52178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153175,'amu*angstrom^2'), symmetry=1, barrier=(3.52178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153175,'amu*angstrom^2'), symmetry=1, barrier=(3.52178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153175,'amu*angstrom^2'), symmetry=1, barrier=(3.52178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153175,'amu*angstrom^2'), symmetry=1, barrier=(3.52178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153175,'amu*angstrom^2'), symmetry=1, barrier=(3.52178,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09891,0.141077,-0.000190524,1.44115e-07,-4.39375e-11,15333.6,43.1766], Tmin=(100,'K'), Tmax=(812.922,'K')), NASAPolynomial(coeffs=[15.9492,0.0510153,-2.20252e-05,4.03109e-09,-2.72786e-13,12440.7,-39.902], Tmin=(812.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)(C=C[C]=C)C(=O)CC(29790)',
    structure = SMILES('[CH2]C(O)([CH]C=C=C)C(=O)CC'),
    E0 = (61.9284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.46296,0.147029,-0.00018258,1.19892e-07,-3.14303e-11,7676.85,40.2103], Tmin=(100,'K'), Tmax=(930.749,'K')), NASAPolynomial(coeffs=[20.9227,0.0465271,-2.06116e-05,3.87892e-09,-2.69408e-13,3323.59,-70.9338], Tmin=(930.749,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.9284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C([O])(CC=C=C)C(=O)CC(29944)',
    structure = SMILES('[CH2]C([O])(CC=C=C)C(=O)CC'),
    E0 = (186.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,1047.53,1180.56,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.166444,'amu*angstrom^2'), symmetry=1, barrier=(3.82689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166444,'amu*angstrom^2'), symmetry=1, barrier=(3.82689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166444,'amu*angstrom^2'), symmetry=1, barrier=(3.82689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166444,'amu*angstrom^2'), symmetry=1, barrier=(3.82689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166444,'amu*angstrom^2'), symmetry=1, barrier=(3.82689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166444,'amu*angstrom^2'), symmetry=1, barrier=(3.82689,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.11207,0.143155,-0.000198808,1.56233e-07,-4.97531e-11,22702.1,44.1557], Tmin=(100,'K'), Tmax=(793.612,'K')), NASAPolynomial(coeffs=[14.8365,0.054377,-2.46707e-05,4.62648e-09,-3.17495e-13,20117.6,-33.0281], Tmin=(793.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJC(C)(C=O)O) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH2]C(O)(CC=C=C)C(=O)[CH]C(29945)',
    structure = SMILES('[CH2]C(O)(CC=C=C)C([O])=CC'),
    E0 = (83.3159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.84635,0.125851,-0.000133214,7.46678e-08,-1.67195e-11,10233.5,44.9231], Tmin=(100,'K'), Tmax=(1084.98,'K')), NASAPolynomial(coeffs=[21.039,0.0414793,-1.65683e-05,2.99415e-09,-2.044e-13,5267.53,-67.3519], Tmin=(1084.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.3159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(C)(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)CC=C=C(29946)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)CC=C=C'),
    E0 = (156.601,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,508.97,638.925,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157617,'amu*angstrom^2'), symmetry=1, barrier=(3.62393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157617,'amu*angstrom^2'), symmetry=1, barrier=(3.62393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157617,'amu*angstrom^2'), symmetry=1, barrier=(3.62393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157617,'amu*angstrom^2'), symmetry=1, barrier=(3.62393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157617,'amu*angstrom^2'), symmetry=1, barrier=(3.62393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157617,'amu*angstrom^2'), symmetry=1, barrier=(3.62393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157617,'amu*angstrom^2'), symmetry=1, barrier=(3.62393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.50666,0.149592,-0.000199405,1.41765e-07,-4.01484e-11,19063.2,43.7994], Tmin=(100,'K'), Tmax=(864.664,'K')), NASAPolynomial(coeffs=[19.7849,0.0464655,-2.04964e-05,3.81843e-09,-2.62313e-13,15208.4,-60.502], Tmin=(864.664,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJC(C)(C=O)O) + radical(CJCC=O)"""),
)

species(
    label = '[CH2]C1(O)CC=C=CO[C]1CC(29947)',
    structure = SMILES('[CH2]C1(O)CC=C=CO[C]1CC'),
    E0 = (189.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52975,0.126395,-0.000119741,5.68629e-08,-1.05591e-11,22981.8,36.0012], Tmin=(100,'K'), Tmax=(1315.84,'K')), NASAPolynomial(coeffs=[28.7242,0.0313863,-1.14347e-05,1.98945e-09,-1.33514e-13,14756.8,-123.359], Tmin=(1315.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(C2CsJOC(O)) + radical(CJC(C)2O)"""),
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
    E0 = (92.9663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (187.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (239.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (162.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (271.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (108.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (211.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (92.9663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (217.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (271.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (262.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (210.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (431.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (290.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (179.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (122.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (145.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (203.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (293.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (266.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (318.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (279.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (249.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (152.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (489.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (171.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (538.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (250.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (250.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (302.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (101.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (482.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (206.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (407.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (236.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (407.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (466.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (204.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (235.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (225.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (170.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (137.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (287.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (238.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (333.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (223.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C=C(O)C(=O)CC(4626)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C#C[CH]CC1(O)CC1([O])CC(29917)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(94.5285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 90.6 to 94.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C#CC1CC([CH2])(O)C1([O])CC(29918)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(146.79,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_csHCt]
Euclidian distance = 3.0
family: Intra_R_Add_Exocyclic
Ea raised from 142.7 to 146.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['[CH]=C1[CH]CC(O)(C1)C(=O)CC(29911)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(69.2285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;triplebond_intra_H;radadd_intra_cs2H] for rate rule [R6;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#CC=CC([CH2])(O)C(=O)CC(29919)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(2.76144,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2560 used for Cds-CsH_Cds-CtH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CtH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]C=C(4699)', 'C=C(O)C(=O)CC(4626)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H5CO(71)', 'C#C[CH]CC(=C)O(27848)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cds-OsCs_Cds;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0132398,'m^3/(mol*s)'), n=2.333, Ea=(3.72143,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CtH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 3.0 to 3.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)', 'C#C[CH]CC(=C)C(=O)CC(29920)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.52832,'m^3/(mol*s)'), n=2.02802, Ea=(14.4047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;YJ] for rate rule [Cds-COCs_Cds;OJ_pri]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#CC[CH]C([CH2])(O)C(=O)CC(29750)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#C[CH]CC(C)([O])C(=O)CC(29921)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C#C[CH][CH]C(C)(O)C(=O)CC(29922)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]#CCCC([CH2])(O)C(=O)CC(29923)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#CCCC([CH2])([O])C(=O)CC(29924)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(363473,'s^-1'), n=1.92229, Ea=(102.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C#C[CH]CC(C)(O)C(=O)[CH]C(29925)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C#CCCC([CH2])(O)C(=O)[CH]C(29926)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.1016,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_1H;Cs_H_out_1H] for rate rule [R5H_CC(O2d)CC;C_rad_out_H/Ct;Cs_H_out_H/NonDeC]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C#C[CH]CC(C)(O)C(=O)C[CH2](29927)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#CCCC([CH2])(O)C(=O)C[CH2](29928)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(256.687,'s^-1'), n=2.005, Ea=(45.3964,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;C_rad_out_2H;Cs_H_out_H/OneDe] for rate rule [R6H_SSSSS;C_rad_out_2H;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[C]#C[CH]CC(C)(O)C(=O)CC(29929)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.31883e+06,'s^-1'), n=1.02765, Ea=(75.0925,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R6HJ_2;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['[CH2]C(O)(CC1[C]=C1)C(=O)CC(29930)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_T;triplebond_intra_H;radadd_intra_cs] for rate rule [R3_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C#C[CH]CC1(O)CO[C]1CC(29931)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C#CC1CC([CH2])(O)[C](CC)O1(29932)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.8435e+15,'s^-1'), n=-1.17677, Ea=(156.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHDe] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_csHCt]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['CC[C]([O])C1(O)CC=C=CC1(29853)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.33254e+09,'s^-1'), n=0.487896, Ea=(59.5573,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6_linear;multiplebond_intra;radadd_intra_cs2H] + [R6_linear;triplebond_intra_H;radadd_intra] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH3CH2OH(54)', 'C#C[CH]CC([CH2])=C=O(29933)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C#CC=CC(C)(O)C(=O)CC(29934)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(S)(23)', 'C#C[CH]CC([CH2])(O)C(C)=O(29935)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C#C[CH]CC[C](O)C(=O)CC(29689)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C#C[CH]C[C](O)CC(=O)CC(29936)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#CC([CH2])C([CH2])(O)C(=O)CC(29694)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['C#CC1CC(O)(C1)C(=O)CC(29697)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#C[CH]CC([CH2])=C(CC)OO(29937)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#C[CH]CO[C](CC)C(=C)O(29685)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#C[CH]CC([CH2])(O)C(=C)OC(29938)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C#C[CH]CC([CH2])(O)C(O)=CC(29939)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(19)', '[CH]=C=CCC(O)=C([O])CC(29862)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction37',
    reactants = ['C3H2(81)', '[CH2]C([CH2])(O)C(=O)CC(5988)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.26928e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction38',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['[CH]=[C]C1CC(O)(C1)C(=O)CC(29940)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.20551e+07,'s^-1'), n=1.225, Ea=(111.218,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 105.9 to 111.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction39',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['[CH2]C1(O)CC=C=CC1([O])CC(29941)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.73685e+10,'s^-1'), n=0.191667, Ea=(142.744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R7_MMSR;carbonylbond_intra_Nd;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['[CH2]C(O)(C[C]=C=C)C(=O)CC(29942)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=C=[C]CC(C)(O)C(=O)CC(29943)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['[CH2]C(O)(C=C[C]=C)C(=O)CC(29790)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/O)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C([O])(CC=C=C)C(=O)CC(29944)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(3027.78,'s^-1'), n=2.16706, Ea=(100.661,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;O_rad_out;Cd_H_out_singleH] + [R6H;O_rad_out;XH_out] for rate rule [R6H;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['[CH2]C(O)(CC=C=C)C(=O)[CH]C(29945)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R7H;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]CC(=O)C([CH2])(O)CC=C=C(29946)'],
    products = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(36369.3,'s^-1'), n=2.22538, Ea=(176.695,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_2H;Cd_H_out_singleH] + [R8H;C_rad_out_2H;XH_out] for rate rule [R8H;C_rad_out_2H;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C#C[CH]CC([CH2])(O)C(=O)CC(29693)'],
    products = ['[CH2]C1(O)CC=C=CO[C]1CC(29947)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.07e+10,'s^-1'), n=0.124, Ea=(130.708,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R7_linear;carbonyl_intra_Nd;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

network(
    label = '4712',
    isomers = [
        'C#C[CH]CC([CH2])(O)C(=O)CC(29693)',
    ],
    reactants = [
        ('C=C(O)C(=O)CC(4626)', 'CH2CHCCH(26391)'),
        ('[CH2]C(O)=C([O])CC(4557)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4712',
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

