species(
    label = '[CH]=C(C=C)C(=C)[CH]C(26579)',
    structure = SMILES('[CH]=C(C=C)C([CH2])=CC'),
    E0 = (489.773,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,686.945],'cm^-1')),
        HinderedRotor(inertia=(0.754518,'amu*angstrom^2'), symmetry=1, barrier=(19.4107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.0999,'amu*angstrom^2'), symmetry=1, barrier=(95.2803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843649,'amu*angstrom^2'), symmetry=1, barrier=(19.3972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843589,'amu*angstrom^2'), symmetry=1, barrier=(19.3961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.205855,0.0719539,-3.72056e-05,-2.53472e-09,5.99275e-12,59052.6,27.7424], Tmin=(100,'K'), Tmax=(1035.62,'K')), NASAPolynomial(coeffs=[16.0731,0.032823,-1.26176e-05,2.29419e-09,-1.59574e-13,54578,-55.0994], Tmin=(1035.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Allyl_P)"""),
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
    label = '[CH2]C(=CC)C1=CC1[CH2](28132)',
    structure = SMILES('[CH2]C(=CC)C1=CC1[CH2]'),
    E0 = (574.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.268869,0.0849679,-6.80368e-05,2.25765e-08,-7.30641e-13,69228.4,26.8681], Tmin=(100,'K'), Tmax=(982.02,'K')), NASAPolynomial(coeffs=[17.7461,0.0293527,-1.02205e-05,1.74708e-09,-1.17402e-13,64833.6,-64.0781], Tmin=(982.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(574.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C1C(=CC)CC1[CH2](27985)',
    structure = SMILES('[CH]=C1C(=CC)CC1[CH2]'),
    E0 = (570.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.643732,0.0601405,-5.19823e-06,-3.60039e-08,1.83076e-11,68755.5,28.0056], Tmin=(100,'K'), Tmax=(973.561,'K')), NASAPolynomial(coeffs=[15.2841,0.03205,-1.1316e-05,2.01158e-09,-1.40562e-13,64385.4,-50.0367], Tmin=(973.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(570.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(Isobutyl) + radical(Cds_P)"""),
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
    label = 'C#CC([CH2])=CC(24773)',
    structure = SMILES('C#CC([CH2])=CC'),
    E0 = (345.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,350,440,435,1725,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.26681,'amu*angstrom^2'), symmetry=1, barrier=(18.6855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812654,'amu*angstrom^2'), symmetry=1, barrier=(18.6845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00731,'amu*angstrom^2'), symmetry=1, barrier=(46.152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45204,0.0489282,-2.77372e-05,2.38041e-09,2.3133e-12,41678.9,19.6959], Tmin=(100,'K'), Tmax=(1060.46,'K')), NASAPolynomial(coeffs=[11.5674,0.0232035,-8.93177e-06,1.61098e-09,-1.10941e-13,38834.6,-32.9939], Tmin=(1060.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(C=C)C(C)=[C]C(28133)',
    structure = SMILES('[CH]=C(C=C)C(C)=[C]C'),
    E0 = (576.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,214.611],'cm^-1')),
        HinderedRotor(inertia=(0.374619,'amu*angstrom^2'), symmetry=1, barrier=(12.2426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.759371,'amu*angstrom^2'), symmetry=1, barrier=(24.8117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374597,'amu*angstrom^2'), symmetry=1, barrier=(12.2422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374655,'amu*angstrom^2'), symmetry=1, barrier=(12.2424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.463689,0.0756547,-5.87549e-05,2.47638e-08,-4.3739e-12,69419.5,27.0076], Tmin=(100,'K'), Tmax=(1309.48,'K')), NASAPolynomial(coeffs=[12.7437,0.0381437,-1.57861e-05,2.88798e-09,-1.97453e-13,66203.5,-35.5471], Tmin=(1309.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(576.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(=CC)C(=C)[C]=C(28134)',
    structure = SMILES('[CH2]C(=C=C)C([CH2])=CC'),
    E0 = (426.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.318671,0.0834232,-5.90699e-05,1.27605e-08,2.28105e-12,51436.6,27.9246], Tmin=(100,'K'), Tmax=(1020.99,'K')), NASAPolynomial(coeffs=[19.0611,0.0290514,-1.08549e-05,1.95446e-09,-1.35814e-13,46355.9,-71.476], Tmin=(1020.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=[C]C)C(=C)C=C(28135)',
    structure = SMILES('[CH2]C(=[C]C)C(=C)C=C'),
    E0 = (480.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.187578,0.0736657,-4.72414e-05,1.13237e-08,3.36571e-13,57938.8,27.8823], Tmin=(100,'K'), Tmax=(1124.87,'K')), NASAPolynomial(coeffs=[15.2868,0.0341102,-1.33458e-05,2.40745e-09,-1.64816e-13,53647.5,-50.7151], Tmin=(1124.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(480.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C([C]=C)C(C)=CC(28136)',
    structure = SMILES('[CH]C(=C=C)C(C)=CC'),
    E0 = (493.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,540,610,2055,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.06619,'amu*angstrom^2'), symmetry=1, barrier=(47.5058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06635,'amu*angstrom^2'), symmetry=1, barrier=(47.5096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06638,'amu*angstrom^2'), symmetry=1, barrier=(47.5101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06654,'amu*angstrom^2'), symmetry=1, barrier=(47.5139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.152776,0.0863719,-6.87658e-05,2.94112e-08,-5.2065e-12,59565.2,27.877], Tmin=(100,'K'), Tmax=(1322.09,'K')), NASAPolynomial(coeffs=[15.1497,0.0400742,-1.62382e-05,2.9241e-09,-1.97969e-13,55519,-50.2215], Tmin=(1322.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C=C)[C](C)C=C(27281)',
    structure = SMILES('[CH]C(C=C)=C(C)C=C'),
    E0 = (441.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.306227,0.0814076,-4.14132e-05,-7.97032e-09,1.00572e-11,53230.5,28.3964], Tmin=(100,'K'), Tmax=(983.524,'K')), NASAPolynomial(coeffs=[18.7377,0.0330267,-1.19629e-05,2.12025e-09,-1.46785e-13,48078.5,-70.3111], Tmin=(983.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(441.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC(=C)C([CH2])=CC(28137)',
    structure = SMILES('[CH]=CC(=C)C([CH2])=CC'),
    E0 = (489.773,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,686.945],'cm^-1')),
        HinderedRotor(inertia=(0.754518,'amu*angstrom^2'), symmetry=1, barrier=(19.4107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.0999,'amu*angstrom^2'), symmetry=1, barrier=(95.2803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843649,'amu*angstrom^2'), symmetry=1, barrier=(19.3972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843589,'amu*angstrom^2'), symmetry=1, barrier=(19.3961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.205855,0.0719539,-3.72056e-05,-2.53472e-09,5.99275e-12,59052.6,27.7424], Tmin=(100,'K'), Tmax=(1035.62,'K')), NASAPolynomial(coeffs=[16.0731,0.032823,-1.26176e-05,2.29419e-09,-1.59574e-13,54578,-55.0994], Tmin=(1035.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C](C=C)C(=C)C=C(27282)',
    structure = SMILES('[CH2]C(C=C)=C([CH2])C=C'),
    E0 = (306.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.284361,0.0633081,6.07385e-06,-6.12019e-08,3.07278e-11,37028.9,27.5912], Tmin=(100,'K'), Tmax=(936.095,'K')), NASAPolynomial(coeffs=[20.0992,0.0250786,-7.08363e-06,1.16637e-09,-8.26656e-14,31284.5,-77.5635], Tmin=(936.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC(=[CH])C(C)=CC(28138)',
    structure = SMILES('[CH]=CC(=[CH])C(C)=CC'),
    E0 = (585.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3115,3125,620,680,785,800,1600,1700,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.719689,'amu*angstrom^2'), symmetry=1, barrier=(16.5471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.719953,'amu*angstrom^2'), symmetry=1, barrier=(16.5531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.720438,'amu*angstrom^2'), symmetry=1, barrier=(16.5643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.719995,'amu*angstrom^2'), symmetry=1, barrier=(16.5541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.174916,0.0774175,-6.0165e-05,2.47136e-08,-4.15601e-12,70546.8,27.9787], Tmin=(100,'K'), Tmax=(1396.8,'K')), NASAPolynomial(coeffs=[15.3807,0.0338727,-1.34029e-05,2.39491e-09,-1.61392e-13,66299,-50.462], Tmin=(1396.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(585.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)[C]1CC1C(28139)',
    structure = SMILES('[CH]=C(C=C)[C]1CC1C'),
    E0 = (510.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.679863,0.050856,3.83419e-05,-9.31813e-08,4.15084e-11,61573,25.5988], Tmin=(100,'K'), Tmax=(943.223,'K')), NASAPolynomial(coeffs=[20.3096,0.0233242,-6.47569e-06,1.11863e-09,-8.3788e-14,55391.6,-81.0938], Tmin=(943.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_P) + radical(Allyl_T)"""),
)

species(
    label = '[CH2]C=C1[CH]C(C)C1=C(28093)',
    structure = SMILES('[CH2]C=C1[CH]C(C)C1=C'),
    E0 = (377.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19133,0.0406398,5.67745e-05,-1.01144e-07,4.1188e-11,45528.8,24.3069], Tmin=(100,'K'), Tmax=(964.532,'K')), NASAPolynomial(coeffs=[15.9418,0.0318532,-1.10276e-05,2.0281e-09,-1.48125e-13,40246.6,-58.9547], Tmin=(964.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(C=CC=CCJ) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C(=CC)C1[CH]CC=1(28140)',
    structure = SMILES('[CH2]C(=CC)C1[CH]CC=1'),
    E0 = (437.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.628647,0.0531241,3.02417e-05,-8.11036e-08,3.5803e-11,52711.4,23.7185], Tmin=(100,'K'), Tmax=(962.949,'K')), NASAPolynomial(coeffs=[19.4041,0.026709,-8.95222e-06,1.65277e-09,-1.22564e-13,46704.2,-78.5695], Tmin=(962.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Allyl_P) + radical(cyclobutene-allyl)"""),
)

species(
    label = '[CH]=C1[CH]CCC1=CC(28097)',
    structure = SMILES('[CH]C1=CCCC1=CC'),
    E0 = (394.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.811404,0.0529502,2.62969e-05,-6.56743e-08,2.70874e-11,47520.2,26.1971], Tmin=(100,'K'), Tmax=(1000.14,'K')), NASAPolynomial(coeffs=[14.1263,0.0403174,-1.56765e-05,2.91134e-09,-2.0695e-13,42825.3,-48.1979], Tmin=(1000.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(3-Methylenecyclopentene) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]=C(C=C)C([CH2])=C(28141)',
    structure = SMILES('[CH]=C(C=C)C([CH2])=C'),
    E0 = (525.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,312.221],'cm^-1')),
        HinderedRotor(inertia=(0.41955,'amu*angstrom^2'), symmetry=1, barrier=(29.0219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.41958,'amu*angstrom^2'), symmetry=1, barrier=(29.0219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.419543,'amu*angstrom^2'), symmetry=1, barrier=(29.022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916759,0.0556429,-1.41995e-05,-2.39087e-08,1.38846e-11,63360.8,22.6903], Tmin=(100,'K'), Tmax=(978.302,'K')), NASAPolynomial(coeffs=[15.5424,0.0238953,-8.53406e-06,1.5414e-09,-1.09263e-13,59156.7,-54.4104], Tmin=(978.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(C=C)C[C]=CC(26587)',
    structure = SMILES('[CH]=C(C=C)C[C]=CC'),
    E0 = (587.869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.651136,'amu*angstrom^2'), symmetry=1, barrier=(14.9709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.65131,'amu*angstrom^2'), symmetry=1, barrier=(14.9749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.651177,'amu*angstrom^2'), symmetry=1, barrier=(14.9718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.651757,'amu*angstrom^2'), symmetry=1, barrier=(14.9852,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3604.07,'J/mol'), sigma=(6.17057,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.95 K, Pc=34.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0860644,0.0814773,-6.82616e-05,2.98871e-08,-5.27052e-12,70858.3,31.6078], Tmin=(100,'K'), Tmax=(1354.66,'K')), NASAPolynomial(coeffs=[17.188,0.0304707,-1.17821e-05,2.09168e-09,-1.40887e-13,66178.3,-56.9728], Tmin=(1354.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(587.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C=CC1=CCC1=CC(27964)',
    structure = SMILES('C=CC1=CCC1=CC'),
    E0 = (254.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12877,0.0429577,4.71271e-05,-9.06887e-08,3.73749e-11,30739.3,23.662], Tmin=(100,'K'), Tmax=(969.778,'K')), NASAPolynomial(coeffs=[16.1444,0.030729,-1.08402e-05,2.01233e-09,-1.47304e-13,25489.6,-60.3698], Tmin=(969.778,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
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
    label = '[CH]=C([C]=CC)C=C(28142)',
    structure = SMILES('[CH]C(=C=CC)C=C'),
    E0 = (531.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.05591,'amu*angstrom^2'), symmetry=1, barrier=(47.2693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05616,'amu*angstrom^2'), symmetry=1, barrier=(47.2752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05584,'amu*angstrom^2'), symmetry=1, barrier=(47.2679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.434748,0.0693313,-4.8366e-05,1.53366e-08,-1.3155e-12,64068.8,25.1719], Tmin=(100,'K'), Tmax=(1151.6,'K')), NASAPolynomial(coeffs=[14.3858,0.0312023,-1.21551e-05,2.16223e-09,-1.46243e-13,60170.7,-47.0766], Tmin=(1151.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(531.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C1C(=C)C(C)C1[CH2](28143)',
    structure = SMILES('[CH]=C1C(=C)C(C)C1[CH2]'),
    E0 = (574.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.647846,0.0566519,1.27049e-05,-6.0408e-08,2.85086e-11,69272.6,26.7978], Tmin=(100,'K'), Tmax=(949.733,'K')), NASAPolynomial(coeffs=[17.4214,0.0281843,-8.9485e-06,1.5521e-09,-1.10099e-13,64184.3,-63.2742], Tmin=(949.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(574.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(Cds_P) + radical(Isobutyl)"""),
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
    label = '[CH]=C(C=C)C(=C)C=C(27280)',
    structure = SMILES('[CH]=C(C=C)C(=C)C=C'),
    E0 = (482.762,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,2980,3010,3040,3070,3100,1330,1380,1430,900,975,1050,1000,1025,1050,1600,1650,1700,180,609.922],'cm^-1')),
        HinderedRotor(inertia=(0.292572,'amu*angstrom^2'), symmetry=1, barrier=(80.19,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06398,'amu*angstrom^2'), symmetry=1, barrier=(24.4629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26372,'amu*angstrom^2'), symmetry=1, barrier=(24.4517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.879132,0.0561998,-6.38828e-06,-2.83999e-08,1.39286e-11,58186.1,25.6759], Tmin=(100,'K'), Tmax=(1016.31,'K')), NASAPolynomial(coeffs=[13.7013,0.0331465,-1.28214e-05,2.35909e-09,-1.6596e-13,54164.2,-43.3557], Tmin=(1016.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)C(=C)C[CH2](28144)',
    structure = SMILES('[CH]=C(C=C)C(=C)C[CH2]'),
    E0 = (556.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,238.13,1209.93],'cm^-1')),
        HinderedRotor(inertia=(0.0157985,'amu*angstrom^2'), symmetry=1, barrier=(16.4521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.717379,'amu*angstrom^2'), symmetry=1, barrier=(16.494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29498,'amu*angstrom^2'), symmetry=1, barrier=(29.774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.71707,'amu*angstrom^2'), symmetry=1, barrier=(16.4868,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.241776,0.0714697,-3.80533e-05,-9.33633e-10,5.31529e-12,67124.6,30.2259], Tmin=(100,'K'), Tmax=(1040.32,'K')), NASAPolynomial(coeffs=[15.8607,0.0325911,-1.25287e-05,2.27599e-09,-1.58096e-13,62728.9,-51.2514], Tmin=(1040.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C(C=C)C(=[CH])CC(28145)',
    structure = SMILES('[CH]=C(C=C)C(=[CH])CC'),
    E0 = (598.75,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.514949,'amu*angstrom^2'), symmetry=1, barrier=(11.8397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.060968,'amu*angstrom^2'), symmetry=1, barrier=(20.0124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.871974,'amu*angstrom^2'), symmetry=1, barrier=(20.0484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.870753,'amu*angstrom^2'), symmetry=1, barrier=(20.0203,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.123022,0.0743523,-4.45945e-05,4.64185e-09,3.64158e-12,72162,29.2686], Tmin=(100,'K'), Tmax=(1042.92,'K')), NASAPolynomial(coeffs=[16.3492,0.0319816,-1.22223e-05,2.21024e-09,-1.53044e-13,67697.2,-54.8744], Tmin=(1042.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(598.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]C)C(=C)C=C(28146)',
    structure = SMILES('[CH]C(=CC)C(=C)C=C'),
    E0 = (461.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.252046,0.0711351,-2.94892e-05,-6.49563e-09,6.06275e-12,55693.7,28.0042], Tmin=(100,'K'), Tmax=(1078.81,'K')), NASAPolynomial(coeffs=[13.8386,0.041221,-1.6347e-05,2.96461e-09,-2.03856e-13,51571.5,-44.0925], Tmin=(1078.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([C]=C)C(=C)CC(28147)',
    structure = SMILES('[CH]C(=C=C)C(=C)CC'),
    E0 = (507.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.378778,0.0852776,-5.96901e-05,1.71783e-08,-4.98406e-13,61188,29.7968], Tmin=(100,'K'), Tmax=(1089.83,'K')), NASAPolynomial(coeffs=[17.2171,0.0364176,-1.40807e-05,2.51566e-09,-1.71489e-13,56419.1,-60.8905], Tmin=(1089.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC(=[CH])C(=C)CC(28148)',
    structure = SMILES('[CH]=CC(=[CH])C(=C)CC'),
    E0 = (598.75,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.514949,'amu*angstrom^2'), symmetry=1, barrier=(11.8397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.060968,'amu*angstrom^2'), symmetry=1, barrier=(20.0124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.871974,'amu*angstrom^2'), symmetry=1, barrier=(20.0484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.870753,'amu*angstrom^2'), symmetry=1, barrier=(20.0203,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.123022,0.0743523,-4.45945e-05,4.64185e-09,3.64158e-12,72162,29.2686], Tmin=(100,'K'), Tmax=(1042.92,'K')), NASAPolynomial(coeffs=[16.3492,0.0319816,-1.22223e-05,2.21024e-09,-1.53044e-13,67697.2,-54.8744], Tmin=(1042.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(598.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C1[CH]CC1=CC(27969)',
    structure = SMILES('[CH2]C=C1[CH]CC1=CC'),
    E0 = (373.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19287,0.0440504,3.92023e-05,-7.72648e-08,3.1259e-11,45011.5,24.8019], Tmin=(100,'K'), Tmax=(985.528,'K')), NASAPolynomial(coeffs=[13.8269,0.0356839,-1.33762e-05,2.48332e-09,-1.7825e-13,40437.3,-46.5385], Tmin=(985.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(12methylenecyclobutane) + radical(Allyl_S) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=C1[CH]CC(C)C1=C(28113)',
    structure = SMILES('[CH]C1=CCC(C)C1=C'),
    E0 = (398.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.814635,0.0494908,4.39989e-05,-8.96499e-08,3.70207e-11,48037.3,24.9914], Tmin=(100,'K'), Tmax=(974.327,'K')), NASAPolynomial(coeffs=[16.1792,0.0365891,-1.33857e-05,2.46954e-09,-1.77925e-13,42661.7,-60.9562], Tmin=(974.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(3-Methylenecyclopentene) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC(=C)C(=C)C=C(27294)',
    structure = SMILES('C=CC(=C)C(=C)C=C'),
    E0 = (235.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.959099,0.0526053,1.05118e-05,-4.60808e-08,1.99632e-11,28466.1,24.2925], Tmin=(100,'K'), Tmax=(1006.91,'K')), NASAPolynomial(coeffs=[13.4171,0.0360526,-1.38972e-05,2.56762e-09,-1.81632e-13,24287.6,-44.1873], Tmin=(1006.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH]=C(C=C)C(C)[C]=C(26583)',
    structure = SMILES('[CH]=C(C=C)C(C)[C]=C'),
    E0 = (590.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.709347,'amu*angstrom^2'), symmetry=1, barrier=(16.3093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709002,'amu*angstrom^2'), symmetry=1, barrier=(16.3013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709488,'amu*angstrom^2'), symmetry=1, barrier=(16.3125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.708522,'amu*angstrom^2'), symmetry=1, barrier=(16.2903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3545.26,'J/mol'), sigma=(6.132,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.76 K, Pc=34.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0558102,0.0865993,-8.11176e-05,4.08577e-08,-8.3685e-12,71147.9,29.2759], Tmin=(100,'K'), Tmax=(1168.98,'K')), NASAPolynomial(coeffs=[15.2455,0.0342415,-1.39336e-05,2.54264e-09,-1.74374e-13,67570.5,-46.9329], Tmin=(1168.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1=CC(C)C1=C(28088)',
    structure = SMILES('C=CC1=CC(C)C1=C'),
    E0 = (257.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865079,0.0514728,2.27911e-05,-6.54477e-08,2.84718e-11,31041.7,21.0026], Tmin=(100,'K'), Tmax=(975.992,'K')), NASAPolynomial(coeffs=[15.9895,0.0315931,-1.13688e-05,2.08898e-09,-1.503e-13,26084,-61.8704], Tmin=(975.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
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
    label = '[CH]=C([C]=C)C=C(27450)',
    structure = SMILES('[CH]C(=C=C)C=C'),
    E0 = (567.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.05661,'amu*angstrom^2'), symmetry=1, barrier=(47.2856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04718,'amu*angstrom^2'), symmetry=1, barrier=(47.0687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21707,0.052238,-2.29274e-05,-8.71584e-09,7.49875e-12,68373.9,20.5529], Tmin=(100,'K'), Tmax=(990.155,'K')), NASAPolynomial(coeffs=[13.1856,0.0233589,-8.67502e-06,1.5483e-09,-1.07222e-13,65049.3,-41.8898], Tmin=(990.155,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
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
    E0 = (489.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (634.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (576.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (570.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (629.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (662.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (652.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (738.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (712.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (893.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (557.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (639.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (893.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (522.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (618.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (812.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (583.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (627.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (627.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (604.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (945.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (757.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (498.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (913.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (574.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (694.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (670.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (743.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (893.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (551.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (636.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (627.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (547.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (514.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (684.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (497.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (911.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['CH3CHCCH2(18175)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH2]C1([CH]C)C=C1C=C(28131)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(144.391,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_D;doublebond_intra_HNd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH2]C(=CC)C1=CC1[CH2](28132)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(86.3586,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 14 used for R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH
Exact match found for rate rule [R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH]=C1C(=CC)CC1[CH2](27985)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.07e+06,'s^-1'), n=1.46, Ea=(80.7863,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 343 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 79.2 to 80.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=[C]C=C(4699)', 'CH3CHCCH2(18175)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C][CH]C(18176)', 'CH2CHCCH(26391)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0345031,'m^3/(mol*s)'), n=2.33733, Ea=(26.7584,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-Cd_Ct-H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H3(30)', 'C#CC([CH2])=CC(24773)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(25300,'cm^3/(mol*s)'), n=2.41, Ea=(20.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2412 used for Ct-Cd_Ct-H;CdsJ-H
Exact match found for rate rule [Ct-Cd_Ct-H;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(C=C)C(C)=[C]C(28133)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH2]C(=CC)C(=C)[C]=C(28134)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH2]C(=[C]C)C(=C)C=C(28135)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_single] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C([C]=C)C(C)=CC(28136)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(60051,'s^-1'), n=2.135, Ea=(63.3667,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SS(Cd)S;Y_rad_out;Cs_H_out_2H] + [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH]=C(C=C)[C](C)C=C(27281)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC(=C)C([CH2])=CC(28137)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.27529e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_singleH] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH2][C](C=C)C(=C)C=C(27282)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 286 used for R5H_DSMS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R5H_DSMS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC(=[CH])C(C)=CC(28138)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C]C=C(4699)', 'C=[C][CH]C(18176)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH]=C(C=C)[C]1CC1C(28139)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.68393e+12,'s^-1'), n=-0.105173, Ea=(93.5715,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe_HNd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH2]C=C1[CH]C(C)C1=C(28093)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH2]C(=CC)C1[CH]CC=1(28140)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R4_D_D;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH]=C1[CH]CCC1=CC(28097)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.12e+09,'s^-1'), n=0.63, Ea=(114.642,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(S)(23)', '[CH]=C(C=C)C([CH2])=C(28141)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(C=C)C[C]=CC(26587)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['C=CC1=CCC1=CC(27964)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(19)', '[CH]=C([C]=CC)C=C(28142)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH]=C1C(=C)C(C)C1[CH2](28143)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.73e+06,'s^-1'), n=1.31, Ea=(85.0598,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 344 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_csHNd
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 83.2 to 85.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH]=C(C=C)C(=C)C=C(27280)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C(C=C)C(=C)C[CH2](28144)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C(C=C)C(=[CH])CC(28145)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH]=C([CH]C)C(=C)C=C(28146)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.27529e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_singleH] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C([C]=C)C(=C)CC(28147)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=CC(=[CH])C(=C)CC(28148)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(315594,'s^-1'), n=1.73223, Ea=(38.2227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;Cd_rad_out_singleH;Cs_H_out] + [R5H_RSSR;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH2]C=C1[CH]CC1=CC(27969)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['[CH]=C1[CH]CC(C)C1=C(28113)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['C=CC(=C)C(=C)C=C(27294)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(C=C)C(C)[C]=C(26583)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    products = ['C=CC1=CC(C)C1=C(28088)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CHCH3(T)(95)', '[CH]=C([C]=C)C=C(27450)'],
    products = ['[CH]=C(C=C)C(=C)[CH]C(26579)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4551',
    isomers = [
        '[CH]=C(C=C)C(=C)[CH]C(26579)',
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
    label = '4551',
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

