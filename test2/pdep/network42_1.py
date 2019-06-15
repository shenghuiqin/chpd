species(
    label = 'C=C[C]1C[CH]CC1(219)',
    structure = SMILES('[CH2]C=C1C[CH]CC1'),
    E0 = (260.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3771.23,'J/mol'), sigma=(6.58809,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=589.06 K, Pc=29.93 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92025,0.0283869,6.30027e-05,-9.23789e-08,3.45631e-11,31370.7,25.2893], Tmin=(100,'K'), Tmax=(996.931,'K')), NASAPolynomial(coeffs=[10.7646,0.0362073,-1.39241e-05,2.63733e-09,-1.9107e-13,27455.2,-28.146], Tmin=(996.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(cyclopentane)"""),
)

species(
    label = '[CH2][CH]C12CCC1C2(251)',
    structure = SMILES('[CH2][CH]C12CCC1C2'),
    E0 = (481.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99619,0.0237657,7.91129e-05,-1.11606e-07,4.18238e-11,58057.2,25.3595], Tmin=(100,'K'), Tmax=(989.49,'K')), NASAPolynomial(coeffs=[12.2926,0.0332437,-1.27209e-05,2.45896e-09,-1.82068e-13,53517.9,-36.8468], Tmin=(989.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(481.965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'C=C[C]1[CH]CCC1(230)',
    structure = SMILES('[CH2]C=C1[CH]CCC1'),
    E0 = (215.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85503,0.0228302,9.52286e-05,-1.35835e-07,5.21811e-11,25995.3,21.8514], Tmin=(100,'K'), Tmax=(970.078,'K')), NASAPolynomial(coeffs=[15.0515,0.0297948,-1.04482e-05,2.01252e-09,-1.5253e-13,20546.9,-56.2986], Tmin=(970.078,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=C1C[CH]CC1(254)',
    structure = SMILES('C[C]=C1C[CH]CC1'),
    E0 = (346.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86334,0.0357018,2.92094e-05,-4.9661e-08,1.77844e-11,41751.6,25.6909], Tmin=(100,'K'), Tmax=(1068.69,'K')), NASAPolynomial(coeffs=[7.88536,0.040773,-1.66631e-05,3.13086e-09,-2.20715e-13,38887.8,-11.1389], Tmin=(1068.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C1CCCC1(255)',
    structure = SMILES('[CH2][C]=C1CCCC1'),
    E0 = (312.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,1685,370,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61778,0.0327225,6.12692e-05,-9.78266e-08,3.81492e-11,37633.3,23.5904], Tmin=(100,'K'), Tmax=(982.698,'K')), NASAPolynomial(coeffs=[13.8569,0.0317177,-1.17067e-05,2.22815e-09,-1.64244e-13,32870.9,-47.2348], Tmin=(982.698,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'CC=C1[CH]C[CH]C1(256)',
    structure = SMILES('CC=C1[CH]C[CH]C1'),
    E0 = (249.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11506,0.0256835,6.3377e-05,-8.75602e-08,3.1586e-11,30113,23.2041], Tmin=(100,'K'), Tmax=(1015.48,'K')), NASAPolynomial(coeffs=[8.74811,0.0393881,-1.57043e-05,2.98421e-09,-2.14612e-13,26712.1,-19.0106], Tmin=(1015.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(cyclopentane)"""),
)

species(
    label = 'CC=C1[CH][CH]CC1(257)',
    structure = SMILES('C[CH]C1=C[CH]CC1'),
    E0 = (199.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08757,0.0187961,9.97016e-05,-1.35474e-07,5.08064e-11,24058.8,18.7453], Tmin=(100,'K'), Tmax=(977.834,'K')), NASAPolynomial(coeffs=[13.0692,0.0327702,-1.20816e-05,2.33389e-09,-1.74742e-13,19095.4,-48.386], Tmin=(977.834,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(cyclopentene-allyl)"""),
)

species(
    label = 'CC=C1C[CH][CH]C1(258)',
    structure = SMILES('CC=C1C[CH][CH]C1'),
    E0 = (294.463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17002,0.0312599,3.15958e-05,-4.54891e-08,1.49374e-11,35489,25.992], Tmin=(100,'K'), Tmax=(1121.59,'K')), NASAPolynomial(coeffs=[5.15631,0.044689,-1.85674e-05,3.46904e-09,-2.41835e-13,33304.5,4.49047], Tmin=(1121.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(cyclopentane)"""),
)

species(
    label = '[CH]1CCC2([CH]C2)C1(259)',
    structure = SMILES('[CH]1CCC2([CH]C2)C1'),
    E0 = (392.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47404,0.0163282,8.50696e-05,-1.07271e-07,3.7907e-11,47225.3,24.0126], Tmin=(100,'K'), Tmax=(1009.56,'K')), NASAPolynomial(coeffs=[7.78671,0.0397718,-1.5871e-05,3.04036e-09,-2.20354e-13,43885.2,-12.8981], Tmin=(1009.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane) + radical(cyclopentane)"""),
)

species(
    label = '[CH2]C1[C]2CCC1C2(260)',
    structure = SMILES('[CH2]C1[C]2CCC1C2'),
    E0 = (488.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16731,0.0212614,8.16259e-05,-1.12746e-07,4.24311e-11,58861,22.9195], Tmin=(100,'K'), Tmax=(972.902,'K')), NASAPolynomial(coeffs=[10.7567,0.0340005,-1.2103e-05,2.24842e-09,-1.63577e-13,54915.4,-29.9714], Tmin=(972.902,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Isobutyl) + radical(bicyclo[2.1.1]hexane-C1)"""),
)

species(
    label = '[CH2]C=[C]CCC=C(261)',
    structure = SMILES('[CH2]C=[C]CCC=C'),
    E0 = (417.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,180,1115.41,1115.53],'cm^-1')),
        HinderedRotor(inertia=(0.280392,'amu*angstrom^2'), symmetry=1, barrier=(6.44676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.783496,'amu*angstrom^2'), symmetry=1, barrier=(18.0141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.783507,'amu*angstrom^2'), symmetry=1, barrier=(18.0144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.783479,'amu*angstrom^2'), symmetry=1, barrier=(18.0137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.717417,0.0627895,-3.41695e-05,4.32614e-09,1.58581e-12,50281.3,29.3918], Tmin=(100,'K'), Tmax=(1164.46,'K')), NASAPolynomial(coeffs=[13.3854,0.0329203,-1.32719e-05,2.42588e-09,-1.66831e-13,46405.9,-37.6255], Tmin=(1164.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(417.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C([CH2])CC=C(262)',
    structure = SMILES('[CH2]C=C([CH2])CC=C'),
    E0 = (316.668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,385.168],'cm^-1')),
        HinderedRotor(inertia=(0.105475,'amu*angstrom^2'), symmetry=1, barrier=(24.6595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07256,'amu*angstrom^2'), symmetry=1, barrier=(24.6602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105477,'amu*angstrom^2'), symmetry=1, barrier=(24.66,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0138396,'amu*angstrom^2'), symmetry=1, barrier=(3.2374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.721808,0.0589118,-1.09021e-05,-2.72676e-08,1.43784e-11,38215.9,27.7148], Tmin=(100,'K'), Tmax=(1005.16,'K')), NASAPolynomial(coeffs=[15.3787,0.0299208,-1.14165e-05,2.10891e-09,-1.49648e-13,33787.5,-50.4433], Tmin=(1005.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'CC=C1C=CCC1(263)',
    structure = SMILES('CC=C1C=CCC1'),
    E0 = (39.7287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74978,0.0307256,6.19802e-05,-9.65052e-08,3.73475e-11,4876.01,20.5542], Tmin=(100,'K'), Tmax=(981.005,'K')), NASAPolynomial(coeffs=[12.7113,0.0326214,-1.21579e-05,2.28974e-09,-1.66989e-13,483.475,-43.5449], Tmin=(981.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.7287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C=C=C1CCCC1(264)',
    structure = SMILES('C=C=C1CCCC1'),
    E0 = (99.2443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76479,0.0290945,6.91117e-05,-1.05636e-07,4.10159e-11,12034.8,20.5437], Tmin=(100,'K'), Tmax=(974.496,'K')), NASAPolynomial(coeffs=[13.4166,0.0313273,-1.13797e-05,2.14347e-09,-1.57591e-13,7386.92,-47.564], Tmin=(974.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.2443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'CC=C1CC=CC1(265)',
    structure = SMILES('CC=C1CC=CC1'),
    E0 = (39.5987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71615,0.0314147,6.08356e-05,-9.58982e-08,3.7277e-11,4861.65,20.5879], Tmin=(100,'K'), Tmax=(979.959,'K')), NASAPolynomial(coeffs=[12.9331,0.0323111,-1.1991e-05,2.25624e-09,-1.64619e-13,421.755,-44.7367], Tmin=(979.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.5987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH2]C=C1CC([CH2])C1(266)',
    structure = SMILES('[CH2]C1C[C](C=C)C1'),
    E0 = (349.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84451,0.0282419,7.0431e-05,-1.06927e-07,4.19642e-11,42125,22.4555], Tmin=(100,'K'), Tmax=(957.817,'K')), NASAPolynomial(coeffs=[12.5261,0.0317807,-1.05117e-05,1.89189e-09,-1.36796e-13,37870.3,-40.1456], Tmin=(957.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_T) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=C1CCC1[CH2](267)',
    structure = SMILES('[CH2]C1CC[C]1C=C'),
    E0 = (349.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8393,0.0286393,6.88026e-05,-1.04857e-07,4.11281e-11,42173.8,22.0985], Tmin=(100,'K'), Tmax=(959.034,'K')), NASAPolynomial(coeffs=[12.3773,0.0320754,-1.06909e-05,1.92635e-09,-1.39061e-13,37973.3,-39.6622], Tmin=(959.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Isobutyl) + radical(Allyl_T)"""),
)

species(
    label = 'C1CC2CCC=1C2(268)',
    structure = SMILES('C1CC2CCC=1C2'),
    E0 = (285.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43316,0.00278501,0.000154356,-2.00842e-07,7.63953e-11,34436.6,17.4738], Tmin=(100,'K'), Tmax=(952.089,'K')), NASAPolynomial(coeffs=[16.5014,0.0244859,-7.14086e-06,1.38321e-09,-1.12321e-13,28095.3,-68.9402], Tmin=(952.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'CH2(T)(22)',
    structure = SMILES('[CH2]'),
    E0 = (381.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([971.045,2816.03,3444.23],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71758,0.00127391,2.17347e-06,-3.48858e-09,1.65209e-12,45872.4,1.75298], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.14632,0.00303671,-9.96474e-07,1.50484e-10,-8.57336e-15,46041.3,4.72342], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(381.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=C1C[CH]CC1(269)',
    structure = SMILES('[CH]=C1C[CH]CC1'),
    E0 = (398.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51782,0.0186668,5.76023e-05,-7.93979e-08,2.92332e-11,47994.3,21.4312], Tmin=(100,'K'), Tmax=(995.235,'K')), NASAPolynomial(coeffs=[8.66658,0.0290998,-1.10932e-05,2.10147e-09,-1.5248e-13,45029.8,-16.9482], Tmin=(995.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1[CH]C[CH]C1(221)',
    structure = SMILES('C=CC1[CH]C[CH]C1'),
    E0 = (341.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23266,0.0235714,6.49863e-05,-8.67847e-08,3.07003e-11,41201,26.1785], Tmin=(100,'K'), Tmax=(1026.87,'K')), NASAPolynomial(coeffs=[8.03455,0.0401068,-1.63354e-05,3.12573e-09,-2.24967e-13,37946,-12.013], Tmin=(1026.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(cyclopentane)"""),
)

species(
    label = '[CH2][CH]C1C=CCC1(124)',
    structure = SMILES('[CH2][CH]C1C=CCC1'),
    E0 = (347.851,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,991.626,991.626,991.626,991.626,991.626,991.626,991.626,991.626,991.626,991.626,991.626,991.626,991.626,991.626,991.626,991.626,991.626,1148.73,3707.18],'cm^-1')),
        HinderedRotor(inertia=(0.0614462,'amu*angstrom^2'), symmetry=1, barrier=(6.19144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0614462,'amu*angstrom^2'), symmetry=1, barrier=(6.19144,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9669,0.0260017,7.01338e-05,-1.0027e-07,3.72998e-11,41926.5,25.3892], Tmin=(100,'K'), Tmax=(999.607,'K')), NASAPolynomial(coeffs=[11.4792,0.0349033,-1.36999e-05,2.64332e-09,-1.93939e-13,37678.4,-32.2351], Tmin=(999.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(Cs_S)"""),
)

species(
    label = 'C=[C]C1C[CH]CC1(224)',
    structure = SMILES('C=[C]C1C[CH]CC1'),
    E0 = (385.225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00247,0.0294826,5.11786e-05,-7.50358e-08,2.72869e-11,46416.3,24.921], Tmin=(100,'K'), Tmax=(1022.49,'K')), NASAPolynomial(coeffs=[8.81932,0.0388755,-1.55023e-05,2.93233e-09,-2.09786e-13,43137.2,-17.3358], Tmin=(1022.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(Cds_S)"""),
)

species(
    label = 'C=CC1C[CH][CH]C1(229)',
    structure = SMILES('C=CC1C[CH][CH]C1'),
    E0 = (333.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30517,0.0251212,5.31212e-05,-7.00495e-08,2.39883e-11,40153.8,25.9274], Tmin=(100,'K'), Tmax=(1051.65,'K')), NASAPolynomial(coeffs=[5.86148,0.0431494,-1.76005e-05,3.31424e-09,-2.34405e-13,37660.8,0.295044], Tmin=(1051.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(cyclopentane)"""),
)

species(
    label = '[CH]=CC1C[CH]CC1(232)',
    structure = SMILES('[CH]=CC1C[CH]CC1'),
    E0 = (394.479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96395,0.028396,5.9255e-05,-8.67127e-08,3.21795e-11,47532.6,24.9877], Tmin=(100,'K'), Tmax=(1003.82,'K')), NASAPolynomial(coeffs=[10.1315,0.0367347,-1.42982e-05,2.70942e-09,-1.9562e-13,43832.9,-24.7072], Tmin=(1003.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[C]1CCCC1(271)',
    structure = SMILES('[CH]C=C1CCCC1'),
    E0 = (293.385,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65124,0.030533,7.79509e-05,-1.14447e-07,4.34504e-11,35389.5,23.8251], Tmin=(100,'K'), Tmax=(983.486,'K')), NASAPolynomial(coeffs=[12.6961,0.0383638,-1.44495e-05,2.72595e-09,-1.9846e-13,30665.8,-42.2465], Tmin=(983.486,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]1CC2CC[C]1C2(272)',
    structure = SMILES('[CH]1CC2CC[C]1C2'),
    E0 = (379.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74453,0.00615605,0.000116902,-1.41981e-07,5.04792e-11,45661.8,20.6671], Tmin=(100,'K'), Tmax=(992.774,'K')), NASAPolynomial(coeffs=[8.74423,0.037626,-1.47191e-05,2.86182e-09,-2.1173e-13,41728.5,-22.0448], Tmin=(992.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(379.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(2-norbornyl) + radical(bridgehead_norbornyl)"""),
)

species(
    label = 'C=CC1C=CCC1(157)',
    structure = SMILES('C=CC1C=CCC1'),
    E0 = (76.6777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3658.67,'J/mol'), sigma=(6.25269,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.48 K, Pc=33.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85069,0.0251389,8.29943e-05,-1.2173e-07,4.71151e-11,9319.49,21.0387], Tmin=(100,'K'), Tmax=(968.097,'K')), NASAPolynomial(coeffs=[14.1361,0.0299511,-1.05691e-05,1.99737e-09,-1.48775e-13,4336.59,-51.2831], Tmin=(968.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.6777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C=CC1=CCCC1(273)',
    structure = SMILES('C=CC1=CCCC1'),
    E0 = (56.8914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81023,0.0276567,7.3231e-05,-1.0984e-07,4.24736e-11,6939.67,20.326], Tmin=(100,'K'), Tmax=(973.738,'K')), NASAPolynomial(coeffs=[13.3973,0.0313183,-1.13732e-05,2.14591e-09,-1.58105e-13,2252.98,-47.7448], Tmin=(973.738,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.8914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C=CC1CC=CC1(274)',
    structure = SMILES('C=CC1CC=CC1'),
    E0 = (74.5396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85758,0.0248011,8.42166e-05,-1.23179e-07,4.76743e-11,9062.27,20.8331], Tmin=(100,'K'), Tmax=(967.465,'K')), NASAPolynomial(coeffs=[14.2085,0.0298069,-1.04792e-05,1.98052e-09,-1.47708e-13,4048.37,-51.906], Tmin=(967.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.5396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH2]C1(C=C)C[CH]C1(275)',
    structure = SMILES('[CH2]C1(C=C)C[CH]C1'),
    E0 = (408.005,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63981,0.0374278,3.51044e-05,-6.36095e-08,2.4513e-11,49169.2,25.3773], Tmin=(100,'K'), Tmax=(1014.67,'K')), NASAPolynomial(coeffs=[10.996,0.0355697,-1.39269e-05,2.62503e-09,-1.88105e-13,45467.5,-28.7822], Tmin=(1014.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Neopentyl) + radical(cyclobutane)"""),
)

species(
    label = 'C=CC12CCC1C2(276)',
    structure = SMILES('C=CC12CCC1C2'),
    E0 = (210.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91457,0.0216721,9.57284e-05,-1.37187e-07,5.31838e-11,25417,19.8773], Tmin=(100,'K'), Tmax=(961.3,'K')), NASAPolynomial(coeffs=[15.0626,0.0280046,-9.40186e-06,1.77686e-09,-1.34554e-13,20068.7,-57.7053], Tmin=(961.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C=CCC(=C)C=C(277)',
    structure = SMILES('C=CCC(=C)C=C'),
    E0 = (138.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759909,0.05678,-3.27195e-06,-3.81897e-08,1.92154e-11,16842.2,25.8832], Tmin=(100,'K'), Tmax=(977.001,'K')), NASAPolynomial(coeffs=[16.3897,0.0268222,-9.52903e-06,1.73432e-09,-1.24e-13,12163.9,-57.4706], Tmin=(977.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(138.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61263,-0.00100893,2.49899e-06,-1.43376e-09,2.58637e-13,-1051.1,2.6527], Tmin=(100,'K'), Tmax=(1817.03,'K')), NASAPolynomial(coeffs=[2.97589,0.00164143,-7.1973e-07,1.25379e-10,-7.91538e-15,-1025.84,5.53764], Tmin=(1817.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.69489,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: BurkeH2O2"""),
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
    E0 = (481.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (366.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (467.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (386.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (528.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (419.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (383.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (486.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (488.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (460.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (431.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (294.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (268.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (294.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (509.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (509.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (285.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (812.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (452.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (510.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (567.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (398.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (530.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (438.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (379.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (283.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (283.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (323.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (565.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (267.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (354.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction155',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['[CH2][CH]C12CCC1C2(251)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.01129e+09,'s^-1'), n=0.843411, Ea=(221.884,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 221.3 to 221.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['C=C[C]1[CH]CCC1(230)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.81031e+09,'s^-1'), n=0.886125, Ea=(106.103,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] + [R2H_S_cy5;C_rad_out_1H;Cs_H_out_H/Cd] for rate rule [R2H_S_cy5;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction160',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['C[C]=C1C[CH]CC1(254)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction161',
    reactants = ['[CH2][C]=C1CCCC1(255)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.508e+06,'s^-1'), n=1.63, Ea=(74.8936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R4H_DSS;Cd_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction162',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['CC=C1[CH]C[CH]C1(256)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.15542e+07,'s^-1'), n=1, Ea=(267.985,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SDS;C_rad_out_single;Cs_H_out_H/NonDeC] + [R4H_SDS;C_rad_out_2H;Cs_H_out_1H] for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction163',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['CC=C1[CH][CH]CC1(257)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(987669,'s^-1'), n=1.76395, Ea=(159.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction164',
    reactants = ['CC=C1C[CH][CH]C1(258)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.90777e+06,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_SSMS;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction165',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['[CH]1CCC2([CH]C2)C1(259)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_NdNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction166',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['[CH2]C1[C]2CCC1C2(260)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.71e+08,'s^-1'), n=0.99, Ea=(228.626,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_Cs_HH_D;doublebond_intra_secNd;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 227.6 to 228.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction167',
    reactants = ['[CH2]C=[C]CCC=C(261)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.47e+07,'s^-1'), n=0.85, Ea=(43.5136,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 7 used for R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction168',
    reactants = ['[CH2]C=C([CH2])CC=C(262)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.12e+09,'s^-1'), n=0.63, Ea=(114.642,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction169',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['CC=C1C=CCC1(263)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction170',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['C=C=C1CCCC1(264)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction171',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['CC=C1CC=CC1(265)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction172',
    reactants = ['[CH2]C=C1CC([CH2])C1(266)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.31121e+11,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction173',
    reactants = ['[CH2]C=C1CCC1[CH2](267)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction174',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['C1CC2CCC=1C2(268)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.24579e+11,'s^-1'), n=0.1555, Ea=(25.532,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Cpri_rad_out_2H] + [R5;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R5_SSDS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
family: Birad_recombination
Ea raised from 22.4 to 25.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction175',
    reactants = ['CH2(T)(22)', '[CH]=C1C[CH]CC1(269)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.06864e+06,'m^3/(mol*s)'), n=0.369106, Ea=(32.7067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction177',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['C=CC1[CH]C[CH]C1(221)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.26e+10,'s^-1'), n=0.77, Ea=(192.464,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 164 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction178',
    reactants = ['[CH2][CH]C1C=CCC1(124)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.52488e+09,'s^-1'), n=1.21745, Ea=(162.572,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_Cd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction179',
    reactants = ['C=[C]C1C[CH]CC1(224)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_Cs2]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction180',
    reactants = ['C=CC1C[CH][CH]C1(229)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.60241e-19,'s^-1'), n=9.03667, Ea=(64.9775,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;Cs_H_out_OneDe] for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_Cd]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction181',
    reactants = ['[CH]=CC1C[CH]CC1(232)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.04e+10,'s^-1'), n=0.59, Ea=(135.98,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 196 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_Cs2
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction182',
    reactants = ['[CH]=C[C]1CCCC1(271)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.692e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R5HJ_2;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction183',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['[CH]1CC2CC[C]1C2(272)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(119.039,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_pri_2H;radadd_intra_csHCs] for rate rule [R5_cyc5_beta;doublebond_intra_pri_2H;radadd_intra_csHCs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 115.9 to 119.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction157',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['C=CC1C=CCC1(157)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['C=CC1=CCCC1(273)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction185',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['C=CC1CC=CC1(274)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction186',
    reactants = ['[CH2]C1(C=C)C[CH]C1(275)'],
    products = ['C=C[C]1C[CH]CC1(219)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction187',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['C=CC12CCC1C2(276)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.16053e+11,'s^-1'), n=0.0244333, Ea=(7.84361,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_single;Cpri_rad_out_single] for rate rule [R3_SS;C_rad_out_OneDe/Cs;Cpri_rad_out_H/NonDeC]
Euclidian distance = 4.12310562562
family: Birad_recombination"""),
)

reaction(
    label = 'reaction188',
    reactants = ['C=C[C]1C[CH]CC1(219)'],
    products = ['C=CCC(=C)C=C(277)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.28917e+12,'s^-1'), n=0.239384, Ea=(94.3194,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

network(
    label = '42',
    isomers = [
        'C=C[C]1C[CH]CC1(219)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '42',
    Tmin = (300,'K'),
    Tmax = (3000,'K'),
    Tcount = 8,
    Tlist = ([302.617,324.619,374.997,470.374,649.057,1000.02,1706.11,2761.25],'K'),
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

