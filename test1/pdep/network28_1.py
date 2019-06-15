species(
    label = 'C=CC1C=CCC1(113)',
    structure = SMILES('C=CC1C=CCC1'),
    E0 = (76.6775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3658.67,'J/mol'), sigma=(6.25269,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.48 K, Pc=33.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85063,0.02514,8.29919e-05,-1.21728e-07,4.71145e-11,9319.48,21.0388], Tmin=(100,'K'), Tmax=(968.098,'K')), NASAPolynomial(coeffs=[14.1364,0.0299506,-1.05687e-05,1.99729e-09,-1.48769e-13,4336.5,-51.2849], Tmin=(968.098,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.6775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C2H3(32)',
    structure = SMILES('[CH]=C'),
    E0 = (286.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,677.08,1086.68,3788.01],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]1C=CCC1(159)',
    structure = SMILES('[CH]1C=CCC1'),
    E0 = (125.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22816,-0.0047173,0.000117611,-1.45673e-07,5.42071e-11,15155.9,11.3916], Tmin=(100,'K'), Tmax=(955.988,'K')), NASAPolynomial(coeffs=[11.0452,0.0181893,-5.59294e-06,1.09697e-09,-8.81254e-14,11120,-39.2607], Tmin=(955.988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-allyl)"""),
)

species(
    label = 'C=C[C]1C[CH]CC1(167)',
    structure = SMILES('[CH2]C=C1C[CH]CC1'),
    E0 = (260.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92032,0.0283861,6.30056e-05,-9.23827e-08,3.45647e-11,31370.7,25.2891], Tmin=(100,'K'), Tmax=(996.92,'K')), NASAPolynomial(coeffs=[10.7644,0.0362078,-1.39244e-05,2.6374e-09,-1.91075e-13,27455.3,-28.1444], Tmin=(996.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(cyclopentane)"""),
)

species(
    label = '[CH2]C[C]1C=CCC1(168)',
    structure = SMILES('[CH2]CC1=C[CH]CC1'),
    E0 = (263.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92326,0.0247926,8.01745e-05,-1.14601e-07,4.33528e-11,31776,22.2283], Tmin=(100,'K'), Tmax=(983.383,'K')), NASAPolynomial(coeffs=[12.8072,0.0329833,-1.23421e-05,2.36913e-09,-1.7544e-13,27098.7,-42.9955], Tmin=(983.383,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-allyl) + radical(RCCJ)"""),
)

species(
    label = 'C=CC1[CH]C[CH]C1(169)',
    structure = SMILES('C=CC1[CH]C[CH]C1'),
    E0 = (341.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23266,0.0235713,6.49855e-05,-8.67829e-08,3.06994e-11,41201,26.178], Tmin=(100,'K'), Tmax=(1026.88,'K')), NASAPolynomial(coeffs=[8.03445,0.0401071,-1.63356e-05,3.12577e-09,-2.24969e-13,37946,-12.013], Tmin=(1026.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(Cs_S)"""),
)

species(
    label = '[CH2]CC1[CH]CC=C1(170)',
    structure = SMILES('[CH2]CC1[CH]CC=C1'),
    E0 = (340.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93987,0.0278243,6.32132e-05,-9.20454e-08,3.42118e-11,41019.3,25.2801], Tmin=(100,'K'), Tmax=(1003.69,'K')), NASAPolynomial(coeffs=[10.8994,0.0358355,-1.40944e-05,2.69993e-09,-1.96548e-13,37018.7,-28.9474], Tmin=(1003.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(340.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(cyclopentene-4)"""),
)

species(
    label = '[CH2]CC1[C]=CCC1(171)',
    structure = SMILES('[CH2]CC1[C]=CCC1'),
    E0 = (410.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73942,0.0318838,5.64504e-05,-8.87183e-08,3.39863e-11,49506.7,24.125], Tmin=(100,'K'), Tmax=(994.077,'K')), NASAPolynomial(coeffs=[12.2707,0.0336632,-1.28625e-05,2.44893e-09,-1.78677e-13,45231.2,-37.5933], Tmin=(994.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C=[C]C1C[CH]CC1(172)',
    structure = SMILES('C=[C]C1C[CH]CC1'),
    E0 = (385.225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,1685,370,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00247,0.0294825,5.1178e-05,-7.50344e-08,2.72861e-11,46416.3,24.9206], Tmin=(100,'K'), Tmax=(1022.5,'K')), NASAPolynomial(coeffs=[8.81919,0.0388759,-1.55025e-05,2.93237e-09,-2.09789e-13,43137.2,-17.3357], Tmin=(1022.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(cyclopentane)"""),
)

species(
    label = 'C=CC1[CH]CC[CH]1(173)',
    structure = SMILES('C=CC1[CH]CC[CH]1'),
    E0 = (350.586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15767,0.022042,7.68224e-05,-1.03551e-07,3.74584e-11,42248.3,26.4384], Tmin=(100,'K'), Tmax=(1011.49,'K')), NASAPolynomial(coeffs=[10.273,0.0369589,-1.50121e-05,2.92385e-09,-2.14443e-13,38201.7,-24.6937], Tmin=(1011.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]C1[CH]CC=C1(174)',
    structure = SMILES('C[CH]C1[CH]CC=C1'),
    E0 = (329.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09315,0.0258438,6.25329e-05,-8.6701e-08,3.11223e-11,39725,25.5016], Tmin=(100,'K'), Tmax=(1025.35,'K')), NASAPolynomial(coeffs=[9.27055,0.038632,-1.58446e-05,3.05498e-09,-2.21212e-13,36109,-19.7603], Tmin=(1025.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(cyclopentene-4)"""),
)

species(
    label = 'C[CH]C1[C]=CCC1(175)',
    structure = SMILES('C[CH]C1[C]=CCC1'),
    E0 = (400.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89106,0.0299284,5.5653e-05,-8.31754e-08,3.07899e-11,48212.4,24.3521], Tmin=(100,'K'), Tmax=(1014.6,'K')), NASAPolynomial(coeffs=[10.6168,0.0364993,-1.46344e-05,2.80891e-09,-2.03739e-13,44332.9,-28.2638], Tmin=(1014.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C=[C]C1[CH]CCC1(176)',
    structure = SMILES('C=[C]C1[CH]CCC1'),
    E0 = (393.886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,1685,370,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92877,0.0279363,6.30834e-05,-9.19057e-08,3.40956e-11,47463.5,25.1765], Tmin=(100,'K'), Tmax=(1006.81,'K')), NASAPolynomial(coeffs=[11.0626,0.0357204,-1.41751e-05,2.72961e-09,-1.99197e-13,43390.6,-30.044], Tmin=(1006.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(Cs_S)"""),
)

species(
    label = 'C=CC1C[CH][CH]C1(177)',
    structure = SMILES('C=CC1C[CH][CH]C1'),
    E0 = (333.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30518,0.0251211,5.31209e-05,-7.00483e-08,2.39877e-11,40153.8,25.9269], Tmin=(100,'K'), Tmax=(1051.66,'K')), NASAPolynomial(coeffs=[5.86137,0.0431497,-1.76007e-05,3.31428e-09,-2.34408e-13,37660.8,0.29507], Tmin=(1051.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(cyclopentane)"""),
)

species(
    label = '[CH2]CC1[CH]C=CC1(97)',
    structure = SMILES('[CH2]CC1[CH]C=CC1'),
    E0 = (279.478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98122,0.0218598,9.08508e-05,-1.27193e-07,4.81716e-11,33706.3,21.4378], Tmin=(100,'K'), Tmax=(977.251,'K')), NASAPolynomial(coeffs=[13.4494,0.0317689,-1.16184e-05,2.23757e-09,-1.67385e-13,28750.2,-47.5148], Tmin=(977.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(cyclopentene-allyl)"""),
)

species(
    label = 'C=C[C]1[CH]CCC1(178)',
    structure = SMILES('[CH2]C=C1[CH]CCC1'),
    E0 = (215.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85517,0.0228285,9.52346e-05,-1.35843e-07,5.21847e-11,25995.3,21.851], Tmin=(100,'K'), Tmax=(970.064,'K')), NASAPolynomial(coeffs=[15.051,0.0297956,-1.04487e-05,2.01263e-09,-1.5254e-13,20547.1,-56.2957], Tmin=(970.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_P) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]CC1C=[C]CC1(179)',
    structure = SMILES('[CH2]CC1C=[C]CC1'),
    E0 = (410.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73942,0.0318838,5.64504e-05,-8.87183e-08,3.39863e-11,49506.7,24.125], Tmin=(100,'K'), Tmax=(994.077,'K')), NASAPolynomial(coeffs=[12.2707,0.0336632,-1.28625e-05,2.44893e-09,-1.78677e-13,45231.2,-37.5933], Tmin=(994.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-vinyl) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CC1C[CH]CC1(180)',
    structure = SMILES('[CH]=CC1C[CH]CC1'),
    E0 = (394.479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96392,0.0283963,5.92531e-05,-8.67094e-08,3.21779e-11,47532.6,24.9873], Tmin=(100,'K'), Tmax=(1003.84,'K')), NASAPolynomial(coeffs=[10.1315,0.0367348,-1.42983e-05,2.70944e-09,-1.95621e-13,43832.9,-24.7077], Tmin=(1003.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentane) + radical(Cds_P)"""),
)

species(
    label = 'C[CH]C1C=C[CH]C1(181)',
    structure = SMILES('C[CH]C1[CH]C=CC1'),
    E0 = (268.774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13488,0.0198907,9.00477e-05,-1.21553e-07,4.48894e-11,32411.9,21.657], Tmin=(100,'K'), Tmax=(990.951,'K')), NASAPolynomial(coeffs=[11.7306,0.0347109,-1.34495e-05,2.61121e-09,-1.93559e-13,27880.7,-37.8168], Tmin=(990.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-allyl) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]C1C=[C]CC1(182)',
    structure = SMILES('C[CH]C1C=[C]CC1'),
    E0 = (400.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89106,0.0299284,5.5653e-05,-8.31754e-08,3.07899e-11,48212.4,24.3521], Tmin=(100,'K'), Tmax=(1014.6,'K')), NASAPolynomial(coeffs=[10.6168,0.0364993,-1.46344e-05,2.80891e-09,-2.03739e-13,44332.9,-28.2638], Tmin=(1014.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopentene-vinyl) + radical(Cs_S)"""),
)

species(
    label = '[CH]=CC1[CH]CCC1(73)',
    structure = SMILES('[CH]=CC1[CH]CCC1'),
    E0 = (403.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88815,0.0268695,7.1117e-05,-1.03571e-07,3.90046e-11,48579.9,25.251], Tmin=(100,'K'), Tmax=(993.5,'K')), NASAPolynomial(coeffs=[12.4137,0.0335161,-1.29355e-05,2.49852e-09,-1.84363e-13,44069.1,-37.6363], Tmin=(993.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_P) + radical(Cs_S)"""),
)

species(
    label = '[CH]=CC(C=C)C[CH2](183)',
    structure = SMILES('[CH]=CC(C=C)C[CH2]'),
    E0 = (483.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,1188.92],'cm^-1')),
        HinderedRotor(inertia=(0.588323,'amu*angstrom^2'), symmetry=1, barrier=(13.5267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127996,'amu*angstrom^2'), symmetry=1, barrier=(2.94287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587045,'amu*angstrom^2'), symmetry=1, barrier=(13.4973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587685,'amu*angstrom^2'), symmetry=1, barrier=(13.512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.887814,0.0683421,-5.20658e-05,2.16265e-08,-3.81158e-12,58292.7,28.0129], Tmin=(100,'K'), Tmax=(1295.93,'K')), NASAPolynomial(coeffs=[11.0732,0.0369036,-1.56761e-05,2.90627e-09,-2.0018e-13,55652.8,-23.7659], Tmin=(1295.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CCC[CH]C=C(184)',
    structure = SMILES('[CH]=CCCC=C[CH2]'),
    E0 = (426.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.00315625,'amu*angstrom^2'), symmetry=1, barrier=(2.24479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805887,'amu*angstrom^2'), symmetry=1, barrier=(18.5289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.8062,'amu*angstrom^2'), symmetry=1, barrier=(18.5361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805567,'amu*angstrom^2'), symmetry=1, barrier=(18.5216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.726582,0.0612088,-2.47002e-05,-8.65034e-09,6.80657e-12,51395.5,29.2831], Tmin=(100,'K'), Tmax=(1059.96,'K')), NASAPolynomial(coeffs=[14.0379,0.0318345,-1.26495e-05,2.33585e-09,-1.63408e-13,47401.8,-41.2388], Tmin=(1059.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=CC([CH2])C=C(185)',
    structure = SMILES('[CH2]C=CC([CH2])C=C'),
    E0 = (375.736,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,3851.24],'cm^-1')),
        HinderedRotor(inertia=(1.07702,'amu*angstrom^2'), symmetry=1, barrier=(24.7628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.89484,'amu*angstrom^2'), symmetry=1, barrier=(89.55,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.530729,'amu*angstrom^2'), symmetry=1, barrier=(12.2025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743502,'amu*angstrom^2'), symmetry=1, barrier=(89.5551,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.547793,0.0685847,-5.08023e-05,2.03371e-08,-3.35663e-12,45320.9,28.6747], Tmin=(100,'K'), Tmax=(1419.14,'K')), NASAPolynomial(coeffs=[13.6134,0.0317576,-1.18765e-05,2.05084e-09,-1.3524e-13,41612.6,-38.9326], Tmin=(1419.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC1[C]CCC1(186)',
    structure = SMILES('C=CC1[C]CCC1'),
    E0 = (336.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64708,0.0310185,6.76139e-05,-1.06291e-07,4.16937e-11,40553.4,22.2396], Tmin=(100,'K'), Tmax=(973.119,'K')), NASAPolynomial(coeffs=[14.3654,0.0303717,-1.09764e-05,2.07413e-09,-1.53306e-13,35633.4,-51.3336], Tmin=(973.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C=CC1C[C]CC1(187)',
    structure = SMILES('C=CC1C[C]CC1'),
    E0 = (325.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.628,0.031315,6.74292e-05,-1.06583e-07,4.19212e-11,39303.8,22.3406], Tmin=(100,'K'), Tmax=(972.075,'K')), NASAPolynomial(coeffs=[14.5506,0.0300787,-1.08097e-05,2.04096e-09,-1.51007e-13,34337.5,-52.2596], Tmin=(972.075,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C[C]C1C=CCC1(188)',
    structure = SMILES('C[C]C1C=CCC1'),
    E0 = (322.947,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64397,0.032184,6.20534e-05,-9.90897e-08,3.88124e-11,38944,21.488], Tmin=(100,'K'), Tmax=(977.051,'K')), NASAPolynomial(coeffs=[13.7363,0.031405,-1.15572e-05,2.17908e-09,-1.59818e-13,34255.2,-48.4719], Tmin=(977.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]CC1C=CCC1(189)',
    structure = SMILES('[CH]CC1C=CCC1'),
    E0 = (397.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57445,0.0330059,6.27492e-05,-1.01638e-07,4.01203e-11,47857.1,21.9008], Tmin=(100,'K'), Tmax=(974.19,'K')), NASAPolynomial(coeffs=[14.4977,0.0304445,-1.10658e-05,2.08867e-09,-1.53981e-13,42942.7,-52.4081], Tmin=(974.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(397.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]1CC2[CH]C1CC2(190)',
    structure = SMILES('[CH]1CC2[CH]C1CC2'),
    E0 = (343.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.82119,-0.000636564,0.000145277,-1.77133e-07,6.42308e-11,41415.4,21.6034], Tmin=(100,'K'), Tmax=(975.933,'K')), NASAPolynomial(coeffs=[11.2195,0.0334231,-1.23273e-05,2.42807e-09,-1.84964e-13,36514.9,-35.4173], Tmin=(975.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(2-norbornyl) + radical(7-norbornyl)"""),
)

species(
    label = '[CH]1CCC2[CH]CC12(70)',
    structure = SMILES('[CH]1CCC2[CH]CC12'),
    E0 = (319.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66787,0.00592453,0.000123207,-1.51788e-07,5.46322e-11,38531.1,22.9388], Tmin=(100,'K'), Tmax=(986.711,'K')), NASAPolynomial(coeffs=[10.1908,0.0357286,-1.37713e-05,2.69645e-09,-2.01726e-13,34111.1,-28.1294], Tmin=(986.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-2) + radical(bicyclo[3.2.0]heptane-C5-6)"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: FFCM1(-)"""),
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
    E0 = (418.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (283.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (285.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (365.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (403.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (499.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (484.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (413.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (337.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (439.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (433.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (358.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (318.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (254.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (419.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (428.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (307.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (408.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (428.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (490.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (433.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (382.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (427.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (417.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (414.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (453.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (438.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (319.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction139',
    reactants = ['C2H3(32)', '[CH]1C=CCC1(159)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.32655e+09,'m^3/(mol*s)'), n=-0.416379, Ea=(6.36255,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;C_rad/H/CdCs] for rate rule [Cd_pri_rad;C_rad/H/CdCs]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction147',
    reactants = ['C=C[C]1C[CH]CC1(167)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction148',
    reactants = ['[CH2]C[C]1C=CCC1(168)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=CC1[CH]C[CH]C1(169)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction150',
    reactants = ['[CH2]CC1[CH]CC=C1(170)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction151',
    reactants = ['[CH2]CC1[C]=CCC1(171)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction152',
    reactants = ['C=[C]C1C[CH]CC1(172)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.94203e+09,'s^-1'), n=0.37, Ea=(99.6901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad_NDe] + [R3radExo;Y_rad;XH_Rrad_NDe] for rate rule [R3radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction153',
    reactants = ['C=CC1[CH]CC[CH]1(173)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction154',
    reactants = ['C[CH]C1[CH]CC=C1(174)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction155',
    reactants = ['C[CH]C1[C]=CCC1(175)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.34494e+09,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction156',
    reactants = ['C=[C]C1[CH]CCC1(176)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction157',
    reactants = ['C=CC1C[CH][CH]C1(177)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.05689e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction158',
    reactants = ['[CH2]CC1[CH]C=CC1(97)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction159',
    reactants = ['C=C[C]1[CH]CCC1(178)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction160',
    reactants = ['[CH2]CC1C=[C]CC1(179)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction161',
    reactants = ['[CH]=CC1C[CH]CC1(180)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction162',
    reactants = ['C[CH]C1C=C[CH]C1(181)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction163',
    reactants = ['C[CH]C1C=[C]CC1(182)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction164',
    reactants = ['[CH]=CC1[CH]CCC1(73)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction165',
    reactants = ['[CH]=CC(C=C)C[CH2](183)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Ypri_rad_out] for rate rule [R5_SSSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction166',
    reactants = ['[CH]=CCC[CH]C=C(184)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R5_SSSD;C_rad_out_H/OneDe;CdsinglepriH_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction167',
    reactants = ['[CH2]C=CC([CH2])C=C(185)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R5_SSDS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction168',
    reactants = ['C=CC1[C]CCC1(186)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.83662e+12,'s^-1'), n=0.345439, Ea=(91.3359,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction169',
    reactants = ['C=CC1C[C]CC1(187)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.83662e+12,'s^-1'), n=0.345439, Ea=(91.3359,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction170',
    reactants = ['C[C]C1C=CCC1(188)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.75493e+12,'s^-1'), n=0.345439, Ea=(91.3359,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction171',
    reactants = ['[CH]CC1C=CCC1(189)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.57909e+12,'s^-1'), n=0.325393, Ea=(56.3425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;singletcarbene;CH2(C)] + [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH2(C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction172',
    reactants = ['[CH]1CC2[CH]C1CC2(190)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.31474e+12,'s^-1'), n=0.239323, Ea=(94.3216,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction173',
    reactants = ['[CH]1CCC2[CH]CC12(70)'],
    products = ['C=CC1C=CCC1(113)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R7JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

network(
    label = '28',
    isomers = [
        'C=CC1C=CCC1(113)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '28',
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

