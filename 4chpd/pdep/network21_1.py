species(
    label = '[O]OC1C=C[CH]CC=C1(192)',
    structure = SMILES('[O]OC1[CH]C=CCC=C1'),
    E0 = (211.195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32235,-0.008698,0.00028297,-4.0762e-07,1.69525e-10,25560.5,31.5361], Tmin=(100,'K'), Tmax=(909.554,'K')), NASAPolynomial(coeffs=[45.5373,-0.0271949,2.33041e-05,-4.61277e-09,2.96105e-13,10239.3,-217.593], Tmin=(909.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(ROOJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1483.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53764,-0.00122828,5.36759e-06,-4.93128e-09,1.45955e-12,-1037.99,4.6718], Tmin=(100,'K'), Tmax=(1087.71,'K')), NASAPolynomial(coeffs=[3.16427,0.00169454,-8.00335e-07,1.5903e-10,-1.14891e-14,-1048.45,6.08303], Tmin=(1087.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C7H8(11)',
    structure = SMILES('C1=CC=CCC=C1'),
    E0 = (165.503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3872.77,'J/mol'), sigma=(6.33517,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=604.92 K, Pc=34.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.976321,0.0470646,2.43918e-05,-7.34711e-08,3.37246e-11,20032,17.1406], Tmin=(100,'K'), Tmax=(948.557,'K')), NASAPolynomial(coeffs=[19.6256,0.016921,-4.63468e-06,8.3158e-10,-6.48659e-14,14312.2,-83.3476], Tmin=(948.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'C7H8O(2)(14)',
    structure = SMILES('C1=CC2OC2C=CC1'),
    E0 = (77.1493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3994.42,'J/mol'), sigma=(6.56534,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=623.92 K, Pc=32.03 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61108,-0.0046343,0.000176491,-2.20656e-07,8.14198e-11,9360.54,20.2194], Tmin=(100,'K'), Tmax=(974.25,'K')), NASAPolynomial(coeffs=[17.1087,0.0261758,-1.00264e-05,2.14765e-09,-1.74993e-13,2248.65,-71.3466], Tmin=(974.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.1493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide)"""),
)

species(
    label = 'O(20)',
    structure = SMILES('[O]'),
    E0 = (243.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29226.7,5.11107], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29226.7,5.11107], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.005,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH]1CC=CC2[CH]C1OO2(897)',
    structure = SMILES('[CH]1CC=CC2[CH]C1OO2'),
    E0 = (283.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.422866,0.0894033,-7.44221e-05,3.02192e-08,-4.88966e-12,34252.9,2.66695], Tmin=(100,'K'), Tmax=(1464.61,'K')), NASAPolynomial(coeffs=[21.2162,0.0303048,-1.38958e-05,2.66867e-09,-1.86961e-13,27914.3,-109.986], Tmin=(1464.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(12dioxolane) + ring(Cycloheptene) + ring(12dioxolane) + radical(CCJCOOH) + radical(CCJCOOH)"""),
)

species(
    label = 'H(19)',
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
    label = '[O]OC1C=CCC=CC=1(898)',
    structure = SMILES('[O]OC1C=CCC=CC=1'),
    E0 = (228.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,492.5,1135,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.437655,0.0854544,-8.07193e-05,3.83373e-08,-7.11734e-12,27605.3,26.9615], Tmin=(100,'K'), Tmax=(1317.31,'K')), NASAPolynomial(coeffs=[20.7361,0.0211594,-7.50646e-06,1.28508e-09,-8.54296e-14,22026.9,-81.0249], Tmin=(1317.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1C=C=CCC=C1(899)',
    structure = SMILES('[O]OC1C=C=CCC=C1'),
    E0 = (361.305,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,492.5,1135,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0022,0.0539597,-8.21476e-06,-2.18542e-08,1.01855e-11,43572.9,25.1592], Tmin=(100,'K'), Tmax=(1094.41,'K')), NASAPolynomial(coeffs=[13.409,0.0335594,-1.44451e-05,2.76889e-09,-1.9708e-13,39363.4,-42.6407], Tmin=(1094.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(ROOJ)"""),
)

species(
    label = '[O]O[C]1C=CCC=CC1(900)',
    structure = SMILES('[O]OC1=C[CH]CC=CC1'),
    E0 = (293.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,492.5,1135,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34674,0.00279873,0.000224235,-3.31651e-07,1.39272e-10,35435.7,30.308], Tmin=(100,'K'), Tmax=(906.949,'K')), NASAPolynomial(coeffs=[38.8684,-0.0175929,1.79907e-05,-3.65376e-09,2.3691e-13,22662.3,-179.946], Tmin=(906.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(ROOJ) + radical(Allyl_S)"""),
)

species(
    label = '[O]OC1C=CCC=[C]C1(901)',
    structure = SMILES('[O]OC1C=CCC=[C]C1'),
    E0 = (378.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,492.5,1135,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04772,0.0099179,0.000208745,-3.18671e-07,1.35322e-10,45728.2,32.1971], Tmin=(100,'K'), Tmax=(906.825,'K')), NASAPolynomial(coeffs=[40.2276,-0.0189687,1.84402e-05,-3.73224e-09,2.42436e-13,32704.2,-185.622], Tmin=(906.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = 'OO[C]1[CH]C=CCC=C1(902)',
    structure = SMILES('OO[C]1[CH]C=CCC=C1'),
    E0 = (246.089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1310,387.5,850,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18026,-0.00240349,0.000261962,-3.81253e-07,1.58517e-10,29759.4,32.2737], Tmin=(100,'K'), Tmax=(913.048,'K')), NASAPolynomial(coeffs=[43.9751,-0.0225463,2.01419e-05,-3.9588e-09,2.50116e-13,14969.5,-208.491], Tmin=(913.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C2CsJOOH) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[O]OC1[C]=CCC=CC1(903)',
    structure = SMILES('[O]OC1[C]=CCC=CC1'),
    E0 = (378.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,492.5,1135,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04772,0.0099179,0.000208745,-3.18671e-07,1.35322e-10,45728.2,32.1971], Tmin=(100,'K'), Tmax=(906.825,'K')), NASAPolynomial(coeffs=[40.2276,-0.0189687,1.84402e-05,-3.73224e-09,2.42436e-13,32704.2,-185.622], Tmin=(906.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]OC1C=CC[C]=CC1(904)',
    structure = SMILES('[O]OC1C=CC[C]=CC1'),
    E0 = (378.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04772,0.0099179,0.000208745,-3.18671e-07,1.35322e-10,45728.2,32.1971], Tmin=(100,'K'), Tmax=(906.825,'K')), NASAPolynomial(coeffs=[40.2276,-0.0189687,1.84402e-05,-3.73224e-09,2.42436e-13,32704.2,-185.622], Tmin=(906.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'OOC1[C]=CCC=C[CH]1(905)',
    structure = SMILES('OOC1[C]=CCC=C[CH]1'),
    E0 = (297.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01247,-0.000346504,0.000262281,-3.86774e-07,1.61878e-10,35894,32.3627], Tmin=(100,'K'), Tmax=(911.582,'K')), NASAPolynomial(coeffs=[46.1024,-0.0263234,2.22048e-05,-4.36354e-09,2.78041e-13,20532.1,-220.166], Tmin=(911.582,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[O]OC1C=[C]CC=CC1(906)',
    structure = SMILES('[O]OC1C=[C]CC=CC1'),
    E0 = (378.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04772,0.0099179,0.000208745,-3.18671e-07,1.35322e-10,45728.2,32.1971], Tmin=(100,'K'), Tmax=(906.825,'K')), NASAPolynomial(coeffs=[40.2276,-0.0189687,1.84402e-05,-3.73224e-09,2.42436e-13,32704.2,-185.622], Tmin=(906.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1C=CC=C[CH]C1(209)',
    structure = SMILES('[O]OC1C=C[CH]C=CC1'),
    E0 = (242.776,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,492.5,1135,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38868,-0.00695603,0.000270838,-3.92395e-07,1.64094e-10,29353.7,30.1708], Tmin=(100,'K'), Tmax=(905.169,'K')), NASAPolynomial(coeffs=[43.6922,-0.0256907,2.31394e-05,-4.6643e-09,3.04355e-13,14804.5,-207.768], Tmin=(905.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(ROOJ)"""),
)

species(
    label = 'OOC1[CH]C=CC[C]=C1(907)',
    structure = SMILES('OOC1[CH]C=CC[C]=C1'),
    E0 = (297.032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01247,-0.000346504,0.000262281,-3.86774e-07,1.61878e-10,35894,32.3627], Tmin=(100,'K'), Tmax=(911.582,'K')), NASAPolynomial(coeffs=[46.1024,-0.0263234,2.22048e-05,-4.36354e-09,2.78041e-13,20532.1,-220.166], Tmin=(911.582,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
)

species(
    label = 'OOC1[CH][C]=CCC=C1(908)',
    structure = SMILES('OOC1[CH][C]=CCC=C1'),
    E0 = (297.032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1310,387.5,850,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01247,-0.000346504,0.000262281,-3.86774e-07,1.61878e-10,35894,32.3627], Tmin=(100,'K'), Tmax=(911.582,'K')), NASAPolynomial(coeffs=[46.1024,-0.0263234,2.22048e-05,-4.36354e-09,2.78041e-13,20532.1,-220.166], Tmin=(911.582,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'OOC1[CH]C=C[CH]C=C1(909)',
    structure = SMILES('OOC1[CH]C=C[CH]C=C1'),
    E0 = (160.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35415,-0.0172283,0.000324396,-4.60517e-07,1.90653e-10,19519.5,30.3338], Tmin=(100,'K'), Tmax=(909.538,'K')), NASAPolynomial(coeffs=[49.5602,-0.0330339,2.68974e-05,-5.29405e-09,3.39831e-13,2635.13,-242.274], Tmin=(909.538,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'OOC1[CH]C=[C]CC=C1(910)',
    structure = SMILES('OOC1[CH]C=[C]CC=C1'),
    E0 = (297.032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1310,387.5,850,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01247,-0.000346504,0.000262281,-3.86774e-07,1.61878e-10,35894,32.3627], Tmin=(100,'K'), Tmax=(911.582,'K')), NASAPolynomial(coeffs=[46.1024,-0.0263234,2.22048e-05,-4.36354e-09,2.78041e-13,20532.1,-220.166], Tmin=(911.582,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
)

species(
    label = '[CH]1C=C[CH]CC=C1(125)',
    structure = SMILES('[CH]1C=C[CH]CC=C1'),
    E0 = (333.186,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.70798,-0.040239,0.000339088,-4.54205e-07,1.84685e-10,40184,19.5412], Tmin=(100,'K'), Tmax=(904.636,'K')), NASAPolynomial(coeffs=[40.627,-0.029508,2.54911e-05,-5.11056e-09,3.33198e-13,26023.7,-199.942], Tmin=(904.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = '[O]OC1C=CCC2[CH]C21(911)',
    structure = SMILES('[O]OC1C=CCC2[CH]C21'),
    E0 = (383.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08592,0.0505844,6.88401e-06,-3.73558e-08,1.53616e-11,46212.5,26.8996], Tmin=(100,'K'), Tmax=(1071.09,'K')), NASAPolynomial(coeffs=[13.0077,0.0363321,-1.55475e-05,2.99083e-09,-2.14008e-13,41922.3,-39.5403], Tmin=(1071.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_2) + radical(ROOJ) + radical(cyclopropane)"""),
)

species(
    label = '[O]OC1[CH]C2CC=CC21(912)',
    structure = SMILES('[O]OC1[CH]C2CC=CC21'),
    E0 = (376.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04582,0.0497291,1.46013e-05,-4.88471e-08,2.02651e-11,45374,26.862], Tmin=(100,'K'), Tmax=(1032.87,'K')), NASAPolynomial(coeffs=[14.0702,0.0344465,-1.42618e-05,2.73802e-09,-1.97406e-13,40808.2,-45.4724], Tmin=(1032.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(CCJCOOH) + radical(ROOJ)"""),
)

species(
    label = '[CH]1C=CCC2[CH]C1OO2(913)',
    structure = SMILES('[CH]1C=CCC2[CH]C1OO2'),
    E0 = (200.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.814211,0.0956653,-8.23658e-05,3.45025e-08,-5.7213e-12,24267,-1.94009], Tmin=(100,'K'), Tmax=(1441.84,'K')), NASAPolynomial(coeffs=[23.1478,0.0291888,-1.32075e-05,2.52553e-09,-1.76806e-13,17357.1,-126.311], Tmin=(1441.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(12dioxolane) + ring(Cycloheptene) + ring(12dioxolane) + radical(CCJCOOH) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]1[CH]C2CC=CC1OO2(692)',
    structure = SMILES('[CH]1[CH]C2CC=CC1OO2'),
    E0 = (287.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.422866,0.0894033,-7.44221e-05,3.02192e-08,-4.88966e-12,34756.1,2.81792], Tmin=(100,'K'), Tmax=(1464.61,'K')), NASAPolynomial(coeffs=[21.2162,0.0303048,-1.38958e-05,2.66867e-09,-1.86961e-13,28417.5,-109.835], Tmin=(1464.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(12dioxane) + ring(Cycloheptane) + ring(12dioxane) + radical(CCJCOOH) + radical(CCJCOOH)"""),
)

species(
    label = '[O]OC1C2[CH]CC=CC21(914)',
    structure = SMILES('[O]OC1C2[CH]CC=CC21'),
    E0 = (351.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10533,0.0496157,1.05692e-05,-4.20312e-08,1.72533e-11,42396.6,26.4639], Tmin=(100,'K'), Tmax=(1054.32,'K')), NASAPolynomial(coeffs=[13.1792,0.0357238,-1.50731e-05,2.89441e-09,-2.07556e-13,38076.9,-40.8364], Tmin=(1054.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + radical(ROOJ) + radical(Cs_S)"""),
)

species(
    label = '[CH]1C=CC2OOC2[CH]C1(343)',
    structure = SMILES('[CH]1C=CC[CH]C2OOC12'),
    E0 = (200.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.814211,0.0956653,-8.23658e-05,3.45025e-08,-5.7213e-12,24267,-1.94009], Tmin=(100,'K'), Tmax=(1441.84,'K')), NASAPolynomial(coeffs=[23.1478,0.0291888,-1.32075e-05,2.52553e-09,-1.76806e-13,17357.1,-126.311], Tmin=(1441.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(12dioxetane) + ring(Cycloheptene) + ring(12dioxetane) + radical(CCJCOOH) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=CCC=CC=CO[O](915)',
    structure = SMILES('[CH]=CCC=CC=CO[O]'),
    E0 = (455.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3120,650,792.5,1650,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.686623,'amu*angstrom^2'), symmetry=1, barrier=(15.7868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.685361,'amu*angstrom^2'), symmetry=1, barrier=(15.7578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.687487,'amu*angstrom^2'), symmetry=1, barrier=(15.8067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.6885,'amu*angstrom^2'), symmetry=1, barrier=(15.83,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.296345,0.0847425,-7.32033e-05,2.76186e-08,-2.55169e-12,54926.1,35.9432], Tmin=(100,'K'), Tmax=(1018.32,'K')), NASAPolynomial(coeffs=[19.3429,0.0243029,-8.78026e-06,1.55095e-09,-1.06657e-13,50060.2,-63.414], Tmin=(1018.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(455.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(ROOJ)"""),
)

species(
    label = 'OOC1C=CCC=CC=1(916)',
    structure = SMILES('OOC1C=CCC=CC=1'),
    E0 = (76.1053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.49989,0.0867669,-6.85138e-05,1.99819e-08,2.10294e-13,9325.96,26.4197], Tmin=(100,'K'), Tmax=(1039.16,'K')), NASAPolynomial(coeffs=[21.1127,0.0247488,-9.55696e-06,1.76708e-09,-1.25016e-13,3690.89,-84.1798], Tmin=(1039.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.1053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'OOC1C=C=CCC=C1(917)',
    structure = SMILES('OOC1C=C=CCC=C1'),
    E0 = (209.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.725273,0.0576975,-3.97488e-06,-3.06438e-08,1.37682e-11,25303,25.3945], Tmin=(100,'K'), Tmax=(1077.19,'K')), NASAPolynomial(coeffs=[15.1953,0.0348732,-1.52321e-05,2.9607e-09,-2.13105e-13,20392.5,-53.8146], Tmin=(1077.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = 'O2(S)(918)',
    structure = SMILES('O=O'),
    E0 = (85.6848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121571,5.31618e-06,-4.89443e-09,1.45845e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69971e-07,1.51275e-10,-1.08782e-14,10302.3,6.16754], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]O[CH]C1C=CCC=C1(919)',
    structure = SMILES('[O]O[CH]C1C=CCC=C1'),
    E0 = (312.513,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,492.5,1135,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.60754,0.0712799,-5.03035e-05,1.83101e-08,-2.76422e-12,37711,28.0961], Tmin=(100,'K'), Tmax=(1507.09,'K')), NASAPolynomial(coeffs=[13.8887,0.0360309,-1.52211e-05,2.79163e-09,-1.90028e-13,33707.7,-41.4256], Tmin=(1507.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cyclohexadiene) + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = 'HO2(26)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263985,1.5223e-05,-1.71671e-08,6.26738e-12,322.677,4.84428], Tmin=(100,'K'), Tmax=(923.913,'K')), NASAPolynomial(coeffs=[4.15133,0.00191146,-4.11274e-07,6.34957e-11,-4.86385e-15,83.4208,3.09341], Tmin=(923.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[C]1=CCC=CC=C1(243)',
    structure = SMILES('[C]1=CCC=CC=C1'),
    E0 = (364.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.898556,0.0536037,-4.67814e-06,-4.06382e-08,2.19826e-11,43964.1,17.7097], Tmin=(100,'K'), Tmax=(939.78,'K')), NASAPolynomial(coeffs=[18.1007,0.0169825,-4.63849e-06,7.70412e-10,-5.5963e-14,39114.8,-72.8103], Tmin=(939.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(364.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'C1=CC2OOC2C=CC1(334)',
    structure = SMILES('C1=CC2OOC2C=CC1'),
    E0 = (-4.40773,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.166765,0.072764,-1.82297e-05,-3.20942e-08,1.75173e-11,-363.68,-1.76976], Tmin=(100,'K'), Tmax=(1032.16,'K')), NASAPolynomial(coeffs=[22.1628,0.0273212,-1.19075e-05,2.39399e-09,-1.78563e-13,-7162.12,-120.807], Tmin=(1032.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.40773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(12dioxetane) + ring(1,4-Cycloheptadiene) + ring(12dioxetane)"""),
)

species(
    label = '[O]C1[CH]C=CCC=C1(920)',
    structure = SMILES('[O]C1[CH]C=CCC=C1'),
    E0 = (218.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91248,-0.0262131,0.000325392,-4.51224e-07,1.8528e-10,26367.5,27.7054], Tmin=(100,'K'), Tmax=(909.363,'K')), NASAPolynomial(coeffs=[45.8741,-0.0321613,2.60451e-05,-5.12021e-09,3.28674e-13,10622.6,-222.817], Tmin=(909.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(CC(C)OJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[O]OC1C=CC=CC=C1(921)',
    structure = SMILES('[O]OC1C=CC=CC=C1'),
    E0 = (211.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,492.5,1135,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.187597,0.0780105,-4.80084e-05,-2.42041e-09,8.6881e-12,25562.4,25.9165], Tmin=(100,'K'), Tmax=(991.15,'K')), NASAPolynomial(coeffs=[21.5287,0.0208399,-7.59969e-06,1.41639e-09,-1.03026e-13,19760.9,-86.2089], Tmin=(991.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1C=[C]CCC=C1(194)',
    structure = SMILES('[O]OC1C=[C]CCC=C1'),
    E0 = (373.867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,492.5,1135,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04566,0.0123472,0.000197102,-3.02548e-07,1.28419e-10,45121.4,32.0142], Tmin=(100,'K'), Tmax=(908.439,'K')), NASAPolynomial(coeffs=[38.6692,-0.015786,1.64692e-05,-3.33849e-09,2.15415e-13,32610.8,-177.119], Tmin=(908.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1[C]=CCCC=C1(195)',
    structure = SMILES('[O]OC1[C]=CCCC=C1'),
    E0 = (373.867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04549,0.0123494,0.000197094,-3.02536e-07,1.28413e-10,45121.4,32.0148], Tmin=(100,'K'), Tmax=(908.445,'K')), NASAPolynomial(coeffs=[38.6697,-0.0157868,1.64697e-05,-3.33861e-09,2.15425e-13,32610.6,-177.122], Tmin=(908.445,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'OOC1[C]=C[CH]CC=C1(922)',
    structure = SMILES('OOC1[C]=C[CH]CC=C1'),
    E0 = (362.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.962452,0.0109263,0.000209985,-3.19267e-07,1.34648e-10,43817.3,30.4469], Tmin=(100,'K'), Tmax=(912.453,'K')), NASAPolynomial(coeffs=[40.5121,-0.0169769,1.6708e-05,-3.32357e-09,2.10605e-13,30544,-189.918], Tmin=(912.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[O]O[C]1C=CCCC=C1(193)',
    structure = SMILES('[O]OC1C=CCC[CH]C=1'),
    E0 = (276.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0677156,0.0819674,-7.05957e-05,3.16582e-08,-5.76107e-12,33408.9,13.241], Tmin=(100,'K'), Tmax=(1303.75,'K')), NASAPolynomial(coeffs=[16.2102,0.0324415,-1.36152e-05,2.52166e-09,-1.74048e-13,29199.7,-68.9189], Tmin=(1303.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(ROOJ) + radical(Allyl_S)"""),
)

species(
    label = '[O]OC1C=CC2C[CH]C21(923)',
    structure = SMILES('[O]OC1C=CC2C[CH]C21'),
    E0 = (363.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23643,0.0450367,2.48254e-05,-5.69997e-08,2.24538e-11,43814.5,27.1304], Tmin=(100,'K'), Tmax=(1039.8,'K')), NASAPolynomial(coeffs=[13.3556,0.0355518,-1.50639e-05,2.92288e-09,-2.1168e-13,39286.6,-41.4642], Tmin=(1039.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(ROOJ) + radical(cyclobutane)"""),
)

species(
    label = '[O]OC1[CH]C2CC2C=C1(924)',
    structure = SMILES('[O]OC1[CH]C2CC2C=C1'),
    E0 = (357.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.990295,0.0496753,1.80861e-05,-5.52172e-08,2.31664e-11,43069.1,27.5317], Tmin=(100,'K'), Tmax=(1017.51,'K')), NASAPolynomial(coeffs=[15.2383,0.032266,-1.31559e-05,2.53725e-09,-1.84596e-13,38171.3,-51.2736], Tmin=(1017.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + radical(CCJCOOH) + radical(ROOJ)"""),
)

species(
    label = '[CH]1[CH]C2OOC2C=CC1(342)',
    structure = SMILES('[CH]1[CH]C2OOC2C=CC1'),
    E0 = (277.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.263012,0.0826633,-6.51586e-05,2.57748e-08,-4.20238e-12,33504.6,0.449824], Tmin=(100,'K'), Tmax=(1406.01,'K')), NASAPolynomial(coeffs=[15.712,0.0387126,-1.82705e-05,3.54291e-09,-2.49418e-13,29160.2,-79.3474], Tmin=(1406.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(12dioxetane) + ring(Cycloheptene) + ring(12dioxetane) + radical(RCCJCC) + radical(CCJCOOH)"""),
)

species(
    label = '[CH]=CC(C=CC=C)O[O](925)',
    structure = SMILES('[CH]=CC(C=CC=C)O[O]'),
    E0 = (433.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3120,650,792.5,1650,492.5,1135,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.758606,'amu*angstrom^2'), symmetry=1, barrier=(17.4418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.757488,'amu*angstrom^2'), symmetry=1, barrier=(17.4161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.757185,'amu*angstrom^2'), symmetry=1, barrier=(17.4092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.758719,'amu*angstrom^2'), symmetry=1, barrier=(17.4444,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.232837,0.0862036,-7.98425e-05,3.77287e-08,-7.09648e-12,52271.4,35.6562], Tmin=(100,'K'), Tmax=(1281.98,'K')), NASAPolynomial(coeffs=[18.5486,0.0276016,-1.12738e-05,2.07066e-09,-1.42713e-13,47455.9,-59.6188], Tmin=(1281.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(ROOJ)"""),
)

species(
    label = 'OOC1C=CC=CC=C1(926)',
    structure = SMILES('OOC1C=CC=CC=C1'),
    E0 = (59.1752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.473805,0.0818503,-4.40869e-05,-1.08562e-08,1.21467e-11,7292.88,26.1856], Tmin=(100,'K'), Tmax=(995.761,'K')), NASAPolynomial(coeffs=[23.4119,0.021997,-8.29958e-06,1.58817e-09,-1.17423e-13,746.456,-97.9333], Tmin=(995.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.1752,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'C1=CC2C=CC(C1)O2(927)',
    structure = SMILES('C1=CC2C=CC(C1)O2'),
    E0 = (40.8973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14773,0.0150209,0.000108289,-1.437e-07,5.27488e-11,5008.23,19.6194], Tmin=(100,'K'), Tmax=(1000.6,'K')), NASAPolynomial(coeffs=[14.7235,0.0303699,-1.3094e-05,2.71657e-09,-2.09132e-13,-793.477,-57.4743], Tmin=(1000.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.8973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5)"""),
)

species(
    label = '[CH2]C1C=CC(C=C1)O[O](928)',
    structure = SMILES('[CH2]C1C=CC(C=C1)O[O]'),
    E0 = (311.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,492.5,1135,1000,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.168045,0.072834,-5.2547e-05,1.94651e-08,-2.91307e-12,37553,31.0085], Tmin=(100,'K'), Tmax=(1572.56,'K')), NASAPolynomial(coeffs=[17.1996,0.0295121,-1.12241e-05,1.94673e-09,-1.28075e-13,32196.3,-58.8691], Tmin=(1572.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cyclohexadiene) + radical(Isobutyl) + radical(ROOJ)"""),
)

species(
    label = '[C]1=CC=CCC=C1(244)',
    structure = SMILES('[C]1=CC=CCC=C1'),
    E0 = (364.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.898556,0.0536037,-4.67814e-06,-4.06382e-08,2.19826e-11,43964.1,17.7097], Tmin=(100,'K'), Tmax=(939.78,'K')), NASAPolynomial(coeffs=[18.1007,0.0169825,-4.63849e-06,7.70412e-10,-5.5963e-14,39114.8,-72.8103], Tmin=(939.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(364.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'C1=CC2C=CC(C1)OO2(684)',
    structure = SMILES('C1=CC2C=CC(C1)OO2'),
    E0 = (-3.51958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.301026,0.0794147,-4.45322e-05,1.93307e-09,3.22262e-12,-256.1,-1.721], Tmin=(100,'K'), Tmax=(1275.7,'K')), NASAPolynomial(coeffs=[24.1894,0.0290365,-1.63526e-05,3.4366e-09,-2.52625e-13,-8653.76,-134.26], Tmin=(1275.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.51958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(12dioxane) + ring(Cycloheptane) + ring(36dihydro12dioxin)"""),
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
    E0 = (211.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (295.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (449.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (588.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (416.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (536.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (381.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (520.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (528.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (356.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (423.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (661.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (330.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (448.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (292.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (474.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (324.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (456.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (382.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (489.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (304.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (357.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (260.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (489.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (289.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (236.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (251.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (340.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (472.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (389.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (219.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (461.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (433.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (531.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (523.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (365.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (694.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (368.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (363.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (289.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (474.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (228.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (291.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (470.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (389.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (218.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)', 'C7H8(11)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.35384,'m^3/(mol*s)'), n=2.24862, Ea=(54.3131,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-CdH;YJ] for rate rule [Cds-CdH_Cds-CdH;O2b]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond
Ea raised from 40.3 to 54.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['[CH]1CC=CC2[CH]C1OO2(897)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.2, Ea=(84.4482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R6;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', '[O]OC1C=CCC=CC=1(898)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(87.5179,'m^3/(mol*s)'), n=1.66467, Ea=(9.4546,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;HJ] for rate rule [Cds-CdOs_Cds;HJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', '[O]OC1C=C=CCC=C1(899)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]O[C]1C=CCC=CC1(900)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(11288.6,'s^-1'), n=2.7035, Ea=(123.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_OneDe/O;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC1C=CCC=[C]C1(901)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OO[C]1[CH]C=CCC=C1(902)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC1[C]=CCC=CC1(903)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC1C=CC[C]=CC1(904)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 200 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['OOC1[C]=CCC=C[CH]1(905)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OC1C=[C]CC=CC1(906)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC1C=CC=C[CH]C1(209)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_1H;Cs_H_out_H/NonDeC] for rate rule [R4H_SDS;C_rad_out_H/Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['OOC1[CH]C=CC[C]=C1(907)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['OOC1[CH][C]=CCC=C1(908)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['OOC1[CH]C=C[CH]C=C1(909)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(34428.8,'s^-1'), n=1.95377, Ea=(81.3736,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R6H_RSSMS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['OOC1[CH]C=[C]CC=C1(910)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.44313e+09,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;XH_out] for rate rule [R6HJ_2;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(2)', '[CH]1C=C[CH]CC=C1(125)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.17528e+06,'m^3/(mol*s)'), n=0.19382, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination
Ea raised from -25.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['[O]OC1C=CCC2[CH]C21(911)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_cs] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['[O]OC1[CH]C2CC=CC21(912)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(171.238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn0c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['[CH]1C=CCC2[CH]C1OO2(913)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.11e+10,'s^-1'), n=0.18, Ea=(278.654,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['[CH]1[CH]C2CC=CC1OO2(692)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.42658e+11,'s^-1'), n=0.0609759, Ea=(92.8718,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['[O]OC1C2[CH]CC=CC21(914)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(146.393,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['[CH]1C=CC2OOC2[CH]C1(343)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=CCC=CC=CO[O](915)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.25555e+08,'s^-1'), n=0.67, Ea=(34.1205,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;doublebond_intra_pri;radadd_intra_cdsingleH] + [R7_linear;doublebond_intra_pri_HNd;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_HNd_O;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['OOC1C=CCC=CC=1(916)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['OOC1C=C=CCC=C1(917)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['C7H8(11)', 'O2(S)(918)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(39.9935,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 40.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['C7H8O(2)(14)', 'O(20)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.27e+09,'s^-1'), n=1.06, Ea=(129.704,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2OO_S;C_sec_rad_intra;OOJ] for rate rule [R2OO_S;C_rad/H/OneDe_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]O[CH]C1C=CCC=C1(919)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.31121e+11,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['HO2(26)', '[C]1=CCC=CC=C1(243)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['C1=CC2OOC2C=CC1(334)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(20)', '[O]C1[CH]C=CCC=C1(920)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(19)', '[O]OC1C=CC=CC=C1(921)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(9.22566,'m^3/(mol*s)'), n=2.04274, Ea=(10.5229,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 101 used for Cds-CdH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CdH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]OC1C=[C]CCC=C1(194)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]OC1[C]=CCCC=C1(195)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 200 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['OOC1[C]=C[CH]CC=C1(922)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(154.108,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]O[C]1C=CCCC=C1(193)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R4H_SDS;C_rad_out_OneDe/O;Cs_H_out_H/NonDeC]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['[O]OC1C=CC2C[CH]C21(923)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(157.636,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn0c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['[O]OC1[CH]C2CC2C=C1(924)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(152.079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['[CH]1[CH]C2OOC2C=CC1(342)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(78.0467,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=CC(C=CC=C)O[O](925)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['OOC1C=CC=CC=C1(926)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_De] for rate rule [R6radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['O(20)', 'C1=CC2C=CC(C1)O2(927)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(3.63e+10,'s^-1'), n=0, Ea=(80.1772,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;C_rad/H/NonDeC_intra;OO] for rate rule [R4OO_SDS;C_rad/H/NonDeC_intra;OOJ]
Euclidian distance = 1.41421356237
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C1C=CC(C=C1)O[O](928)'],
    products = ['[O]OC1C=C[CH]CC=C1(192)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.31121e+11,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['HO2(26)', '[C]1=CC=CCC=C1(244)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[O]OC1C=C[CH]CC=C1(192)'],
    products = ['C1=CC2C=CC(C1)OO2(684)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '21',
    isomers = [
        '[O]OC1C=C[CH]CC=C1(192)',
    ],
    reactants = [
        ('O2(2)', 'C7H8(11)'),
        ('C7H8O(2)(14)', 'O(20)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '21',
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

