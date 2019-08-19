species(
    label = 'CC([O])C1[CH]C=CC=CC1(2281)',
    structure = SMILES('CC([O])C1C=C[CH]C=CC1'),
    E0 = (194.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.557254,0.00663242,0.000280611,-4.13634e-07,1.72925e-10,23573.3,34.7764], Tmin=(100,'K'), Tmax=(908.77,'K')), NASAPolynomial(coeffs=[46.4148,-0.015159,1.93868e-05,-3.98434e-09,2.55989e-13,7803.55,-222.98], Tmin=(908.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(CC(C)OJ)"""),
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
    label = 'CH3CHO(52)',
    structure = SMILES('CC=O'),
    E0 = (-178.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1305.64,1305.66,1305.67,3976.84],'cm^-1')),
        HinderedRotor(inertia=(0.136163,'amu*angstrom^2'), symmetry=1, barrier=(3.13064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-178.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'CC1OC2[CH]CC1[CH]C=C2(2595)',
    structure = SMILES('CC1OC2[CH]C=CC1C[CH]2'),
    E0 = (100.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56892,0.0721578,9.19282e-05,-2.07521e-07,9.38906e-11,12303.7,6.51982], Tmin=(100,'K'), Tmax=(937.878,'K')), NASAPolynomial(coeffs=[45.9696,-0.00248834,6.43139e-06,-1.11264e-09,4.99978e-14,-2247.44,-249.814], Tmin=(937.878,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(Cycloheptane) + ring(Oxane) + radical(CCJCO) + radical(C=CCJCO)"""),
)

species(
    label = 'CC1OC2[CH]C=CCC1[CH]2(2614)',
    structure = SMILES('CC1OC2[CH]C=CCC1[CH]2'),
    E0 = (92.4589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46279,0.102451,-5.94579e-05,1.86015e-09,6.63979e-12,11332.1,2.40772], Tmin=(100,'K'), Tmax=(1069.46,'K')), NASAPolynomial(coeffs=[23.9357,0.0395968,-1.63813e-05,3.10983e-09,-2.21639e-13,4061.45,-130.424], Tmin=(1069.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.4589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Tetrahydrofuran) + ring(Cycloheptene) + ring(Tetrahydrofuran) + radical(C=CCJCO) + radical(CCJCO)"""),
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
    label = 'CC(=O)C1C=C[CH]C=CC1(2881)',
    structure = SMILES('CC(=O)C1C=C[CH]C=CC1'),
    E0 = (23.7275,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,375,552.5,462.5,1710,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.777186,0.00972387,0.000246449,-3.64897e-07,1.52478e-10,3027.01,32.6388], Tmin=(100,'K'), Tmax=(909.749,'K')), NASAPolynomial(coeffs=[41.3728,-0.00949206,1.55175e-05,-3.22502e-09,2.05945e-13,-10950.5,-195.597], Tmin=(909.749,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.7275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC([O])C1C=CC=C=CC1(2882)',
    structure = SMILES('CC([O])C1C=CC=C=CC1'),
    E0 = (298.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0906556,0.0697548,8.59214e-06,-6.02648e-08,2.75603e-11,36019.3,28.2752], Tmin=(100,'K'), Tmax=(995.337,'K')), NASAPolynomial(coeffs=[19.529,0.0379903,-1.44915e-05,2.72039e-09,-1.96412e-13,29781.5,-78.0021], Tmin=(995.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C1C=C=CC=CC1(2883)',
    structure = SMILES('CC([O])C1C=C=CC=CC1'),
    E0 = (298.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0906556,0.0697548,8.59214e-06,-6.02648e-08,2.75603e-11,36019.3,28.2752], Tmin=(100,'K'), Tmax=(995.337,'K')), NASAPolynomial(coeffs=[19.529,0.0379903,-1.44915e-05,2.72039e-09,-1.96412e-13,29781.5,-78.0021], Tmin=(995.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CC(C)OJ)"""),
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
    label = 'CH3(34)',
    structure = SMILES('[CH3]'),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=CC1C=C[CH]C=CC1(1470)',
    structure = SMILES('O=CC1C=C[CH]C=CC1'),
    E0 = (78.5983,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.156,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49025,-0.00625321,0.000265446,-3.7964e-07,1.5747e-10,9601.15,28.4883], Tmin=(100,'K'), Tmax=(907.009,'K')), NASAPolynomial(coeffs=[40.3214,-0.0173795,1.90368e-05,-3.88592e-09,2.51515e-13,-4029.28,-191.368], Tmin=(907.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.5983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C)"""),
)

species(
    label = 'C[C](O)C1C=C[CH]C=CC1(2884)',
    structure = SMILES('C[C](O)C1C=C[CH]C=CC1'),
    E0 = (140.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.252766,0.0192084,0.000238206,-3.67117e-07,1.56045e-10,17116.3,35.4859], Tmin=(100,'K'), Tmax=(906.276,'K')), NASAPolynomial(coeffs=[45.1038,-0.0132896,1.81384e-05,-3.7814e-09,2.46008e-13,2191.91,-213.968], Tmin=(906.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(C2CsJOH)"""),
)

species(
    label = 'CC([O])C1C=CC[C]=CC1(2845)',
    structure = SMILES('CC([O])C1C=CC[C]=CC1'),
    E0 = (330.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.216137,0.0235074,0.00021852,-3.39923e-07,1.44164e-10,39947.8,36.8033], Tmin=(100,'K'), Tmax=(910.883,'K')), NASAPolynomial(coeffs=[42.9545,-0.00844422,1.46917e-05,-3.05325e-09,1.9415e-13,25701.5,-200.858], Tmin=(910.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'CC([O])C1C=[C]CC=CC1(2847)',
    structure = SMILES('CC([O])C1C=[C]CC=CC1'),
    E0 = (330.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.216137,0.0235074,0.00021852,-3.39923e-07,1.44164e-10,39947.8,36.8033], Tmin=(100,'K'), Tmax=(910.883,'K')), NASAPolynomial(coeffs=[42.9545,-0.00844422,1.46917e-05,-3.05325e-09,1.9415e-13,25701.5,-200.858], Tmin=(910.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'CC(O)[C]1C=C[CH]C=CC1(2885)',
    structure = SMILES('CC(O)C1=C[CH]C=C[CH]C1'),
    E0 = (98.6787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.540007,0.0059734,0.000283908,-4.18471e-07,1.74922e-10,12058.5,34.951], Tmin=(100,'K'), Tmax=(909.363,'K')), NASAPolynomial(coeffs=[47.2419,-0.0167548,2.00352e-05,-4.08904e-09,2.6209e-13,-3989.34,-227.456], Tmin=(909.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.6787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH2]C(O)C1C=C[CH]C=CC1(2886)',
    structure = SMILES('[CH2]C(O)C1C=C[CH]C=CC1'),
    E0 = (175.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.246624,0.015752,0.000255237,-3.89918e-07,1.6543e-10,21324.7,36.7176], Tmin=(100,'K'), Tmax=(906.144,'K')), NASAPolynomial(coeffs=[47.306,-0.0171673,2.03471e-05,-4.19961e-09,2.73377e-13,5619.21,-225.281], Tmin=(906.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(CJCO) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC([O])C1C=CCC=[C]C1(2842)',
    structure = SMILES('CC([O])C1C=CCC=[C]C1'),
    E0 = (330.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.216137,0.0235074,0.00021852,-3.39923e-07,1.44164e-10,39947.8,36.8033], Tmin=(100,'K'), Tmax=(910.883,'K')), NASAPolynomial(coeffs=[42.9545,-0.00844422,1.46917e-05,-3.05325e-09,1.9415e-13,25701.5,-200.858], Tmin=(910.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C1[C]=CCC=CC1(2844)',
    structure = SMILES('CC([O])C1[C]=CCC=CC1'),
    E0 = (330.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.216137,0.0235074,0.00021852,-3.39923e-07,1.44164e-10,39947.8,36.8033], Tmin=(100,'K'), Tmax=(910.883,'K')), NASAPolynomial(coeffs=[42.9545,-0.00844422,1.46917e-05,-3.05325e-09,1.9415e-13,25701.5,-200.858], Tmin=(910.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(O)C1C=C[CH][CH]C=C1(2832)',
    structure = SMILES('CC(O)C1[CH]C=C[CH]C=C1'),
    E0 = (105.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.598442,0.00532904,0.000283933,-4.17482e-07,1.74484e-10,12838.1,33.7141], Tmin=(100,'K'), Tmax=(908.288,'K')), NASAPolynomial(coeffs=[46.5503,-0.0158981,1.98439e-05,-4.07807e-09,2.62592e-13,-2981.25,-224.689], Tmin=(908.288,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC(O)C1[C]=C[CH]C=CC1(2887)',
    structure = SMILES('CC(O)C1[C]=C[CH]C=CC1'),
    E0 = (201.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3615,1277.5,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.372554,0.015091,0.000250405,-3.79981e-07,1.6064e-10,24475.7,36.1052], Tmin=(100,'K'), Tmax=(906.84,'K')), NASAPolynomial(coeffs=[45.2724,-0.0138354,1.85056e-05,-3.8437e-09,2.49333e-13,9378.28,-214.462], Tmin=(906.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC([O])[C]1C=CCC=CC1(2841)',
    structure = SMILES('CC([O])C1=C[CH]CC=CC1'),
    E0 = (227.315,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.383512,0.0143893,0.000252032,-3.78439e-07,1.58466e-10,27530.7,35.6495], Tmin=(100,'K'), Tmax=(913.302,'K')), NASAPolynomial(coeffs=[44.9297,-0.0113731,1.62267e-05,-3.29988e-09,2.07013e-13,12331.5,-213.884], Tmin=(913.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C1[CH]C=CCC=C1(2282)',
    structure = SMILES('CC([O])C1C=C[CH]CC=C1'),
    E0 = (232.001,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4625.32,'J/mol'), sigma=(7.7745,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=722.47 K, Pc=22.33 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202721,0.0254186,0.00021092,-3.29116e-07,1.39412e-10,28094.3,32.858], Tmin=(100,'K'), Tmax=(912.381,'K')), NASAPolynomial(coeffs=[41.7601,-0.00547265,1.2959e-05,-2.7103e-09,1.70539e-13,14213.6,-198.333], Tmin=(912.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(O)C1C=C[CH]C=[C]C1(2888)',
    structure = SMILES('CC(O)C1C=C[CH]C=[C]C1'),
    E0 = (201.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.372554,0.015091,0.000250405,-3.79981e-07,1.6064e-10,24475.7,36.1052], Tmin=(100,'K'), Tmax=(906.84,'K')), NASAPolynomial(coeffs=[45.2724,-0.0138354,1.85056e-05,-3.8437e-09,2.49333e-13,9378.28,-214.462], Tmin=(906.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(Cds_S)"""),
)

species(
    label = 'CC(O)C1C=[C][CH]C=CC1(2889)',
    structure = SMILES('CC(O)C1C=[C]C=C[CH]C1'),
    E0 = (184.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27093,0.102224,-7.19941e-05,1.74879e-08,1.51162e-12,22446.7,20.1988], Tmin=(100,'K'), Tmax=(1040.44,'K')), NASAPolynomial(coeffs=[21.4225,0.038696,-1.46005e-05,2.62252e-09,-1.81028e-13,16440.7,-96.353], Tmin=(1040.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Allyl_S) + radical(C=CJC=C)"""),
)

species(
    label = 'C[C]([O])C1C=CCC=CC1(2843)',
    structure = SMILES('C[C]([O])C1C=CCC=CC1'),
    E0 = (269.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,360,370,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0963564,0.0276249,0.000206319,-3.27052e-07,1.39565e-10,32588.4,36.184], Tmin=(100,'K'), Tmax=(910.362,'K')), NASAPolynomial(coeffs=[42.7847,-0.00789646,1.43233e-05,-2.99069e-09,1.90804e-13,18515.6,-200.358], Tmin=(910.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C([O])C1C=CCC=CC1(2846)',
    structure = SMILES('[CH2]C([O])C1C=CCC=CC1'),
    E0 = (304.295,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0903086,0.0241675,0.000223353,-3.49856e-07,1.48951e-10,36796.9,37.4153], Tmin=(100,'K'), Tmax=(909.976,'K')), NASAPolynomial(coeffs=[44.9862,-0.0117728,1.65312e-05,-3.40871e-09,2.18157e-13,21943.2,-211.666], Tmin=(909.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(O)C1C=C[CH][C]=CC1(2890)',
    structure = SMILES('CC(O)C1[CH]C=C[C]=CC1'),
    E0 = (183.886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28167,0.102621,-7.38492e-05,2.0066e-08,4.2066e-13,22318.2,20.2033], Tmin=(100,'K'), Tmax=(1053.71,'K')), NASAPolynomial(coeffs=[21.3145,0.0389943,-1.4807e-05,2.66148e-09,-1.83434e-13,16326.5,-95.8276], Tmin=(1053.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(Allyl_S)"""),
)

species(
    label = 'C[CH][O](2420)',
    structure = SMILES('C[CH][O]'),
    E0 = (157.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1642.51],'cm^-1')),
        HinderedRotor(inertia=(0.123965,'amu*angstrom^2'), symmetry=1, barrier=(2.85019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65562,0.0114444,2.34936e-06,-4.83164e-09,1.17966e-12,18963.9,10.3625], Tmin=(100,'K'), Tmax=(1718.65,'K')), NASAPolynomial(coeffs=[6.06294,0.0136322,-6.35953e-06,1.18407e-09,-7.90642e-14,16985.9,-5.90233], Tmin=(1718.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'CC([O])C1C=CC2[CH]C2C1(2891)',
    structure = SMILES('CC([O])C1C=CC2[CH]C2C1'),
    E0 = (334.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281761,0.0619958,2.50431e-05,-6.8901e-08,2.8208e-11,40349.5,31.74], Tmin=(100,'K'), Tmax=(1023.36,'K')), NASAPolynomial(coeffs=[16.3227,0.0454722,-1.84202e-05,3.50539e-09,-2.51933e-13,34648.5,-57.8324], Tmin=(1023.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.229,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + radical(cyclopropane) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C1CC=CC2[CH]C21(2892)',
    structure = SMILES('CC([O])C1CC=CC2[CH]C21'),
    E0 = (333.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.282366,0.0622629,2.36238e-05,-6.68539e-08,2.73266e-11,40220.5,31.7036], Tmin=(100,'K'), Tmax=(1027.6,'K')), NASAPolynomial(coeffs=[16.1426,0.0458897,-1.8694e-05,3.56001e-09,-2.55623e-13,34565.8,-56.8984], Tmin=(1027.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + radical(cyclopropane) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC1OC2C=C[CH][CH]C1C2(2549)',
    structure = SMILES('CC1OC2[CH]C=C[CH]C1C2'),
    E0 = (33.6695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46066,0.102427,-5.89257e-05,2.3438e-09,6.10986e-12,4261.14,-1.20114], Tmin=(100,'K'), Tmax=(1084.4,'K')), NASAPolynomial(coeffs=[23.6696,0.0409916,-1.71885e-05,3.26986e-09,-2.32679e-13,-3027.19,-132.951], Tmin=(1084.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.6695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Tetrahydrofuran) + ring(Cycloheptene) + ring(Tetrahydrofuran) + radical(Allyl_S) + radical(C=CCJCO)"""),
)

species(
    label = 'CC1OC2[CH][CH]C=CCC12(2893)',
    structure = SMILES('CC1OC2[CH]C=C[CH]CC12'),
    E0 = (33.6695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46066,0.102427,-5.89257e-05,2.3438e-09,6.10986e-12,4261.14,-1.20114], Tmin=(100,'K'), Tmax=(1084.4,'K')), NASAPolynomial(coeffs=[23.6696,0.0409916,-1.71885e-05,3.26986e-09,-2.32679e-13,-3027.19,-132.951], Tmin=(1084.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.6695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(Cycloheptene) + ring(Oxetane) + radical(Allyl_S) + radical(C=CCJCO)"""),
)

species(
    label = 'CC(=O)C1C=CCC=CC1(2853)',
    structure = SMILES('CC(=O)C1C=CCC=CC1'),
    E0 = (-77.9975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.462153,0.0220796,0.000208875,-3.20135e-07,1.34565e-10,-9201.55,34.0984], Tmin=(100,'K'), Tmax=(913.726,'K')), NASAPolynomial(coeffs=[39.125,-0.00231636,1.11225e-05,-2.34822e-09,1.45871e-13,-22314,-182.029], Tmin=(913.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.9975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene)"""),
)

species(
    label = 'CC(O)C1C=C=CC=CC1(2894)',
    structure = SMILES('CC(O)C1C=C=CC=CC1'),
    E0 = (67.7472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.25449,0.0737621,2.63357e-06,-5.51552e-08,2.59208e-11,8318.8,29.055], Tmin=(100,'K'), Tmax=(996.028,'K')), NASAPolynomial(coeffs=[19.6004,0.0397718,-1.50704e-05,2.80618e-09,-2.01252e-13,2094.43,-78.0454], Tmin=(996.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.7472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = 'CC(O)C1C=CC=C=CC1(2895)',
    structure = SMILES('CC(O)C1C=CC=C=CC1'),
    E0 = (67.7472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.25449,0.0737621,2.63357e-06,-5.51552e-08,2.59208e-11,8318.8,29.055], Tmin=(100,'K'), Tmax=(996.028,'K')), NASAPolynomial(coeffs=[19.6004,0.0397718,-1.50704e-05,2.80618e-09,-2.01252e-13,2094.43,-78.0454], Tmin=(996.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.7472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = 'CH2(S)(40)',
    structure = SMILES('[CH2]'),
    E0 = (418.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1358.21,2621.43,3089.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19331,-0.00233105,8.15676e-06,-6.62986e-09,1.93233e-12,50366.2,-0.746734], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.13502,0.00289594,-8.16668e-07,1.13573e-10,-6.36263e-15,50504.1,4.06031], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(418.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]CC1[CH]C=CC=CC1(1355)',
    structure = SMILES('[O]CC1C=C[CH]C=CC1'),
    E0 = (226.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4471.47,'J/mol'), sigma=(7.43333,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=698.43 K, Pc=24.7 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46299,-0.00542411,0.000267902,-3.80522e-07,1.5702e-10,27354.1,30.6786], Tmin=(100,'K'), Tmax=(908.68,'K')), NASAPolynomial(coeffs=[39.2648,-0.0125822,1.6845e-05,-3.46892e-09,2.22522e-13,13909.7,-184.249], Tmin=(908.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(CCOJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC1OC2C=CCC1C=C2(2287)',
    structure = SMILES('CC1OC2C=CCC1C=C2'),
    E0 = (-102.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46215,0.0723789,7.84919e-05,-1.84997e-07,8.30475e-11,-12064.3,6.75717], Tmin=(100,'K'), Tmax=(952.792,'K')), NASAPolynomial(coeffs=[44.1351,0.000651887,2.96938e-06,-2.99816e-10,-1.11233e-14,-26186.4,-239.53], Tmin=(952.792,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(Cycloheptane) + ring(36dihydro2hpyran)"""),
)

species(
    label = 'C=CC1[CH]C1C=CC(C)[O](2896)',
    structure = SMILES('C=CC1[CH]C1C=CC(C)[O]'),
    E0 = (411.368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1380,1390,370,380,2900,435,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.187375,0.0709554,9.98704e-06,-6.35361e-08,2.90033e-11,49645.6,40.1337], Tmin=(100,'K'), Tmax=(994.632,'K')), NASAPolynomial(coeffs=[20.3624,0.0377057,-1.43584e-05,2.70932e-09,-1.96656e-13,43114.6,-71.1782], Tmin=(994.632,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(cyclopropane)"""),
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
    label = 'C[CH]C1C=C[CH]C=CC1(623)',
    structure = SMILES('C[CH]C1C=C[CH]C=CC1'),
    E0 = (331.083,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32356,-0.00827351,0.000298118,-4.19497e-07,1.72249e-10,39979.2,32.4153], Tmin=(100,'K'), Tmax=(910.239,'K')), NASAPolynomial(coeffs=[41.9376,-0.0124003,1.76048e-05,-3.61553e-09,2.29943e-13,25362.8,-199.379], Tmin=(910.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cs_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC([O])C1=CC=CC=CC1(2897)',
    structure = SMILES('CC([O])C1=CC=CC=CC1'),
    E0 = (161.357,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06884,0.0897393,-2.63698e-05,-4.10749e-08,2.4982e-11,19608.8,31.9183], Tmin=(100,'K'), Tmax=(969.224,'K')), NASAPolynomial(coeffs=[26.4534,0.0274109,-9.23432e-06,1.70191e-09,-1.25715e-13,11866.3,-112.42], Tmin=(969.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])[C]1CC=CC=CC1(2898)',
    structure = SMILES('CC([O])[C]1CC=CC=CC1'),
    E0 = (226.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.915863,0.0900587,-3.4404e-05,-2.16467e-08,1.49716e-11,27460,20.4035], Tmin=(100,'K'), Tmax=(1015.92,'K')), NASAPolynomial(coeffs=[21.9472,0.0384509,-1.49189e-05,2.78347e-09,-1.98535e-13,20832.4,-100.014], Tmin=(1015.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CCJ(C)CO) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C1C[C]=CC=CC1(2899)',
    structure = SMILES('CC([O])C1C[C]=CC=CC1'),
    E0 = (311.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29731,0.102153,-7.31785e-05,2.02421e-08,-8.22361e-14,37725.5,22.0544], Tmin=(100,'K'), Tmax=(1092.24,'K')), NASAPolynomial(coeffs=[21.9319,0.0384956,-1.51613e-05,2.77869e-09,-1.93219e-13,31373.9,-97.9094], Tmin=(1092.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C[C]([O])C1CC=CC=CC1(2900)',
    structure = SMILES('C[C]([O])C1CC=CC=CC1'),
    E0 = (250.767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,360,370,350,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58116,0.108039,-9.07531e-05,3.88939e-08,-6.63478e-12,30373.4,22.0345], Tmin=(100,'K'), Tmax=(1408.92,'K')), NASAPolynomial(coeffs=[23.9226,0.0356324,-1.36663e-05,2.4184e-09,-1.62548e-13,23186.9,-109.75], Tmin=(1408.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'CC([O])C1CC=[C]C=CC1(2901)',
    structure = SMILES('CC([O])C1CC=[C]C=CC1'),
    E0 = (273.134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34962,0.104216,-7.806e-05,2.48444e-08,-1.42373e-12,33054.5,21.3631], Tmin=(100,'K'), Tmax=(1079.01,'K')), NASAPolynomial(coeffs=[21.4331,0.0393168,-1.50301e-05,2.70099e-09,-1.85606e-13,26999.4,-95.559], Tmin=(1079.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([O])C1CC=CC=CC1(2902)',
    structure = SMILES('[CH2]C([O])C1CC=CC=CC1'),
    E0 = (285.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3644,0.102186,-6.64959e-05,8.49461e-09,5.18769e-12,34572,22.4518], Tmin=(100,'K'), Tmax=(1030.96,'K')), NASAPolynomial(coeffs=[23.3338,0.0361978,-1.39e-05,2.55707e-09,-1.80142e-13,27893.7,-105.146], Tmin=(1030.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'CC(O)C1[CH]C=[C]C=CC1(2903)',
    structure = SMILES('CC(O)C1[CH]C=[C]C=CC1'),
    E0 = (183.886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3615,1277.5,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28167,0.102621,-7.38492e-05,2.0066e-08,4.2066e-13,22318.2,20.2033], Tmin=(100,'K'), Tmax=(1053.71,'K')), NASAPolynomial(coeffs=[21.3145,0.0389943,-1.4807e-05,2.66148e-09,-1.83434e-13,16326.5,-95.8276], Tmin=(1053.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(Allyl_S)"""),
)

species(
    label = 'CC([O])C1C[CH]C2C=CC21(2904)',
    structure = SMILES('CC([O])C1C[CH]C2C=CC21'),
    E0 = (349.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.438106,0.0573418,3.74712e-05,-8.09184e-08,3.2144e-11,42188.5,30.9532], Tmin=(100,'K'), Tmax=(1021.22,'K')), NASAPolynomial(coeffs=[16.1751,0.0456217,-1.86373e-05,3.57629e-09,-2.58632e-13,36371.3,-58.0438], Tmin=(1021.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(Cs_S)"""),
)

species(
    label = 'CC([O])C1CC2[CH]C=CC21(2905)',
    structure = SMILES('CC([O])C1CC2[CH]C=CC21'),
    E0 = (253.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.437953,0.0531391,5.99662e-05,-1.10149e-07,4.3818e-11,30628.3,28.1267], Tmin=(100,'K'), Tmax=(995.145,'K')), NASAPolynomial(coeffs=[18.3303,0.0426159,-1.67142e-05,3.21659e-09,-2.36216e-13,24027.2,-73.3807], Tmin=(995.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(cyclopentene-allyl) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC1OC2[CH][CH]C1CC=C2(2615)',
    structure = SMILES('CC1OC2[CH][CH]C1CC=C2'),
    E0 = (176.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11406,0.0640996,9.89211e-05,-2.05512e-07,9.08954e-11,21452.5,11.1184], Tmin=(100,'K'), Tmax=(942.925,'K')), NASAPolynomial(coeffs=[42.6611,0.00101358,4.22498e-06,-6.5388e-10,1.73857e-14,7746.33,-226.403], Tmin=(942.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(Cycloheptane) + ring(Oxane) + radical(CCJCO) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C=CC=CC=CC(C)[O](2906)',
    structure = SMILES('[CH2]C=CC=CC=CC(C)[O]'),
    E0 = (226.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3001,3007,3013,3019,3025,975,980,985,990,995,1000,1300,1315,1330,1345,1360,1375,400,420,440,460,480,500,1630,1640,1650,1660,1670,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.87887,0.0869412,-2.07027e-05,-4.27473e-08,2.46748e-11,27448.4,39.2989], Tmin=(100,'K'), Tmax=(967.85,'K')), NASAPolynomial(coeffs=[23.8767,0.0327091,-1.11671e-05,2.01116e-09,-1.44438e-13,20404.6,-90.9567], Tmin=(967.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(O)C1=CC=CC=CC1(2907)',
    structure = SMILES('CC(O)C1=CC=CC=CC1'),
    E0 = (-69.0038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.23026,0.0937175,-3.22246e-05,-3.61022e-08,2.34016e-11,-8091.76,32.6895], Tmin=(100,'K'), Tmax=(968.867,'K')), NASAPolynomial(coeffs=[26.5151,0.029209,-9.82283e-06,1.78996e-09,-1.30743e-13,-15816.6,-112.408], Tmin=(968.867,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.0038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'CC(=O)C1CC=CC=CC1(2908)',
    structure = SMILES('CC(=O)C1CC=CC=CC1'),
    E0 = (-94.3605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07071,0.0975742,-6.92227e-05,2.19806e-08,-1.96585e-12,-11154.8,20.1306], Tmin=(100,'K'), Tmax=(1194.52,'K')), NASAPolynomial(coeffs=[20.6195,0.0400926,-1.60666e-05,2.93236e-09,-2.01573e-13,-17417.6,-92.892], Tmin=(1194.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.3605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = 'CC([O])C1[CH]CC=CC=C1(2280)',
    structure = SMILES('CC([O])C1[CH]CC=CC=C1'),
    E0 = (269.751,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4625.32,'J/mol'), sigma=(7.7745,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=722.47 K, Pc=22.33 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04707,0.0957413,-5.7185e-05,5.53456e-09,4.56276e-12,32638.2,23.2735], Tmin=(100,'K'), Tmax=(1073.77,'K')), NASAPolynomial(coeffs=[21.1771,0.0395547,-1.58582e-05,2.94928e-09,-2.07305e-13,26331.9,-92.6686], Tmin=(1073.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(Cs_S)"""),
)

species(
    label = 'CC([O])[CH]C1C=CC=CC1(2878)',
    structure = SMILES('CC([O])[CH]C1C=CC=CC1'),
    E0 = (266.852,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0378215,0.0656404,2.83352e-05,-8.36346e-08,3.64056e-11,32261,34.9142], Tmin=(100,'K'), Tmax=(982.445,'K')), NASAPolynomial(coeffs=[20.4535,0.0377711,-1.39438e-05,2.61886e-09,-1.91103e-13,25553.3,-77.2278], Tmin=(982.445,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = 'CC1OC2C=CC=CCC12(2288)',
    structure = SMILES('CC1OC2C=CC=CCC12'),
    E0 = (-128.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99782,0.162149,-0.000302372,3.14488e-07,-1.22904e-10,-15311,-13.7435], Tmin=(100,'K'), Tmax=(832.435,'K')), NASAPolynomial(coeffs=[-2.29853,0.0927313,-4.95955e-05,9.7865e-09,-6.83055e-13,-12805.8,2.39916], Tmin=(832.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-128.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(1,3-Cycloheptadiene) + ring(Oxetane)"""),
)

species(
    label = 'CC([O])C1C=CC=CC=C1(2819)',
    structure = SMILES('CC([O])C1C=CC=CC=C1'),
    E0 = (166.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29987,0.101346,-6.94014e-05,1.05337e-08,5.06903e-12,20174.6,29.3078], Tmin=(100,'K'), Tmax=(1010.44,'K')), NASAPolynomial(coeffs=[23.6217,0.0327445,-1.21783e-05,2.21553e-09,-1.55921e-13,13604,-98.7753], Tmin=(1010.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C1C=CC=[C]CC1(2862)',
    structure = SMILES('CC([O])C1C=CC=[C]CC1'),
    E0 = (313.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28476,0.101742,-7.13092e-05,1.76974e-08,9.71955e-13,37853.9,22.0429], Tmin=(100,'K'), Tmax=(1075.6,'K')), NASAPolynomial(coeffs=[21.9879,0.0382793,-1.49996e-05,2.74986e-09,-1.91626e-13,31512.1,-98.1378], Tmin=(1075.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'CC([O])[C]1C=CC=CCC1(2857)',
    structure = SMILES('CC([O])C1=C[CH]C=CCC1'),
    E0 = (187.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.498708,0.00727809,0.000280581,-4.14618e-07,1.73361e-10,22793.7,36.0137], Tmin=(100,'K'), Tmax=(909.852,'K')), NASAPolynomial(coeffs=[47.107,-0.0160167,1.95787e-05,-3.99547e-09,2.55499e-13,6795.21,-225.751], Tmin=(909.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(CC(C)OJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC([O])C1C=C[C]=CCC1(2866)',
    structure = SMILES('CC([O])C1C=C[C]=CCC1'),
    E0 = (274.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33379,0.103764,-7.60391e-05,2.20912e-08,-2.76144e-13,33182.8,21.34], Tmin=(100,'K'), Tmax=(1060.68,'K')), NASAPolynomial(coeffs=[21.4888,0.0391029,-1.48704e-05,2.67277e-09,-1.84071e-13,27137.1,-95.7868], Tmin=(1060.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'CC(=O)C1C=C[CH][CH]CC1(2796)',
    structure = SMILES('CC(=O)C1C=C[CH][CH]CC1'),
    E0 = (144.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,375,552.5,462.5,1710,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.52199,0.028412,0.000175159,-2.73397e-07,1.15004e-10,17517.4,34.4582], Tmin=(100,'K'), Tmax=(915.565,'K')), NASAPolynomial(coeffs=[34.0182,0.00685558,6.03697e-06,-1.38925e-09,8.30154e-14,6153.69,-152.749], Tmin=(915.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(RCCJCC) + radical(Allyl_S)"""),
)

species(
    label = 'CC([O])C1[C]=CC=CCC1(2859)',
    structure = SMILES('CC([O])C1[C]=CC=CCC1'),
    E0 = (313.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28476,0.101742,-7.13092e-05,1.76974e-08,9.71955e-13,37853.9,22.0429], Tmin=(100,'K'), Tmax=(1075.6,'K')), NASAPolynomial(coeffs=[21.9879,0.0382793,-1.49996e-05,2.74986e-09,-1.91626e-13,31512.1,-98.1378], Tmin=(1075.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'CC([O])C1C=[C]C=CCC1(2865)',
    structure = SMILES('CC([O])C1C=[C]C=CCC1'),
    E0 = (274.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33379,0.103764,-7.60391e-05,2.20912e-08,-2.76144e-13,33182.8,21.34], Tmin=(100,'K'), Tmax=(1060.68,'K')), NASAPolynomial(coeffs=[21.4888,0.0391029,-1.48704e-05,2.67277e-09,-1.84071e-13,27137.1,-95.7868], Tmin=(1060.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([O])C1C=CC=CCC1(2864)',
    structure = SMILES('[CH2]C([O])C1C=CC=CCC1'),
    E0 = (286.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35978,0.101857,-6.48629e-05,6.17465e-09,6.18283e-12,34700.7,22.4695], Tmin=(100,'K'), Tmax=(1023.29,'K')), NASAPolynomial(coeffs=[23.4889,0.0358223,-1.36502e-05,2.50808e-09,-1.76916e-13,27987.1,-105.939], Tmin=(1023.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'CC([O])C1[CH]C2C=CC2C1(2909)',
    structure = SMILES('CC([O])C1[CH]C2C=CC2C1'),
    E0 = (349.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.438082,0.0573421,3.74702e-05,-8.09171e-08,3.21434e-11,42188.5,30.9533], Tmin=(100,'K'), Tmax=(1021.23,'K')), NASAPolynomial(coeffs=[16.1752,0.0456216,-1.86372e-05,3.57627e-09,-2.5863e-13,36371.2,-58.0445], Tmin=(1021.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + radical(Cs_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(=O)C1C=CC=CCC1(2798)',
    structure = SMILES('CC(=O)C1C=CC=CCC1'),
    E0 = (-95.4953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.959691,0.0992419,-7.6406e-05,3.02067e-08,-4.84251e-12,-11299,19.062], Tmin=(100,'K'), Tmax=(1465.75,'K')), NASAPolynomial(coeffs=[20.4474,0.0408224,-1.66212e-05,3.01472e-09,-2.04596e-13,-17574.4,-92.3999], Tmin=(1465.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.4953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = 'CC(O)C1C=CC=CC=C1(2838)',
    structure = SMILES('CC(O)C1C=CC=CC=C1'),
    E0 = (-64.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46677,0.105389,-7.54844e-05,1.58014e-08,3.36415e-12,-7525.72,30.0985], Tmin=(100,'K'), Tmax=(1013.04,'K')), NASAPolynomial(coeffs=[23.7087,0.0344998,-1.27422e-05,2.2978e-09,-1.60471e-13,-14089.7,-98.9067], Tmin=(1013.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'CC([O])C1=CC=CC[CH]C1(2910)',
    structure = SMILES('CC([O])C1=CC=CC[CH]C1'),
    E0 = (263.162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30379,0.100811,-7.76698e-05,3.03909e-08,-4.7559e-12,31855.3,25.2524], Tmin=(100,'K'), Tmax=(1520.87,'K')), NASAPolynomial(coeffs=[23.2659,0.0361905,-1.3935e-05,2.45274e-09,-1.63388e-13,24381.9,-103.584], Tmin=(1520.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(RCCJCC) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1C=CC=CC1C(C)[O](2879)',
    structure = SMILES('[CH2]C1C=CC=CC1C(C)[O]'),
    E0 = (264.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.156544,0.069205,2.06622e-05,-7.87182e-08,3.58921e-11,31928.8,33.6865], Tmin=(100,'K'), Tmax=(964.043,'K')), NASAPolynomial(coeffs=[20.7266,0.0365216,-1.24504e-05,2.24552e-09,-1.61535e-13,25394.7,-79.3039], Tmin=(964.043,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Isobutyl) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC1OC2C=CC=CC1C2(2580)',
    structure = SMILES('CC1OC2C=CC=CC1C2'),
    E0 = (-127.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00096,0.161901,-0.000300967,3.12362e-07,-1.21938e-10,-15182,-13.6975], Tmin=(100,'K'), Tmax=(832.208,'K')), NASAPolynomial(coeffs=[-2.05641,0.0922118,-4.92642e-05,9.71853e-09,-6.78273e-13,-12750.2,1.11438], Tmin=(832.208,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Tetrahydrofuran) + ring(1,3-Cycloheptadiene) + ring(Tetrahydrofuran)"""),
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
    E0 = (194.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (215.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (316.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (306.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (526.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (526.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (194.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (246.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (360.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (487.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (487.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (277.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (333.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (479.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (479.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (297.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (303.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (342.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (423.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (238.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (229.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (331.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (383.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (246.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (496.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (439.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (439.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (237.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (243.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (233.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (212.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (212.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (645.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (201.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (558.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (574.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (378.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (324.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (352.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (512.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (375.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (420.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (366.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (410.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (421.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (373.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (220.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (253.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (272.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (283.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (315.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (426.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (202.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (377.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (513.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (250.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (421.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (300.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (357.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (535.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (317.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (421.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (233.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (212.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (379.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (423.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (201.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['C7H8(11)', 'CH3CHO(52)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC1OC2[CH]CC1[CH]C=C2(2595)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(414000,'s^-1'), n=1.14, Ea=(20.92,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R8;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC1OC2[CH]C=CCC1[CH]2(2614)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.291e+11,'s^-1'), n=0.234, Ea=(121.888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra;radadd_intra_O] for rate rule [R9;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'CC(=O)C1C=C[CH]C=CC1(2881)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0366254,'m^3/(mol*s)'), n=1.743, Ea=(71.4418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsCs_O;YJ] for rate rule [CO-CsCs_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(19)', 'CC([O])C1C=CC=C=CC1(2882)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(19)', 'CC([O])C1C=C=CC=CC1(2883)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', 'CH3CHO(52)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.000108643,'m^3/(mol*s)'), n=3.00879, Ea=(40.0106,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-OneDe] for rate rule [CO-CsH_O;CsJ-CdCsH]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 36.4 to 40.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3(34)', 'O=CC1C=C[CH]C=CC1(1470)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.61258,'m^3/(mol*s)'), n=1.485, Ea=(32.0285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-HHH] for rate rule [CO-CsH_O;CsJ-HHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['C[C](O)C1C=C[CH]C=CC1(2884)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC([O])C1C=CC[C]=CC1(2845)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CC([O])C1C=[C]CC=CC1(2847)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC(O)[C]1C=C[CH]C=CC1(2885)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7e-08,'s^-1'), n=6.3, Ea=(82.634,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['[CH2]C(O)C1C=C[CH]C=CC1(2886)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CC([O])C1C=CCC=[C]C1(2842)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CC([O])C1[C]=CCC=CC1(2844)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC(O)C1C=C[CH][CH]C=C1(2832)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(657459,'s^-1'), n=1.87237, Ea=(102.584,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC(O)C1[C]=C[CH]C=CC1(2887)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC([O])[C]1C=CCC=CC1(2841)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.30755e+06,'s^-1'), n=1.805, Ea=(115.478,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_SDS;C_rad_out_Cs2;Cs_H_out_H/Cd]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['CC([O])C1[CH]C=CCC=C1(2282)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(285077,'s^-1'), n=1.88833, Ea=(191.836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] + [R4H_SDS;C_rad_out_single;Cs_H_out_H/Cd] + [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CC(O)C1C=C[CH]C=[C]C1(2888)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(116983,'s^-1'), n=1.68042, Ea=(36.5782,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;Y_rad_out;XH_out] for rate rule [R5H_CCC;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC(O)C1C=[C][CH]C=CC1(2889)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C[C]([O])C1C=CCC=CC1(2843)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(60205.5,'s^-1'), n=1.86417, Ea=(61.9987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5H_SSMS;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([O])C1C=CCC=CC1(2846)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(42204,'s^-1'), n=1.945, Ea=(78.973,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;C_rad_out_2H;Cs_H_out_H/Cd] for rate rule [R6H_RSSMS;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CC(O)C1C=C[CH][C]=CC1(2890)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(13085.7,'s^-1'), n=1.78967, Ea=(62.3074,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;XH_out] for rate rule [R6H_DSSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', 'C[CH][O](2420)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.42928e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [C_rad/H/CdCs;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC([O])C1C=CC2[CH]C2C1(2891)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_csHDe] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC([O])C1CC=CC2[CH]C21(2892)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_csHDe] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC1OC2C=C[CH][CH]C1C2(2549)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.95208e+10,'s^-1'), n=0.228823, Ea=(43.0818,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_beta_long;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC1OC2[CH][CH]C=CCC12(2893)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC(=O)C1C=CCC=CC1(2853)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC(O)C1C=C=CC=CC1(2894)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC(O)C1C=CC=C=CC1(2895)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_De] for rate rule [R6radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2(S)(40)', '[O]CC1[CH]C=CC=CC1(1355)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC1OC2C=CCC1C=C2(2287)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=CC1[CH]C1C=CC(C)[O](2896)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(20)', 'C[CH]C1C=C[CH]C=CC1(623)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(19)', 'CC([O])C1=CC=CC=CC1(2897)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(7.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(4.93712,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2570 used for Cds-CsCs_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C7H8(11)', 'C[CH][O](2420)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.11106,'m^3/(mol*s)'), n=1.92886, Ea=(1.83261,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds-OneDeH;CJ] + [Cds-CsH_Cds-(Cd-Cd-Cd)H;YJ] for rate rule [Cds-CsH_Cds-(Cd-Cd-Cd)H;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CC([O])[C]1CC=CC=CC1(2898)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(5.88e+08,'s^-1'), n=1.27, Ea=(125.938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_Cs2;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_Cs2;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CC([O])C1C[C]=CC=CC1(2899)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C[C]([O])C1CC=CC=CC1(2900)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.59423e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CC([O])C1CC=[C]C=CC1(2901)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C([O])C1CC=CC=CC1(2902)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(775195,'s^-1'), n=1.84, Ea=(81.1696,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_single;Cs_H_out_H/(Cd-Cd-Cd)] + [R4H_SSS;C_rad_out_2H;Cs_H_out_H/Cd] for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['CC(O)C1[CH]C=[C]C=CC1(2903)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.61756e+09,'s^-1'), n=1.29078, Ea=(226.405,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC([O])C1C[CH]C2C=CC21(2904)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(5.65845e+18,'s^-1'), n=-1.57151, Ea=(227.049,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 131 used for Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd
Exact match found for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC([O])C1CC2[CH]C=CC21(2905)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(6.33004e+14,'s^-1'), n=-0.792922, Ea=(179.405,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 112 used for Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC1OC2[CH][CH]C1CC=C2(2615)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(3.125e+09,'s^-1'), n=0.76, Ea=(25.9408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri_HCd;radadd_intra] for rate rule [Rn2c7_gamma_short;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C=CC=CC=CC(C)[O](2906)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H] for rate rule [R7_SDSD_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction49',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC(O)C1=CC=CC=CC1(2907)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction50',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC(=O)C1CC=CC=CC1(2908)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['CC([O])C1[CH]CC=CC=C1(2280)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(4.25e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;CdHC_2] for rate rule [1_3_pentadiene;CH(CJ)_1;CdHC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction52',
    reactants = ['CC([O])[CH]C1C=CC=CC1(2878)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction53',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC1OC2C=CC=CCC12(2288)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction54',
    reactants = ['H(19)', 'CC([O])C1C=CC=CC=C1(2819)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(187.105,'m^3/(mol*s)'), n=1.3603, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 20 used for Cds-CsH_Cds-(Cd-Cd-Cd)H;HJ
Exact match found for rate rule [Cds-CsH_Cds-(Cd-Cd-Cd)H;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -1.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction55',
    reactants = ['CC([O])C1C=CC=[C]CC1(2862)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC([O])[C]1C=CC=CCC1(2857)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(2.76e-23,'s^-1'), n=10.17, Ea=(56.5677,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_1H;Cs_H_out_OneDe] for rate rule [R3H_SS_Cs;C_rad_out_H/(Cd-Cd-Cd);Cs_H_out_(CdCdCd)]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['CC([O])C1C=C[C]=CCC1(2866)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC(=O)C1C=C[CH][CH]CC1(2796)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.78e+06,'s^-1'), n=1.75, Ea=(105.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/OneDe;XH_out] for rate rule [R4H_SSS;C_rad_out_H/(Cd-Cd-Cd);XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['CC([O])C1[C]=CC=CCC1(2859)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['CC([O])C1C=[C]C=CCC1(2865)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(2.10713e+09,'s^-1'), n=0.59575, Ea=(261.003,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SDS;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] + [R4H_SDS;Cd_rad_out_Cd;Cs_H_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2]C([O])C1C=CC=CCC1(2864)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(15400,'s^-1'), n=1.87, Ea=(30.5432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_H/Cd] for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC([O])C1[CH]C2C=CC2C1(2909)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(5.65845e+18,'s^-1'), n=-1.57151, Ea=(227.049,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 131 used for Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd
Exact match found for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction63',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC(=O)C1C=CC=CCC1(2798)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction64',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC(O)C1C=CC=CC=C1(2838)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction65',
    reactants = ['CC([O])C1=CC=CC[CH]C1(2910)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(3.94565e+09,'s^-1'), n=0.909333, Ea=(116.555,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH(CJ)_1;unsaturated_end] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH(CJ)_1;Cd(C)C_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2]C1C=CC=CC1C(C)[O](2879)'],
    products = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction67',
    reactants = ['CC([O])C1[CH]C=CC=CC1(2281)'],
    products = ['CC1OC2C=CC=CC1C2(2580)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R5_SSSS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '535',
    isomers = [
        'CC([O])C1[CH]C=CC=CC1(2281)',
    ],
    reactants = [
        ('C7H8(11)', 'CH3CHO(52)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '535',
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

