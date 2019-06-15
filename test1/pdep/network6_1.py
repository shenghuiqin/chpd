species(
    label = '[O]OC1[CH]C=CCCC1(62)',
    structure = SMILES('[O]OC1[CH]C=CCCC1'),
    E0 = (120.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.546607,0.0561385,2.66407e-05,-7.24519e-08,3.08589e-11,14629,25.0356], Tmin=(100,'K'), Tmax=(991.402,'K')), NASAPolynomial(coeffs=[17.4754,0.0351587,-1.3216e-05,2.49659e-09,-1.81817e-13,8946.73,-68.2187], Tmin=(991.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CCJCO) + radical(ROOJ)"""),
)

species(
    label = 'C7H10(1)',
    structure = SMILES('C1C=CCCCC=1'),
    E0 = (68.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3862.82,'J/mol'), sigma=(6.55596,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=603.36 K, Pc=31.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98716,0.0214517,9.2224e-05,-1.30105e-07,4.97519e-11,8284.97,18.5915], Tmin=(100,'K'), Tmax=(968.934,'K')), NASAPolynomial(coeffs=[13.7239,0.0305102,-1.08313e-05,2.05835e-09,-1.53811e-13,3310.9,-51.592], Tmin=(968.934,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'O2(9)',
    structure = SMILES('[O][O]'),
    E0 = (-8.64439,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1498.97],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78246,-0.00299673,9.8473e-06,-9.6813e-09,3.24373e-12,-1063.94,3.65768], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.66096,0.000656366,-1.4115e-07,2.05798e-11,-1.29913e-15,-1215.98,3.41536], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64439,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]1CCCC2[CH]C1OO2(662)',
    structure = SMILES('[CH]1CCCC2[CH]C1OO2'),
    E0 = (218.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36859,0.0304845,9.74286e-05,-1.44278e-07,5.57298e-11,26455.9,28.152], Tmin=(100,'K'), Tmax=(977.384,'K')), NASAPolynomial(coeffs=[17.5046,0.0331589,-1.21287e-05,2.37926e-09,-1.8111e-13,20019.7,-66.1152], Tmin=(977.384,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(CCJCOOH) + radical(CCJCOOH)"""),
)

species(
    label = 'H(10)',
    structure = SMILES('[H]'),
    E0 = (211.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(211.8,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]OC1=CC=CCCC1(663)',
    structure = SMILES('[O]OC1=CC=CCCC1'),
    E0 = (105.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.800388,0.055667,6.0236e-06,-4.34354e-08,1.93827e-11,12781,26.948], Tmin=(100,'K'), Tmax=(1009.48,'K')), NASAPolynomial(coeffs=[14.4737,0.0351559,-1.35266e-05,2.51413e-09,-1.78814e-13,8304.92,-47.6439], Tmin=(1009.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(ROOJ)"""),
)

species(
    label = '[O]OC1C=C=CCCC1(664)',
    structure = SMILES('[O]OC1C=C=CCCC1'),
    E0 = (235.574,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773237,0.0550315,1.12729e-05,-5.13422e-08,2.27882e-11,28463.3,26.6705], Tmin=(100,'K'), Tmax=(994.993,'K')), NASAPolynomial(coeffs=[15.3911,0.03357,-1.26105e-05,2.34066e-09,-1.6753e-13,23707.8,-53.058], Tmin=(994.993,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(ROOJ)"""),
)

species(
    label = '[O]O[C]1CC=CCCC1(665)',
    structure = SMILES('[O]O[C]1CC=CCCC1'),
    E0 = (190.435,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.870082,0.0550154,8.59402e-06,-4.11597e-08,1.70556e-11,23028.6,27.666], Tmin=(100,'K'), Tmax=(1049.13,'K')), NASAPolynomial(coeffs=[12.5239,0.041641,-1.66893e-05,3.12375e-09,-2.20734e-13,18874,-37.2619], Tmin=(1049.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(ROOJ) + radical(C2CsJOOH)"""),
)

species(
    label = '[O]OC1C[C]=CCCC1(666)',
    structure = SMILES('[O]OC1C[C]=CCCC1'),
    E0 = (241.378,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716786,0.0569463,9.11191e-06,-4.65214e-08,2.01352e-11,29162.5,27.7003], Tmin=(100,'K'), Tmax=(1019.28,'K')), NASAPolynomial(coeffs=[14.3458,0.0383644,-1.49073e-05,2.78408e-09,-1.98126e-13,24571.1,-47.2063], Tmin=(1019.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'OO[C]1[CH]C=CCCC1(667)',
    structure = SMILES('OO[C]1[CH]C=CCCC1'),
    E0 = (155.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.376305,0.062697,5.11344e-06,-4.61317e-08,2.02508e-11,18829,25.8779], Tmin=(100,'K'), Tmax=(1036.01,'K')), NASAPolynomial(coeffs=[16.4435,0.0389323,-1.58885e-05,3.03898e-09,-2.18757e-13,13446.1,-62.1175], Tmin=(1036.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C2CsJOOH) + radical(C=CCJCO)"""),
)

species(
    label = '[O]OC1[CH]CCC=CC1(668)',
    structure = SMILES('[O]OC1[CH]CCC=CC1'),
    E0 = (203.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.779627,0.0543222,1.81484e-05,-5.67967e-08,2.39507e-11,24660,29.3058], Tmin=(100,'K'), Tmax=(1008.98,'K')), NASAPolynomial(coeffs=[14.7162,0.0376864,-1.45262e-05,2.72294e-09,-1.9504e-13,19882.1,-47.795], Tmin=(1008.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(CCJCOOH) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1CC=[C]CCC1(669)',
    structure = SMILES('[O]OC1CC=[C]CCC1'),
    E0 = (241.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716767,0.0569465,9.11111e-06,-4.65204e-08,2.01347e-11,29162.5,27.7004], Tmin=(100,'K'), Tmax=(1019.28,'K')), NASAPolynomial(coeffs=[14.3459,0.0383643,-1.49073e-05,2.78406e-09,-1.98125e-13,24571,-47.2068], Tmin=(1019.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1C[CH]CC=CC1(670)',
    structure = SMILES('[O]OC1C[CH]CC=CC1'),
    E0 = (197.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00915,0.0527075,1.06241e-05,-4.09782e-08,1.66009e-11,23932,28.7437], Tmin=(100,'K'), Tmax=(1051.91,'K')), NASAPolynomial(coeffs=[11.4311,0.0425657,-1.6964e-05,3.1562e-09,-2.21936e-13,20108,-29.8182], Tmin=(1051.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC) + radical(ROOJ)"""),
)

species(
    label = 'OOC1[CH]C=CCC[CH]1(671)',
    structure = SMILES('OOC1[CH]C=CCC[CH]1'),
    E0 = (168.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281704,0.0620467,1.45492e-05,-6.16654e-08,2.71259e-11,20460.7,27.5329], Tmin=(100,'K'), Tmax=(1005.47,'K')), NASAPolynomial(coeffs=[18.6915,0.0348872,-1.36749e-05,2.62653e-09,-1.92116e-13,14429.3,-72.967], Tmin=(1005.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CCJCO) + radical(CCJCOOH)"""),
)

species(
    label = '[O]OC1CC=C[CH]CC1(672)',
    structure = SMILES('[O]OC1CC=C[CH]CC1'),
    E0 = (144.648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.968011,0.0469029,4.35243e-05,-8.49853e-08,3.42936e-11,17523.9,25.2173], Tmin=(100,'K'), Tmax=(990.286,'K')), NASAPolynomial(coeffs=[15.3972,0.036678,-1.37824e-05,2.5995e-09,-1.88958e-13,12309.6,-56.1523], Tmin=(990.286,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(ROOJ)"""),
)

species(
    label = 'OOC1[CH]C=CC[CH]C1(673)',
    structure = SMILES('OOC1[CH]C=CC[CH]C1'),
    E0 = (162.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.515821,0.0603844,7.15708e-06,-4.59632e-08,1.97994e-11,19732.5,26.954], Tmin=(100,'K'), Tmax=(1037.54,'K')), NASAPolynomial(coeffs=[15.345,0.0398661,-1.61683e-05,3.0726e-09,-2.20053e-13,14682.5,-54.6416], Tmin=(1037.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CCJCO) + radical(RCCJCC)"""),
)

species(
    label = 'OOC1[CH][C]=CCCC1(674)',
    structure = SMILES('OOC1[CH][C]=CCCC1'),
    E0 = (206.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.220296,0.0646547,5.56376e-06,-5.14479e-08,2.33307e-11,24963.1,25.9223], Tmin=(100,'K'), Tmax=(1013.36,'K')), NASAPolynomial(coeffs=[18.3091,0.035585,-1.40673e-05,2.69027e-09,-1.95415e-13,19123.5,-72.3101], Tmin=(1013.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = 'OOC1[CH]C=C[CH]CC1(675)',
    structure = SMILES('OOC1[CH]C=C[CH]CC1'),
    E0 = (109.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.468385,0.0546469,3.98598e-05,-8.97754e-08,3.74384e-11,13324.6,23.4507], Tmin=(100,'K'), Tmax=(990.154,'K')), NASAPolynomial(coeffs=[19.384,0.0338597,-1.29204e-05,2.50057e-09,-1.85827e-13,6851.87,-81.3889], Tmin=(990.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CCJCO) + radical(Allyl_S)"""),
)

species(
    label = 'OOC1[CH]C=[C]CCC1(676)',
    structure = SMILES('OOC1[CH]C=[C]CCC1'),
    E0 = (206.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.220296,0.0646547,5.56376e-06,-5.14479e-08,2.33307e-11,24963.1,25.9223], Tmin=(100,'K'), Tmax=(1013.36,'K')), NASAPolynomial(coeffs=[18.3091,0.035585,-1.40673e-05,2.69027e-09,-1.95415e-13,19123.5,-72.3101], Tmin=(1013.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]1C=C[CH]CCC1(51)',
    structure = SMILES('[CH]1C=C[CH]CCC1'),
    E0 = (233.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33755,0.00875945,0.000132706,-1.72324e-07,6.44554e-11,28131.9,15.3206], Tmin=(100,'K'), Tmax=(967.659,'K')), NASAPolynomial(coeffs=[14.2911,0.0302458,-1.05031e-05,2.05604e-09,-1.58625e-13,22499.1,-59.1067], Tmin=(967.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Allyl_S)"""),
)

species(
    label = '[O]OC1CCCC2[CH]C21(677)',
    structure = SMILES('[O]OC1CCCC2[CH]C21'),
    E0 = (276.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17104,0.0441178,4.27603e-05,-7.88652e-08,3.08632e-11,33410.4,28.3331], Tmin=(100,'K'), Tmax=(1006.25,'K')), NASAPolynomial(coeffs=[13.4358,0.03921,-1.52856e-05,2.89594e-09,-2.09056e-13,28722.3,-41.9438], Tmin=(1006.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(ROOJ) + radical(bicyclo[4.1.0]heptane-C3-7)"""),
)

species(
    label = '[CH]1[CH]C2CCCC1OO2(678)',
    structure = SMILES('[CH]1[CH]C2CCCC1OO2'),
    E0 = (223.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35431,0.0316878,9.24132e-05,-1.37987e-07,5.32343e-11,27027.1,26.7693], Tmin=(100,'K'), Tmax=(980.402,'K')), NASAPolynomial(coeffs=[17.0604,0.0339692,-1.26094e-05,2.46928e-09,-1.86855e-13,20758.1,-64.9578], Tmin=(980.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(CCJCOOH) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]CCC=CC=CO[O](679)',
    structure = SMILES('[CH2]CCC=CC=CO[O]'),
    E0 = (284.423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.623348,0.0912953,-8.16512e-05,3.79422e-08,-7.01613e-12,34383.4,38.0045], Tmin=(100,'K'), Tmax=(1307.87,'K')), NASAPolynomial(coeffs=[19.5339,0.0296461,-1.09454e-05,1.90087e-09,-1.26794e-13,29110.8,-64.6526], Tmin=(1307.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(RCCJ)"""),
)

species(
    label = 'OOC1=CC=CCCC1(680)',
    structure = SMILES('OOC1=CC=CCCC1'),
    E0 = (-46.7992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.521751,0.0593845,1.04347e-05,-5.25623e-08,2.31551e-11,-5488.39,27.1941], Tmin=(100,'K'), Tmax=(1008.11,'K')), NASAPolynomial(coeffs=[16.2574,0.0364785,-1.43009e-05,2.69202e-09,-1.93254e-13,-10669.7,-58.8111], Tmin=(1008.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.7992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'OOC1C=C=CCCC1(681)',
    structure = SMILES('OOC1C=C=CCCC1'),
    E0 = (83.5723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493365,0.0587636,1.56321e-05,-6.04015e-08,2.6532e-11,10193.9,26.9211], Tmin=(100,'K'), Tmax=(995.898,'K')), NASAPolynomial(coeffs=[17.1802,0.0348834,-1.33795e-05,2.5173e-09,-1.81868e-13,4730.81,-64.2557], Tmin=(995.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.5723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'O2(S)(682)',
    structure = SMILES('O=O'),
    E0 = (85.6848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.0012157,5.31615e-06,-4.8944e-09,1.45844e-12,10304.5,4.68369], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15383,0.00167803,-7.69968e-07,1.51274e-10,-1.08781e-14,10302.3,6.16752], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O(11)',
    structure = SMILES('[O]'),
    E0 = (242.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.16827,-0.00327932,6.64306e-06,-6.12807e-09,2.11266e-12,29122.3,2.05193], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.54364,-2.73162e-05,-4.1903e-09,4.95482e-12,-4.79554e-16,29226,4.92229], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(242.977,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C1=CC2OC2CCC1(683)',
    structure = SMILES('C1=CC2OC2CCC1'),
    E0 = (-17.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65264,0.0264344,9.44822e-05,-1.37705e-07,5.31834e-11,-2037.78,21.0339], Tmin=(100,'K'), Tmax=(970.462,'K')), NASAPolynomial(coeffs=[15.7639,0.0315699,-1.12933e-05,2.1686e-09,-1.63344e-13,-7757.39,-61.9789], Tmin=(970.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[O]O[CH]C1C=CCCC1(684)',
    structure = SMILES('[O]O[CH]C1C=CCCC1'),
    E0 = (183.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.70662,0.0569659,1.0356e-05,-4.78139e-08,2.05108e-11,22215.8,29.8483], Tmin=(100,'K'), Tmax=(1021.2,'K')), NASAPolynomial(coeffs=[14.3269,0.0390132,-1.52676e-05,2.85669e-09,-2.03349e-13,17588.3,-45.1843], Tmin=(1021.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(494.711,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(CCsJOOH) + radical(ROOJ)"""),
)

species(
    label = 'HO2(17)',
    structure = SMILES('[O]O'),
    E0 = (2.49012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.3018,-0.00474912,2.11583e-05,-2.42764e-08,9.29225e-12,264.018,3.71666], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.17229,0.00188118,-3.46277e-07,1.94658e-11,1.76257e-16,31.0207,2.95768], Tmin=(1000,'K'), Tmax=(5000,'K'))], Tmin=(200,'K'), Tmax=(5000,'K'), E0=(2.49012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]1C=CC=CCC1(43)',
    structure = SMILES('[CH]1C=CCCC=C1'),
    E0 = (167.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2863,0.00959575,0.000128536,-1.73375e-07,6.6947e-11,20259.5,17.7361], Tmin=(100,'K'), Tmax=(947.927,'K')), NASAPolynomial(coeffs=[15.9797,0.0233649,-6.47559e-06,1.20576e-09,-9.641e-14,14448.7,-64.5505], Tmin=(947.927,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CCJC=C)"""),
)

species(
    label = 'C1=CC2OOC2CCC1(685)',
    structure = SMILES('C1=CC2OOC2CCC1'),
    E0 = (-5.69596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19083,0.0355319,8.4852e-05,-1.33057e-07,5.23263e-11,-560.367,23.419], Tmin=(100,'K'), Tmax=(973.533,'K')), NASAPolynomial(coeffs=[17.7725,0.0328749,-1.19335e-05,2.30181e-09,-1.73383e-13,-6891.58,-72.0683], Tmin=(973.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.69596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[O]C1[CH]C=CCCC1(686)',
    structure = SMILES('[O]C1[CH]C=CCCC1'),
    E0 = (127.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13134,0.038037,7.38268e-05,-1.20787e-07,4.7798e-11,15503.7,20.0508], Tmin=(100,'K'), Tmax=(979.355,'K')), NASAPolynomial(coeffs=[17.9206,0.0316677,-1.16895e-05,2.27943e-09,-1.72244e-13,9232.1,-75.8269], Tmin=(979.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O]OC1C=CC=CCC1(687)',
    structure = SMILES('[O]OC1C=CC=CCC1'),
    E0 = (121.738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.842545,0.0533214,1.51516e-05,-5.4567e-08,2.37224e-11,14769.6,26.6027], Tmin=(100,'K'), Tmax=(995.853,'K')), NASAPolynomial(coeffs=[15.0932,0.0339889,-1.28265e-05,2.38646e-09,-1.70991e-13,10051.6,-51.5266], Tmin=(995.853,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(ROOJ)"""),
)

species(
    label = '[O]OC1C=CC[CH]CC1(688)',
    structure = SMILES('[O]OC1C=CC[CH]CC1'),
    E0 = (198.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02151,0.0523808,1.14447e-05,-4.17397e-08,1.68465e-11,23999.2,28.7396], Tmin=(100,'K'), Tmax=(1051.21,'K')), NASAPolynomial(coeffs=[11.3932,0.0426097,-1.69843e-05,3.16101e-09,-2.22344e-13,20177.9,-29.6195], Tmin=(1051.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1C=[C]CCCC1(689)',
    structure = SMILES('[O]OC1C=[C]CCCC1'),
    E0 = (241.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.72889,0.0566222,9.92444e-06,-4.7275e-08,2.03786e-11,29229.7,27.6971], Tmin=(100,'K'), Tmax=(1018.96,'K')), NASAPolynomial(coeffs=[14.311,0.0384035,-1.49249e-05,2.78825e-09,-1.98483e-13,24639.7,-47.0249], Tmin=(1018.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]OC1C=CCC[CH]C1(690)',
    structure = SMILES('[O]OC1C=CCC[CH]C1'),
    E0 = (198.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02151,0.0523808,1.14447e-05,-4.17397e-08,1.68465e-11,23999.2,28.7396], Tmin=(100,'K'), Tmax=(1051.21,'K')), NASAPolynomial(coeffs=[11.3932,0.0426097,-1.69843e-05,3.16101e-09,-2.22344e-13,20177.9,-29.6195], Tmin=(1051.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJCC) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1[C]=CCCCC1(691)',
    structure = SMILES('[O]OC1[C]=CCCCC1'),
    E0 = (241.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728931,0.0566217,9.92616e-06,-4.72772e-08,2.03795e-11,29229.7,27.697], Tmin=(100,'K'), Tmax=(1018.95,'K')), NASAPolynomial(coeffs=[14.3108,0.0384038,-1.49251e-05,2.78829e-09,-1.98487e-13,24639.7,-47.0238], Tmin=(1018.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1[CH]CCCC=C1(63)',
    structure = SMILES('[O]OC1[CH]CCCC=C1'),
    E0 = (204.513,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.791654,0.0539991,1.89583e-05,-5.75473e-08,2.41931e-11,24727.2,29.3028], Tmin=(100,'K'), Tmax=(1008.76,'K')), NASAPolynomial(coeffs=[14.682,0.0377243,-1.45431e-05,2.72696e-09,-1.95385e-13,19950.4,-47.6175], Tmin=(1008.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(ROOJ) + radical(CCJCOOH)"""),
)

species(
    label = 'OOC1[C]=C[CH]CCC1(692)',
    structure = SMILES('OOC1[C]=C[CH]CCC1'),
    E0 = (231.047,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.655702,0.0550733,2.33326e-05,-6.48226e-08,2.70427e-11,27925.1,26.094], Tmin=(100,'K'), Tmax=(1011.1,'K')), NASAPolynomial(coeffs=[16.1834,0.0371646,-1.46632e-05,2.80016e-09,-2.03144e-13,22560.5,-59.9905], Tmin=(1011.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[O]O[C]1C=CCCCC1(693)',
    structure = SMILES('[O]OC1=C[CH]CCCC1'),
    E0 = (130.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.982444,0.0473415,4.02916e-05,-8.03088e-08,3.2323e-11,15827.3,25.4812], Tmin=(100,'K'), Tmax=(995.249,'K')), NASAPolynomial(coeffs=[14.8739,0.0375245,-1.42635e-05,2.68952e-09,-1.94748e-13,10783.3,-52.92], Tmin=(995.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(ROOJ) + radical(Allyl_S)"""),
)

species(
    label = '[CH]1[CH]C2OOC2CCC1(694)',
    structure = SMILES('[CH]1[CH]C2OOC2CCC1'),
    E0 = (245.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36106,0.0352578,7.43725e-05,-1.135e-07,4.31106e-11,29620.5,27.8237], Tmin=(100,'K'), Tmax=(999.305,'K')), NASAPolynomial(coeffs=[14.7815,0.0382071,-1.51163e-05,2.94816e-09,-2.18229e-13,24108.8,-51.0701], Tmin=(999.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cycloheptane) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]CC(C=CC=C)O[O](695)',
    structure = SMILES('[CH2]CC(C=CC=C)O[O]'),
    E0 = (267.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,492.5,1135,1000,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.596758,0.0908189,-7.85888e-05,3.49024e-08,-6.18352e-12,32337.1,38.0388], Tmin=(100,'K'), Tmax=(1356.89,'K')), NASAPolynomial(coeffs=[19.878,0.030461,-1.18649e-05,2.11963e-09,-1.43464e-13,26780.7,-66.9888], Tmin=(1356.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(RCCJ)"""),
)

species(
    label = 'OOC1C=CC=CCC1(696)',
    structure = SMILES('OOC1C=CC=CCC1'),
    E0 = (-30.2638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.56275,0.0570526,1.95143e-05,-6.36311e-08,2.74682e-11,-3499.75,26.853], Tmin=(100,'K'), Tmax=(996.619,'K')), NASAPolynomial(coeffs=[16.882,0.0353027,-1.35958e-05,2.56316e-09,-1.85332e-13,-8925.24,-62.7231], Tmin=(996.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.2638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'C1=CC2CCCC1O2(697)',
    structure = SMILES('C1=CC2CCCC1O2'),
    E0 = (-102.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94555,0.0142151,0.000135137,-1.84073e-07,7.0679e-11,-12206.7,19.1362], Tmin=(100,'K'), Tmax=(958.191,'K')), NASAPolynomial(coeffs=[17.5827,0.0275552,-8.81796e-06,1.71342e-09,-1.35216e-13,-18812.4,-74.4687], Tmin=(958.191,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH2]C1C=CC(CC1)O[O](698)',
    structure = SMILES('[CH2]C1C=CC(CC1)O[O]'),
    E0 = (197.001,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.625357,0.0597423,3.00781e-06,-4.29005e-08,1.99223e-11,23828,30.0141], Tmin=(100,'K'), Tmax=(989.98,'K')), NASAPolynomial(coeffs=[14.5391,0.0370737,-1.34792e-05,2.43447e-09,-1.70842e-13,19429.1,-45.2751], Tmin=(989.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(494.711,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Isobutyl) + radical(ROOJ)"""),
)

species(
    label = '[C]1=CCCCC=C1(45)',
    structure = SMILES('[C]1=CCCCC=C1'),
    E0 = (267.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89666,0.0281884,6.2347e-05,-9.60905e-08,3.74486e-11,32217.1,19.2007], Tmin=(100,'K'), Tmax=(970.591,'K')), NASAPolynomial(coeffs=[12.303,0.0303972,-1.07588e-05,1.99276e-09,-1.4509e-13,28073,-41.6352], Tmin=(970.591,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(C=CJC=C)"""),
)

species(
    label = 'C1=CC2CCCC1OO2(699)',
    structure = SMILES('C1=CC2CCCC1OO2'),
    E0 = (-57.1003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38895,0.0283663,0.000106983,-1.57402e-07,6.13316e-11,-6747.37,21.3203], Tmin=(100,'K'), Tmax=(966.536,'K')), NASAPolynomial(coeffs=[18.3669,0.0311926,-1.08334e-05,2.09959e-09,-1.61053e-13,-13443.3,-77.6714], Tmin=(966.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.1003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
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
    E0 = (120.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (218.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (328.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (463.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (312.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (441.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (290.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (331.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (378.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (239.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (252.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (563.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (196.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (357.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (174.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (383.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (224.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (367.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (223.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (311.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (198.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (145.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (153.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (250.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (343.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (245.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (128.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (370.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (335.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (323.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (442.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (326.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (378.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (246.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (332.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (548.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (245.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (294.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (138.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (174.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (356.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (299.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (127.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C7H10(1)', 'O2(9)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.971438,'m^3/(mol*s)'), n=2.13093, Ea=(60.9849,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-CdH;YJ] for rate rule [Cds-CsH_Cds-CdH;O2b]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond
Ea raised from -1.7 to 61.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['[CH]1CCCC2[CH]C1OO2(662)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.2, Ea=(98.5227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R6;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 95.0 to 98.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(10)', '[O]OC1=CC=CCCC1(663)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(231.051,'m^3/(mol*s)'), n=1.52235, Ea=(11.3887,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(10)', '[O]OC1C=C=CCCC1(664)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]O[C]1CC=CCCC1(665)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(20108.5,'s^-1'), n=2.606, Ea=(121.63,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_NDMustO;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC1C[C]=CCCC1(666)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OO[C]1[CH]C=CCCC1(667)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC1[CH]CCC=CC1(668)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 219 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['[O]OC1CC=[C]CCC1(669)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.96e+10,'s^-1'), n=0.83, Ea=(257.734,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 201 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC1C[CH]CC=CC1(670)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.502,'s^-1'), n=3.86, Ea=(41.6308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 332 used for R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OOC1[CH]C=CCC[CH]1(671)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 323 used for R4H_SSS;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC1CC=C[CH]CC1(672)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 334 used for R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['OOC1[CH]C=CC[CH]C1(673)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.788e+07,'s^-1'), n=1.26, Ea=(76.0233,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 260 used for R5H_SSSS_OCC;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R5H_SSSS_OCC;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['OOC1[CH][C]=CCCC1(674)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['OOC1[CH]C=C[CH]CC1(675)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(120.747,'s^-1'), n=2.0375, Ea=(54.2874,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R6H_SSSSS_OO;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['OOC1[CH]C=[C]CCC1(676)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.44313e+09,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;XH_out] for rate rule [R6HJ_2;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(9)', '[CH]1C=C[CH]CCC1(51)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.16777e+06,'m^3/(mol*s)'), n=0.193752, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination
Ea raised from -25.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['[O]OC1CCCC2[CH]C21(677)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.18204e+14,'s^-1'), n=-0.478585, Ea=(247.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_cs] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['[CH]1[CH]C2CCCC1OO2(678)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.84441e+10,'s^-1'), n=0.419705, Ea=(103.275,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 99.9 to 103.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CCC=CC=CO[O](679)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_HNd;radadd_intra_cs2H] for rate rule [R7_linear;doublebond_intra_pri_HNd_O;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['OOC1=CC=CCCC1(680)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['OOC1C=C=CCCC1(681)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['C7H10(1)', 'O2(S)(682)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(33.3443,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 33.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['O(11)', 'C1=CC2OC2CCC1(683)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.27e+09,'s^-1'), n=1.06, Ea=(129.704,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2OO_S;C_sec_rad_intra;OOJ] for rate rule [R2OO_S;C_rad/H/OneDe_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]O[CH]C1C=CCCC1(684)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['HO2(17)', '[CH]1C=CC=CCC1(43)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(8.00406e+10,'s^-1'), n=0.563333, Ea=(124.683,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_HNd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['C1=CC2OOC2CCC1(685)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O(11)', '[O]C1[CH]C=CCCC1(686)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(9.16281e+06,'m^3/(mol*s)'), n=0.183438, Ea=(0.123953,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(10)', '[O]OC1C=CC=CCC1(687)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.35e+08,'cm^3/(mol*s)'), n=1.64, Ea=(1.58992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2557 used for Cds-CsH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]OC1C=CC[CH]CC1(688)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]OC1C=[C]CCCC1(689)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]OC1C=CCC[CH]C1(690)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 219 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['[O]OC1[C]=CCCCC1(691)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.96e+10,'s^-1'), n=0.83, Ea=(257.734,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 201 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]OC1[CH]CCCC=C1(63)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(0.502,'s^-1'), n=3.86, Ea=(41.6308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 332 used for R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['OOC1[C]=C[CH]CCC1(692)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]O[C]1C=CCCCC1(693)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R4H_SDS;C_rad_out_NDMustO;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['[CH]1[CH]C2OOC2CCC1(694)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.01908e+09,'s^-1'), n=0.467551, Ea=(124.869,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic
Ea raised from 121.8 to 124.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]CC(C=CC=C)O[O](695)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.2, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 46 used for R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['OOC1C=CC=CCC1(696)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_De] for rate rule [R6radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['O(11)', 'C1=CC2CCCC1O2(697)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.63e+10,'s^-1'), n=0, Ea=(54.392,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;C_rad/H/NonDeC_intra;OO] for rate rule [R4OO_SDS;C_rad/H/NonDeC_intra;OOJ]
Euclidian distance = 1.41421356237
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C1C=CC(CC1)O[O](698)'],
    products = ['[O]OC1[CH]C=CCCC1(62)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['HO2(17)', '[C]1=CCCCC=C1(45)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O]OC1[CH]C=CCCC1(62)'],
    products = ['C1=CC2CCCC1OO2(699)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '6',
    isomers = [
        '[O]OC1[CH]C=CCCC1(62)',
    ],
    reactants = [
        ('C7H10(1)', 'O2(9)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '6',
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

