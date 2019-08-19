species(
    label = 'S(214)(213)',
    structure = SMILES('C=CC=CCOO'),
    E0 = (-23.7352,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.19487,'amu*angstrom^2'), symmetry=1, barrier=(27.4724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.801598,'amu*angstrom^2'), symmetry=1, barrier=(18.4303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.801555,'amu*angstrom^2'), symmetry=1, barrier=(18.4293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.801522,'amu*angstrom^2'), symmetry=1, barrier=(18.4286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3857.61,'J/mol'), sigma=(6.30207,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=602.55 K, Pc=34.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.560968,0.066755,-5.38307e-05,2.21667e-08,-3.65487e-12,-2723.52,26.6722], Tmin=(100,'K'), Tmax=(1446.91,'K')), NASAPolynomial(coeffs=[16.092,0.0238198,-9.3209e-06,1.65903e-09,-1.1157e-13,-7218,-53.9942], Tmin=(1446.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.7352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'OH(5)(5)',
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
    label = 'C5H7O(216)(215)',
    structure = SMILES('C=CC=CC[O]'),
    E0 = (130.465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,324.17,324.798,326.428],'cm^-1')),
        HinderedRotor(inertia=(0.0834456,'amu*angstrom^2'), symmetry=1, barrier=(6.17438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147528,'amu*angstrom^2'), symmetry=1, barrier=(11.0638,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3714.19,'J/mol'), sigma=(6.10713,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=580.15 K, Pc=37 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80098,0.0446979,-2.65813e-05,7.81281e-09,-9.36109e-13,15773.2,22.3196], Tmin=(100,'K'), Tmax=(1862.48,'K')), NASAPolynomial(coeffs=[11.7735,0.0232802,-9.33191e-06,1.63842e-09,-1.07321e-13,12058.5,-31.9938], Tmin=(1862.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = 'HO2(8)(9)',
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
    label = 'C5H7(211)(210)',
    structure = SMILES('[CH2]C=CC=C'),
    E0 = (175.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.52628,'amu*angstrom^2'), symmetry=1, barrier=(35.0921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5297,'amu*angstrom^2'), symmetry=1, barrier=(35.1709,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3140.68,'J/mol'), sigma=(5.4037,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=490.57 K, Pc=45.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2987,0.0238036,3.93739e-05,-6.93177e-08,2.88419e-11,21235.8,16.7723], Tmin=(100,'K'), Tmax=(942.053,'K')), NASAPolynomial(coeffs=[11.9106,0.0173605,-5.0926e-06,8.77895e-10,-6.40305e-14,17899.7,-37.1203], Tmin=(942.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'H(3)(3)',
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
    label = 'S(213)(212)',
    structure = SMILES('C=CC=CCO[O]'),
    E0 = (128.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.617223,'amu*angstrom^2'), symmetry=1, barrier=(14.1912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.616803,'amu*angstrom^2'), symmetry=1, barrier=(14.1815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.616486,'amu*angstrom^2'), symmetry=1, barrier=(14.1742,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3857.61,'J/mol'), sigma=(6.30207,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=602.55 K, Pc=34.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03674,0.0608182,-5.10545e-05,2.27959e-08,-4.14589e-12,15537.5,25.7137], Tmin=(100,'K'), Tmax=(1305.45,'K')), NASAPolynomial(coeffs=[12.7278,0.0249956,-9.89297e-06,1.77541e-09,-1.20324e-13,12485.1,-33.805], Tmin=(1305.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'CH2OOH(33)(34)',
    structure = SMILES('[CH2]OO'),
    E0 = (54.8878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.878704,'amu*angstrom^2'), symmetry=1, barrier=(20.2031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0333,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3467.13,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85644,0.0201433,-7.52661e-06,-9.66998e-09,6.63679e-12,6647.28,10.9503], Tmin=(100,'K'), Tmax=(928.333,'K')), NASAPolynomial(coeffs=[9.91697,0.00263661,-1.08458e-07,-1.03744e-11,-4.90104e-16,4779.83,-25.5851], Tmin=(928.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.8878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CsJOOH)"""),
)

species(
    label = 'CH2CHCHCH(1116)',
    structure = SMILES('[CH]=CC=C'),
    E0 = (346.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.31937,'amu*angstrom^2'), symmetry=1, barrier=(30.3349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64255,0.0163337,3.86225e-05,-6.71377e-08,2.83603e-11,41729.6,13.282], Tmin=(100,'K'), Tmax=(937.724,'K')), NASAPolynomial(coeffs=[12.9705,0.00669127,-1.0007e-06,1.67602e-10,-1.71436e-14,38279.7,-43.9476], Tmin=(937.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH2CHCHCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=CC=C[CH]OO(1526)',
    structure = SMILES('[CH2]C=CC=COO'),
    E0 = (91.6536,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.20561,'amu*angstrom^2'), symmetry=1, barrier=(27.7194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2057,'amu*angstrom^2'), symmetry=1, barrier=(27.7213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20579,'amu*angstrom^2'), symmetry=1, barrier=(27.7235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20561,'amu*angstrom^2'), symmetry=1, barrier=(27.7194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.899917,0.0557286,-1.7294e-05,-2.40356e-08,1.51562e-11,11146.3,26.3877], Tmin=(100,'K'), Tmax=(954.006,'K')), NASAPolynomial(coeffs=[16.991,0.0183508,-5.83461e-06,1.0174e-09,-7.27006e-14,6706.82,-57.6612], Tmin=(954.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.6536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=[C]COO(1527)',
    structure = SMILES('C=CC=[C]COO'),
    E0 = (214.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1685,370,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.713623,'amu*angstrom^2'), symmetry=1, barrier=(16.4076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714834,'amu*angstrom^2'), symmetry=1, barrier=(16.4354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31608,'amu*angstrom^2'), symmetry=1, barrier=(30.2593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714278,'amu*angstrom^2'), symmetry=1, barrier=(16.4227,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.961748,0.0663851,-6.19497e-05,3.08678e-08,-6.30013e-12,25860.8,25.6989], Tmin=(100,'K'), Tmax=(1164.68,'K')), NASAPolynomial(coeffs=[12.266,0.0275615,-1.19484e-05,2.24691e-09,-1.56602e-13,23227.7,-30.5608], Tmin=(1164.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C2H3(28)(29)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=CCOO(2784)',
    structure = SMILES('[CH]=CCOO'),
    E0 = (177.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350],'cm^-1')),
        HinderedRotor(inertia=(0.5528,'amu*angstrom^2'), symmetry=1, barrier=(12.71,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.552819,'amu*angstrom^2'), symmetry=1, barrier=(12.7104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.552853,'amu*angstrom^2'), symmetry=1, barrier=(12.7112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71936,0.0552392,-8.46289e-05,7.36628e-08,-2.53867e-11,21398.1,18.8788], Tmin=(100,'K'), Tmax=(832.126,'K')), NASAPolynomial(coeffs=[6.68181,0.0228008,-1.06812e-05,2.0218e-09,-1.38691e-13,20869.4,-2.36441], Tmin=(832.126,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_P)"""),
)

species(
    label = 'C=C[C]=CCOO(1528)',
    structure = SMILES('C=C[C]=CCOO'),
    E0 = (175.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1685,370,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.888053,'amu*angstrom^2'), symmetry=1, barrier=(20.4181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.887069,'amu*angstrom^2'), symmetry=1, barrier=(20.3955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.47831,'amu*angstrom^2'), symmetry=1, barrier=(56.9812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12501,'amu*angstrom^2'), symmetry=1, barrier=(25.8662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.863667,0.0689364,-6.83298e-05,3.71636e-08,-8.27514e-12,21192,25.1756], Tmin=(100,'K'), Tmax=(1076.53,'K')), NASAPolynomial(coeffs=[11.7776,0.0283842,-1.18263e-05,2.17281e-09,-1.4938e-13,18842.1,-28.283], Tmin=(1076.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=CCOO(1529)',
    structure = SMILES('C=C=C[CH]COO'),
    E0 = (152.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.48424,'amu*angstrom^2'), symmetry=1, barrier=(34.1257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48397,'amu*angstrom^2'), symmetry=1, barrier=(34.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48438,'amu*angstrom^2'), symmetry=1, barrier=(34.1287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48424,'amu*angstrom^2'), symmetry=1, barrier=(34.1257,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.757171,0.0721266,-6.86944e-05,3.43735e-08,-7.04756e-12,18475.8,23.2457], Tmin=(100,'K'), Tmax=(1156.69,'K')), NASAPolynomial(coeffs=[13.0408,0.0296478,-1.36075e-05,2.62365e-09,-1.85305e-13,15634.1,-37.8038], Tmin=(1156.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=CC=CCOO(1530)',
    structure = SMILES('[CH]=CC=CCOO'),
    E0 = (223.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.904379,'amu*angstrom^2'), symmetry=1, barrier=(20.7934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.904711,'amu*angstrom^2'), symmetry=1, barrier=(20.8011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.904858,'amu*angstrom^2'), symmetry=1, barrier=(20.8045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.904853,'amu*angstrom^2'), symmetry=1, barrier=(20.8044,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.726593,0.067598,-6.17808e-05,2.91805e-08,-5.52996e-12,26985.6,26.4717], Tmin=(100,'K'), Tmax=(1266.37,'K')), NASAPolynomial(coeffs=[14.6304,0.0236808,-9.76104e-06,1.79517e-09,-1.23669e-13,23464.1,-43.8895], Tmin=(1266.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=C[CH]C[CH]OO(2785)',
    structure = SMILES('[CH2]C=CC[CH]OO'),
    E0 = (199.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,312.97,2071.53],'cm^-1')),
        HinderedRotor(inertia=(0.139574,'amu*angstrom^2'), symmetry=1, barrier=(9.69805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.338808,'amu*angstrom^2'), symmetry=1, barrier=(23.5421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.952172,'amu*angstrom^2'), symmetry=1, barrier=(66.1771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.614117,'amu*angstrom^2'), symmetry=1, barrier=(42.6815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.338702,'amu*angstrom^2'), symmetry=1, barrier=(23.5426,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.951538,0.0686429,-5.97536e-05,2.81807e-08,-5.58191e-12,24068,26.7243], Tmin=(100,'K'), Tmax=(1175.14,'K')), NASAPolynomial(coeffs=[11.2089,0.0337284,-1.51872e-05,2.89785e-09,-2.03242e-13,21657.2,-24.4171], Tmin=(1175.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCsJOOH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C[C]=CCOO(2786)',
    structure = SMILES('[CH2]C[C]=CCOO'),
    E0 = (308.91,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,341.888,3156.17],'cm^-1')),
        HinderedRotor(inertia=(0.108703,'amu*angstrom^2'), symmetry=1, barrier=(9.01591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470691,'amu*angstrom^2'), symmetry=1, barrier=(39.0407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108698,'amu*angstrom^2'), symmetry=1, barrier=(9.01572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1087,'amu*angstrom^2'), symmetry=1, barrier=(9.01596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470678,'amu*angstrom^2'), symmetry=1, barrier=(39.0399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12989,0.0683192,-7.04511e-05,4.48538e-08,-1.2515e-11,37252.1,28.9814], Tmin=(100,'K'), Tmax=(842.221,'K')), NASAPolynomial(coeffs=[7.30164,0.0390077,-1.82475e-05,3.53205e-09,-2.49421e-13,36212.5,0.265983], Tmin=(842.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C[CH]COO(2787)',
    structure = SMILES('C=[C]C[CH]COO'),
    E0 = (309.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1685,370,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82257,0.0754843,-8.93805e-05,6.45146e-08,-1.98474e-11,37349.5,29.9861], Tmin=(100,'K'), Tmax=(778.04,'K')), NASAPolynomial(coeffs=[8.00498,0.0385587,-1.81911e-05,3.51595e-09,-2.47337e-13,36231.8,-2.86227], Tmin=(778.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCOOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC=[C]COO(2788)',
    structure = SMILES('[CH2]CC=[C]COO'),
    E0 = (308.91,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,341.888,3156.17],'cm^-1')),
        HinderedRotor(inertia=(0.108703,'amu*angstrom^2'), symmetry=1, barrier=(9.01591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470691,'amu*angstrom^2'), symmetry=1, barrier=(39.0407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108698,'amu*angstrom^2'), symmetry=1, barrier=(9.01572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1087,'amu*angstrom^2'), symmetry=1, barrier=(9.01596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470678,'amu*angstrom^2'), symmetry=1, barrier=(39.0399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12989,0.0683192,-7.04511e-05,4.48538e-08,-1.2515e-11,37252.1,28.9814], Tmin=(100,'K'), Tmax=(842.221,'K')), NASAPolynomial(coeffs=[7.30164,0.0390077,-1.82475e-05,3.53205e-09,-2.49421e-13,36212.5,0.265983], Tmin=(842.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC[CH]COO(2789)',
    structure = SMILES('[CH]=CC[CH]COO'),
    E0 = (318.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,313.231,313.239],'cm^-1')),
        HinderedRotor(inertia=(0.00171814,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110609,'amu*angstrom^2'), symmetry=1, barrier=(7.70137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284612,'amu*angstrom^2'), symmetry=1, barrier=(19.8166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110606,'amu*angstrom^2'), symmetry=1, barrier=(7.70139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.683839,'amu*angstrom^2'), symmetry=1, barrier=(47.6119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.887198,0.0731356,-7.66235e-05,4.62265e-08,-1.18111e-11,38461.3,29.6855], Tmin=(100,'K'), Tmax=(927.594,'K')), NASAPolynomial(coeffs=[9.46151,0.0361616,-1.68338e-05,3.25578e-09,-2.30007e-13,36870.6,-11.0363], Tmin=(927.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCOOH) + radical(Cds_P)"""),
)

species(
    label = 'C[CH]C=[C]COO(2790)',
    structure = SMILES('C[CH]C=[C]COO'),
    E0 = (244.776,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1685,370,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,338.607,338.608],'cm^-1')),
        HinderedRotor(inertia=(0.517019,'amu*angstrom^2'), symmetry=1, barrier=(42.0662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00239765,'amu*angstrom^2'), symmetry=1, barrier=(10.651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216349,'amu*angstrom^2'), symmetry=1, barrier=(17.6034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0558685,'amu*angstrom^2'), symmetry=1, barrier=(4.54573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.51701,'amu*angstrom^2'), symmetry=1, barrier=(42.0662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40048,0.0612656,-4.80273e-05,2.1081e-08,-4.07764e-12,29529.9,25.1031], Tmin=(100,'K'), Tmax=(1155.47,'K')), NASAPolynomial(coeffs=[8.2536,0.0375414,-1.72289e-05,3.3113e-09,-2.32922e-13,27946.2,-8.94956], Tmin=(1155.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]CCOO(2791)',
    structure = SMILES('[CH]C=CCCOO'),
    E0 = (229.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.985892,0.0679546,-5.00701e-05,1.98777e-08,-3.39956e-12,27747.1,27.2735], Tmin=(100,'K'), Tmax=(1307.33,'K')), NASAPolynomial(coeffs=[10.2515,0.0396049,-1.75424e-05,3.29044e-09,-2.27601e-13,25324.5,-19.911], Tmin=(1307.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C=CCCO[O](1853)',
    structure = SMILES('[CH2]C=CCCO[O]'),
    E0 = (162.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,237.338,237.416],'cm^-1')),
        HinderedRotor(inertia=(0.215518,'amu*angstrom^2'), symmetry=1, barrier=(8.52564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213208,'amu*angstrom^2'), symmetry=1, barrier=(8.52501,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.799865,'amu*angstrom^2'), symmetry=1, barrier=(31.3474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00294395,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11993,0.063438,-5.17591e-05,2.31361e-08,-4.3601e-12,19663.9,26.7264], Tmin=(100,'K'), Tmax=(1228.64,'K')), NASAPolynomial(coeffs=[10.6475,0.0324195,-1.38897e-05,2.58791e-09,-1.78997e-13,17322.7,-21.2005], Tmin=(1228.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC=C[CH]OO(2792)',
    structure = SMILES('[CH2]CC=C[CH]OO'),
    E0 = (188.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,1077.12,1077.19],'cm^-1')),
        HinderedRotor(inertia=(1.02757,'amu*angstrom^2'), symmetry=1, barrier=(23.6259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225645,'amu*angstrom^2'), symmetry=1, barrier=(5.18802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02797,'amu*angstrom^2'), symmetry=1, barrier=(23.6351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02835,'amu*angstrom^2'), symmetry=1, barrier=(23.6438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0499988,'amu*angstrom^2'), symmetry=1, barrier=(41.2508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.823166,0.0593004,-2.96073e-05,-3.22558e-09,4.99561e-12,22778.3,30.1917], Tmin=(100,'K'), Tmax=(1076.77,'K')), NASAPolynomial(coeffs=[14.9554,0.0263459,-1.09258e-05,2.0646e-09,-1.46336e-13,18601.8,-44.2946], Tmin=(1076.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJO) + radical(RCCJ)"""),
)

species(
    label = 'C=CC[CH]CO[O](2793)',
    structure = SMILES('C=CC[CH]CO[O]'),
    E0 = (223.795,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,308.387,308.387,308.388],'cm^-1')),
        HinderedRotor(inertia=(0.00177262,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10416,'amu*angstrom^2'), symmetry=1, barrier=(7.02903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104155,'amu*angstrom^2'), symmetry=1, barrier=(7.02905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104153,'amu*angstrom^2'), symmetry=1, barrier=(7.02905,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16407,0.0667565,-6.73373e-05,4.1775e-08,-1.12862e-11,27014.7,29.0463], Tmin=(100,'K'), Tmax=(871.071,'K')), NASAPolynomial(coeffs=[7.53177,0.0375165,-1.69868e-05,3.24072e-09,-2.27037e-13,25905.3,-0.795429], Tmin=(871.071,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.795,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCJCOOH)"""),
)

species(
    label = 'C[CH]C=C[CH]OO(2794)',
    structure = SMILES('C[CH]C=C[CH]OO'),
    E0 = (124.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3050,390,425,1340,1360,335,370,462.364,1194.9],'cm^-1')),
        HinderedRotor(inertia=(0.0244086,'amu*angstrom^2'), symmetry=1, barrier=(24.736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07575,'amu*angstrom^2'), symmetry=1, barrier=(24.7336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0244131,'amu*angstrom^2'), symmetry=1, barrier=(24.735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07551,'amu*angstrom^2'), symmetry=1, barrier=(24.7282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.69265,'amu*angstrom^2'), symmetry=1, barrier=(84.9013,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01131,0.0530507,-9.33707e-06,-2.48315e-08,1.26504e-11,15060,26.6216], Tmin=(100,'K'), Tmax=(1031.53,'K')), NASAPolynomial(coeffs=[14.9457,0.0265758,-1.09131e-05,2.08655e-09,-1.50299e-13,10719,-48.1433], Tmin=(1031.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]CC=CCO[O](2795)',
    structure = SMILES('[CH2]CC=CCO[O]'),
    E0 = (223.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,312.477,1766.85],'cm^-1')),
        HinderedRotor(inertia=(0.123721,'amu*angstrom^2'), symmetry=1, barrier=(8.53813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124474,'amu*angstrom^2'), symmetry=1, barrier=(8.54487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123412,'amu*angstrom^2'), symmetry=1, barrier=(8.5404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124626,'amu*angstrom^2'), symmetry=1, barrier=(8.54893,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44294,0.0600444,-5.05484e-05,2.56907e-08,-5.84317e-12,26918.4,28.1359], Tmin=(100,'K'), Tmax=(1001.1,'K')), NASAPolynomial(coeffs=[7.08743,0.0374909,-1.67547e-05,3.18602e-09,-2.23079e-13,25788.2,0.898368], Tmin=(1001.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ) + radical(RCCJ)"""),
)

species(
    label = 'CC=C[CH]CO[O](1834)',
    structure = SMILES('CC=C[CH]CO[O]'),
    E0 = (128.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,208.544,208.564],'cm^-1')),
        HinderedRotor(inertia=(0.366342,'amu*angstrom^2'), symmetry=1, barrier=(11.3051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.28311,'amu*angstrom^2'), symmetry=1, barrier=(70.4728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.799706,'amu*angstrom^2'), symmetry=1, barrier=(24.6833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.366262,'amu*angstrom^2'), symmetry=1, barrier=(11.3051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04295,0.0681206,-6.15e-05,3.14634e-08,-6.88003e-12,15505,23.9297], Tmin=(100,'K'), Tmax=(1064.62,'K')), NASAPolynomial(coeffs=[9.60267,0.0359597,-1.61863e-05,3.08761e-09,-2.16608e-13,13682.4,-17.9018], Tmin=(1064.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJCO) + radical(ROOJ)"""),
)

species(
    label = 'CC=CC=COO(2796)',
    structure = SMILES('CC=CC=COO'),
    E0 = (-26.4022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.587381,0.0650119,-4.18743e-05,3.63527e-09,4.39987e-12,-3043.74,25.9518], Tmin=(100,'K'), Tmax=(1004.13,'K')), NASAPolynomial(coeffs=[16.515,0.0219846,-8.10458e-06,1.46826e-09,-1.03169e-13,-7271.93,-56.0818], Tmin=(1004.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.4022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH)"""),
)

species(
    label = 'C=CC[C]COO(2797)',
    structure = SMILES('C=CC[C]COO'),
    E0 = (298.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31811,0.0701921,-4.39538e-05,-4.75583e-08,6.93882e-11,36035,34.2204], Tmin=(100,'K'), Tmax=(485.358,'K')), NASAPolynomial(coeffs=[7.02746,0.0428677,-2.04781e-05,3.94276e-09,-2.75355e-13,35248.5,8.40941], Tmin=(485.358,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=C[C]CCOO(2798)',
    structure = SMILES('C=C[C]CCOO'),
    E0 = (298.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31811,0.0701921,-4.39538e-05,-4.75583e-08,6.93882e-11,36035,34.2204], Tmin=(100,'K'), Tmax=(485.358,'K')), NASAPolynomial(coeffs=[7.02746,0.0428677,-2.04781e-05,3.94276e-09,-2.75355e-13,35248.5,8.40941], Tmin=(485.358,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C[C]C=CCOO(2799)',
    structure = SMILES('C[C]C=CCOO'),
    E0 = (293.351,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.755959,0.0777737,-0.00010017,8.43739e-08,-3.02104e-11,35392.6,36.6764], Tmin=(100,'K'), Tmax=(767.711,'K')), NASAPolynomial(coeffs=[5.91867,0.0440572,-2.09724e-05,4.03381e-09,-2.81562e-13,34800.8,14.4425], Tmin=(767.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]CC=CCOO(2800)',
    structure = SMILES('[CH]CC=CCOO'),
    E0 = (315.787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51102,0.0650393,-4.00786e-05,-4.09396e-08,6.01307e-11,38059.5,24.9648], Tmin=(100,'K'), Tmax=(484.807,'K')), NASAPolynomial(coeffs=[6.41033,0.0416056,-2.01384e-05,3.92197e-09,-2.76398e-13,37384.9,2.81654], Tmin=(484.807,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'OOCC1C=CC1(2801)',
    structure = SMILES('OOCC1C=CC1'),
    E0 = (24.3929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24744,0.0465434,7.58328e-06,-4.22633e-08,1.89796e-11,3045.23,23.9525], Tmin=(100,'K'), Tmax=(1001.66,'K')), NASAPolynomial(coeffs=[14.5872,0.0258037,-1.00748e-05,1.91309e-09,-1.388e-13,-1259.09,-48.5726], Tmin=(1001.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.3929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
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
    label = 'Ar(8)',
    structure = SMILES('[Ar]'),
    E0 = (-6.19426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,-745,4.3663], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,-745,4.3663], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.19426,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: BurkeH2O2"""),
)

transitionState(
    label = 'TS1',
    E0 = (158.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (178.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (340.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (402.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (306.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (425.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (463.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (391.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (368.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (435.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (225.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (331.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (332.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (372.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (382.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (253.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (254.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (180.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (196.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (248.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (132.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (248.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (153.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (139.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (316.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (320.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (319.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (333.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (98.6778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction2',
    reactants = ['OH(5)(5)', 'C5H7O(216)(215)'],
    products = ['S(214)(213)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 103 used for O_pri_rad;O_rad/NonDe
Exact match found for rate rule [O_pri_rad;O_rad/NonDe]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction1',
    reactants = ['HO2(8)(9)', 'C5H7(211)(210)'],
    products = ['S(214)(213)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.74095e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cs_rad;O_rad/NonDe] for rate rule [C_rad/H2/Cd;O_rad/NonDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)(3)', 'S(213)(212)'],
    products = ['S(214)(213)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2OOH(33)(34)', 'CH2CHCHCH(1116)'],
    products = ['S(214)(213)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;Cd_rad] for rate rule [C_rad/H2/O;Cd_pri_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)(3)', 'C=CC=C[CH]OO(1526)'],
    products = ['S(214)(213)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;H_rad] for rate rule [C_rad/H/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)(3)', 'C=CC=[C]COO(1527)'],
    products = ['S(214)(213)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H3(28)(29)', '[CH]=CCOO(2784)'],
    products = ['S(214)(213)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 89 used for Cd_pri_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)(3)', 'C=C[C]=CCOO(1528)'],
    products = ['S(214)(213)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)(3)', 'C=[C]C=CCOO(1529)'],
    products = ['S(214)(213)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)(3)', '[CH]=CC=CCOO(1530)'],
    products = ['S(214)(213)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C[CH]C[CH]OO(2785)'],
    products = ['S(214)(213)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.64935e+10,'s^-1'), n=0.2551, Ea=(25.8154,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C[C]=CCOO(2786)'],
    products = ['S(214)(213)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]C[CH]COO(2787)'],
    products = ['S(214)(213)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CC=[C]COO(2788)'],
    products = ['S(214)(213)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC[CH]COO(2789)'],
    products = ['S(214)(213)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[CH]C=[C]COO(2790)'],
    products = ['S(214)(213)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C[CH]CCOO(2791)'],
    products = ['S(214)(213)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C=CCCO[O](1853)'],
    products = ['S(214)(213)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CC=C[CH]OO(2792)'],
    products = ['S(214)(213)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=CC[CH]CO[O](2793)'],
    products = ['S(214)(213)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C[CH]C=C[CH]OO(2794)'],
    products = ['S(214)(213)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]CC=CCO[O](2795)'],
    products = ['S(214)(213)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC=C[CH]CO[O](1834)'],
    products = ['S(214)(213)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['S(214)(213)'],
    products = ['CC=CC=COO(2796)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.8582e+08,'s^-1'), n=1.23767, Ea=(163.714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH_end;CdH2_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=CC[C]COO(2797)'],
    products = ['S(214)(213)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;CsJ2C;CH2(C=C)] for rate rule [CsJ2-C;CsJ2C;CH2(C=C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[C]CCOO(2798)'],
    products = ['S(214)(213)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.65514e+12,'s^-1'), n=-0.428961, Ea=(21.7426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2(C=C);CH2(C)] + [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C[C]C=CCOO(2799)'],
    products = ['S(214)(213)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.00341e+13,'s^-1'), n=-1.27142, Ea=(26.2576,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]CC=CCOO(2800)'],
    products = ['S(214)(213)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;singletcarbene;CH2(C=C)] for rate rule [CsJ2-C;CsJ2H;CH2(C=C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['S(214)(213)'],
    products = ['OOCC1C=CC1(2801)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH2_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

network(
    label = '106',
    isomers = [
        'S(214)(213)',
    ],
    reactants = [
        ('OH(5)(5)', 'C5H7O(216)(215)'),
        ('HO2(8)(9)', 'C5H7(211)(210)'),
        ('H(3)(3)', 'S(213)(212)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '106',
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

