species(
    label = 'C7H6O(490)(489)',
    structure = SMILES('C1=CC2=CC=C(C1)O2'),
    E0 = (41.2019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4154,'J/mol'), sigma=(6.51034,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=648.85 K, Pc=34.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20282,0.0266771,3.98414e-05,-5.77521e-08,1.86939e-11,5030.4,17.6585], Tmin=(100,'K'), Tmax=(1194.95,'K')), NASAPolynomial(coeffs=[12.8261,0.0324323,-1.92456e-05,4.14704e-09,-3.09659e-13,-458.206,-47.8272], Tmin=(1194.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.2019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(Furan)"""),
)

species(
    label = 'C1=CC23CC2C=C1O3(5820)',
    structure = SMILES('C1=CC23CC2C=C1O3'),
    E0 = (404.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04665,0.0511015,-6.05254e-06,-2.6159e-08,1.19741e-11,48755.5,15.5777], Tmin=(100,'K'), Tmax=(1099.77,'K')), NASAPolynomial(coeffs=[16.2131,0.0259106,-1.25724e-05,2.57316e-09,-1.90234e-13,43607,-67.2744], Tmin=(1099.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s3_5_5_diene_0_4) + polycyclic(s1_3_5_ene_1) + polycyclic(s2_3_5_ene_1) - ring(2,3-Dihydrofuran) - ring(Cyclopropane) - ring(25dihydrofuran)"""),
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
    label = 'C1=C[C]2C=CC(=C1)O2(5776)',
    structure = SMILES('[CH]1C=CC2=CC=C1O2'),
    E0 = (158.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.114,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97064,0.0308161,3.1354e-05,-5.28177e-08,1.77798e-11,19101.5,14.9964], Tmin=(100,'K'), Tmax=(1187.38,'K')), NASAPolynomial(coeffs=[14.8064,0.0286413,-1.77762e-05,3.89431e-09,-2.93431e-13,13158.4,-61.3232], Tmin=(1187.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(Furan) + radical(C=CCJCO)"""),
)

species(
    label = '[C]1=CC2=CC=C(C1)O2(5821)',
    structure = SMILES('[C]1=CC2=CC=C(C1)O2'),
    E0 = (279.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.114,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26644,0.02998,2.03047e-05,-3.62658e-08,1.14361e-11,33629.8,17.9149], Tmin=(100,'K'), Tmax=(1287.56,'K')), NASAPolynomial(coeffs=[13.1449,0.0296638,-1.83303e-05,3.93332e-09,-2.90346e-13,28053.3,-48.094], Tmin=(1287.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(Furan) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CCC2=CC=C1O2(5822)',
    structure = SMILES('[C]1=CCC2=CC=C1O2'),
    E0 = (240.197,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.114,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19905,0.0322078,1.49001e-05,-3.106e-08,9.86646e-12,28959.5,17.2785], Tmin=(100,'K'), Tmax=(1304.15,'K')), NASAPolynomial(coeffs=[12.9019,0.0300912,-1.79878e-05,3.8083e-09,-2.78965e-13,23556.3,-47.2116], Tmin=(1304.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(Furan) + radical(C=CJC=C)"""),
)

species(
    label = '[C]1C=C2C=CCC=1O2(5823)',
    structure = SMILES('[C]1C=C2C=CCC=1O2'),
    E0 = (240.197,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.114,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19905,0.0322078,1.49001e-05,-3.106e-08,9.86646e-12,28959.5,17.2785], Tmin=(100,'K'), Tmax=(1304.15,'K')), NASAPolynomial(coeffs=[12.9019,0.0300912,-1.79878e-05,3.8083e-09,-2.78965e-13,23556.3,-47.2116], Tmin=(1304.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(Furan) + radical(C=CJC=C)"""),
)

species(
    label = '[C]1C=C2CC=CC=1O2(5824)',
    structure = SMILES('[C]1C=C2CC=CC=1O2'),
    E0 = (240.197,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.114,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19905,0.0322078,1.49001e-05,-3.106e-08,9.86646e-12,28959.5,17.2785], Tmin=(100,'K'), Tmax=(1304.15,'K')), NASAPolynomial(coeffs=[12.9019,0.0300912,-1.79878e-05,3.8083e-09,-2.78965e-13,23556.3,-47.2116], Tmin=(1304.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(Furan) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]1C=CC2=C[CH]C1O2(5750)',
    structure = SMILES('[CH]1C=C[C]2C=CC1O2'),
    E0 = (291.745,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31636,0.0114208,0.000107973,-1.42957e-07,5.27226e-11,35172.1,20.4515], Tmin=(100,'K'), Tmax=(997.225,'K')), NASAPolynomial(coeffs=[15.1774,0.0251985,-1.1072e-05,2.35617e-09,-1.84923e-13,29356.9,-57.856], Tmin=(997.225,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(C=CCJC(O)C=C) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH]1CC=C2C=C[C]1O2(5802)',
    structure = SMILES('[CH]1CC=C2[CH]C=C1O2'),
    E0 = (332.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,303.741,808.702,808.702,808.702,808.702,808.702,808.702,808.702,808.702,808.702,808.702,808.702,1621.26,1621.26,1621.26,1621.26,1621.26,1621.26,1621.26,1621.26,1621.26,1621.26,1621.26,1621.26],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56426,0.0179482,5.76215e-05,-7.24768e-08,2.34538e-11,40022.1,18.8423], Tmin=(100,'K'), Tmax=(1139.15,'K')), NASAPolynomial(coeffs=[9.62729,0.0344441,-1.84782e-05,3.88277e-09,-2.88255e-13,35733.4,-27.9139], Tmin=(1139.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(2,3-Dihydrofuran) + radical(C=CCJCO) + radical(CCJCO)"""),
)

species(
    label = '[C]1=CC[C]2C=CC1O2(5790)',
    structure = SMILES('[C]1=CC[C]2C=CC1O2'),
    E0 = (459.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02119,0.0301439,3.40852e-05,-5.56259e-08,1.97924e-11,55340.7,21.1943], Tmin=(100,'K'), Tmax=(1089.73,'K')), NASAPolynomial(coeffs=[10.7697,0.0319612,-1.51202e-05,3.04885e-09,-2.23309e-13,51419.4,-31.0076], Tmin=(1089.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(Cds_S) + radical(C2CsJOCs)"""),
)

species(
    label = '[C]1C[C]2C=CCC=1O2(5825)',
    structure = SMILES('[C]1CC2=C[CH]CC=1O2'),
    E0 = (394.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,2471.73,4000,4000,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85425,0.0172023,4.66258e-05,-5.47403e-08,1.62628e-11,47479.8,18.1788], Tmin=(100,'K'), Tmax=(1232.08,'K')), NASAPolynomial(coeffs=[7.20116,0.0376011,-2.02245e-05,4.16595e-09,-3.02783e-13,43789.2,-14.3299], Tmin=(1232.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(2,3-Dihydrofuran) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[C]1C[C]2CC=CC=1O2(5826)',
    structure = SMILES('[C]1C[C]2CC=CC=1O2'),
    E0 = (400.898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35586,0.0444142,1.98477e-06,-3.08761e-08,1.27994e-11,48323.6,17.5008], Tmin=(100,'K'), Tmax=(1123.72,'K')), NASAPolynomial(coeffs=[15.9387,0.0239283,-1.2615e-05,2.67047e-09,-2.00232e-13,43062.2,-63.3817], Tmin=(1123.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(2,3-Dihydrofuran) + radical(C2CsJOC(O)) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CC2=C[CH]C(C1)O2(5827)',
    structure = SMILES('[C]1=C[C]2C=CC(C1)O2'),
    E0 = (459.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02119,0.0301439,3.40852e-05,-5.56259e-08,1.97924e-11,55340.7,21.1943], Tmin=(100,'K'), Tmax=(1089.73,'K')), NASAPolynomial(coeffs=[10.7697,0.0319612,-1.51202e-05,3.04885e-09,-2.23309e-13,51419.4,-31.0076], Tmin=(1089.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(Cds_S) + radical(C2CsJOCs)"""),
)

species(
    label = '[C]1C=C2C[CH]CC=1O2(5828)',
    structure = SMILES('[C]1C=C2C[CH]CC=1O2'),
    E0 = (331.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,304.188,812.235,812.235,812.235,812.235,812.235,812.235,812.235,812.235,812.235,812.235,812.235,1622.91,1622.91,1622.91,1622.91,1622.91,1622.91,1622.91,1622.91,1622.91,1622.91,1622.91,1622.91],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93353,0.0220537,3.1861e-05,-3.78759e-08,1.0315e-11,39933.2,20.9867], Tmin=(100,'K'), Tmax=(1377.85,'K')), NASAPolynomial(coeffs=[7.24358,0.0404703,-2.1859e-05,4.40759e-09,-3.12976e-13,35809.6,-11.8422], Tmin=(1377.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(Furan) + radical(C=CJC=C) + radical(RCCJCC)"""),
)

species(
    label = '[C]1=CC2C=C[C](C1)O2(5786)',
    structure = SMILES('[C]1=CC2C=C[C](C1)O2'),
    E0 = (459.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02119,0.0301439,3.40852e-05,-5.56259e-08,1.97924e-11,55340.7,21.1943], Tmin=(100,'K'), Tmax=(1089.73,'K')), NASAPolynomial(coeffs=[10.7697,0.0319612,-1.51202e-05,3.04885e-09,-2.23309e-13,51419.4,-31.0076], Tmin=(1089.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(C2CsJOCs) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CC2=CC[C](C1)O2(5829)',
    structure = SMILES('[C]1=CC2=CC[C](C1)O2'),
    E0 = (400.898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35586,0.0444142,1.98477e-06,-3.08761e-08,1.27994e-11,48323.6,17.5008], Tmin=(100,'K'), Tmax=(1123.72,'K')), NASAPolynomial(coeffs=[15.9387,0.0239283,-1.2615e-05,2.67047e-09,-2.00232e-13,43062.2,-63.3817], Tmin=(1123.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(2,3-Dihydrofuran) + radical(C2CsJOC(O)) + radical(Cds_S)"""),
)

species(
    label = '[C]1C=C2CC[CH]C=1O2(5830)',
    structure = SMILES('[C]1C=C2CC[CH]C=1O2'),
    E0 = (337.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.57405,0.0202395,4.93329e-05,-6.17642e-08,1.92928e-11,40609.2,22.1866], Tmin=(100,'K'), Tmax=(1183.88,'K')), NASAPolynomial(coeffs=[8.71811,0.0371843,-1.99086e-05,4.12863e-09,-3.02646e-13,36512.2,-19.6513], Tmin=(1183.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(Furan) + radical(CCJCO) + radical(C=CJC=C)"""),
)

species(
    label = '[C]1=C[C]2CC=C(C1)O2(5831)',
    structure = SMILES('[C]1=C[C]2CC=C(C1)O2'),
    E0 = (418.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76102,0.0322328,3.47751e-05,-6.29838e-08,2.34122e-11,50429.6,20.3515], Tmin=(100,'K'), Tmax=(1080.71,'K')), NASAPolynomial(coeffs=[15.3988,0.0241972,-1.29798e-05,2.81443e-09,-2.15232e-13,45003.4,-57.9685], Tmin=(1080.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(2,3-Dihydrofuran) + radical(C2CsJOC(O)) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CC2C=CC[C]1O2(5803)',
    structure = SMILES('[C]1=CC2C=CC[C]1O2'),
    E0 = (459.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02119,0.0301439,3.40852e-05,-5.56259e-08,1.97924e-11,55340.7,21.1943], Tmin=(100,'K'), Tmax=(1089.73,'K')), NASAPolynomial(coeffs=[10.7697,0.0319612,-1.51202e-05,3.04885e-09,-2.23309e-13,51419.4,-31.0076], Tmin=(1089.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(C2CsJOCs) + radical(Cds_S)"""),
)

species(
    label = '[C]1C=C2[CH]CCC=1O2(5832)',
    structure = SMILES('[C]1C=C2[CH]CCC=1O2'),
    E0 = (337.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.57405,0.0202395,4.93329e-05,-6.17642e-08,1.92928e-11,40609.2,22.1866], Tmin=(100,'K'), Tmax=(1183.88,'K')), NASAPolynomial(coeffs=[8.71811,0.0371843,-1.99086e-05,4.12863e-09,-3.02646e-13,36512.2,-19.6513], Tmin=(1183.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(Furan) + radical(C=CJC=C) + radical(CCJCO)"""),
)

species(
    label = '[CH]1C[C]2C=CC=C1O2(5782)',
    structure = SMILES('[CH]1C=CC2=CC[C]1O2'),
    E0 = (279.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16823,0.0438885,1.81541e-05,-5.44821e-08,2.22806e-11,33790.7,14.8937], Tmin=(100,'K'), Tmax=(1061.32,'K')), NASAPolynomial(coeffs=[18.6187,0.021426,-1.1305e-05,2.46894e-09,-1.908e-13,27647.6,-81.8235], Tmin=(1061.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(2,3-Dihydrofuran) + radical(C2CsJOC(O)) + radical(C=CCJCO)"""),
)

species(
    label = '[C]1[CH]CC2C=CC=1O2(5800)',
    structure = SMILES('[C]1=CCC2C=C[C]1O2'),
    E0 = (459.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02119,0.0301439,3.40852e-05,-5.56259e-08,1.97924e-11,55340.7,21.1943], Tmin=(100,'K'), Tmax=(1089.73,'K')), NASAPolynomial(coeffs=[10.7697,0.0319612,-1.51202e-05,3.04885e-09,-2.23309e-13,51419.4,-31.0076], Tmin=(1089.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(C2CsJOCs) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CC2C[CH]C=C1O2(5797)',
    structure = SMILES('[C]1=CC2CC=C[C]1O2'),
    E0 = (459.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02119,0.0301439,3.40852e-05,-5.56259e-08,1.97924e-11,55340.7,21.1943], Tmin=(100,'K'), Tmax=(1089.73,'K')), NASAPolynomial(coeffs=[10.7697,0.0319612,-1.51202e-05,3.04885e-09,-2.23309e-13,51419.4,-31.0076], Tmin=(1089.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(Cds_S) + radical(C2CsJOCs)"""),
)

species(
    label = '[C]1=CC[C]2CC=C1O2(5833)',
    structure = SMILES('[C]1=CC[C]2CC=C1O2'),
    E0 = (362.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.287,0.0466684,-3.54915e-06,-2.54603e-08,1.11273e-11,43653.3,16.869], Tmin=(100,'K'), Tmax=(1122.26,'K')), NASAPolynomial(coeffs=[15.5222,0.0246128,-1.24063e-05,2.57469e-09,-1.91134e-13,38652,-61.4966], Tmin=(1122.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(2,3-Dihydrofuran) + radical(C2CsJOC(O)) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]1C=C2C=C[C](C1)O2(5783)',
    structure = SMILES('C1=C[C]2C=C[C](C1)O2'),
    E0 = (402.305,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97461,0.0358572,1.03846e-05,-2.75166e-08,9.44083e-12,48466,21.5381], Tmin=(100,'K'), Tmax=(1209.37,'K')), NASAPolynomial(coeffs=[9.6968,0.0338035,-1.62004e-05,3.19765e-09,-2.28154e-13,44880.6,-24.2863], Tmin=(1209.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(C2CsJOCs) + radical(C2CsJOCs)"""),
)

species(
    label = '[C]1=CCC2=CC[C]1O2(5834)',
    structure = SMILES('[C]1[CH]CC2=CCC=1O2'),
    E0 = (394.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,1102.18,2471.73,4000,4000,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85425,0.0172023,4.66258e-05,-5.47403e-08,1.62628e-11,47479.8,18.1788], Tmin=(100,'K'), Tmax=(1232.08,'K')), NASAPolynomial(coeffs=[7.20116,0.0376011,-2.02245e-05,4.16595e-09,-3.02783e-13,43789.2,-14.3299], Tmin=(1232.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(2,3-Dihydrofuran) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'C1=CC2=CCCC=1O2(5835)',
    structure = SMILES('C1=CC2=CCCC=1O2'),
    E0 = (416.685,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13465,0.0356104,1.00714e-05,-2.65699e-08,8.30809e-12,50186.7,16.6997], Tmin=(100,'K'), Tmax=(1364.08,'K')), NASAPolynomial(coeffs=[13.832,0.0317127,-1.90751e-05,4.01426e-09,-2.91762e-13,44166.9,-53.7331], Tmin=(1364.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(1,2-Cyclopentadiene)"""),
)

species(
    label = 'C1=CC2=CCC(=C1)O2(5836)',
    structure = SMILES('C1=CC2=CCC(=C1)O2'),
    E0 = (96.3021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0718,0.0290701,3.36211e-05,-5.3567e-08,1.79983e-11,11662.7,16.3737], Tmin=(100,'K'), Tmax=(1170.54,'K')), NASAPolynomial(coeffs=[13.1437,0.0300046,-1.72583e-05,3.70649e-09,-2.77322e-13,6414.65,-50.1304], Tmin=(1170.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.3021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(1,3-Cyclohexadiene) + ring(2,3-Dihydrofuran)"""),
)

species(
    label = 'C7H6O(489)(488)',
    structure = SMILES('C1=CC2C=CC(=C1)O2'),
    E0 = (65.5427,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4079.36,'J/mol'), sigma=(6.42564,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=637.19 K, Pc=34.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37721,0.0309341,6.52369e-05,-1.06142e-07,3.97479e-11,8000.39,16.533], Tmin=(100,'K'), Tmax=(1058.94,'K')), NASAPolynomial(coeffs=[22.6207,0.0182415,-1.24709e-05,3.01959e-09,-2.45304e-13,-286.211,-105.054], Tmin=(1058.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.5427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(1,3-Cyclohexadiene) + ring(25dihydrofuran)"""),
)

species(
    label = 'C1=CCC2C=CC=1O2(5837)',
    structure = SMILES('C1=CCC2C=CC=1O2'),
    E0 = (149.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.107632,0.0945118,-8.57857e-05,3.38133e-08,-5.23806e-12,18155.6,5.75732], Tmin=(100,'K'), Tmax=(1506.45,'K')), NASAPolynomial(coeffs=[25.3045,0.0270361,-1.85988e-05,4.08029e-09,-3.03763e-13,10499.2,-127.254], Tmin=(1506.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(six-inringtwodouble-12) + ring(25dihydrofuran)"""),
)

species(
    label = '[CH]=C1C=CCC(=[CH])O1(5838)',
    structure = SMILES('[CH]=C1C=CCC(=[CH])O1'),
    E0 = (510.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3115,3125,620,680,785,800,1600,1700,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34878,0.0393886,3.54817e-05,-7.36014e-08,2.99665e-11,61472.3,21.8419], Tmin=(100,'K'), Tmax=(1004.08,'K')), NASAPolynomial(coeffs=[16.599,0.024317,-1.02467e-05,2.07175e-09,-1.56793e-13,56107.1,-63.2607], Tmin=(1004.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[O]C1=CC=[C]C=CC1(5766)',
    structure = SMILES('[O]C1=CC=[C]C=CC1'),
    E0 = (287.745,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.269515,0.0682813,-3.44407e-05,-1.92435e-08,1.69793e-11,34754.8,21.5267], Tmin=(100,'K'), Tmax=(918.739,'K')), NASAPolynomial(coeffs=[21.4963,0.0121483,-2.03359e-06,2.27115e-10,-1.59173e-14,29323.1,-87.4149], Tmin=(918.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]C1C=CC[C]=CC=1(2905)',
    structure = SMILES('[O]C1C=CC[C]=CC=1'),
    E0 = (327.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.246637,0.070671,-4.70103e-05,-1.6633e-09,9.24222e-12,39505.9,21.6702], Tmin=(100,'K'), Tmax=(943.141,'K')), NASAPolynomial(coeffs=[20.3521,0.014511,-3.98979e-06,6.53317e-10,-4.65847e-14,34418.8,-81.0139], Tmin=(943.141,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC1=CC=C([CH2])O1(5839)',
    structure = SMILES('[CH]=CC1=CC=C([CH2])O1'),
    E0 = (363.694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3120,650,792.5,1650,2750,3150,900,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14062,0.0531449,-2.41922e-05,-5.14945e-09,5.12553e-12,43853.6,25.3761], Tmin=(100,'K'), Tmax=(1080.88,'K')), NASAPolynomial(coeffs=[13.7232,0.0248694,-1.0333e-05,1.95664e-09,-1.38793e-13,40065.2,-41.2478], Tmin=(1080.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Furan) + radical(Cds_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH]=CCC1=CC=[C]O1(5840)',
    structure = SMILES('[CH]=CCC1=CC=[C]O1'),
    E0 = (480.476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2750,3150,900,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65957,0.0420678,-3.6208e-06,-1.86215e-08,7.95532e-12,57880,29.5971], Tmin=(100,'K'), Tmax=(1130.61,'K')), NASAPolynomial(coeffs=[11.0346,0.0289431,-1.27997e-05,2.47076e-09,-1.7575e-13,54479.1,-22.4479], Tmin=(1130.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(480.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Furan) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=CC1=CC=[C]O1(5841)',
    structure = SMILES('[CH2]C=CC1=CC=[C]O1'),
    E0 = (321.247,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73641,0.0386022,1.12985e-05,-3.62551e-08,1.47393e-11,38728.3,26.6324], Tmin=(100,'K'), Tmax=(1042.59,'K')), NASAPolynomial(coeffs=[10.8725,0.0292036,-1.20866e-05,2.2976e-09,-1.64002e-13,35429,-24.5114], Tmin=(1042.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Furan) + radical(C=CC=CCJ) + radical(C=CJO)"""),
)

species(
    label = '[C]1CC2=CC=C(C1)O2(5844)',
    structure = SMILES('[C]1CC2=CC=C(C1)O2'),
    E0 = (364.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61103,0.0297682,2.07277e-05,-3.09964e-08,8.62986e-12,43906.2,27.2607], Tmin=(100,'K'), Tmax=(1423.07,'K')), NASAPolynomial(coeffs=[10.4795,0.0389679,-2.19789e-05,4.47451e-09,-3.18224e-13,38495.7,-24.6176], Tmin=(1423.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(364.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CsJ2_singlet-CsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Cyclohexane) - ring(Tetrahydrofuran) + ring(Cyclohexane) + ring(Furan)"""),
)

species(
    label = '[C]1CCC2=CC=C1O2(5845)',
    structure = SMILES('[C]1CCC2=CC=C1O2'),
    E0 = (293.612,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43916,0.0254141,3.57849e-05,-4.87463e-08,1.48995e-11,35376.4,21.483], Tmin=(100,'K'), Tmax=(1249.34,'K')), NASAPolynomial(coeffs=[10.2136,0.0364112,-2.0507e-05,4.28459e-09,-3.13207e-13,30633.1,-28.9638], Tmin=(1249.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CsJ2_singlet-(Cds-Cds-Cds-Cds)Cs_6_ring) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Cyclohexane) - ring(Tetrahydrofuran) + ring(Cyclohexane) + ring(Furan)"""),
)

species(
    label = '[C]1C=C2C=CCC1O2(5846)',
    structure = SMILES('[C]1C=C2C=CCC1O2'),
    E0 = (324.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4941,0.029894,6.81868e-05,-1.14393e-07,4.60139e-11,39192.2,15.7101], Tmin=(100,'K'), Tmax=(980.698,'K')), NASAPolynomial(coeffs=[20.2987,0.0158791,-6.25373e-06,1.38705e-09,-1.15363e-13,32489.5,-90.0135], Tmin=(980.698,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(CsJ2_singlet-(Cds-Cds-Cds-Cds)Cs_6_ring) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Cyclopentane) + ring(36dihydro2hpyran) + ring(Cyclopentane)"""),
)

species(
    label = '[C]1C=C2CC=CC1O2(5847)',
    structure = SMILES('[C]1C=C2CC=CC1O2'),
    E0 = (413.896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72995,0.0288937,6.02243e-05,-9.56982e-08,3.66538e-11,49880.4,25.334], Tmin=(100,'K'), Tmax=(1011.43,'K')), NASAPolynomial(coeffs=[15.9416,0.0248499,-1.11356e-05,2.32651e-09,-1.78785e-13,44337.6,-56.5797], Tmin=(1011.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Cyclopentane) + ring(36dihydro2hpyran) + ring(Cyclopentane)"""),
)

species(
    label = 'C1=C2CC3C=C(O2)C13(5848)',
    structure = SMILES('C1=C2CC3C=C(O2)C13'),
    E0 = (521.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2168,0.047027,-1.22272e-06,-2.91585e-08,1.23683e-11,62829.9,21.5014], Tmin=(100,'K'), Tmax=(1132.45,'K')), NASAPolynomial(coeffs=[16.9264,0.023673,-1.28534e-05,2.74574e-09,-2.06551e-13,57211.3,-65.3403], Tmin=(1132.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(521.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_5_ene_m) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Tetrahydrofuran) - ring(Cyclobutane) + ring(2,3-Dihydrofuran) + ring(Cyclobutene) + polycyclic(s2_4_5_diene_1_5) - ring(Cyclobutene) - ring(2,3-Dihydrofuran) - ring(Cyclopentene)"""),
)

species(
    label = 'C1=CC23C=CC2(C1)O3(5849)',
    structure = SMILES('C1=CC23C=CC2(C1)O3'),
    E0 = (284.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37516,0.0365909,3.95867e-05,-7.60886e-08,2.93589e-11,34280.2,21.6958], Tmin=(100,'K'), Tmax=(1057.61,'K')), NASAPolynomial(coeffs=[19.1796,0.0210667,-1.18834e-05,2.67929e-09,-2.10567e-13,27616.3,-78.8973], Tmin=(1057.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + polycyclic(s2_3_5_ene_1) + polycyclic(s2_3_4_ene_1) - ring(Cyclopentene) - ring(Ethylene_oxide) - ring(Cyclobutene)"""),
)

species(
    label = '[CH]1C[C]2OC3=CC2C13(5850)',
    structure = SMILES('[CH]1C[C]2OC3=CC2C13'),
    E0 = (696.238,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4652,0.0512791,-1.99311e-05,-3.58161e-09,2.33583e-12,83832.7,15.7767], Tmin=(100,'K'), Tmax=(1506.22,'K')), NASAPolynomial(coeffs=[17.8555,0.0273934,-1.57042e-05,3.20501e-09,-2.27568e-13,76667.3,-77.4064], Tmin=(1506.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(696.238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_5_5_ene_1) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + polycyclic(s3_4_5_ene_3) - ring(Cyclobutene) - ring(Cyclopentane) - ring(2,3-Dihydrofuran) + radical(C2CsJOC(O)) + radical(Cs_S)"""),
)

species(
    label = '[CH]1C2=CC3[C](CC13)O2(5851)',
    structure = SMILES('[CH]1C2=CC3C[C](O2)C13'),
    E0 = (670.473,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44747,0.0390495,2.29138e-05,-5.38148e-08,2.06514e-11,80745.6,20.1636], Tmin=(100,'K'), Tmax=(1094.68,'K')), NASAPolynomial(coeffs=[16.9839,0.0238511,-1.32256e-05,2.88652e-09,-2.20835e-13,74853.2,-67.5734], Tmin=(1094.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(670.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s3_5_5_ene_0) + polycyclic(s2_4_5_ene_1) + polycyclic(s2_4_5_ane) - ring(Tetrahydrofuran) - ring(Cyclopentene) - ring(Cyclobutane) + radical(C2CsJOC(O)) + radical(CCJCO)"""),
)

species(
    label = '[CH]1CC2=CC3[C](O2)C13(5852)',
    structure = SMILES('[CH]1CC2=CC3[C](O2)C13'),
    E0 = (442.992,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29512,0.033728,1.75653e-05,-3.39396e-08,1.01473e-11,53343.4,18.0839], Tmin=(100,'K'), Tmax=(1393.96,'K')), NASAPolynomial(coeffs=[16.3778,0.0318534,-2.18849e-05,4.75962e-09,-3.5001e-13,45673.2,-67.964], Tmin=(1393.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(2,3-Dihydrofuran) + polycyclic(s2_3_6_ane) + polycyclic(s2_3_5_ene_1) - ring(Oxane) - ring(Cyclopropane) - ring(2,3-Dihydrofuran) + radical(Cs_S) + radical(C2CsJOC(O))"""),
)

species(
    label = 'O=C1CC=CC2=CC12(5842)',
    structure = SMILES('O=C1CC=CC2=CC12'),
    E0 = (247.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52566,0.0408833,1.1268e-05,-4.19555e-08,1.79496e-11,29840.7,18.1543], Tmin=(100,'K'), Tmax=(1022.22,'K')), NASAPolynomial(coeffs=[13.8954,0.0239074,-9.9374e-06,1.94942e-09,-1.434e-13,25669.8,-49.8256], Tmin=(1022.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_3_6_ane) - ring(Cyclohexane) - ring(Cyclopropane) + ring(Cyclohexene) + ring(Cyclopropene)"""),
)

species(
    label = 'O=C1C=CCC2=CC12(5843)',
    structure = SMILES('O=C1C=CCC2=CC12'),
    E0 = (246.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0322,0.0305474,3.17137e-05,-5.35763e-08,1.95226e-11,29770,19.7581], Tmin=(100,'K'), Tmax=(1065.16,'K')), NASAPolynomial(coeffs=[10.1062,0.0316433,-1.40707e-05,2.76902e-09,-2.00884e-13,26267.7,-28.0698], Tmin=(1065.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + polycyclic(s2_3_6_diene_2_5)"""),
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
    E0 = (408.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (369.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (490.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (456.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (456.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (456.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (319.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (360.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (487.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (417.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (423.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (552.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (431.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (522.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (461.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (376.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (426.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (484.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (354.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (340.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (477.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (484.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (422.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (427.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (428.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (601.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (211.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (111.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (313.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (517.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (294.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (334.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (372.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (488.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (328.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (381.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (315.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (351.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (440.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (524.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (285.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (846.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (670.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (442.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (269.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (355.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction150',
    reactants = ['C7H6O(490)(489)'],
    products = ['C1=CC23CC2C=C1O3(5820)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.67327e+12,'s^-1'), n=-0.101958, Ea=(366.829,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_5_unsaturated_hexane] for rate rule [linear_1_3_5_hexatriene]
Euclidian distance = 1.0
family: Intra_Diels_alder_monocyclic"""),
)

reaction(
    label = 'reaction151',
    reactants = ['H(3)(3)', 'C1=C[C]2C=CC(=C1)O2(5776)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.28e-13,'cm^3/(molecule*s)'), n=0.611, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 53 used for C_rad/H/CdCd;H_rad
Exact match found for rate rule [C_rad/H/CdCd;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -1.8 to -1.8 kJ/mol.
Ea raised from -1.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction152',
    reactants = ['H(3)(3)', '[C]1=CC2=CC=C(C1)O2(5821)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction153',
    reactants = ['H(3)(3)', '[C]1=CCC2=CC=C1O2(5822)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction154',
    reactants = ['H(3)(3)', '[C]1C=C2C=CCC=1O2(5823)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [H_rad;Cd_rad/Cd]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction155',
    reactants = ['H(3)(3)', '[C]1C=C2CC=CC=1O2(5824)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]1C=CC2=C[CH]C1O2(5750)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction157',
    reactants = ['[CH]1CC=C2C=C[C]1O2(5802)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.89319e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction158',
    reactants = ['[C]1=CC[C]2C=CC1O2(5790)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction159',
    reactants = ['[C]1C[C]2C=CCC=1O2(5825)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction160',
    reactants = ['[C]1C[C]2CC=CC=1O2(5826)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction161',
    reactants = ['[C]1=CC2=C[CH]C(C1)O2(5827)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]1C=C2C[CH]CC=1O2(5828)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.94203e+09,'s^-1'), n=0.37, Ea=(99.6901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad_NDe] + [R3radExo;Y_rad;XH_Rrad_NDe] for rate rule [R3radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction163',
    reactants = ['[C]1=CC2C=C[C](C1)O2(5786)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction164',
    reactants = ['[C]1=CC2=CC[C](C1)O2(5829)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(60.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction165',
    reactants = ['[C]1C=C2CC[CH]C=1O2(5830)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction166',
    reactants = ['[C]1=C[C]2CC=C(C1)O2(5831)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction167',
    reactants = ['[C]1=CC2C=CC[C]1O2(5803)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction168',
    reactants = ['[C]1C=C2[CH]CCC=1O2(5832)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_De] for rate rule [R4radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]1C[C]2C=CC=C1O2(5782)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(60.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction170',
    reactants = ['[C]1[CH]CC2C=CC=1O2(5800)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_De] for rate rule [R4radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction171',
    reactants = ['[C]1=CC2C[CH]C=C1O2(5797)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction172',
    reactants = ['[C]1=CC[C]2CC=C1O2(5833)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(60.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction173',
    reactants = ['[CH]1C=C2C=C[C](C1)O2(5783)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction174',
    reactants = ['[C]1=CCC2=CC[C]1O2(5834)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_NDe] for rate rule [R6radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction175',
    reactants = ['C1=CC2=CCCC=1O2(5835)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(9.66e+08,'s^-1'), n=1.162, Ea=(185.28,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH2(C)_1;unsaturated_end] for rate rule [1_3_pentadiene;CH2(C)_1;CddC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction176',
    reactants = ['C1=CC2=CCC(=C1)O2(5836)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.0262e+09,'s^-1'), n=1.0765, Ea=(115.616,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_pentadiene;CH_end;CdHC_2] + [1_3_pentadiene;CH2(C)_1;unsaturated_end] for rate rule [1_3_pentadiene;CH2(C)_1;CdHC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction132',
    reactants = ['C7H6O(489)(488)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.125e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1_3_pentadiene;CH_end;CdHC_2]
Euclidian distance = 0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction177',
    reactants = ['C1=CCC2C=CC=1O2(5837)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(8.45351e+08,'s^-1'), n=1.02346, Ea=(163.761,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH_end;CddC_2] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH_end;CddC_2]
Euclidian distance = 1.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction178',
    reactants = ['[CH]=C1C=CCC(=[CH])O1(5838)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_DSSD;CdsingleH_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction179',
    reactants = ['[O]C1=CC=[C]C=CC1(5766)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SDSD;O_rad;CdsinglepriDe_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction180',
    reactants = ['[O]C1C=CC[C]=CC=1(2905)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SDSD;O_rad;CdsinglepriND_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction181',
    reactants = ['[CH]=CC1=CC=C([CH2])O1(5839)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;C_rad_out_2H;Ypri_rad_out] for rate rule [R6_SSSSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction182',
    reactants = ['[CH]=CCC1=CC=[C]O1(5840)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6;CdsingleH_rad_out;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction183',
    reactants = ['[CH2]C=CC1=CC=[C]O1(5841)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using an average for rate rule [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction186',
    reactants = ['[C]1CC2=CC=C(C1)O2(5844)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.36993e+12,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;CsJ2C;CH2(C=C)] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C=C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 4.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction187',
    reactants = ['[C]1CCC2=CC=C1O2(5845)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(5.65514e+12,'s^-1'), n=-0.428961, Ea=(21.7426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2(C=C);CH2(C)] + [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction188',
    reactants = ['[C]1C=C2C=CCC1O2(5846)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.33447e+13,'s^-1'), n=-1.27142, Ea=(26.2576,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [CsJ2-C;CsJ2(C=C);CH]
Euclidian distance = 0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction189',
    reactants = ['[C]1C=C2CC=CC1O2(5847)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.33447e+13,'s^-1'), n=-1.27142, Ea=(26.2576,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [CsJ2-C;CsJ2(C=C);CH]
Euclidian distance = 0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction190',
    reactants = ['C7H6O(490)(489)'],
    products = ['C1=C2CC3C=C(O2)C13(5848)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(483.075,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction191',
    reactants = ['C7H6O(490)(489)'],
    products = ['C1=CC23C=CC2(C1)O3(5849)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(244.528,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1,3-butadiene_backbone;C=C_1;C=C_2]
Euclidian distance = 0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction192',
    reactants = ['[CH]1C[C]2OC3=CC2C13(5850)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction193',
    reactants = ['[CH]1C2=CC3[C](CC13)O2(5851)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction194',
    reactants = ['[CH]1CC2=CC3[C](O2)C13(5852)'],
    products = ['C7H6O(490)(489)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R7JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction184',
    reactants = ['C7H6O(490)(489)'],
    products = ['O=C1CC=CC2=CC12(5842)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(5460.43,'s^-1'), n=2.73, Ea=(228.457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.2360679775
family: ketoenol"""),
)

reaction(
    label = 'reaction185',
    reactants = ['C7H6O(490)(489)'],
    products = ['O=C1C=CCC2=CC12(5843)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

network(
    label = '373',
    isomers = [
        'C7H6O(490)(489)',
    ],
    reactants = [
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '373',
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

