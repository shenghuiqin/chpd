species(
    label = 'S(227)(226)',
    structure = SMILES('C=CC=CC(O)O[O]'),
    E0 = (-78.6532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.60723,'amu*angstrom^2'), symmetry=1, barrier=(13.9614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.606212,'amu*angstrom^2'), symmetry=1, barrier=(13.938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.606969,'amu*angstrom^2'), symmetry=1, barrier=(13.9554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.606708,'amu*angstrom^2'), symmetry=1, barrier=(13.9494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4426.56,'J/mol'), sigma=(7.01777,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.42 K, Pc=29.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.475664,0.0774889,-8.0698e-05,4.44867e-08,-9.89095e-12,-9332.79,29.8388], Tmin=(100,'K'), Tmax=(1084.84,'K')), NASAPolynomial(coeffs=[13.9574,0.0277788,-1.19634e-05,2.24672e-09,-1.56684e-13,-12257.9,-36.3003], Tmin=(1084.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.6532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'O2(2)(2)',
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
    label = 'C5H7O(224)(223)',
    structure = SMILES('[CH2]C=CC=CO'),
    E0 = (-32.8413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.31979,'amu*angstrom^2'), symmetry=1, barrier=(30.3447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32385,'amu*angstrom^2'), symmetry=1, barrier=(30.4378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31983,'amu*angstrom^2'), symmetry=1, barrier=(30.3455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3794.25,'J/mol'), sigma=(6.17562,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.65 K, Pc=36.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25554,0.0396093,3.48643e-05,-9.07575e-08,4.31029e-11,-3831.66,21.0706], Tmin=(100,'K'), Tmax=(908.805,'K')), NASAPolynomial(coeffs=[21.7542,0.00426328,2.6289e-06,-6.6847e-10,4.32897e-14,-9823.72,-88.3316], Tmin=(908.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.8413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ)"""),
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
    label = 'C5H6O(217)(216)',
    structure = SMILES('C=CC=CC=O'),
    E0 = (-29.5668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.9508,'amu*angstrom^2'), symmetry=1, barrier=(21.8608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.95303,'amu*angstrom^2'), symmetry=1, barrier=(21.912,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3519.98,'J/mol'), sigma=(5.63814,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.81 K, Pc=44.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58677,0.0449815,-2.33922e-05,-1.17435e-10,2.4605e-12,-3462.46,19.7432], Tmin=(100,'K'), Tmax=(1159.43,'K')), NASAPolynomial(coeffs=[12.6828,0.020301,-9.05764e-06,1.75766e-09,-1.25367e-13,-6949.62,-39.3724], Tmin=(1159.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.5668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=C[CH]C1OOC1O(7350)',
    structure = SMILES('C=C[CH]C1OOC1O'),
    E0 = (-55.2269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.499151,0.0708451,-5.71247e-05,2.30532e-08,-3.74359e-12,-6511.35,24.1797], Tmin=(100,'K'), Tmax=(1452.32,'K')), NASAPolynomial(coeffs=[16.4836,0.0268203,-1.16546e-05,2.18071e-09,-1.50634e-13,-11154.3,-58.9009], Tmin=(1452.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.2269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C1C=CC(O)OO1(7351)',
    structure = SMILES('[CH2]C1C=CC(O)OO1'),
    E0 = (-66.8177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.940765,0.0473136,2.28984e-05,-6.99474e-08,3.13536e-11,-7908.21,25.92], Tmin=(100,'K'), Tmax=(973.244,'K')), NASAPolynomial(coeffs=[19.9696,0.0176509,-6.20437e-06,1.23907e-09,-9.73149e-14,-13911.3,-77.1789], Tmin=(973.244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.8177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro12dioxin) + radical(CJCOOH)"""),
)

species(
    label = 'C=CC=C[C](O)OO(7352)',
    structure = SMILES('[CH2]C=CC=C(O)OO'),
    E0 = (-65.4831,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3615,1277.5,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.802397,'amu*angstrom^2'), symmetry=1, barrier=(18.4487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.98572,'amu*angstrom^2'), symmetry=1, barrier=(68.6475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.802445,'amu*angstrom^2'), symmetry=1, barrier=(18.4498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.974965,'amu*angstrom^2'), symmetry=1, barrier=(22.4164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.802195,'amu*angstrom^2'), symmetry=1, barrier=(18.444,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.1195,0.0829667,-8.5835e-05,4.50809e-08,-9.25714e-12,-7720.98,29.8344], Tmin=(100,'K'), Tmax=(1194.49,'K')), NASAPolynomial(coeffs=[18.5577,0.0204223,-7.29391e-06,1.24576e-09,-8.2715e-14,-12182.9,-63.5917], Tmin=(1194.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.4831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=CC([O])OO(7353)',
    structure = SMILES('C=CC=CC([O])OO'),
    E0 = (-4.95273,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,281.056,281.058,281.061],'cm^-1')),
        HinderedRotor(inertia=(0.208929,'amu*angstrom^2'), symmetry=1, barrier=(11.7113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208917,'amu*angstrom^2'), symmetry=1, barrier=(11.7113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.390435,'amu*angstrom^2'), symmetry=1, barrier=(21.8859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608441,'amu*angstrom^2'), symmetry=1, barrier=(34.106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.825757,0.0751621,-7.62569e-05,4.33477e-08,-1.04449e-11,-485.824,28.698], Tmin=(100,'K'), Tmax=(976.692,'K')), NASAPolynomial(coeffs=[10.1013,0.0371744,-1.79153e-05,3.52489e-09,-2.51554e-13,-2297.68,-15.8322], Tmin=(976.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.95273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = 'C=CC=[C]C(O)OO(7354)',
    structure = SMILES('C=CC=[C]C(O)OO'),
    E0 = (7.18391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.289035,0.0842883,-9.55147e-05,5.71196e-08,-1.37804e-11,995.492,30.2301], Tmin=(100,'K'), Tmax=(1002.14,'K')), NASAPolynomial(coeffs=[14.0008,0.0295585,-1.35954e-05,2.62336e-09,-1.85416e-13,-1752.73,-35.9505], Tmin=(1002.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.18391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=C[C]=CC(O)OO(7355)',
    structure = SMILES('C=C[C]=CC(O)OO'),
    E0 = (-31.6624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.833653,'amu*angstrom^2'), symmetry=1, barrier=(19.1673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833514,'amu*angstrom^2'), symmetry=1, barrier=(19.1641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.28261,'amu*angstrom^2'), symmetry=1, barrier=(52.4818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833667,'amu*angstrom^2'), symmetry=1, barrier=(19.1676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833644,'amu*angstrom^2'), symmetry=1, barrier=(19.1671,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.148032,0.0873692,-0.000103822,6.59637e-08,-1.68434e-11,-3671.62,29.8587], Tmin=(100,'K'), Tmax=(952.08,'K')), NASAPolynomial(coeffs=[13.8466,0.0298173,-1.31498e-05,2.47313e-09,-1.71897e-13,-6280.05,-35.5563], Tmin=(952.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.6624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=CC(O)OO(7356)',
    structure = SMILES('C=C=C[CH]C(O)OO'),
    E0 = (-47.5799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21998,'amu*angstrom^2'), symmetry=1, barrier=(28.0497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22004,'amu*angstrom^2'), symmetry=1, barrier=(28.051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21993,'amu*angstrom^2'), symmetry=1, barrier=(28.0486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22001,'amu*angstrom^2'), symmetry=1, barrier=(28.0503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21994,'amu*angstrom^2'), symmetry=1, barrier=(28.0489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.188145,0.0857514,-9.10099e-05,4.99889e-08,-1.10947e-11,-5586.8,28.452], Tmin=(100,'K'), Tmax=(1082.36,'K')), NASAPolynomial(coeffs=[15.1394,0.030498,-1.44378e-05,2.82619e-09,-2.01471e-13,-8823.4,-44.863], Tmin=(1082.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.5799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=CC=CC(O)OO(7357)',
    structure = SMILES('[CH]=CC=CC(O)OO'),
    E0 = (16.4383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3120,650,792.5,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.822094,'amu*angstrom^2'), symmetry=1, barrier=(18.9016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.822887,'amu*angstrom^2'), symmetry=1, barrier=(18.9198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.822338,'amu*angstrom^2'), symmetry=1, barrier=(18.9072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.822125,'amu*angstrom^2'), symmetry=1, barrier=(18.9023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.822649,'amu*angstrom^2'), symmetry=1, barrier=(18.9143,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.156604,0.0843671,-9.17345e-05,5.12235e-08,-1.14038e-11,2115.7,30.6292], Tmin=(100,'K'), Tmax=(1088.84,'K')), NASAPolynomial(coeffs=[15.9618,0.0263037,-1.17441e-05,2.24669e-09,-1.58439e-13,-1326.1,-46.9665], Tmin=(1088.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.4383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
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
    label = 'C=CC=C[CH]O[O](2605)',
    structure = SMILES('[CH2]C=CC=CO[O]'),
    E0 = (243.658,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.920477,'amu*angstrom^2'), symmetry=1, barrier=(21.1636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92045,'amu*angstrom^2'), symmetry=1, barrier=(21.163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.920435,'amu*angstrom^2'), symmetry=1, barrier=(21.1626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19261,0.0517996,-2.08373e-05,-1.61991e-08,1.20077e-11,29415.5,26.0962], Tmin=(100,'K'), Tmax=(937.741,'K')), NASAPolynomial(coeffs=[15.1311,0.0171574,-5.11526e-06,8.4127e-10,-5.79576e-14,25710.4,-46.0702], Tmin=(937.741,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CC=CCJ)"""),
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
    label = 'C=CC=CC([O])O[O](7130)',
    structure = SMILES('C=CC=CC([O])O[O]'),
    E0 = (147.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,238.978,239.718,239.956],'cm^-1')),
        HinderedRotor(inertia=(0.267201,'amu*angstrom^2'), symmetry=1, barrier=(10.852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265643,'amu*angstrom^2'), symmetry=1, barrier=(10.8528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26736,'amu*angstrom^2'), symmetry=1, barrier=(10.852,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.887179,0.0739964,-8.9668e-05,6.44041e-08,-1.94919e-11,17793.4,29.2337], Tmin=(100,'K'), Tmax=(793.237,'K')), NASAPolynomial(coeffs=[8.51287,0.0355434,-1.69548e-05,3.2939e-09,-2.324e-13,16583.6,-5.78956], Tmin=(793.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = 'CH2CHCHCH(913)',
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
    label = '[O]O[CH]O(7358)',
    structure = SMILES('[O]O[CH]O'),
    E0 = (9.28844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,492.5,1135,1000],'cm^-1')),
        HinderedRotor(inertia=(0.156407,'amu*angstrom^2'), symmetry=1, barrier=(3.59612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156804,'amu*angstrom^2'), symmetry=1, barrier=(3.60523,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58857,0.0452524,-0.00011197,1.22556e-07,-4.65581e-11,1154.51,14.2974], Tmin=(100,'K'), Tmax=(891.504,'K')), NASAPolynomial(coeffs=[-0.0459,0.0212294,-1.12418e-05,2.13324e-09,-1.41938e-13,3048.61,34.6933], Tmin=(891.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.28844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = 'C=CC=C[C](O)O[O](7359)',
    structure = SMILES('[CH2]C=CC=C(O)O[O]'),
    E0 = (86.5216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3615,1277.5,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.599637,'amu*angstrom^2'), symmetry=1, barrier=(13.7868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.600111,'amu*angstrom^2'), symmetry=1, barrier=(13.7977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13495,'amu*angstrom^2'), symmetry=1, barrier=(26.0948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.599616,'amu*angstrom^2'), symmetry=1, barrier=(13.7864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.249628,0.0781669,-8.65039e-05,4.94876e-08,-1.10922e-11,10544.9,29.2669], Tmin=(100,'K'), Tmax=(1096.33,'K')), NASAPolynomial(coeffs=[16.1613,0.0201118,-7.07173e-06,1.18503e-09,-7.74178e-14,7056.06,-48.961], Tmin=(1096.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.5216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=[C]C(O)O[O](7360)',
    structure = SMILES('C=CC=[C]C(O)O[O]'),
    E0 = (159.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.569127,'amu*angstrom^2'), symmetry=1, barrier=(13.0854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569458,'amu*angstrom^2'), symmetry=1, barrier=(13.093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569151,'amu*angstrom^2'), symmetry=1, barrier=(13.0859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569176,'amu*angstrom^2'), symmetry=1, barrier=(13.0865,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.459432,0.0817952,-0.000104066,7.1487e-08,-1.97592e-11,19270,30.3777], Tmin=(100,'K'), Tmax=(881.376,'K')), NASAPolynomial(coeffs=[12.2814,0.0281428,-1.27554e-05,2.42014e-09,-1.68517e-13,17186.1,-25.1635], Tmin=(881.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(ROOJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH]=CC(O)O[O](7361)',
    structure = SMILES('[CH]=CC(O)O[O]'),
    E0 = (116.671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,492.5,1135,1000],'cm^-1')),
        HinderedRotor(inertia=(0.459796,'amu*angstrom^2'), symmetry=1, barrier=(10.5716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459704,'amu*angstrom^2'), symmetry=1, barrier=(10.5695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459674,'amu*angstrom^2'), symmetry=1, barrier=(10.5688,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80913,0.0543703,-6.66677e-05,3.28693e-08,1.44191e-12,14105.2,22.8668], Tmin=(100,'K'), Tmax=(576,'K')), NASAPolynomial(coeffs=[8.244,0.0211563,-1.00494e-05,1.91807e-09,-1.32744e-13,13173.6,-6.27994], Tmin=(576,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(ROOJ)"""),
)

species(
    label = 'C=C[C]=CC(O)O[O](7362)',
    structure = SMILES('C=C[C]=CC(O)O[O]'),
    E0 = (120.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.711232,'amu*angstrom^2'), symmetry=1, barrier=(16.3526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712078,'amu*angstrom^2'), symmetry=1, barrier=(16.3721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.711266,'amu*angstrom^2'), symmetry=1, barrier=(16.3534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.711587,'amu*angstrom^2'), symmetry=1, barrier=(16.3608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.289161,0.0852885,-0.000114123,8.30173e-08,-2.41438e-11,14604.1,30.1067], Tmin=(100,'K'), Tmax=(841.914,'K')), NASAPolynomial(coeffs=[12.3456,0.028009,-1.20742e-05,2.2126e-09,-1.50139e-13,12573.9,-25.9845], Tmin=(841.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(ROOJ)"""),
)

species(
    label = 'C=[C]C=CC(O)O[O](7363)',
    structure = SMILES('C=C=C[CH]C(O)O[O]'),
    E0 = (104.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.880655,'amu*angstrom^2'), symmetry=1, barrier=(20.248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.880537,'amu*angstrom^2'), symmetry=1, barrier=(20.2453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.880741,'amu*angstrom^2'), symmetry=1, barrier=(20.25,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6583,'amu*angstrom^2'), symmetry=1, barrier=(38.1275,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.418222,0.0824927,-9.66286e-05,6.02512e-08,-1.5212e-11,12685.3,28.3901], Tmin=(100,'K'), Tmax=(957.724,'K')), NASAPolynomial(coeffs=[13.0544,0.0297163,-1.39687e-05,2.71149e-09,-1.91973e-13,10264.9,-32.0263], Tmin=(957.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(ROOJ)"""),
)

species(
    label = '[CH]=CC=CC(O)O[O](7364)',
    structure = SMILES('[CH]=CC=CC(O)O[O]'),
    E0 = (168.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,492.5,1135,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.642974,'amu*angstrom^2'), symmetry=1, barrier=(14.7832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.642923,'amu*angstrom^2'), symmetry=1, barrier=(14.7821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.642757,'amu*angstrom^2'), symmetry=1, barrier=(14.7782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64296,'amu*angstrom^2'), symmetry=1, barrier=(14.7829,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.414976,0.0808008,-9.63884e-05,6.03653e-08,-1.5088e-11,20386.5,30.464], Tmin=(100,'K'), Tmax=(974.617,'K')), NASAPolynomial(coeffs=[13.8935,0.0254812,-1.12461e-05,2.1242e-09,-1.48232e-13,17759.3,-34.2155], Tmin=(974.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(ROOJ)"""),
)

species(
    label = 'C=CC1[CH]C(O)OO1(7323)',
    structure = SMILES('C=CC1[CH]C(O)OO1'),
    E0 = (-56.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03075,0.047502,2.30603e-05,-7.57959e-08,3.67907e-11,-6731.23,29.2158], Tmin=(100,'K'), Tmax=(906.734,'K')), NASAPolynomial(coeffs=[19.0555,0.0152985,-1.93271e-06,1.24922e-10,-7.80083e-15,-11944.9,-66.7032], Tmin=(906.734,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxolane) + radical(CCJCOOH)"""),
)

species(
    label = 'OC1[CH]C=CCOO1(7305)',
    structure = SMILES('OC1C=C[CH]COO1'),
    E0 = (-134.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.89044,0.0160801,0.000182109,-2.88263e-07,1.23453e-10,-15965.4,25.0312], Tmin=(100,'K'), Tmax=(911.314,'K')), NASAPolynomial(coeffs=[40.3276,-0.0204779,1.75379e-05,-3.46134e-09,2.21254e-13,-28823.2,-192.676], Tmin=(911.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptane) + radical(C=CCJCO)"""),
)

species(
    label = 'CC=CC=C(O)O[O](7365)',
    structure = SMILES('CC=CC=C(O)O[O]'),
    E0 = (-31.5342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.594335,'amu*angstrom^2'), symmetry=1, barrier=(13.6649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.595335,'amu*angstrom^2'), symmetry=1, barrier=(13.6879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594405,'amu*angstrom^2'), symmetry=1, barrier=(13.6665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594252,'amu*angstrom^2'), symmetry=1, barrier=(13.663,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.182477,0.084482,-0.000100396,6.29376e-08,-1.56534e-11,-3655.61,27.9562], Tmin=(100,'K'), Tmax=(983.323,'K')), NASAPolynomial(coeffs=[14.7002,0.0254279,-1.03149e-05,1.86629e-09,-1.27029e-13,-6510.79,-41.8394], Tmin=(983.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.5342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + radical(ROOJ)"""),
)

species(
    label = 'C=CC=C=CO(6203)',
    structure = SMILES('C=CC=C=CO'),
    E0 = (25.7073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.24386,'amu*angstrom^2'), symmetry=1, barrier=(28.5988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24416,'amu*angstrom^2'), symmetry=1, barrier=(28.6057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02364,0.048058,2.85443e-06,-5.66469e-08,3.08127e-11,3215.32,19.3983], Tmin=(100,'K'), Tmax=(909.916,'K')), NASAPolynomial(coeffs=[22.142,0.00160986,2.95307e-06,-6.91045e-10,4.50216e-14,-2548.21,-91.0442], Tmin=(909.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.7073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[O]OC(O)C1C=CC1(7366)',
    structure = SMILES('[O]OC(O)C1C=CC1'),
    E0 = (-23.8334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.864061,0.0575622,-2.31389e-05,-1.22537e-08,8.90783e-12,-2743.64,29.2466], Tmin=(100,'K'), Tmax=(1022.21,'K')), NASAPolynomial(coeffs=[15.4596,0.0243787,-9.55986e-06,1.79125e-09,-1.27901e-13,-6977.82,-47.6042], Tmin=(1022.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.8334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(ROOJ)"""),
)

species(
    label = 'O(4)(4)',
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
    label = 'C=CC=CC([O])O(7367)',
    structure = SMILES('C=CC=CC([O])O'),
    E0 = (-76.4577,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.830092,'amu*angstrom^2'), symmetry=1, barrier=(19.0855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.831697,'amu*angstrom^2'), symmetry=1, barrier=(19.1224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.831731,'amu*angstrom^2'), symmetry=1, barrier=(19.1231,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37943,0.0598699,-5.1679e-05,2.45689e-08,-4.95755e-12,-9103.3,25.9347], Tmin=(100,'K'), Tmax=(1147.37,'K')), NASAPolynomial(coeffs=[9.64425,0.0310565,-1.40099e-05,2.68154e-09,-1.8848e-13,-10999.8,-15.0744], Tmin=(1147.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.4577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ)"""),
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
    E0 = (-38.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (43.3898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-62.1961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (109.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (100.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (66.9498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1.37779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (71.9725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (167.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (273.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (358.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (355.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (298.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (370.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (403.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (336.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (320.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (380.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (10.0679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (0.273836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (132.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (44.5658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (100.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (43.7598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (166.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)(2)', 'C5H7O(224)(223)'],
    products = ['S(227)(226)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.19343e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;Y_rad] for rate rule [C_rad/H/OneDeO;O2_birad]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction2',
    reactants = ['S(227)(226)'],
    products = ['C=C[CH]C1OOC1O(7350)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_O] for rate rule [R5_SS_D;doublebond_intra_HCd_pri;radadd_intra_O]
Euclidian distance = 3.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['S(227)(226)'],
    products = ['[CH2]C1C=CC(O)OO1(7351)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(17945.7,'s^-1'), n=1.45333, Ea=(16.4571,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SSSR;doublebond_intra_2H_pri;radadd_intra] for rate rule [R7_SSSM_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=CC=C[C](O)OO(7352)'],
    products = ['S(227)(226)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.47099e+07,'s^-1'), n=1.47622, Ea=(175.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_O;C_rad_out_OneDe/O;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=CC=CC([O])OO(7353)'],
    products = ['S(227)(226)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(94.0113,'s^-1'), n=2.81534, Ea=(105.641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['S(227)(226)'],
    products = ['C=CC=[C]C(O)OO(7354)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C[C]=CC(O)OO(7355)'],
    products = ['S(227)(226)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]C=CC(O)OO(7356)'],
    products = ['S(227)(226)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62236e+06,'s^-1'), n=1.59783, Ea=(119.552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R6H_SMSSR;Y_rad_out;XH_out] for rate rule [R6H_SMSSR;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CC=CC(O)OO(7357)'],
    products = ['S(227)(226)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7H;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OH(5)(5)', 'C=CC=C[CH]O[O](2605)'],
    products = ['S(227)(226)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.81233e+06,'m^3/(mol*s)'), n=0.00568686, Ea=(1.481,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;C_rad/H/OneDe] + [O_rad;C_sec_rad] for rate rule [O_pri_rad;C_rad/H/OneDeO]
Euclidian distance = 2.2360679775
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)(3)', 'C=CC=CC([O])O[O](7130)'],
    products = ['S(227)(226)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['CH2CHCHCH(913)', '[O]O[CH]O(7358)'],
    products = ['S(227)(226)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.18e+06,'m^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/NonDeO;Y_rad] for rate rule [C_rad/H/O2;Cd_pri_rad]
Euclidian distance = 2.2360679775
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)(3)', 'C=CC=C[C](O)O[O](7359)'],
    products = ['S(227)(226)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""Estimated using template [C_rad/OneDe;H_rad] for rate rule [C_rad/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)(3)', 'C=CC=[C]C(O)O[O](7360)'],
    products = ['S(227)(226)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C2H3(28)(29)', '[CH]=CC(O)O[O](7361)'],
    products = ['S(227)(226)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 89 used for Cd_pri_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)(3)', 'C=C[C]=CC(O)O[O](7362)'],
    products = ['S(227)(226)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)(3)', 'C=[C]C=CC(O)O[O](7363)'],
    products = ['S(227)(226)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)(3)', '[CH]=CC=CC(O)O[O](7364)'],
    products = ['S(227)(226)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['S(227)(226)'],
    products = ['C=CC1[CH]C(O)OO1(7323)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.81589e+09,'s^-1'), n=0.403324, Ea=(88.7211,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_linear;doublebond_intra_pri_HCd;radadd_intra] + [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['S(227)(226)'],
    products = ['OC1[CH]C=CCOO1(7305)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.76476e+07,'s^-1'), n=0.815689, Ea=(78.927,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC=CC=C(O)O[O](7365)'],
    products = ['S(227)(226)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.02873e+09,'s^-1'), n=1.23767, Ea=(163.714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH3_1;unsaturated_end]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction22',
    reactants = ['S(227)(226)'],
    products = ['HO2(8)(9)', 'C5H6O(217)(216)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.813e+10,'s^-1','*|/',10), n=0.493, Ea=(123.219,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 13 used for R2OO_O_HNd
Exact match found for rate rule [R2OO_O_HNd]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction23',
    reactants = ['S(227)(226)'],
    products = ['HO2(8)(9)', 'C=CC=C=CO(6203)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction24',
    reactants = ['S(227)(226)'],
    products = ['[O]OC(O)C1C=CC1(7366)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH2_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(4)(4)', 'C=CC=CC([O])O(7367)'],
    products = ['S(227)(226)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '602',
    isomers = [
        'S(227)(226)',
    ],
    reactants = [
        ('O2(2)(2)', 'C5H7O(224)(223)'),
        ('HO2(8)(9)', 'C5H6O(217)(216)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '602',
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

