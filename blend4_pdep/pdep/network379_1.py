species(
    label = 'S(256)(255)',
    structure = SMILES('CC(=O)C(C=C[O])O[O]'),
    E0 = (-169.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,492.5,1135,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,303.487,303.491,303.494],'cm^-1')),
        HinderedRotor(inertia=(0.00183026,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24462,'amu*angstrom^2'), symmetry=1, barrier=(15.9894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244631,'amu*angstrom^2'), symmetry=1, barrier=(15.9894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244631,'amu*angstrom^2'), symmetry=1, barrier=(15.9894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4827.56,'J/mol'), sigma=(7.44261,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=754.05 K, Pc=26.57 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0203165,0.0783122,-6.92035e-05,2.73987e-08,-3.35436e-12,-20176,35.0699], Tmin=(100,'K'), Tmax=(1078.58,'K')), NASAPolynomial(coeffs=[19.3026,0.0204884,-8.03024e-06,1.48192e-09,-1.04115e-13,-25149.2,-63.3449], Tmin=(1078.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ)"""),
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
    label = 'S(247)(246)',
    structure = SMILES('CC(=O)C=CC=O'),
    E0 = (-249.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,338.299],'cm^-1')),
        HinderedRotor(inertia=(0.102894,'amu*angstrom^2'), symmetry=1, barrier=(8.35729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1029,'amu*angstrom^2'), symmetry=1, barrier=(8.35732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102901,'amu*angstrom^2'), symmetry=1, barrier=(8.35737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3947.01,'J/mol'), sigma=(6.07703,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=616.51 K, Pc=39.91 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72359,0.0430426,-2.26393e-05,4.45245e-09,-2.82809e-13,-30034.5,19.0288], Tmin=(100,'K'), Tmax=(2443.83,'K')), NASAPolynomial(coeffs=[29.6812,0.00674437,-5.16283e-06,9.95166e-10,-6.31701e-14,-45547.2,-139.895], Tmin=(2443.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'O2(S)(162)(161)',
    structure = SMILES('O=O'),
    E0 = (85.6848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1857.18,'J/mol'), sigma=(4.34667,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=290.09 K, Pc=51.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121571,5.31618e-06,-4.89443e-09,1.45845e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69971e-07,1.51275e-10,-1.08782e-14,10302.3,6.16754], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CC1([O])OOC1C=C[O](6413)',
    structure = SMILES('CC1([O])OOC1C=C[O]'),
    E0 = (-32.3986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0894349,0.0683932,-2.22035e-05,-3.13979e-08,1.9588e-11,-3740,28.2159], Tmin=(100,'K'), Tmax=(971.184,'K')), NASAPolynomial(coeffs=[22.8649,0.0160565,-5.41706e-06,1.04447e-09,-8.06832e-14,-10119.5,-91.0648], Tmin=(971.184,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.3986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxetane) + radical(C=COJ) + radical(CC(C)(O)OJ)"""),
)

species(
    label = 'CC(=O)C1OOC1[CH][O](6414)',
    structure = SMILES('CC(=O)C1OOC1[CH][O]'),
    E0 = (22.5367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.858411,0.0758373,-8.40328e-05,5.74119e-08,-1.71763e-11,2817.67,28.678], Tmin=(100,'K'), Tmax=(788.763,'K')), NASAPolynomial(coeffs=[7.39162,0.0427055,-2.10249e-05,4.15676e-09,-2.96782e-13,1787.05,-1.29061], Tmin=(788.763,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.5367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + ring(12dioxetane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'CC1([O])OC=CC1O[O](6415)',
    structure = SMILES('CC1([O])OC=CC1O[O]'),
    E0 = (-94.4207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.391626,0.0551486,2.59655e-05,-9.47701e-08,4.7006e-11,-11203.6,25.6907], Tmin=(100,'K'), Tmax=(911.135,'K')), NASAPolynomial(coeffs=[26.3412,0.00522773,2.78547e-06,-7.15557e-10,4.57102e-14,-18588.9,-111.664], Tmin=(911.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.4207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(ROOJ) + radical(CC(C)(O)OJ)"""),
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
    label = 'CC(=O)C(C=C=O)O[O](6416)',
    structure = SMILES('CC(=O)C(C=C=O)O[O]'),
    E0 = (-194.901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,552.5,462.5,1710,2120,512.5,787.5,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,183.087,183.088],'cm^-1')),
        HinderedRotor(inertia=(0.360068,'amu*angstrom^2'), symmetry=1, barrier=(8.56503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36007,'amu*angstrom^2'), symmetry=1, barrier=(8.56504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.360068,'amu*angstrom^2'), symmetry=1, barrier=(8.56503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.360068,'amu*angstrom^2'), symmetry=1, barrier=(8.56503,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0445393,0.0962449,-0.000155122,1.35636e-07,-4.6224e-11,-23307.6,31.3565], Tmin=(100,'K'), Tmax=(850.195,'K')), NASAPolynomial(coeffs=[9.34359,0.035157,-1.67554e-05,3.15246e-09,-2.14307e-13,-24262.2,-8.31187], Tmin=(850.195,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-194.901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + radical(ROOJ)"""),
)

species(
    label = 'CC(=O)C(C=[C]O)O[O](5953)',
    structure = SMILES('CC(=O)C(C=[C]O)O[O]'),
    E0 = (-70.7471,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,552.5,462.5,1710,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,385.743,385.806],'cm^-1')),
        HinderedRotor(inertia=(0.109175,'amu*angstrom^2'), symmetry=1, barrier=(11.5093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108973,'amu*angstrom^2'), symmetry=1, barrier=(11.5093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108786,'amu*angstrom^2'), symmetry=1, barrier=(11.5055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109043,'amu*angstrom^2'), symmetry=1, barrier=(11.5081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108999,'amu*angstrom^2'), symmetry=1, barrier=(11.5084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0386541,0.0849002,-9.2117e-05,5.05791e-08,-1.09062e-11,-8360.09,36.8533], Tmin=(100,'K'), Tmax=(1134.05,'K')), NASAPolynomial(coeffs=[17.8149,0.0219263,-8.82032e-06,1.61104e-09,-1.11096e-13,-12409.4,-51.5253], Tmin=(1134.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-70.7471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = 'CC(=O)[C](C=C[O])OO(6417)',
    structure = SMILES('CC([O])=C(C=C[O])OO'),
    E0 = (-175.586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,325,375,415,465,420,450,1700,1750,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.772913,'amu*angstrom^2'), symmetry=1, barrier=(17.7708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.772023,'amu*angstrom^2'), symmetry=1, barrier=(17.7503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.773237,'amu*angstrom^2'), symmetry=1, barrier=(17.7782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.773365,'amu*angstrom^2'), symmetry=1, barrier=(17.7812,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.712899,0.0959549,-0.00010832,5.81129e-08,-1.15007e-11,-20941.7,32.8293], Tmin=(100,'K'), Tmax=(1001.74,'K')), NASAPolynomial(coeffs=[22.2631,0.0156446,-5.1852e-06,8.69811e-10,-5.84036e-14,-26118.5,-80.92], Tmin=(1001.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-175.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC(=O)C([C]=CO)O[O](5952)',
    structure = SMILES('CC(=O)C([C]=CO)O[O]'),
    E0 = (-72.6495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,552.5,462.5,1710,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,284.855,284.861],'cm^-1')),
        HinderedRotor(inertia=(0.245674,'amu*angstrom^2'), symmetry=1, barrier=(14.1528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246053,'amu*angstrom^2'), symmetry=1, barrier=(14.1543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246345,'amu*angstrom^2'), symmetry=1, barrier=(14.1531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246057,'amu*angstrom^2'), symmetry=1, barrier=(14.1549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245506,'amu*angstrom^2'), symmetry=1, barrier=(14.152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.52352,0.0906086,-0.000100108,5.42437e-08,-1.13281e-11,-8567.25,35.9374], Tmin=(100,'K'), Tmax=(1182.98,'K')), NASAPolynomial(coeffs=[21.5392,0.016008,-5.51527e-06,9.36207e-10,-6.25306e-14,-13787.2,-74.21], Tmin=(1182.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-72.6495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)C([C]=C[O])OO(6418)',
    structure = SMILES('CC(=O)C([C]=C[O])OO'),
    E0 = (-83.1916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.232895,0.0854749,-8.54778e-05,4.20311e-08,-8.08818e-12,-9846.84,35.5494], Tmin=(100,'K'), Tmax=(1266.17,'K')), NASAPolynomial(coeffs=[20.2581,0.0207416,-8.79058e-06,1.65398e-09,-1.15954e-13,-15035.9,-68.1438], Tmin=(1266.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.1916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'CC(=O)[C](C=CO)O[O](5949)',
    structure = SMILES('CC([O])=C(C=CO)O[O]'),
    E0 = (-165.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.663511,'amu*angstrom^2'), symmetry=1, barrier=(15.2554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.663521,'amu*angstrom^2'), symmetry=1, barrier=(15.2557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.663217,'amu*angstrom^2'), symmetry=1, barrier=(15.2487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.66288,'amu*angstrom^2'), symmetry=1, barrier=(15.2409,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.75027,0.0979818,-0.000111451,5.43589e-08,-7.40179e-12,-19673,32.316], Tmin=(100,'K'), Tmax=(888.795,'K')), NASAPolynomial(coeffs=[23.4403,0.0111392,-2.06174e-06,1.91268e-10,-8.43846e-15,-24843.1,-86.4317], Tmin=(888.795,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(=O)C(C=C[O])OO(6419)',
    structure = SMILES('C=C([O])C(C=C[O])OO'),
    E0 = (-161.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,274.161,274.266,274.301,274.355],'cm^-1')),
        HinderedRotor(inertia=(0.401359,'amu*angstrom^2'), symmetry=1, barrier=(21.4151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401013,'amu*angstrom^2'), symmetry=1, barrier=(21.4153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401029,'amu*angstrom^2'), symmetry=1, barrier=(21.415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400919,'amu*angstrom^2'), symmetry=1, barrier=(21.416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.47839,0.0864741,-7.62318e-05,2.36588e-08,5.5462e-13,-19271.5,35.5645], Tmin=(100,'K'), Tmax=(988.481,'K')), NASAPolynomial(coeffs=[23.2244,0.0151412,-5.28994e-06,9.72516e-10,-7.0883e-14,-25158.5,-84.5882], Tmin=(988.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'CC(=O)C(C=[C][O])OO(6420)',
    structure = SMILES('CC(=O)C([CH][C]=O)OO'),
    E0 = (-103.288,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3025,407.5,1350,352.5,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.534813,0.079626,-8.42345e-05,4.74629e-08,-1.0963e-11,-12300.6,36.6002], Tmin=(100,'K'), Tmax=(1034.46,'K')), NASAPolynomial(coeffs=[12.8462,0.0320211,-1.52062e-05,2.9773e-09,-2.12101e-13,-14847.8,-23.2125], Tmin=(1034.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C(=O)C(C=CO)O[O](5951)',
    structure = SMILES('C=C([O])C(C=CO)O[O]'),
    E0 = (-151.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,492.5,1135,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.748245,'amu*angstrom^2'), symmetry=1, barrier=(17.2036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.748669,'amu*angstrom^2'), symmetry=1, barrier=(17.2134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.749542,'amu*angstrom^2'), symmetry=1, barrier=(17.2334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.750113,'amu*angstrom^2'), symmetry=1, barrier=(17.2466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19421,0.09652,-0.000107527,5.66813e-08,-1.11897e-11,-17973.4,37.4841], Tmin=(100,'K'), Tmax=(1385.61,'K')), NASAPolynomial(coeffs=[25.9458,0.00796529,-6.12001e-07,-7.52479e-11,9.81156e-15,-24514.7,-98.7667], Tmin=(1385.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-151.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ)"""),
)

species(
    label = 'CC([O])=CC=C[O](5355)',
    structure = SMILES('CC([O])=CC=C[O]'),
    E0 = (-86.1882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.912938,'amu*angstrom^2'), symmetry=1, barrier=(20.9902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.915391,'amu*angstrom^2'), symmetry=1, barrier=(21.0466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.756235,0.0562742,-1.51561e-05,-3.6247e-08,2.26938e-11,-10235.2,24.2724], Tmin=(100,'K'), Tmax=(917.521,'K')), NASAPolynomial(coeffs=[21.0957,0.00729401,3.01004e-08,-1.33448e-10,7.27657e-15,-15638.3,-81.2071], Tmin=(917.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.1882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CHCHO(45)(45)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (150.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.3215,0.0099565,8.23489e-06,-1.92986e-08,8.37464e-12,18150.6,9.47663], Tmin=(100,'K'), Tmax=(969.366,'K')), NASAPolynomial(coeffs=[7.73867,0.0039342,-1.33182e-06,2.69015e-10,-2.15546e-14,16720.8,-14.6539], Tmin=(969.366,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ketene(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CC(=O)[CH]O[O](5948)',
    structure = SMILES('CC(=O)[CH]O[O]'),
    E0 = (-19.5836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,180],'cm^-1')),
        HinderedRotor(inertia=(0.208148,'amu*angstrom^2'), symmetry=1, barrier=(4.78572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208468,'amu*angstrom^2'), symmetry=1, barrier=(4.7931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.815334,'amu*angstrom^2'), symmetry=1, barrier=(18.7461,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59074,0.0559261,-7.39611e-05,5.36886e-08,-1.56601e-11,-2271.27,20.4408], Tmin=(100,'K'), Tmax=(837.501,'K')), NASAPolynomial(coeffs=[9.23765,0.0194001,-8.53516e-06,1.60338e-09,-1.10847e-13,-3552,-15.0944], Tmin=(837.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.5836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(OCJC=O)"""),
)

species(
    label = 'CC(=O)C(O[O])C1[CH]O1(6421)',
    structure = SMILES('CC(=O)C(O[O])C1[CH]O1'),
    E0 = (-27.5304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.495918,0.0847456,-8.79714e-05,4.66108e-08,-9.3918e-12,-3137.16,35.2951], Tmin=(100,'K'), Tmax=(1374.42,'K')), NASAPolynomial(coeffs=[19.3399,0.017087,-3.29369e-06,2.80951e-10,-8.49027e-15,-7651.83,-63.2979], Tmin=(1374.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.5304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + ring(Ethylene_oxide) + radical(CCsJO) + radical(ROOJ)"""),
)

species(
    label = 'C[C]1OOOC1C=C[O](6422)',
    structure = SMILES('C[C]1OOOC1C=C[O]'),
    E0 = (100.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.381547,0.0604188,5.2015e-06,-6.88092e-08,3.71914e-11,12268.3,28.9704], Tmin=(100,'K'), Tmax=(897.528,'K')), NASAPolynomial(coeffs=[23.2637,0.0108036,6.0924e-07,-3.95851e-10,2.9406e-14,6051.75,-90.6981], Tmin=(897.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.771,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(123trioxolane) + radical(C=COJ) + radical(C2CsJOO)"""),
)

species(
    label = 'CC(=O)C1[CH]C([O])OO1(6423)',
    structure = SMILES('CC(=O)C1[CH]C([O])OO1'),
    E0 = (-70.0359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.999164,0.0537831,-5.52914e-06,-3.46641e-08,1.85982e-11,-8304.18,32.6398], Tmin=(100,'K'), Tmax=(940.556,'K')), NASAPolynomial(coeffs=[14.935,0.0248836,-7.86929e-06,1.32119e-09,-9.07448e-14,-12268.9,-40.8791], Tmin=(940.556,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-70.0359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + ring(12dioxolane) + radical(CCJCOOH) + radical(CCOJ)"""),
)

species(
    label = 'C[C]1OOC=CC1O[O](6424)',
    structure = SMILES('C[C]1OOC=CC1O[O]'),
    E0 = (162.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.537034,0.0669556,-4.28927e-05,5.26469e-09,3.68199e-12,19679.3,27.3924], Tmin=(100,'K'), Tmax=(997.326,'K')), NASAPolynomial(coeffs=[15.3688,0.0261075,-9.4881e-06,1.67331e-09,-1.14837e-13,15794,-48.7696], Tmin=(997.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(34dihydro12dioxin) + radical(C2CsJOOC) + radical(ROOJ)"""),
)

species(
    label = 'CC(=O)C(C=C=O)OO(6425)',
    structure = SMILES('CC(=O)C(C=C=O)OO'),
    E0 = (-346.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.171403,0.0992483,-0.0001483,1.23626e-07,-4.1407e-11,-41580.1,31.3745], Tmin=(100,'K'), Tmax=(799.267,'K')), NASAPolynomial(coeffs=[10.7971,0.0370576,-1.78896e-05,3.4278e-09,-2.37356e-13,-43100.3,-17.6266], Tmin=(799.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-346.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = 'CO(10)(11)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(762.44,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35308e-07,1.51269e-10,-9.88872e-15,-14292.7,6.51157], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'CC(C=C[O])O[O](6426)',
    structure = SMILES('CC(C=C[O])O[O]'),
    E0 = (-33.1488,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.775067,'amu*angstrom^2'), symmetry=1, barrier=(17.8203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766247,'amu*angstrom^2'), symmetry=1, barrier=(17.6175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.77262,'amu*angstrom^2'), symmetry=1, barrier=(17.764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.880669,0.0573995,-3.05022e-05,-9.21118e-09,9.69881e-12,-3864.49,27.3311], Tmin=(100,'K'), Tmax=(967.888,'K')), NASAPolynomial(coeffs=[17.4221,0.0149823,-4.97193e-06,8.97615e-10,-6.53378e-14,-8281.74,-58.2093], Tmin=(967.888,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.1488,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=COJ)"""),
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
    label = 'CC(=O)C1C=COO1(6427)',
    structure = SMILES('CC(=O)C1C=COO1'),
    E0 = (-172.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35083,0.0404807,2.6044e-05,-6.58895e-08,2.84389e-11,-20641.1,27.2469], Tmin=(100,'K'), Tmax=(981.359,'K')), NASAPolynomial(coeffs=[17.1687,0.0191068,-7.16295e-06,1.42109e-09,-1.08926e-13,-25821.1,-59.3413], Tmin=(981.359,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene)"""),
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
    label = 'CC(=O)C=[C]C=O(6097)',
    structure = SMILES('CC(=O)C=C=C[O]'),
    E0 = (-53.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.863481,'amu*angstrom^2'), symmetry=1, barrier=(19.8531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.863734,'amu*angstrom^2'), symmetry=1, barrier=(19.8589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18871,0.0546892,-4.44353e-05,1.77663e-08,-2.81993e-12,-6296.83,23.6411], Tmin=(100,'K'), Tmax=(1501.54,'K')), NASAPolynomial(coeffs=[15.0289,0.0178197,-7.60341e-06,1.41325e-09,-9.72104e-14,-10453.1,-48.7556], Tmin=(1501.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'CC(=O)C1C=COOO1(6428)',
    structure = SMILES('CC(=O)C1C=COOO1'),
    E0 = (-140.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.573585,0.0600637,-1.47703e-05,-2.85468e-08,1.60254e-11,-16819.4,27.8498], Tmin=(100,'K'), Tmax=(1000.55,'K')), NASAPolynomial(coeffs=[18.8993,0.0212895,-8.34492e-06,1.62206e-09,-1.20391e-13,-22212.9,-69.1982], Tmin=(1000.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-140.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(123trioxene)"""),
)

species(
    label = 'C[C](C=CC=O)OO[O](6362)',
    structure = SMILES('CC(=CC=C[O])OO[O]'),
    E0 = (88.4836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.918305,0.0861327,-8.6381e-05,4.14718e-08,-7.4608e-12,10838.1,37.7141], Tmin=(100,'K'), Tmax=(1550.53,'K')), NASAPolynomial(coeffs=[23.9907,0.0101357,-1.50552e-06,9.64815e-11,-2.45366e-15,4524.61,-88.8325], Tmin=(1550.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.4836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ)"""),
)

species(
    label = 'CC(=CO[O])OC=C[O](6429)',
    structure = SMILES('CC(=CO[O])OC=C[O]'),
    E0 = (-71.6009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.714069,'amu*angstrom^2'), symmetry=1, barrier=(16.4179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.715183,'amu*angstrom^2'), symmetry=1, barrier=(16.4435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714522,'amu*angstrom^2'), symmetry=1, barrier=(16.4283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714981,'amu*angstrom^2'), symmetry=1, barrier=(16.4388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.850937,0.0910875,-9.73515e-05,5.02545e-08,-9.83582e-12,-8423.62,36.3889], Tmin=(100,'K'), Tmax=(1364.13,'K')), NASAPolynomial(coeffs=[23.9162,0.0116959,-2.61064e-06,3.16659e-10,-1.73341e-14,-14551.1,-88.4809], Tmin=(1364.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.6009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ)"""),
)

species(
    label = 'CC(O)=C(C=C[O])O[O](6430)',
    structure = SMILES('CC(O)=C(C=C[O])O[O]'),
    E0 = (-161.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,492.5,1135,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.741633,'amu*angstrom^2'), symmetry=1, barrier=(17.0516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.741369,'amu*angstrom^2'), symmetry=1, barrier=(17.0455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.740083,'amu*angstrom^2'), symmetry=1, barrier=(17.016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.740658,'amu*angstrom^2'), symmetry=1, barrier=(17.0292,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.723453,0.0952342,-0.000101557,4.50655e-08,-4.91781e-12,-19232.2,32.1267], Tmin=(100,'K'), Tmax=(924.57,'K')), NASAPolynomial(coeffs=[23.8773,0.0115426,-2.66933e-06,3.63264e-10,-2.33164e-14,-24753.2,-89.8843], Tmin=(924.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ)"""),
)

species(
    label = 'C=C(O)C(C=C[O])O[O](6431)',
    structure = SMILES('C=C(O)C(C=C[O])O[O]'),
    E0 = (-147.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,492.5,1135,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.855138,'amu*angstrom^2'), symmetry=1, barrier=(19.6613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855416,'amu*angstrom^2'), symmetry=1, barrier=(19.6677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855859,'amu*angstrom^2'), symmetry=1, barrier=(19.6779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.856221,'amu*angstrom^2'), symmetry=1, barrier=(19.6862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28387,0.0949558,-0.000100884,5.03266e-08,-9.38944e-12,-17527.2,37.7255], Tmin=(100,'K'), Tmax=(1475.33,'K')), NASAPolynomial(coeffs=[26.8453,0.00759499,-7.81537e-07,-4.98659e-12,3.24867e-15,-24619.7,-104.828], Tmin=(1475.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-147.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(ROOJ)"""),
)

species(
    label = 'CC(=O)C([O])C=C[O](6432)',
    structure = SMILES('CC(=O)C([O])C=C[O]'),
    E0 = (-148.805,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,434.582,434.588,434.597,434.6],'cm^-1')),
        HinderedRotor(inertia=(0.000892584,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0880741,'amu*angstrom^2'), symmetry=1, barrier=(11.8043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0880783,'amu*angstrom^2'), symmetry=1, barrier=(11.8044,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.703049,0.0612643,-3.46594e-05,-4.90676e-09,7.77238e-12,-17768.4,31.7307], Tmin=(100,'K'), Tmax=(990.205,'K')), NASAPolynomial(coeffs=[17.4217,0.0182443,-6.62883e-06,1.22474e-09,-8.83507e-14,-22281.3,-54.8315], Tmin=(990.205,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=OCOJ)"""),
)

species(
    label = '[CH]=CC(O[O])C(C)=O(5946)',
    structure = SMILES('[CH]=CC(O[O])C(C)=O'),
    E0 = (145.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3120,650,792.5,1650,304.835],'cm^-1')),
        HinderedRotor(inertia=(0.00181409,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15717,'amu*angstrom^2'), symmetry=1, barrier=(10.3656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157176,'amu*angstrom^2'), symmetry=1, barrier=(10.3658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157191,'amu*angstrom^2'), symmetry=1, barrier=(10.3656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01925,0.0672561,-6.43522e-05,3.30204e-08,-7.00154e-12,17593.2,30.4758], Tmin=(100,'K'), Tmax=(1115.92,'K')), NASAPolynomial(coeffs=[11.5936,0.0293522,-1.3402e-05,2.58172e-09,-1.82305e-13,15233.2,-21.699], Tmin=(1115.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'CC1([O])C(C=O)C1O[O](6433)',
    structure = SMILES('CC1([O])C(C=O)C1O[O]'),
    E0 = (-4.54661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.231627,0.079698,-7.88388e-05,4.03022e-08,-8.23393e-12,-408.455,29.5031], Tmin=(100,'K'), Tmax=(1181.91,'K')), NASAPolynomial(coeffs=[16.021,0.0262606,-1.10188e-05,2.04732e-09,-1.42079e-13,-4140.73,-49.3102], Tmin=(1181.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.54661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(ROOJ)"""),
)

species(
    label = 'CC(=O)C(=CC=O)O[O](6434)',
    structure = SMILES('CC(=O)C(=CC=O)O[O]'),
    E0 = (-180.039,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,552.5,462.5,1710,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.265885,'amu*angstrom^2'), symmetry=1, barrier=(6.11321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265613,'amu*angstrom^2'), symmetry=1, barrier=(6.10698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265759,'amu*angstrom^2'), symmetry=1, barrier=(6.11032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265761,'amu*angstrom^2'), symmetry=1, barrier=(6.11037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.484235,0.0877272,-0.00014253,1.35714e-07,-5.16599e-11,-21537.2,30.568], Tmin=(100,'K'), Tmax=(780.707,'K')), NASAPolynomial(coeffs=[5.6344,0.0435126,-2.33263e-05,4.67403e-09,-3.31755e-13,-21798.1,10.4758], Tmin=(780.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-180.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ)"""),
)

species(
    label = 'CH3CO(55)(54)',
    structure = SMILES('C[C]=O'),
    E0 = (-22.2282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,2865.56],'cm^-1')),
        HinderedRotor(inertia=(0.0188671,'amu*angstrom^2'), symmetry=1, barrier=(18.7749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03587,0.000877295,3.071e-05,-3.92476e-08,1.52969e-11,-2682.07,7.86177], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.31372,0.00917378,-3.32204e-06,5.39475e-10,-3.24524e-14,-3645.04,-1.67576], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-22.2282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH3CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]OC=CC=O(6435)',
    structure = SMILES('[O]OC=CC=O'),
    E0 = (-13.6319,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.63794,'amu*angstrom^2'), symmetry=1, barrier=(14.6675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.639209,'amu*angstrom^2'), symmetry=1, barrier=(14.6967,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20762,0.0433675,-4.7744e-05,2.88714e-08,-7.35003e-12,-1578.3,19.5449], Tmin=(100,'K'), Tmax=(929.546,'K')), NASAPolynomial(coeffs=[7.5802,0.0202477,-1.04347e-05,2.11246e-09,-1.53051e-13,-2577.08,-5.98205], Tmin=(929.546,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.6319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ)"""),
)

species(
    label = 'CC(=O)[C](CC=O)O[O](6436)',
    structure = SMILES('CC([O])=C(CC=O)O[O]'),
    E0 = (-143.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,195.03,4000],'cm^-1')),
        HinderedRotor(inertia=(0.438924,'amu*angstrom^2'), symmetry=1, barrier=(11.878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.44056,'amu*angstrom^2'), symmetry=1, barrier=(11.878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0044289,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.44195,'amu*angstrom^2'), symmetry=1, barrier=(11.8767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.88453,0.0744854,-7.44146e-05,4.03769e-08,-9.25182e-12,-17209.2,31.0882], Tmin=(100,'K'), Tmax=(1022.11,'K')), NASAPolynomial(coeffs=[10.6973,0.0360832,-1.80572e-05,3.61787e-09,-2.60813e-13,-19215.1,-16.4673], Tmin=(1022.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(ROOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC(=O)C(C[C]=O)O[O](6437)',
    structure = SMILES('CC(=O)C(C[C]=O)O[O]'),
    E0 = (-151.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,552.5,462.5,1710,2750,2850,1437.5,1250,1305,750,350,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,278.891,3222.13],'cm^-1')),
        HinderedRotor(inertia=(0.178804,'amu*angstrom^2'), symmetry=1, barrier=(9.86709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00216697,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178724,'amu*angstrom^2'), symmetry=1, barrier=(9.86727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294302,'amu*angstrom^2'), symmetry=1, barrier=(16.2459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178726,'amu*angstrom^2'), symmetry=1, barrier=(9.86696,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.509316,0.0823027,-0.000101005,7.05196e-08,-2.03673e-11,-18062.5,34.9586], Tmin=(100,'K'), Tmax=(835.267,'K')), NASAPolynomial(coeffs=[10.313,0.0353547,-1.66962e-05,3.22969e-09,-2.27388e-13,-19700.3,-10.5742], Tmin=(835.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-151.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(=O)C(CC=O)O[O](6438)',
    structure = SMILES('C=C([O])C(CC=O)O[O]'),
    E0 = (-149.152,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,492.5,1135,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,194.615,490.341,3530.53],'cm^-1')),
        HinderedRotor(inertia=(0.578179,'amu*angstrom^2'), symmetry=1, barrier=(13.2935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.578179,'amu*angstrom^2'), symmetry=1, barrier=(13.2935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.578097,'amu*angstrom^2'), symmetry=1, barrier=(13.2916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0779853,'amu*angstrom^2'), symmetry=1, barrier=(13.2955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.259831,0.0861003,-0.00010372,6.77161e-08,-1.78945e-11,-17807.5,34.1387], Tmin=(100,'K'), Tmax=(917.951,'K')), NASAPolynomial(coeffs=[12.8302,0.0313235,-1.42099e-05,2.70756e-09,-1.89422e-13,-20115.2,-25.4299], Tmin=(917.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[C]1OC(C=O)C1O[O](6439)',
    structure = SMILES('C[C]1OC(C=O)C1O[O]'),
    E0 = (-17.6332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.250946,0.0734371,-6.92847e-05,3.58848e-08,-7.34885e-12,-1978.13,32.3431], Tmin=(100,'K'), Tmax=(1279.82,'K')), NASAPolynomial(coeffs=[14.713,0.0242514,-6.96606e-06,9.89467e-10,-5.71058e-14,-5353.52,-39.721], Tmin=(1279.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.6332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(ROOJ) + radical(C2CsJOCs)"""),
)

species(
    label = 'CC(=O)C1[CH][CH]OOO1(6440)',
    structure = SMILES('CC(=O)C1[CH][CH]OOO1'),
    E0 = (120.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.501531,0.0639341,-1.98848e-05,-3.1305e-08,2.07825e-11,14658.6,29.2856], Tmin=(100,'K'), Tmax=(903.965,'K')), NASAPolynomial(coeffs=[18.6773,0.0189905,-4.1869e-06,5.411e-10,-3.39015e-14,9922.8,-64.586], Tmin=(903.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + ring(123trioxane) + radical(CCJCOOH) + radical(CCsJOO)"""),
)

species(
    label = 'CC(=O)C(=CC=O)OO(6441)',
    structure = SMILES('CC(=O)C(=CC=O)OO'),
    E0 = (-332.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.30471,0.0902912,-0.000134211,1.22084e-07,-4.6476e-11,-39811.2,30.4566], Tmin=(100,'K'), Tmax=(718.803,'K')), NASAPolynomial(coeffs=[6.85107,0.0458423,-2.47194e-05,5.01257e-09,-3.60177e-13,-40545.2,2.47673], Tmin=(718.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-332.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'CC(=O)C1OC1C=O(6442)',
    structure = SMILES('CC(=O)C1OC1C=O'),
    E0 = (-355.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44912,0.0380278,3.79886e-05,-8.68094e-08,3.98216e-11,-42605.4,27.689], Tmin=(100,'K'), Tmax=(908.136,'K')), NASAPolynomial(coeffs=[17.6276,0.014412,-1.70036e-06,9.74322e-11,-6.65235e-15,-47508.5,-59.6207], Tmin=(908.136,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-355.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(Ethylene_oxide)"""),
)

species(
    label = 'CC(=O)C([CH]O[O])C=O(6443)',
    structure = SMILES('CC(=O)C([CH]O[O])C=O'),
    E0 = (-117.97,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,552.5,462.5,1710,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180,1859.53],'cm^-1')),
        HinderedRotor(inertia=(0.182859,'amu*angstrom^2'), symmetry=1, barrier=(4.20429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180178,'amu*angstrom^2'), symmetry=1, barrier=(4.14266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0424156,'amu*angstrom^2'), symmetry=1, barrier=(0.975217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182445,'amu*angstrom^2'), symmetry=1, barrier=(4.19477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.54554,'amu*angstrom^2'), symmetry=1, barrier=(58.527,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05992,0.0744412,-6.87441e-05,8.52708e-09,2.41508e-11,-14092.2,31.771], Tmin=(100,'K'), Tmax=(532.483,'K')), NASAPolynomial(coeffs=[7.66136,0.0409355,-1.96676e-05,3.80995e-09,-2.67549e-13,-15023.2,1.94198], Tmin=(532.483,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = 'CC(=O)C1OOC1C=O(6444)',
    structure = SMILES('CC(=O)C1OOC1C=O'),
    E0 = (-303.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.943631,0.0539929,-1.01311e-05,-2.37574e-08,1.18741e-11,-36423.3,29.7279], Tmin=(100,'K'), Tmax=(1060.88,'K')), NASAPolynomial(coeffs=[15.5778,0.0269883,-1.17828e-05,2.31252e-09,-1.68177e-13,-41113.7,-49.2103], Tmin=(1060.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-303.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(12dioxetane)"""),
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
    E0 = (-169.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-32.3986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (22.7658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-89.6851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (49.3854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (92.0321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-0.476503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (77.6959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-23.4256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (20.4458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-128.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-58.9897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-70.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-94.8099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (134.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (44.5282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (100.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-70.0359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (162.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-144.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (284.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (71.3152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (9.62836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-140.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (170.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (126.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-79.3381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (56.7167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (94.1996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (388.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-4.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (41.2074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (58.0708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (14.1665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (5.68845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-67.9826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-17.6332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (120.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (-90.7816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (-164.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (-39.3246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (4.78051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (-160.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)(2)', 'S(247)(246)'],
    products = ['S(256)(255)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.233e+09,'cm^3/(mol*s)'), n=1.533, Ea=(89.5748,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-COH_Cds;YJ] for rate rule [Cds-COH_Cds;O2b]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 85.1 to 89.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['S(256)(255)'],
    products = ['CC1([O])OOC1C=C[O](6413)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(136.63,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 135.8 to 136.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['S(256)(255)'],
    products = ['CC(=O)C1OOC1[CH][O](6414)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(191.794,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['S(256)(255)'],
    products = ['CC1([O])OC=CC1O[O](6415)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(79.3435,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra] for rate rule [R6_SMS_CO;carbonylbond_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)(3)', 'CC(=O)C(C=C=O)O[O](6416)'],
    products = ['S(256)(255)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CC(=O)C(C=[C]O)O[O](5953)'],
    products = ['S(256)(255)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC(=O)[C](C=C[O])OO(6417)'],
    products = ['S(256)(255)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.47099e+07,'s^-1'), n=1.47622, Ea=(175.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_O;C_rad_out_TwoDe;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CC(=O)C([C]=CO)O[O](5952)'],
    products = ['S(256)(255)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['S(256)(255)'],
    products = ['CC(=O)C([C]=C[O])OO(6418)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC(=O)[C](C=CO)O[O](5949)'],
    products = ['S(256)(255)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.21576e+08,'s^-1'), n=0.977901, Ea=(185.49,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;XH_out] for rate rule [R4H_SDS;C_rad_out_OneDe/O;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(=O)C(C=C[O])OO(6419)'],
    products = ['S(256)(255)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.06724e+06,'s^-1'), n=1.18977, Ea=(32.9114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_2H;XH_out] for rate rule [R5H_SSSS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC(=O)C(C=[C][O])OO(6420)'],
    products = ['S(256)(255)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(=O)C(C=CO)O[O](5951)'],
    products = ['S(256)(255)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(12584.6,'s^-1'), n=1.925, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;C_rad_out_2H;XH_out] for rate rule [R6H_RSSMS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O2(2)(2)', 'CC([O])=CC=C[O](5355)'],
    products = ['S(256)(255)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.34662e+06,'m^3/(mol*s)'), n=-0.0521589, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [C_sec_rad;O2_birad] + [C_rad/H/TwoDe;Y_rad] for rate rule [C_rad/H/TwoDe;O2_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -5.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['CHCHO(45)(45)', 'CC(=O)[CH]O[O](5948)'],
    products = ['S(256)(255)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/OneDe] for rate rule [Cd_rad;C_rad/H/OneDeO]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['S(256)(255)'],
    products = ['CC(=O)C(O[O])C1[CH]O1(6421)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.85157e+11,'s^-1'), n=0.224969, Ea=(213.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;doublebond_intra_pri_HNd_Cs;radadd_intra] + [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['S(256)(255)'],
    products = ['C[C]1OOOC1C=C[O](6422)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.81184e+09,'s^-1'), n=0.551229, Ea=(269.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 268.2 to 269.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['S(256)(255)'],
    products = ['CC(=O)C1[CH]C([O])OO1(6423)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(98.9927,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 95.9 to 99.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['S(256)(255)'],
    products = ['C[C]1OOC=CC1O[O](6424)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.39072e+10,'s^-1'), n=0.346137, Ea=(331.548,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra] for rate rule [R6_SMS_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 330.2 to 331.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['S(256)(255)'],
    products = ['CC(=O)C(C=C=O)OO(6425)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CO(10)(11)', 'CC(C=C[O])O[O](6426)'],
    products = ['S(256)(255)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(538,'cm^3/(mol*s)'), n=3.29, Ea=(437.228,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [CO;Cs_Cs]
Euclidian distance = 0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction22',
    reactants = ['S(256)(255)'],
    products = ['O(4)(4)', 'CC(=O)C1C=COO1(6427)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(240.344,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OO] for rate rule [R4OO_SDS;Y_rad_intra;OOJ]
Euclidian distance = 1.41421356237
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['S(256)(255)'],
    products = ['HO2(8)(9)', 'CC(=O)C=[C]C=O(6097)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction24',
    reactants = ['S(256)(255)'],
    products = ['CC(=O)C1C=COOO1(6428)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(28.0463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination
Ea raised from 25.8 to 28.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C[C](C=CC=O)OO[O](6362)'],
    products = ['S(256)(255)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(104,'s^-1'), n=3.21, Ea=(82.0482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CH3;R_O_R]
Euclidian distance = 1.41421356237
family: ketoenol"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CC(=CO[O])OC=C[O](6429)'],
    products = ['S(256)(255)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(855.663,'s^-1'), n=2.935, Ea=(197.924,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CH3;R_O_C]
Euclidian distance = 2.2360679775
family: ketoenol"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CC(O)=C(C=C[O])O[O](6430)'],
    products = ['S(256)(255)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(104,'s^-1'), n=3.21, Ea=(82.0482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C(O)C(C=C[O])O[O](6431)'],
    products = ['S(256)(255)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O(4)(4)', 'CC(=O)C([O])C=C[O](6432)'],
    products = ['S(256)(255)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['O(4)(4)', '[CH]=CC(O[O])C(C)=O(5946)'],
    products = ['S(256)(255)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['S(256)(255)'],
    products = ['CC1([O])C(C=O)C1O[O](6433)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.90568e+10,'s^-1'), n=0.237, Ea=(164.979,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHDe] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_csHDe]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)(3)', 'CC(=O)C(=CC=O)O[O](6434)'],
    products = ['S(256)(255)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(87.5179,'m^3/(mol*s)'), n=1.66467, Ea=(9.4546,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;HJ] for rate rule [Cds-COOs_Cds;HJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH3CO(55)(54)', '[O]OC=CC=O(6435)'],
    products = ['S(256)(255)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cds-OsH_Cds;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CC(=O)[C](CC=O)O[O](6436)'],
    products = ['S(256)(255)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.4044e+07,'s^-1'), n=1.72645, Ea=(158.144,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_noH;Cs_H_out_H/OneDe] + [R2H_S;C_rad_out_OneDe;Cs_H_out_1H] for rate rule [R2H_S;C_rad_out_OneDe/O;Cs_H_out_H/CO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CC(=O)C(C[C]=O)O[O](6437)'],
    products = ['S(256)(255)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(=O)C(CC=O)O[O](6438)'],
    products = ['S(256)(255)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(46663.9,'s^-1'), n=1.985, Ea=(81.1696,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_2H;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['S(256)(255)'],
    products = ['C[C]1OC(C=O)C1O[O](6439)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(6.89861e+07,'s^-1'), n=1.13751, Ea=(151.395,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_csHCO]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic
Ea raised from 150.3 to 151.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['S(256)(255)'],
    products = ['CC(=O)C1[CH][CH]OOO1(6440)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(289.759,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 288.4 to 289.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction39',
    reactants = ['S(256)(255)'],
    products = ['CC(=O)C(=CC=O)OO(6441)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['S(256)(255)'],
    products = ['O2(S)(162)(161)', 'S(247)(246)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(4.73174,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 4.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction41',
    reactants = ['S(256)(255)'],
    products = ['O(4)(4)', 'CC(=O)C1OC1C=O(6442)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.27e+09,'s^-1'), n=1.06, Ea=(129.704,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2OO_S;C_sec_rad_intra;OOJ] for rate rule [R2OO_S;C_rad/H/OneDe_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CC(=O)C([CH]O[O])C=O(6443)'],
    products = ['S(256)(255)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;CO]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction43',
    reactants = ['S(256)(255)'],
    products = ['CC(=O)C1OOC1C=O(6444)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '379',
    isomers = [
        'S(256)(255)',
    ],
    reactants = [
        ('O2(2)(2)', 'S(247)(246)'),
        ('O2(S)(162)(161)', 'S(247)(246)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '379',
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

