species(
    label = 'S(246)(245)',
    structure = SMILES('CC(=O)C=CC(O)O[O]'),
    E0 = (-299.068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,552.5,462.5,1710,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,188.281,4000],'cm^-1')),
        HinderedRotor(inertia=(0.308706,'amu*angstrom^2'), symmetry=1, barrier=(7.76576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308705,'amu*angstrom^2'), symmetry=1, barrier=(7.76575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308706,'amu*angstrom^2'), symmetry=1, barrier=(7.76575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308706,'amu*angstrom^2'), symmetry=1, barrier=(7.76576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01342,'amu*angstrom^2'), symmetry=1, barrier=(25.4935,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4827.56,'J/mol'), sigma=(7.44261,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=754.05 K, Pc=26.57 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.253431,0.0911207,-0.000133361,1.18791e-07,-4.33916e-11,-35843.1,34.0443], Tmin=(100,'K'), Tmax=(775.642,'K')), NASAPolynomial(coeffs=[6.8213,0.0456165,-2.28629e-05,4.48083e-09,-3.15127e-13,-36512,6.2827], Tmin=(775.642,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-299.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + radical(ROOJ)"""),
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
    label = 'S(245)(244)',
    structure = SMILES('CC([O])=CC=CO'),
    E0 = (-227.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.898987,'amu*angstrom^2'), symmetry=1, barrier=(20.6695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.896645,'amu*angstrom^2'), symmetry=1, barrier=(20.6156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.896937,'amu*angstrom^2'), symmetry=1, barrier=(20.6223,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4358.71,'J/mol'), sigma=(6.881,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=680.82 K, Pc=30.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.342212,0.0632525,-1.84618e-05,-4.28291e-08,2.77283e-11,-27232.2,24.3496], Tmin=(100,'K'), Tmax=(902.368,'K')), NASAPolynomial(coeffs=[24.4528,0.00343297,2.75282e-06,-7.11476e-10,4.87923e-14,-33499.4,-100.109], Tmin=(902.368,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-227.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ)"""),
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
    label = 'CC(=O)[CH]C1OOC1O(5962)',
    structure = SMILES('CC([O])=CC1OOC1O'),
    E0 = (-261.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.241474,0.0780092,-6.99856e-05,3.22176e-08,-5.96179e-12,-31309.6,29.9294], Tmin=(100,'K'), Tmax=(1291.37,'K')), NASAPolynomial(coeffs=[16.3929,0.0279804,-1.18741e-05,2.21754e-09,-1.53981e-13,-35481,-52.1216], Tmin=(1291.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-261.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(12dioxetane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC1([O])C=CC(O)OO1(5963)',
    structure = SMILES('CC1([O])C=CC(O)OO1'),
    E0 = (-260.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0253632,0.0666509,-6.37275e-06,-5.5054e-08,3.00848e-11,-31192.8,26.3402], Tmin=(100,'K'), Tmax=(943.67,'K')), NASAPolynomial(coeffs=[25.0625,0.0123691,-2.49914e-06,4.28554e-10,-3.73709e-14,-38226.6,-105.229], Tmin=(943.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro12dioxin) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = 'CC(=O)C=C[C](O)OO(5964)',
    structure = SMILES('CC([O])=CC=C(O)OO'),
    E0 = (-260.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.652161,'amu*angstrom^2'), symmetry=1, barrier=(14.9945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.652552,'amu*angstrom^2'), symmetry=1, barrier=(15.0035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.652435,'amu*angstrom^2'), symmetry=1, barrier=(15.0008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.6524,'amu*angstrom^2'), symmetry=1, barrier=(15,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.652736,'amu*angstrom^2'), symmetry=1, barrier=(15.0077,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.814974,0.104032,-0.000130152,8.14264e-08,-1.97675e-11,-31130.9,32.333], Tmin=(100,'K'), Tmax=(1016.55,'K')), NASAPolynomial(coeffs=[20.191,0.0213766,-8.18827e-06,1.44142e-09,-9.68932e-14,-35401.6,-69.3534], Tmin=(1016.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC(=O)C=CC([O])OO(5965)',
    structure = SMILES('CC(=O)C=CC([O])OO'),
    E0 = (-225.368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.313116,0.0924997,-0.00014328,1.38496e-07,-5.39521e-11,-26983.8,33.9271], Tmin=(100,'K'), Tmax=(782.606,'K')), NASAPolynomial(coeffs=[3.75876,0.0536193,-2.79933e-05,5.56175e-09,-3.93436e-13,-26871.8,22.3098], Tmin=(782.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-225.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + radical(CCOJ)"""),
)

species(
    label = 'CC(=O)C=[C]C(O)OO(5966)',
    structure = SMILES('CC(=O)C=[C]C(O)OO'),
    E0 = (-213.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0625862,0.0995655,-0.000154492,1.40419e-07,-5.14853e-11,-25509.4,34.8917], Tmin=(100,'K'), Tmax=(783.696,'K')), NASAPolynomial(coeffs=[7.36021,0.0465301,-2.39855e-05,4.73536e-09,-3.33621e-13,-26207.6,3.85829], Tmin=(783.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)[C]=CC(O)OO(5967)',
    structure = SMILES('CC([O])=C=CC(O)OO'),
    E0 = (-254.63,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,230.284,230.3],'cm^-1')),
        HinderedRotor(inertia=(0.332386,'amu*angstrom^2'), symmetry=1, barrier=(12.5212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332333,'amu*angstrom^2'), symmetry=1, barrier=(12.5203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.33263,'amu*angstrom^2'), symmetry=1, barrier=(12.5202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.509479,'amu*angstrom^2'), symmetry=1, barrier=(19.1823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332891,'amu*angstrom^2'), symmetry=1, barrier=(12.5197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.357751,0.101734,-0.000136398,9.83002e-08,-2.84643e-11,-30473.3,33.6702], Tmin=(100,'K'), Tmax=(842.642,'K')), NASAPolynomial(coeffs=[13.8874,0.0341162,-1.60369e-05,3.08046e-09,-2.15473e-13,-32874.2,-32.6165], Tmin=(842.642,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-254.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(=O)C=CC(O)OO(5968)',
    structure = SMILES('C=C([O])C=CC(O)OO'),
    E0 = (-306.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.475573,0.0999096,-0.000121689,7.61487e-08,-1.88335e-11,-36733,33.2425], Tmin=(100,'K'), Tmax=(988.863,'K')), NASAPolynomial(coeffs=[17.3811,0.0276789,-1.21238e-05,2.28287e-09,-1.59147e-13,-40264.6,-52.7058], Tmin=(988.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-306.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = 'CC(=O)C=C[CH]O[O](5969)',
    structure = SMILES('CC(=O)C=C[CH]O[O]'),
    E0 = (25.1498,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,256.177,1520.47],'cm^-1')),
        HinderedRotor(inertia=(0.168238,'amu*angstrom^2'), symmetry=1, barrier=(7.83024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168067,'amu*angstrom^2'), symmetry=1, barrier=(7.82951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0143296,'amu*angstrom^2'), symmetry=1, barrier=(23.5085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.504783,'amu*angstrom^2'), symmetry=1, barrier=(23.5085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55801,0.0570581,-4.09152e-05,1.47587e-08,-2.24572e-12,3109.27,27.8442], Tmin=(100,'K'), Tmax=(1447.72,'K')), NASAPolynomial(coeffs=[10.513,0.0323157,-1.52791e-05,2.95335e-09,-2.07091e-13,516.436,-18.6715], Tmin=(1447.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.1498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + radical(C=CCJO) + radical(ROOJ)"""),
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
    label = 'CC(=O)C=CC([O])O[O](5970)',
    structure = SMILES('CC(=O)C=CC([O])O[O]'),
    E0 = (-73.3629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,492.5,1135,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,2390.47],'cm^-1')),
        HinderedRotor(inertia=(0.223752,'amu*angstrom^2'), symmetry=1, barrier=(5.14449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223898,'amu*angstrom^2'), symmetry=1, barrier=(5.14787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224133,'amu*angstrom^2'), symmetry=1, barrier=(5.15325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224211,'amu*angstrom^2'), symmetry=1, barrier=(5.15505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.54446,0.0892722,-0.000149064,1.48648e-07,-5.76608e-11,-8712.02,33.8562], Tmin=(100,'K'), Tmax=(814.784,'K')), NASAPolynomial(coeffs=[2.34675,0.0516466,-2.68169e-05,5.2763e-09,-3.6954e-13,-8050.48,31.3922], Tmin=(814.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.3629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH]=CC(C)=O(4552)',
    structure = SMILES('[CH]=CC(C)=O'),
    E0 = (120.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710],'cm^-1')),
        HinderedRotor(inertia=(0.223718,'amu*angstrom^2'), symmetry=1, barrier=(5.14373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225086,'amu*angstrom^2'), symmetry=1, barrier=(5.17517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72572,0.0297242,-1.46437e-05,2.73924e-09,-1.34393e-13,14549.4,16.4009], Tmin=(100,'K'), Tmax=(2060.72,'K')), NASAPolynomial(coeffs=[13.6528,0.0132567,-6.10919e-06,1.09503e-09,-7.04091e-14,9038.89,-46.66], Tmin=(2060.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[O]O[CH]O(4536)',
    structure = SMILES('[O]O[CH]O'),
    E0 = (9.28844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,3615,1277.5,1000],'cm^-1')),
        HinderedRotor(inertia=(0.156407,'amu*angstrom^2'), symmetry=1, barrier=(3.59612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156804,'amu*angstrom^2'), symmetry=1, barrier=(3.60523,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58857,0.0452524,-0.00011197,1.22556e-07,-4.65581e-11,1154.51,14.2974], Tmin=(100,'K'), Tmax=(891.504,'K')), NASAPolynomial(coeffs=[-0.0459,0.0212294,-1.12418e-05,2.13324e-09,-1.41938e-13,3048.61,34.6933], Tmin=(891.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.28844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OCJO) + radical(ROOJ)"""),
)

species(
    label = 'CC(=O)C=C[C](O)O[O](5971)',
    structure = SMILES('CC([O])=CC=C(O)O[O]'),
    E0 = (-108.288,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,492.5,1135,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.519831,'amu*angstrom^2'), symmetry=1, barrier=(11.9519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519873,'amu*angstrom^2'), symmetry=1, barrier=(11.9529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519779,'amu*angstrom^2'), symmetry=1, barrier=(11.9507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519477,'amu*angstrom^2'), symmetry=1, barrier=(11.9438,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.574307,0.100727,-0.000135931,9.22579e-08,-2.42443e-11,-12859.4,32.2275], Tmin=(100,'K'), Tmax=(941.622,'K')), NASAPolynomial(coeffs=[18.4224,0.0200291,-7.38068e-06,1.24459e-09,-8.04421e-14,-16437,-58.2783], Tmin=(941.622,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CH3(15)(16)',
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
    label = '[O]OC(O)C=C[C]=O(5972)',
    structure = SMILES('[O]OC(O)[CH]C=C=O'),
    E0 = (-121.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2120,512.5,787.5,3025,407.5,1350,352.5,3615,1277.5,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,279.421,279.613],'cm^-1')),
        HinderedRotor(inertia=(0.535153,'amu*angstrom^2'), symmetry=1, barrier=(29.7346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204908,'amu*angstrom^2'), symmetry=1, barrier=(11.3726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.874689,'amu*angstrom^2'), symmetry=1, barrier=(48.5102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204945,'amu*angstrom^2'), symmetry=1, barrier=(11.3726,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.544891,0.0814303,-0.000115448,8.72034e-08,-2.6414e-11,-14471.4,27.8088], Tmin=(100,'K'), Tmax=(806.875,'K')), NASAPolynomial(coeffs=[11.6646,0.0263043,-1.29655e-05,2.52761e-09,-1.77776e-13,-16265.8,-23.4511], Tmin=(806.875,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-121.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + radical(C=CCJCO) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(=O)C=CC(O)O[O](5973)',
    structure = SMILES('C=C([O])C=CC(O)O[O]'),
    E0 = (-154.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,492.5,1135,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.507923,'amu*angstrom^2'), symmetry=1, barrier=(11.6781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.508029,'amu*angstrom^2'), symmetry=1, barrier=(11.6806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507844,'amu*angstrom^2'), symmetry=1, barrier=(11.6763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.50786,'amu*angstrom^2'), symmetry=1, barrier=(11.6767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.282728,0.0971691,-0.000129448,8.95654e-08,-2.44276e-11,-18459.5,33.3084], Tmin=(100,'K'), Tmax=(900.713,'K')), NASAPolynomial(coeffs=[15.7004,0.0261874,-1.12356e-05,2.06745e-09,-1.41184e-13,-21338.6,-42.1289], Tmin=(900.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-154.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ)"""),
)

species(
    label = 'CC(=O)C=[C]C(O)O[O](5974)',
    structure = SMILES('CC(=O)C=[C]C(O)O[O]'),
    E0 = (-61.2263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,552.5,462.5,1710,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,3026.18],'cm^-1')),
        HinderedRotor(inertia=(0.714095,'amu*angstrom^2'), symmetry=1, barrier=(16.4185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31529,'amu*angstrom^2'), symmetry=1, barrier=(7.24913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.317718,'amu*angstrom^2'), symmetry=1, barrier=(7.30495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.314065,'amu*angstrom^2'), symmetry=1, barrier=(7.22098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312699,'amu*angstrom^2'), symmetry=1, barrier=(7.18956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.161536,0.0964385,-0.000160714,1.51307e-07,-5.56075e-11,-7237.25,34.8458], Tmin=(100,'K'), Tmax=(822.596,'K')), NASAPolynomial(coeffs=[5.94712,0.0445586,-2.28095e-05,4.44998e-09,-3.09727e-13,-7385.67,12.947], Tmin=(822.596,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.2263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + radical(ROOJ) + radical(Cds_S)"""),
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
    label = '[CH]=CC(O)O[O](5975)',
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
    label = 'CC(=O)[C]=CC(O)O[O](5976)',
    structure = SMILES('CC([O])=C=CC(O)O[O]'),
    E0 = (-102.626,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,492.5,1135,1000,540,610,2055,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.426532,'amu*angstrom^2'), symmetry=1, barrier=(9.80681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.426513,'amu*angstrom^2'), symmetry=1, barrier=(9.80636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.426494,'amu*angstrom^2'), symmetry=1, barrier=(9.80595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.426487,'amu*angstrom^2'), symmetry=1, barrier=(9.80577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.282599,0.100566,-0.000150559,1.2137e-07,-3.88108e-11,-12194.9,34.1463], Tmin=(100,'K'), Tmax=(823.186,'K')), NASAPolynomial(coeffs=[12.6849,0.0317626,-1.46304e-05,2.73883e-09,-1.868e-13,-14133.6,-24.6994], Tmin=(823.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.626,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(ROOJ)"""),
)

species(
    label = 'CC(=O)C1[CH]C(O)OO1(5955)',
    structure = SMILES('CC(=O)C1[CH]C(O)OO1'),
    E0 = (-295.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.673556,0.0560896,8.33384e-06,-6.18041e-08,3.15781e-11,-35433.8,32.9487], Tmin=(100,'K'), Tmax=(915.435,'K')), NASAPolynomial(coeffs=[19.4755,0.0187251,-3.83432e-06,5.05348e-10,-3.45644e-14,-40752.9,-66.3491], Tmin=(915.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-295.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + ring(12dioxolane) + radical(CCJCOOH)"""),
)

species(
    label = 'C[C]1C=CC(O)OOO1(5977)',
    structure = SMILES('CC1=C[CH]C(O)OOO1'),
    E0 = (-99.8456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.610818,0.0173277,0.000197004,-3.13825e-07,1.34907e-10,-11833.1,30.785], Tmin=(100,'K'), Tmax=(910.275,'K')), NASAPolynomial(coeffs=[44.2782,-0.0250276,2.03943e-05,-4.017e-09,2.58123e-13,-25978,-209.808], Tmin=(910.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.8456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cycloheptane) + radical(C=CCJCO)"""),
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
    label = 'CC=CC(O)O[O](5978)',
    structure = SMILES('CC=CC(O)O[O]'),
    E0 = (-166.451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,492.5,1135,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.37824,'amu*angstrom^2'), symmetry=1, barrier=(8.69647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.37844,'amu*angstrom^2'), symmetry=1, barrier=(8.70108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.378937,'amu*angstrom^2'), symmetry=1, barrier=(8.7125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.378248,'amu*angstrom^2'), symmetry=1, barrier=(8.69666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09057,0.068927,-8.44761e-05,6.35248e-08,-2.02488e-11,-19919.2,26.8073], Tmin=(100,'K'), Tmax=(754.18,'K')), NASAPolynomial(coeffs=[7.57152,0.0345536,-1.61106e-05,3.09248e-09,-2.16403e-13,-20896.8,-2.63111], Tmin=(754.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ)"""),
)

species(
    label = 'CC(=O)C=C=CO(5366)',
    structure = SMILES('CC(=O)C=C=CO'),
    E0 = (-194.708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,180],'cm^-1')),
        HinderedRotor(inertia=(0.857163,'amu*angstrom^2'), symmetry=1, barrier=(19.7079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855569,'amu*angstrom^2'), symmetry=1, barrier=(19.6712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.854091,'amu*angstrom^2'), symmetry=1, barrier=(19.6372,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.801099,0.0615574,-4.83841e-05,1.36454e-08,3.99662e-13,-23295.1,23.6111], Tmin=(100,'K'), Tmax=(1036.66,'K')), NASAPolynomial(coeffs=[16.3184,0.0171946,-6.63719e-06,1.23191e-09,-8.74905e-14,-27345.8,-55.8294], Tmin=(1036.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-194.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C(O)C=CC(O)O[O](5979)',
    structure = SMILES('C=C(O)C=CC(O)O[O]'),
    E0 = (-292.546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,492.5,1135,1000,3580,3650,1210,1345,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.680992,'amu*angstrom^2'), symmetry=1, barrier=(15.6573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681829,'amu*angstrom^2'), symmetry=1, barrier=(15.6766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682303,'amu*angstrom^2'), symmetry=1, barrier=(15.6875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682415,'amu*angstrom^2'), symmetry=1, barrier=(15.6901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.68112,'amu*angstrom^2'), symmetry=1, barrier=(15.6603,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.71585,0.101999,-0.000125205,7.70356e-08,-1.84473e-11,-35013.8,33.357], Tmin=(100,'K'), Tmax=(1028.59,'K')), NASAPolynomial(coeffs=[19.8112,0.0221717,-8.79008e-06,1.58184e-09,-1.07849e-13,-39236.5,-66.2526], Tmin=(1028.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-292.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'CC(O)=C=CC(O)O[O](5980)',
    structure = SMILES('CC(O)=C=CC(O)O[O]'),
    E0 = (-240.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,492.5,1135,1000,540,610,2055,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.600908,'amu*angstrom^2'), symmetry=1, barrier=(13.8161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602491,'amu*angstrom^2'), symmetry=1, barrier=(13.8524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.601225,'amu*angstrom^2'), symmetry=1, barrier=(13.8234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602076,'amu*angstrom^2'), symmetry=1, barrier=(13.8429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60159,'amu*angstrom^2'), symmetry=1, barrier=(13.8317,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.487856,0.102492,-0.000135099,9.2686e-08,-2.51701e-11,-28758.8,33.3916], Tmin=(100,'K'), Tmax=(902.501,'K')), NASAPolynomial(coeffs=[16.0994,0.0289756,-1.29124e-05,2.42842e-09,-1.68207e-13,-31752.8,-44.9307], Tmin=(902.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-240.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
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
    label = 'CC(=O)C=CC([O])O(5981)',
    structure = SMILES('CC(=O)C=CC([O])O'),
    E0 = (-296.873,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3615,1277.5,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,2702.11],'cm^-1')),
        HinderedRotor(inertia=(0.172389,'amu*angstrom^2'), symmetry=1, barrier=(3.96356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174244,'amu*angstrom^2'), symmetry=1, barrier=(4.00622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169445,'amu*angstrom^2'), symmetry=1, barrier=(3.89588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178574,'amu*angstrom^2'), symmetry=1, barrier=(4.10576,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00467,0.0753597,-0.00011121,1.08518e-07,-4.30018e-11,-35606.9,30.6843], Tmin=(100,'K'), Tmax=(784.351,'K')), NASAPolynomial(coeffs=[2.57013,0.0487889,-2.48491e-05,4.90152e-09,-3.45761e-13,-35280.7,27.1569], Tmin=(784.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-296.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + radical(CCOJ)"""),
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
    E0 = (-233.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-177.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-260.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-85.1829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-119.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-153.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-221.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-215.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (55.0027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (138.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (129.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (104.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (14.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (57.0508,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (150.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (94.4429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (113.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-241.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-99.8456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (151.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-175.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-120.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-88.3671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-158.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-53.8679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)(2)', 'S(245)(244)'],
    products = ['S(246)(245)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.19343e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;Y_rad] for rate rule [C_rad/H/OneDeO;O2_birad]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction2',
    reactants = ['S(246)(245)'],
    products = ['CC(=O)[CH]C1OOC1O(5962)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_O] for rate rule [R5_SS_D;doublebond_intra_HDe_pri;radadd_intra_O]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['S(246)(245)'],
    products = ['CC1([O])C=CC(O)OO1(5963)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(24026.2,'s^-1'), n=1.34471, Ea=(38.3686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SSSR;multiplebond_intra;radadd_intra] for rate rule [R7_SSSM_CO;carbonylbond_intra_Nd;radadd_intra_O]
Euclidian distance = 3.0
family: Intra_R_Add_Exocyclic
Ea raised from 35.8 to 38.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CC(=O)C=C[C](O)OO(5964)'],
    products = ['S(246)(245)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.47099e+07,'s^-1'), n=1.47622, Ea=(175.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_O;C_rad_out_OneDe/O;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CC(=O)C=CC([O])OO(5965)'],
    products = ['S(246)(245)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(94.0113,'s^-1'), n=2.81534, Ea=(105.641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['S(246)(245)'],
    products = ['CC(=O)C=[C]C(O)OO(5966)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC(=O)[C]=CC(O)OO(5967)'],
    products = ['S(246)(245)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleDe;O_H_out]
Euclidian distance = 3.31662479036
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['S(246)(245)'],
    products = ['[CH2]C(=O)C=CC(O)OO(5968)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(271800,'s^-1'), n=1.51, Ea=(83.4708,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 266 used for R7H_OOCs4;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R7H_OOCs4;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)(5)', 'CC(=O)C=C[CH]O[O](5969)'],
    products = ['S(246)(245)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.81233e+06,'m^3/(mol*s)'), n=0.00568686, Ea=(1.481,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;C_rad/H/OneDe] + [O_rad;C_sec_rad] for rate rule [O_pri_rad;C_rad/H/OneDeO]
Euclidian distance = 2.2360679775
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)(3)', 'CC(=O)C=CC([O])O[O](5970)'],
    products = ['S(246)(245)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC(C)=O(4552)', '[O]O[CH]O(4536)'],
    products = ['S(246)(245)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.18e+06,'m^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/NonDeO;Y_rad] for rate rule [C_rad/H/O2;Cd_pri_rad]
Euclidian distance = 2.2360679775
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)(3)', 'CC(=O)C=C[C](O)O[O](5971)'],
    products = ['S(246)(245)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""Estimated using template [C_rad/OneDe;H_rad] for rate rule [C_rad/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH3(15)(16)', '[O]OC(O)C=C[C]=O(5972)'],
    products = ['S(246)(245)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.2e+13,'cm^3/(mol*s)','+|-',8.4e+12), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using template [CO_sec_rad;C_methyl] for rate rule [CO_rad/OneDe;C_methyl]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)(3)', '[CH2]C(=O)C=CC(O)O[O](5973)'],
    products = ['S(246)(245)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.83648e+06,'m^3/(mol*s)'), n=0.3308, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;H_rad] for rate rule [C_rad/H2/CO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -1.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)(3)', 'CC(=O)C=[C]C(O)O[O](5974)'],
    products = ['S(246)(245)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CH3CO(55)(54)', '[CH]=CC(O)O[O](5975)'],
    products = ['S(246)(245)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.6647e+07,'m^3/(mol*s)'), n=-0.0666667, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;CO_rad/NonDe] + [Cd_pri_rad;CO_rad] for rate rule [Cd_pri_rad;CO_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)(3)', 'CC(=O)[C]=CC(O)O[O](5976)'],
    products = ['S(246)(245)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['S(246)(245)'],
    products = ['CC(=O)C1[CH]C(O)OO1(5955)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(57.3526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_HCO;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['S(246)(245)'],
    products = ['C[C]1C=CC(O)OOO1(5977)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(199.223,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;carbonyl_intra;radadd_intra] for rate rule [R7_linear;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic
Ea raised from 191.1 to 199.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['CO(10)(11)', 'CC=CC(O)O[O](5978)'],
    products = ['S(246)(245)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(538,'cm^3/(mol*s)'), n=3.29, Ea=(437.228,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO;R_R'] for rate rule [CO;C_methyl_Cd_pri]
Euclidian distance = 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction21',
    reactants = ['S(246)(245)'],
    products = ['HO2(8)(9)', 'S(247)(246)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.813e+10,'s^-1','*|/',10), n=0.493, Ea=(123.219,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 13 used for R2OO_O_HNd
Exact match found for rate rule [R2OO_O_HNd]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction22',
    reactants = ['S(246)(245)'],
    products = ['HO2(8)(9)', 'CC(=O)C=C=CO(5366)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C(O)C=CC(O)O[O](5979)'],
    products = ['S(246)(245)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(204.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CC(O)=C=CC(O)O[O](5980)'],
    products = ['S(246)(245)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(104,'s^-1'), n=3.21, Ea=(82.0482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(4)(4)', 'CC(=O)C=CC([O])O(5981)'],
    products = ['S(246)(245)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

network(
    label = '368',
    isomers = [
        'S(246)(245)',
    ],
    reactants = [
        ('O2(2)(2)', 'S(245)(244)'),
        ('HO2(8)(9)', 'S(247)(246)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '368',
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

