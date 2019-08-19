species(
    label = 'S(659)(658)',
    structure = SMILES('CC(=O)C(C=CO)O[O]'),
    E0 = (-310.491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,552.5,462.5,1710,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,276.687,276.694],'cm^-1')),
        HinderedRotor(inertia=(0.272732,'amu*angstrom^2'), symmetry=1, barrier=(14.8205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272779,'amu*angstrom^2'), symmetry=1, barrier=(14.8205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272788,'amu*angstrom^2'), symmetry=1, barrier=(14.8206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272723,'amu*angstrom^2'), symmetry=1, barrier=(14.8206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272706,'amu*angstrom^2'), symmetry=1, barrier=(14.8206,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4827.56,'J/mol'), sigma=(7.44261,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=754.05 K, Pc=26.57 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.373521,0.0846528,-7.07195e-05,1.92703e-08,1.93585e-12,-37175.7,34.9243], Tmin=(100,'K'), Tmax=(975.29,'K')), NASAPolynomial(coeffs=[21.9864,0.0177386,-5.93452e-06,1.04964e-09,-7.45415e-14,-42716.2,-78.4345], Tmin=(975.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-310.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ)"""),
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
    label = 'CC1([O])OOC1C=CO(5939)',
    structure = SMILES('CC1([O])OOC1C=CO'),
    E0 = (-173.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.32149,0.0753665,-2.56699e-05,-3.74487e-08,2.42204e-11,-20737.1,28.2802], Tmin=(100,'K'), Tmax=(947.582,'K')), NASAPolynomial(coeffs=[26.075,0.0124398,-2.83293e-06,4.9878e-10,-4.18252e-14,-27917.1,-109.136], Tmin=(947.582,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxetane) + radical(CC(C)(O)OJ)"""),
)

species(
    label = 'CC(=O)C1OOC1[CH]O(5940)',
    structure = SMILES('CC(=O)C1OOC1[CH]O'),
    E0 = (-203.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.471693,0.0791428,-7.48527e-05,3.77274e-08,-7.84162e-12,-24309.8,29.1867], Tmin=(100,'K'), Tmax=(1137.93,'K')), NASAPolynomial(coeffs=[13.2655,0.0341709,-1.55718e-05,2.99728e-09,-2.11544e-13,-27221.5,-34.1891], Tmin=(1137.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + ring(12dioxetane) + radical(CCsJOH)"""),
)

species(
    label = 'CC(=O)[C](C=CO)OO(5941)',
    structure = SMILES('CC([O])=C(C=CO)OO'),
    E0 = (-317.049,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.796372,'amu*angstrom^2'), symmetry=1, barrier=(18.3102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.796425,'amu*angstrom^2'), symmetry=1, barrier=(18.3114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.796121,'amu*angstrom^2'), symmetry=1, barrier=(18.3044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.79663,'amu*angstrom^2'), symmetry=1, barrier=(18.3161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.796538,'amu*angstrom^2'), symmetry=1, barrier=(18.314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06309,0.102182,-0.000109035,4.82876e-08,-5.16454e-12,-37941.4,32.6777], Tmin=(100,'K'), Tmax=(921.132,'K')), NASAPolynomial(coeffs=[25.2842,0.012354,-2.79113e-06,3.69359e-10,-2.33188e-14,-43838.2,-97.93], Tmin=(921.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC(=O)C([C]=CO)OO(5942)',
    structure = SMILES('CC(=O)C([C]=CO)OO'),
    E0 = (-224.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.940105,0.0959202,-0.000101018,5.1612e-08,-1.01434e-11,-26831,36.6784], Tmin=(100,'K'), Tmax=(1256.61,'K')), NASAPolynomial(coeffs=[24.1916,0.0159212,-5.52366e-06,9.4904e-10,-6.40164e-14,-33147.1,-90.3078], Tmin=(1256.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-224.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=O)C(C=CO)OO(5943)',
    structure = SMILES('C=C([O])C(C=CO)OO'),
    E0 = (-303.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.64657,0.102202,-0.000109524,5.52182e-08,-1.04181e-11,-36235.4,38.3575], Tmin=(100,'K'), Tmax=(1456.13,'K')), NASAPolynomial(coeffs=[28.3019,0.00832177,-8.54076e-07,-1.06534e-11,4.23164e-15,-43726.2,-113.154], Tmin=(1456.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-303.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC(=O)C(C=[C]O)OO(5944)',
    structure = SMILES('CC(=O)C(C=[C]O)OO'),
    E0 = (-222.752,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,1685,370,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.420669,0.0898533,-9.19774e-05,4.68275e-08,-9.3303e-12,-26625.4,37.4667], Tmin=(100,'K'), Tmax=(1225.83,'K')), NASAPolynomial(coeffs=[20.3915,0.0219414,-8.87643e-06,1.63326e-09,-1.13244e-13,-31727.9,-67.1778], Tmin=(1225.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.752,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = 'CC(=O)C(C=C[O])OO(5945)',
    structure = SMILES('CC(=O)C(C=C[O])OO'),
    E0 = (-321.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.284176,0.0818939,-6.44102e-05,1.78857e-08,5.37413e-13,-38446.5,35.2585], Tmin=(100,'K'), Tmax=(1052.82,'K')), NASAPolynomial(coeffs=[21.0509,0.0218684,-8.85589e-06,1.68294e-09,-1.20909e-13,-44104.6,-74.3051], Tmin=(1052.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-321.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ)"""),
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
    label = '[O]O[CH]C=CO(5947)',
    structure = SMILES('[O]O[CH]C=CO'),
    E0 = (-14.9999,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.12219,'amu*angstrom^2'), symmetry=1, barrier=(25.8013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12188,'amu*angstrom^2'), symmetry=1, barrier=(25.7941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12189,'amu*angstrom^2'), symmetry=1, barrier=(25.7944,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41255,0.0421239,-1.65469e-06,-4.44282e-08,2.4907e-11,-1697.1,22.1858], Tmin=(100,'K'), Tmax=(913.122,'K')), NASAPolynomial(coeffs=[19.9116,0.000280065,2.70015e-06,-6.02097e-10,3.85777e-14,-6709.41,-74.3272], Tmin=(913.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.9999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = 'CHCHOH(49)(49)',
    structure = SMILES('[CH]=CO'),
    E0 = (118.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,443.815,444.347,445.091,446.033],'cm^-1')),
        HinderedRotor(inertia=(0.000109588,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1868.27,'J/mol'), sigma=(4.162,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.0819,0.0104649,2.52045e-05,-4.72208e-08,2.08838e-11,14310.7,10.1558], Tmin=(100,'K'), Tmax=(919.505,'K')), NASAPolynomial(coeffs=[11.3957,-0.000441188,1.78863e-06,-3.67188e-10,2.20229e-14,11713.9,-35.0624], Tmin=(919.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCHOH""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[O]OC([C]=O)C=CO(5950)',
    structure = SMILES('[O]OC([C]=O)C=CO'),
    E0 = (-95.6599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1855,455,950,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.777517,'amu*angstrom^2'), symmetry=1, barrier=(17.8766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.7778,'amu*angstrom^2'), symmetry=1, barrier=(17.8832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.775318,'amu*angstrom^2'), symmetry=1, barrier=(17.8261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.775879,'amu*angstrom^2'), symmetry=1, barrier=(17.839,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.259946,0.0778139,-8.8216e-05,4.64867e-08,-9.09223e-12,-11338.2,33.746], Tmin=(100,'K'), Tmax=(1421.55,'K')), NASAPolynomial(coeffs=[22.4357,0.00373233,7.37487e-07,-2.87279e-10,2.30345e-14,-16758.2,-80.0987], Tmin=(1421.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.6599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(ROOJ)"""),
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
    label = 'C[C]1OOOC1C=CO(5954)',
    structure = SMILES('C[C]1OOOC1C=CO'),
    E0 = (-40.6916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0333047,0.0673988,1.93785e-06,-7.55376e-08,4.2341e-11,-4728.63,29.0511], Tmin=(100,'K'), Tmax=(888.247,'K')), NASAPolynomial(coeffs=[26.6543,0.00688586,3.36446e-06,-9.81527e-10,7.15542e-14,-11823.5,-109.789], Tmin=(888.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.6916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(123trioxolane) + radical(C2CsJOO)"""),
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
    label = 'CC(C=CO)O[O](5956)',
    structure = SMILES('CC(C=CO)O[O]'),
    E0 = (-174.611,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.800067,'amu*angstrom^2'), symmetry=1, barrier=(18.3951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805294,'amu*angstrom^2'), symmetry=1, barrier=(18.5153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.801427,'amu*angstrom^2'), symmetry=1, barrier=(18.4264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800868,'amu*angstrom^2'), symmetry=1, barrier=(18.4135,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.476086,0.0642878,-3.36172e-05,-1.58067e-08,1.46087e-11,-20861.9,27.3733], Tmin=(100,'K'), Tmax=(933.199,'K')), NASAPolynomial(coeffs=[20.6432,0.0113492,-2.37929e-06,3.50083e-10,-2.63369e-14,-26084.7,-76.3438], Tmin=(933.199,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ)"""),
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
    label = 'CC(=O)C(CC=O)O[O](5957)',
    structure = SMILES('CC(=O)C(CC=O)O[O]'),
    E0 = (-311.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.744521,0.0773689,-7.8821e-05,4.65974e-08,-1.18264e-11,-37309.8,33.2997], Tmin=(100,'K'), Tmax=(926.174,'K')), NASAPolynomial(coeffs=[9.24887,0.0406399,-1.93358e-05,3.77944e-09,-2.68635e-13,-38885.1,-7.07658], Tmin=(926.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-311.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(ROOJ)"""),
)

species(
    label = 'CC(=CC=CO)OO[O](5733)',
    structure = SMILES('CC(=CC=CO)OO[O]'),
    E0 = (-52.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58642,0.096167,-0.000100714,4.97998e-08,-9.10344e-12,-6147.84,38.6995], Tmin=(100,'K'), Tmax=(1562.55,'K')), NASAPolynomial(coeffs=[26.4065,0.00758166,5.732e-07,-3.47052e-10,2.89471e-14,-12829.7,-102.232], Tmin=(1562.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'CC(=CO[O])OC=CO(5958)',
    structure = SMILES('CC(=CO[O])OC=CO'),
    E0 = (-213.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.51593,0.101086,-0.000111561,5.84267e-08,-1.14159e-11,-25409.7,37.363], Tmin=(100,'K'), Tmax=(1421.37,'K')), NASAPolynomial(coeffs=[26.8116,0.00845797,-1.86133e-07,-2.00809e-10,1.97448e-14,-32158.5,-104.675], Tmin=(1421.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'CC(O)=C(C=CO)O[O](5959)',
    structure = SMILES('CC(O)=C(C=CO)O[O]'),
    E0 = (-302.849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,492.5,1135,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.771773,'amu*angstrom^2'), symmetry=1, barrier=(17.7446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.770535,'amu*angstrom^2'), symmetry=1, barrier=(17.7161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.771718,'amu*angstrom^2'), symmetry=1, barrier=(17.7433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.771801,'amu*angstrom^2'), symmetry=1, barrier=(17.7452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.772936,'amu*angstrom^2'), symmetry=1, barrier=(17.7713,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75529,0.109618,-0.000131366,7.38157e-08,-1.54243e-11,-36202.5,34.4124], Tmin=(100,'K'), Tmax=(1322.08,'K')), NASAPolynomial(coeffs=[28.152,0.00602195,1.04745e-06,-4.55486e-10,3.85121e-14,-42964.6,-113.89], Tmin=(1322.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-302.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'C=C(O)C(C=CO)O[O](5960)',
    structure = SMILES('C=C(O)C(C=CO)O[O]'),
    E0 = (-288.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,492.5,1135,1000,3580,3650,1210,1345,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.858067,'amu*angstrom^2'), symmetry=1, barrier=(19.7287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860227,'amu*angstrom^2'), symmetry=1, barrier=(19.7783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.856281,'amu*angstrom^2'), symmetry=1, barrier=(19.6876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858561,'amu*angstrom^2'), symmetry=1, barrier=(19.74,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.86028,'amu*angstrom^2'), symmetry=1, barrier=(19.7795,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95221,0.104991,-0.000115214,5.86436e-08,-1.10257e-11,-34513.1,38.7119], Tmin=(100,'K'), Tmax=(1502.57,'K')), NASAPolynomial(coeffs=[29.436,0.00479436,1.42071e-06,-4.74734e-10,3.66504e-14,-42067.5,-119.249], Tmin=(1502.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-288.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
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
    label = 'CC(=O)C([O])C=CO(5961)',
    structure = SMILES('CC(=O)C([O])C=CO'),
    E0 = (-290.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3615,1277.5,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,299.922,299.998,300.046],'cm^-1')),
        HinderedRotor(inertia=(0.00187524,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216907,'amu*angstrom^2'), symmetry=1, barrier=(13.8393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216828,'amu*angstrom^2'), symmetry=1, barrier=(13.8373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216812,'amu*angstrom^2'), symmetry=1, barrier=(13.8377,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.303015,0.0681119,-3.77041e-05,-1.14662e-08,1.25998e-11,-34766,31.7558], Tmin=(100,'K'), Tmax=(947.661,'K')), NASAPolynomial(coeffs=[20.5618,0.0147452,-4.11195e-06,6.94839e-10,-5.07958e-14,-40049.1,-72.5074], Tmin=(947.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=OCOJ)"""),
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
    E0 = (-236.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-173.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-203.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-141.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-164.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-270.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-189.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-231.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (173.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (42.7635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-37.2281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (102.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (46.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (40.5277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (60.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (139.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (141.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-40.6916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-252.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (143.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-131.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-166.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (29.0692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-15.1394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-220.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-84.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-47.2631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)(2)', 'S(245)(244)'],
    products = ['S(659)(658)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(6.34662e+06,'m^3/(mol*s)'), n=-0.0521589, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [C_sec_rad;O2_birad] + [C_rad/H/TwoDe;Y_rad] for rate rule [C_rad/H/TwoDe;O2_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -5.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction2',
    reactants = ['S(659)(658)'],
    products = ['CC1([O])OOC1C=CO(5939)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(136.63,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 135.8 to 136.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['S(659)(658)'],
    products = ['CC(=O)C1OOC1[CH]O(5940)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.65e+06,'s^-1'), n=1.3, Ea=(107.323,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 107.0 to 107.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CC(=O)[C](C=CO)OO(5941)'],
    products = ['S(659)(658)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.47099e+07,'s^-1'), n=1.47622, Ea=(175.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_O;C_rad_out_TwoDe;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['S(659)(658)'],
    products = ['CC(=O)C([C]=CO)OO(5942)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(=O)C(C=CO)OO(5943)'],
    products = ['S(659)(658)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.06724e+06,'s^-1'), n=1.18977, Ea=(32.9114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_2H;XH_out] for rate rule [R5H_SSSS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC(=O)C(C=[C]O)OO(5944)'],
    products = ['S(659)(658)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 3.31662479036
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['S(659)(658)'],
    products = ['CC(=O)C(C=C[O])OO(5945)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(146928,'s^-1'), n=1.47286, Ea=(78.6904,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6H;O_rad_out;XH_out] + [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)(5)', '[CH]=CC(O[O])C(C)=O(5946)'],
    products = ['S(659)(658)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.44289e+07,'m^3/(mol*s)'), n=0.0225, Ea=(0.0523,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;Cd_pri_rad] + [O_pri_rad;Y_rad] for rate rule [O_pri_rad;Cd_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)(3)', 'S(256)(255)'],
    products = ['S(659)(658)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH3CO(55)(54)', '[O]O[CH]C=CO(5947)'],
    products = ['S(659)(658)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.64e+13,'cm^3/(mol*s)','*|/',2), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [C_sec_rad;CO_rad/NonDe] for rate rule [C_rad/H/OneDeO;CO_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CHCHOH(49)(49)', 'CC(=O)[CH]O[O](5948)'],
    products = ['S(659)(658)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/OneDe] for rate rule [Cd_pri_rad;C_rad/H/OneDeO]
Euclidian distance = 2.2360679775
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)(3)', 'CC(=O)[C](C=CO)O[O](5949)'],
    products = ['S(659)(658)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.62e+13,'cm^3/(mol*s)'), n=0.228, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/TwoDe;H_rad] for rate rule [C_rad/TwoDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['CH3(15)(16)', '[O]OC([C]=O)C=CO(5950)'],
    products = ['S(659)(658)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.2e+13,'cm^3/(mol*s)','+|-',8.4e+12), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 72 used for C_methyl;CO_rad/NonDe
Exact match found for rate rule [CO_rad/NonDe;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)(3)', '[CH2]C(=O)C(C=CO)O[O](5951)'],
    products = ['S(659)(658)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.83648e+06,'m^3/(mol*s)'), n=0.3308, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;H_rad] for rate rule [C_rad/H2/CO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -1.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)(3)', 'CC(=O)C([C]=CO)O[O](5952)'],
    products = ['S(659)(658)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)(3)', 'CC(=O)C(C=[C]O)O[O](5953)'],
    products = ['S(659)(658)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['S(659)(658)'],
    products = ['C[C]1OOOC1C=CO(5954)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.81184e+09,'s^-1'), n=0.551229, Ea=(269.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 268.2 to 269.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['S(659)(658)'],
    products = ['CC(=O)C1[CH]C(O)OO1(5955)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.94e+07,'s^-1'), n=0.93, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_HNd;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_HNd_O;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CO(10)(11)', 'CC(C=CO)O[O](5956)'],
    products = ['S(659)(658)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(538,'cm^3/(mol*s)'), n=3.29, Ea=(437.228,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [CO;Cs_Cs]
Euclidian distance = 0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction21',
    reactants = ['S(659)(658)'],
    products = ['HO2(8)(9)', 'CC(=O)C=C=CO(5366)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction22',
    reactants = ['S(659)(658)'],
    products = ['CC(=O)C(CC=O)O[O](5957)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC(=CC=CO)OO[O](5733)'],
    products = ['S(659)(658)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(104,'s^-1'), n=3.21, Ea=(82.0482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CH3;R_O_R]
Euclidian distance = 1.41421356237
family: ketoenol"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CC(=CO[O])OC=CO(5958)'],
    products = ['S(659)(658)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(855.663,'s^-1'), n=2.935, Ea=(197.924,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CH3;R_O_C]
Euclidian distance = 2.2360679775
family: ketoenol"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CC(O)=C(C=CO)O[O](5959)'],
    products = ['S(659)(658)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(104,'s^-1'), n=3.21, Ea=(82.0482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C(O)C(C=CO)O[O](5960)'],
    products = ['S(659)(658)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O(4)(4)', 'CC(=O)C([O])C=CO(5961)'],
    products = ['S(659)(658)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

network(
    label = '367',
    isomers = [
        'S(659)(658)',
    ],
    reactants = [
        ('O2(2)(2)', 'S(245)(244)'),
        ('H(3)(3)', 'S(256)(255)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '367',
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

