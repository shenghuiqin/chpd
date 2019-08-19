species(
    label = 'O=[C]CO[CH]C=O(1105)',
    structure = SMILES('[O]C=COC[C]=O'),
    E0 = (-164.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.20415,'amu*angstrom^2'), symmetry=1, barrier=(27.6858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20305,'amu*angstrom^2'), symmetry=1, barrier=(27.6604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20395,'amu*angstrom^2'), symmetry=1, barrier=(27.6811,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.781893,0.0506223,2.01007e-06,-6.11693e-08,3.3304e-11,-19680.8,24.3926], Tmin=(100,'K'), Tmax=(920.378,'K')), NASAPolynomial(coeffs=[25.3896,-0.00366686,4.67014e-06,-9.34381e-10,5.77036e-14,-26440.8,-104.4], Tmin=(920.378,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CsCJ=O)"""),
)

species(
    label = 'CH2CO(28)',
    structure = SMILES('C=C=O'),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'OCHCHO(48)',
    structure = SMILES('O=CC=O'),
    E0 = (-225.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,180,857.339,857.946,861.419,1717.87],'cm^-1')),
        HinderedRotor(inertia=(0.00964493,'amu*angstrom^2'), symmetry=1, barrier=(5.0669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.68412,0.000478013,4.26391e-05,-5.79018e-08,2.31669e-11,-27198.5,4.51187], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[8.72507,0.00633097,-2.35575e-06,3.89783e-10,-2.37487e-14,-29102.4,-20.3904], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-225.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""OCHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[O][CH]C1OCC1=O(4754)',
    structure = SMILES('[O][CH]C1OCC1=O'),
    E0 = (21.6095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78431,0.0397954,-1.16349e-05,-1.14922e-08,6.25145e-12,2686.32,24.3797], Tmin=(100,'K'), Tmax=(1102.27,'K')), NASAPolynomial(coeffs=[12.2647,0.019817,-9.01513e-06,1.78185e-09,-1.29135e-13,-720.874,-32.1774], Tmin=(1102.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.6095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'H(3)',
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
    label = '[O]C=COC=C=O(4755)',
    structure = SMILES('[O]C=COC=C=O'),
    E0 = (-124.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.27169,'amu*angstrom^2'), symmetry=1, barrier=(29.2387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27825,'amu*angstrom^2'), symmetry=1, barrier=(29.3895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01813,0.0513186,-1.99376e-05,-2.98551e-08,2.06906e-11,-14893.2,23.6763], Tmin=(100,'K'), Tmax=(911.12,'K')), NASAPolynomial(coeffs=[21.6464,-0.00098806,3.19512e-06,-6.98058e-10,4.5656e-14,-20240.1,-82.6368], Tmin=(911.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]COC=C=O(4756)',
    structure = SMILES('O=[C]COC=C=O'),
    E0 = (-124.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3010,987.5,1337.5,450,1655,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03577,'amu*angstrom^2'), symmetry=1, barrier=(23.8145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03712,'amu*angstrom^2'), symmetry=1, barrier=(23.8454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03678,'amu*angstrom^2'), symmetry=1, barrier=(23.8377,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.336903,0.0627665,-6.67288e-05,3.25598e-08,-5.8629e-12,-14855.8,25.3974], Tmin=(100,'K'), Tmax=(1572.81,'K')), NASAPolynomial(coeffs=[19.952,0.00247765,6.90943e-07,-2.22909e-10,1.64318e-14,-19739.2,-74.0263], Tmin=(1572.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C=C[O](2537)',
    structure = SMILES('[O]C=C[O]'),
    E0 = (-18.8461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,714.947,715.113],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7018,-0.0027881,5.18222e-05,-5.96735e-08,2.03775e-11,-2247.83,13.3237], Tmin=(100,'K'), Tmax=(1035.3,'K')), NASAPolynomial(coeffs=[6.65087,0.0113883,-5.76543e-06,1.26602e-09,-9.87778e-14,-4228.83,-7.62441], Tmin=(1035.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.8461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C=CO[CH]C=O(4757)',
    structure = SMILES('[O]C=COC=C[O]'),
    E0 = (-164.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.39708,'amu*angstrom^2'), symmetry=1, barrier=(32.1216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39822,'amu*angstrom^2'), symmetry=1, barrier=(32.1478,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.75379,0.0473164,2.13002e-05,-8.91986e-08,4.5705e-11,-19687,24.5388], Tmin=(100,'K'), Tmax=(907.428,'K')), NASAPolynomial(coeffs=[28.1154,-0.00907459,8.35811e-06,-1.69861e-09,1.11423e-13,-27296.7,-119.376], Tmin=(907.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]COC=[C]O(4758)',
    structure = SMILES('O=[C]COC=[C]O'),
    E0 = (-66.4735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3615,1277.5,1000,1855,455,950,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01195,'amu*angstrom^2'), symmetry=1, barrier=(23.2668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01202,'amu*angstrom^2'), symmetry=1, barrier=(23.2684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01146,'amu*angstrom^2'), symmetry=1, barrier=(23.2554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0111,'amu*angstrom^2'), symmetry=1, barrier=(23.2471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.826517,0.0755393,-8.32937e-05,4.0935e-08,-7.22246e-12,-7794.63,31.9118], Tmin=(100,'K'), Tmax=(1687.69,'K')), NASAPolynomial(coeffs=[23.0497,-0.00193049,4.11914e-06,-9.25482e-10,6.43462e-14,-12880.1,-86.963], Tmin=(1687.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.4735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(C=CJO)"""),
)

species(
    label = 'O=[C]CO[C]=CO(4759)',
    structure = SMILES('O=[C]CO[C]=CO'),
    E0 = (-66.4735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3615,1277.5,1000,1855,455,950,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01181,'amu*angstrom^2'), symmetry=1, barrier=(23.2634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01151,'amu*angstrom^2'), symmetry=1, barrier=(23.2566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01135,'amu*angstrom^2'), symmetry=1, barrier=(23.253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01154,'amu*angstrom^2'), symmetry=1, barrier=(23.2572,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.826496,0.0755391,-8.32932e-05,4.09346e-08,-7.22234e-12,-7794.63,31.9118], Tmin=(100,'K'), Tmax=(1687.71,'K')), NASAPolynomial(coeffs=[23.0487,-0.00192917,4.11849e-06,-9.25346e-10,6.43359e-14,-12879.5,-86.9574], Tmin=(1687.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.4735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C=[C]OCC=O(4760)',
    structure = SMILES('[O]C=[C]OCC=O'),
    E0 = (-79.6818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,205.789,205.951,207.632,208.569],'cm^-1')),
        HinderedRotor(inertia=(0.767041,'amu*angstrom^2'), symmetry=1, barrier=(23.4449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.78147,'amu*angstrom^2'), symmetry=1, barrier=(23.4432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779098,'amu*angstrom^2'), symmetry=1, barrier=(23.4295,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.940589,0.0520301,-1.54337e-05,-3.36728e-08,2.11811e-11,-9459.25,26.2725], Tmin=(100,'K'), Tmax=(929.61,'K')), NASAPolynomial(coeffs=[21.5233,0.0026118,1.14034e-06,-2.59288e-10,1.28692e-14,-14977.5,-80.6228], Tmin=(929.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.6818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[O][C]=COCC=O(4761)',
    structure = SMILES('[O][C]=COCC=O'),
    E0 = (-79.6818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,205.789,205.951,207.632,208.569],'cm^-1')),
        HinderedRotor(inertia=(0.767041,'amu*angstrom^2'), symmetry=1, barrier=(23.4449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.78147,'amu*angstrom^2'), symmetry=1, barrier=(23.4432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779098,'amu*angstrom^2'), symmetry=1, barrier=(23.4295,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.940589,0.0520301,-1.54337e-05,-3.36728e-08,2.11811e-11,-9459.25,26.2725], Tmin=(100,'K'), Tmax=(929.61,'K')), NASAPolynomial(coeffs=[21.5233,0.0026118,1.14034e-06,-2.59288e-10,1.28692e-14,-14977.5,-80.6228], Tmin=(929.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.6818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'O=[C][CH]OC=CO(4762)',
    structure = SMILES('O=[C][CH]OC=CO'),
    E0 = (-112.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13064,'amu*angstrom^2'), symmetry=1, barrier=(25.9956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12849,'amu*angstrom^2'), symmetry=1, barrier=(25.9463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13079,'amu*angstrom^2'), symmetry=1, barrier=(25.9992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13151,'amu*angstrom^2'), symmetry=1, barrier=(26.0157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95756,0.0896085,-0.000105421,5.31741e-08,-9.46183e-12,-13255,32.3237], Tmin=(100,'K'), Tmax=(1704.65,'K')), NASAPolynomial(coeffs=[27.3513,-0.00942452,8.34875e-06,-1.7333e-09,1.17969e-13,-18850.9,-111.812], Tmin=(1704.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(CCsJOC(O))"""),
)

species(
    label = 'C=[C][O](173)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97566e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02476e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CsCJ=O) + radical(CJC=O)"""),
)

species(
    label = '[O]C[C]=O(113)',
    structure = SMILES('[O]C[C]=O'),
    E0 = (65.8932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,298.422],'cm^-1')),
        HinderedRotor(inertia=(0.460283,'amu*angstrom^2'), symmetry=1, barrier=(28.9433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01457,0.0123696,2.06309e-05,-3.93785e-08,1.66633e-11,7969.01,14.0082], Tmin=(100,'K'), Tmax=(965.91,'K')), NASAPolynomial(coeffs=[10.9319,0.00281306,-6.04133e-07,1.77115e-10,-1.9139e-14,5355.84,-29.5236], Tmin=(965.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.8932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CsCJ=O) + radical(C=OCOJ)"""),
)

species(
    label = 'CHCHO(47)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (245.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06864,0.0187233,-1.21319e-05,-3.33727e-10,2.32882e-12,29739.4,14.7866], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.96288,0.00799899,-4.30606e-06,1.11076e-09,-1.11415e-13,28725.6,-5.17392], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(245.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'O=[C]COC1[CH]O1(4763)',
    structure = SMILES('O=[C]COC1[CH]O1'),
    E0 = (-12.2359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.73035,0.0607864,-6.33757e-05,3.35799e-08,-6.72016e-12,-1344.37,26.0072], Tmin=(100,'K'), Tmax=(1406.54,'K')), NASAPolynomial(coeffs=[15.0104,0.0116208,-1.81959e-06,7.93169e-11,2.90178e-15,-4515.22,-44.749], Tmin=(1406.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.2359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CCsJO) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C1[CH]OCC1=O(4764)',
    structure = SMILES('[O]C1[CH]OCC1=O'),
    E0 = (-41.9711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98534,0.0267379,4.04564e-05,-7.96494e-08,3.45482e-11,-4959.31,23.6588], Tmin=(100,'K'), Tmax=(941.124,'K')), NASAPolynomial(coeffs=[16.751,0.00848885,-1.39683e-06,2.49614e-10,-2.49119e-14,-9709.65,-57.1528], Tmin=(941.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.9711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + ring(Cyclopentane) + radical(CCsJOCs) + radical(C=OCOJ)"""),
)

species(
    label = 'O=C=COC=CO(4765)',
    structure = SMILES('O=C=COC=CO'),
    E0 = (-266.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.949793,0.0762712,-8.48029e-05,4.21361e-08,-7.47565e-12,-31821.5,29.3556], Tmin=(100,'K'), Tmax=(1696.54,'K')), NASAPolynomial(coeffs=[22.4404,-0.001969,4.7904e-06,-1.09314e-09,7.68104e-14,-36434.7,-86.0582], Tmin=(1696.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH)"""),
)

species(
    label = 'O=C=COCC=O(4766)',
    structure = SMILES('O=C=COCC=O'),
    E0 = (-279.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.982001,0.0544209,-3.02474e-05,-1.2637e-08,1.22567e-11,-33486,22.3888], Tmin=(100,'K'), Tmax=(942.02,'K')), NASAPolynomial(coeffs=[19.2668,0.0064494,-1.10464e-06,1.7272e-10,-1.57766e-14,-38247.3,-71.72], Tmin=(942.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH)"""),
)

species(
    label = 'O=C1COC=CO1(4767)',
    structure = SMILES('O=C1COC=CO1'),
    E0 = (-465.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83306,0.0143946,0.000114943,-1.73744e-07,7.05104e-11,-55889.5,12.8434], Tmin=(100,'K'), Tmax=(948.402,'K')), NASAPolynomial(coeffs=[24.9015,0.00158407,1.58448e-06,-1.34737e-10,-1.2265e-14,-64064.7,-117.257], Tmin=(948.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-465.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclohexane)"""),
)

species(
    label = 'CO(12)',
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
    label = '[CH2]OC=C[O](2549)',
    structure = SMILES('[CH2]OC=C[O]'),
    E0 = (-23.2817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17869,'amu*angstrom^2'), symmetry=1, barrier=(27.1004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17224,'amu*angstrom^2'), symmetry=1, barrier=(26.9521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87374,0.0271007,3.67127e-05,-8.65645e-08,4.08996e-11,-2705.16,18.2416], Tmin=(100,'K'), Tmax=(908.37,'K')), NASAPolynomial(coeffs=[21.3383,-0.00659169,6.44909e-06,-1.31006e-09,8.52988e-14,-8387.52,-85.6062], Tmin=(908.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.2817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COCJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C(C=O)C[C]=O(1104)',
    structure = SMILES('[O]C(C=O)C[C]=O'),
    E0 = (-76.0907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,769.311,4000],'cm^-1')),
        HinderedRotor(inertia=(0.469049,'amu*angstrom^2'), symmetry=1, barrier=(10.7844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.433925,'amu*angstrom^2'), symmetry=1, barrier=(10.7843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256679,'amu*angstrom^2'), symmetry=1, barrier=(10.7835,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4220.17,'J/mol'), sigma=(6.6134,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=659.18 K, Pc=33.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80208,0.0510391,-5.38874e-05,3.2043e-08,-7.9556e-12,-9074.71,27.9793], Tmin=(100,'K'), Tmax=(959.602,'K')), NASAPolynomial(coeffs=[8.42412,0.0234354,-1.0738e-05,2.0652e-09,-1.45531e-13,-10345.6,-3.695], Tmin=(959.602,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.0907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=OCOJ)"""),
)

species(
    label = '[C]=O(361)',
    structure = SMILES('[C]=O'),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08916,0.00200416,-1.61661e-05,2.55058e-08,-1.16424e-11,52802.7,4.52505], Tmin=(100,'K'), Tmax=(856.11,'K')), NASAPolynomial(coeffs=[0.961625,0.00569045,-3.48044e-06,7.19202e-10,-5.08041e-14,53738.7,21.4663], Tmin=(856.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'O(4)',
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
    label = '[CH]=COC[C]=O(4561)',
    structure = SMILES('[CH]=COC[C]=O'),
    E0 = (149.671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07186,'amu*angstrom^2'), symmetry=1, barrier=(24.6442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07213,'amu*angstrom^2'), symmetry=1, barrier=(24.6504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07107,'amu*angstrom^2'), symmetry=1, barrier=(24.6261,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33743,0.0453203,-1.34801e-05,-2.88325e-08,1.80932e-11,18109.4,21.5317], Tmin=(100,'K'), Tmax=(932.627,'K')), NASAPolynomial(coeffs=[19.13,0.00275851,6.9297e-07,-1.61928e-10,6.53065e-15,13322.8,-70.935], Tmin=(932.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CsCJ=O) + radical(Cds_P)"""),
)

species(
    label = 'O=[C]COC[C]=O(4768)',
    structure = SMILES('O=[C]COC[C]=O'),
    E0 = (-97.2089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1850,1860,440,470,900,1000,398.181,398.279,1459.08],'cm^-1')),
        HinderedRotor(inertia=(0.0761029,'amu*angstrom^2'), symmetry=1, barrier=(8.56027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.076167,'amu*angstrom^2'), symmetry=1, barrier=(8.56064,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203497,'amu*angstrom^2'), symmetry=1, barrier=(22.8876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203471,'amu*angstrom^2'), symmetry=1, barrier=(22.8871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80317,0.0515898,-4.86149e-05,2.44123e-08,-5.15165e-12,-11615.2,23.9609], Tmin=(100,'K'), Tmax=(1104.56,'K')), NASAPolynomial(coeffs=[9.21422,0.0247522,-1.21699e-05,2.41592e-09,-1.73198e-13,-13252.4,-12.5302], Tmin=(1104.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.2089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(CsCJ=O)"""),
)

species(
    label = 'O=C1CO[CH][CH]O1(4769)',
    structure = SMILES('O=C1CO[CH][CH]O1'),
    E0 = (-148.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52242,0.0369636,2.08474e-05,-5.72104e-08,2.43493e-11,-17770.1,18.3386], Tmin=(100,'K'), Tmax=(1014.32,'K')), NASAPolynomial(coeffs=[17.9535,0.0139014,-6.76578e-06,1.50328e-09,-1.2003e-13,-23250.3,-71.7484], Tmin=(1014.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + ring(Cyclohexanone) + radical(CCsJOCs) + radical(CCsJOC(O))"""),
)

species(
    label = 'O=CC1OCC1=O(1107)',
    structure = SMILES('O=CC1OCC1=O'),
    E0 = (-306.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23491,0.0246692,3.43812e-05,-6.0915e-08,2.40891e-11,-36772.2,20.9774], Tmin=(100,'K'), Tmax=(1001.6,'K')), NASAPolynomial(coeffs=[12.5707,0.0180335,-7.55972e-06,1.531e-09,-1.16145e-13,-40580.3,-37.5774], Tmin=(1001.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-306.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
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
    E0 = (-164.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (21.6095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (97.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (119.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-19.5444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-7.94626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (96.3057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (83.8719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-35.3732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-35.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-17.6776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (141.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (311.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (17.0377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-41.9711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-139.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-139.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-157.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-116.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (149.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (415.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (392.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-43.3582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (36.3651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-113.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-156.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=[C]CO[CH]C=O(1105)'],
    products = ['CH2CO(28)', 'OCHCHO(48)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=[C]CO[CH]C=O(1105)'],
    products = ['[O][CH]C1OCC1=O(4754)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.48517e+07,'s^-1'), n=1.03851, Ea=(186.364,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 184.2 to 186.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[O]C=COC=C=O(4755)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.82e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ck;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'O=[C]COC=C=O(4756)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2CO(28)', '[O]C=C[O](2537)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.6e+11,'cm^3/(mol*s)'), n=0, Ea=(60.12,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cd_R;O_rad/OneDe] for rate rule [Cds-HH_Ck;O_rad/OneDe]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C=CO[CH]C=O(4757)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.005257,'s^-1'), n=4.42, Ea=(156.9,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeO;XH_out] for rate rule [R2H_S;C_rad_out_H/NonDeO;CO_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=[C]COC=[C]O(4758)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=[C]CO[C]=CO(4759)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C=[C]OCC=O(4760)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SSS_OCs;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O][C]=COCC=O(4761)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;CO_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=[C][CH]OC=CO(4762)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(722272,'s^-1'), n=1.6737, Ea=(94.6126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C][O](173)', '[O]C=C[O](2537)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.27681e+06,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_rad/OneDe;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C[C]=O(113)', 'CHCHO(47)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.87866e+07,'m^3/(mol*s)'), n=-0.180109, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;O_rad/NonDe] for rate rule [Cd_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=[C]CO[CH]C=O(1105)'],
    products = ['O=[C]COC1[CH]O1(4763)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_O;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=[C]CO[CH]C=O(1105)'],
    products = ['[O]C1[CH]OCC1=O(4764)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(122.784,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 119.1 to 122.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=[C]CO[CH]C=O(1105)'],
    products = ['O=C=COC=CO(4765)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=[C]CO[CH]C=O(1105)'],
    products = ['O=C=COCC=O(4766)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=[C]CO[CH]C=O(1105)'],
    products = ['O=C1COC=CO1(4767)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CO(12)', '[CH2]OC=C[O](2549)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1230.59,'m^3/(mol*s)'), n=1.0523, Ea=(25.6182,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [COm;C_pri_rad] for rate rule [COm;C_rad/H2/O]
Euclidian distance = 1.0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=[C]CO[CH]C=O(1105)'],
    products = ['[O]C(C=O)C[C]=O(1104)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[C]=O(361)', '[CH2]OC=C[O](2549)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(4)', '[CH]=COC[C]=O(4561)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=[C][O](173)', 'OCHCHO(48)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.2635,'m^3/(mol*s)'), n=2.48, Ea=(21.7568,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-DeH;YJ] for rate rule [Od_CO-DeH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=[C]COC[C]=O(4768)'],
    products = ['O=[C]CO[CH]C=O(1105)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00282701,'s^-1'), n=4.515, Ea=(133.574,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=[C]CO[CH]C=O(1105)'],
    products = ['O=C1CO[CH][CH]O1(4769)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(50.8858,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra_CO] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_CO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=[C]CO[CH]C=O(1105)'],
    products = ['O=CC1OCC1=O(1107)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

network(
    label = '684',
    isomers = [
        'O=[C]CO[CH]C=O(1105)',
    ],
    reactants = [
        ('CH2CO(28)', 'OCHCHO(48)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '684',
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

