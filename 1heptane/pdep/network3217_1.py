species(
    label = 'CCCC(CC[CH]OC[CH]O)OO(18034)',
    structure = SMILES('CCCC(CC[CH]OC[CH]O)OO'),
    E0 = (-292.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,180,180,180,214.606,926.085,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.10223,0.189989,-0.000249355,1.85574e-07,-5.64541e-11,-34876.2,56.131], Tmin=(100,'K'), Tmax=(798.894,'K')), NASAPolynomial(coeffs=[18.7093,0.0757771,-3.49197e-05,6.63788e-09,-4.61049e-13,-38521.2,-48.8007], Tmin=(798.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-292.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CCsJOCs)"""),
)

species(
    label = 'CH2CHOH(42)',
    structure = SMILES('C=CO'),
    E0 = (-138.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.72808,'amu*angstrom^2'), symmetry=1, barrier=(39.7321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28758,0.0197013,1.96383e-06,-1.9439e-08,1.02617e-11,-16537.3,14.1333], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[7.49818,0.0103957,-3.66891e-06,5.85206e-10,-3.47374e-14,-18164.3,-13.8388], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-138.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CCCC(CCC=O)OO(17459)',
    structure = SMILES('CCCC(CCC=O)OO'),
    E0 = (-399.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.184,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4428.2,'J/mol'), sigma=(7.53396,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.67 K, Pc=23.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.712285,0.111707,-0.000104517,5.79224e-08,-1.40607e-11,-47875,38.4312], Tmin=(100,'K'), Tmax=(955.014,'K')), NASAPolynomial(coeffs=[10.6432,0.0641444,-2.98101e-05,5.76987e-09,-4.08009e-13,-50043.8,-15.8297], Tmin=(955.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-399.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH)"""),
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
    label = 'CCCC(CC[CH]OCC=O)OO(18053)',
    structure = SMILES('CCCC(CC[CH]OCC=O)OO'),
    E0 = (-379.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.229,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.97413,0.161427,-0.000168879,9.82819e-08,-2.38091e-11,-45383.5,52.044], Tmin=(100,'K'), Tmax=(983.033,'K')), NASAPolynomial(coeffs=[18.8337,0.0726892,-3.34736e-05,6.45292e-09,-4.55426e-13,-49671,-52.7929], Tmin=(983.033,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-379.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(685.944,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(CCsJOCs)"""),
)

species(
    label = 'CCCC(CC=COC[CH]O)OO(18054)',
    structure = SMILES('CCCC(CC=COC[CH]O)OO'),
    E0 = (-396.014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3615,1310,387.5,850,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.229,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.93267,0.170597,-0.000177402,9.59678e-08,-2.07319e-11,-47340.4,53.9865], Tmin=(100,'K'), Tmax=(1121.42,'K')), NASAPolynomial(coeffs=[28.3819,0.0553343,-2.3227e-05,4.31347e-09,-2.99281e-13,-54588.1,-105.616], Tmin=(1121.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-396.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(685.944,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCsJOH)"""),
)

species(
    label = 'CCCC(CC[CH]OC=CO)OO(18055)',
    structure = SMILES('CCCC(CC[CH]OC=CO)OO'),
    E0 = (-420.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3615,1310,387.5,850,1000,180,180,180,250.79,1467.71,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150061,'amu*angstrom^2'), symmetry=1, barrier=(3.45021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150061,'amu*angstrom^2'), symmetry=1, barrier=(3.45021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150061,'amu*angstrom^2'), symmetry=1, barrier=(3.45021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150061,'amu*angstrom^2'), symmetry=1, barrier=(3.45021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150061,'amu*angstrom^2'), symmetry=1, barrier=(3.45021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150061,'amu*angstrom^2'), symmetry=1, barrier=(3.45021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150061,'amu*angstrom^2'), symmetry=1, barrier=(3.45021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150061,'amu*angstrom^2'), symmetry=1, barrier=(3.45021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150061,'amu*angstrom^2'), symmetry=1, barrier=(3.45021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150061,'amu*angstrom^2'), symmetry=1, barrier=(3.45021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150061,'amu*angstrom^2'), symmetry=1, barrier=(3.45021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.229,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-5.74285,0.18588,-0.000194345,9.93315e-08,-1.94563e-11,-50168.1,57.3379], Tmin=(100,'K'), Tmax=(1313.82,'K')), NASAPolynomial(coeffs=[43.3334,0.0298871,-8.73643e-06,1.33825e-09,-8.45328e-14,-62495.9,-190.661], Tmin=(1313.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-420.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(685.944,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCsJOC(O))"""),
)

species(
    label = '[CH2][CH]O(284)',
    structure = SMILES('[CH2][CH]O'),
    E0 = (143.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.217851,'amu*angstrom^2'), symmetry=1, barrier=(5.00882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0197382,'amu*angstrom^2'), symmetry=1, barrier=(14.867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84769,0.0229973,-1.85068e-05,7.77211e-09,-1.30725e-12,17300.5,13.0245], Tmin=(100,'K'), Tmax=(1418.34,'K')), NASAPolynomial(coeffs=[7.97636,0.00853337,-3.21015e-06,5.82141e-10,-3.99268e-14,15845.7,-13.5107], Tmin=(1418.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CJCO) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C(CCC)OO(331)',
    structure = SMILES('[CH2]C(CCC)OO'),
    E0 = (-55.0871,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.14,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.424678,0.0777384,-6.31189e-05,2.74429e-08,-4.97666e-12,-6496.19,29.796], Tmin=(100,'K'), Tmax=(1280.05,'K')), NASAPolynomial(coeffs=[13.2212,0.0377512,-1.6261e-05,3.03877e-09,-2.10445e-13,-9772.24,-35.0991], Tmin=(1280.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.0871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJCOOH)"""),
)

species(
    label = 'C=COC[CH]O(4708)',
    structure = SMILES('C=COC[CH]O'),
    E0 = (-165.037,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3025,407.5,1350,352.5,277.924,277.924,277.924,277.924],'cm^-1')),
        HinderedRotor(inertia=(0.329728,'amu*angstrom^2'), symmetry=1, barrier=(18.0732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329728,'amu*angstrom^2'), symmetry=1, barrier=(18.0732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329727,'amu*angstrom^2'), symmetry=1, barrier=(18.0732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329727,'amu*angstrom^2'), symmetry=1, barrier=(18.0732,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0175152,0.0723657,-7.52207e-05,3.78216e-08,-7.15704e-12,-19692.9,25.2358], Tmin=(100,'K'), Tmax=(1455.71,'K')), NASAPolynomial(coeffs=[19.9968,0.00912081,-1.45197e-06,9.9646e-11,-2.40605e-15,-24625.5,-75.6176], Tmin=(1455.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.037,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCsJOH)"""),
)

species(
    label = 'CCCC(CC[CH]OCC[O])OO(18056)',
    structure = SMILES('CCCC(CC[CH]OCC[O])OO'),
    E0 = (-246.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.60646,0.180296,-0.000237798,1.89897e-07,-6.34258e-11,-29434.3,55.6514], Tmin=(100,'K'), Tmax=(751.256,'K')), NASAPolynomial(coeffs=[13.5504,0.0851186,-4.01192e-05,7.69565e-09,-5.36681e-13,-31904.1,-21.4945], Tmin=(751.256,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-246.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = 'CCCC(C[CH]COC[CH]O)OO(18057)',
    structure = SMILES('CCCC(C[CH]COC[CH]O)OO'),
    E0 = (-272.867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.47855,0.176127,-0.000215271,1.5335e-07,-4.54749e-11,-32559.5,56.9514], Tmin=(100,'K'), Tmax=(813.088,'K')), NASAPolynomial(coeffs=[16.1756,0.0794321,-3.6874e-05,7.06854e-09,-4.94734e-13,-35755.4,-33.8003], Tmin=(813.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CCJCO)"""),
)

species(
    label = 'CCCC(CC[CH]O[CH]CO)OO(18058)',
    structure = SMILES('CCCC(CC[CH]O[CH]CO)OO'),
    E0 = (-292.155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,180,180,180,256.793,881.462,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153723,'amu*angstrom^2'), symmetry=1, barrier=(3.5344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153723,'amu*angstrom^2'), symmetry=1, barrier=(3.5344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153723,'amu*angstrom^2'), symmetry=1, barrier=(3.5344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153723,'amu*angstrom^2'), symmetry=1, barrier=(3.5344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153723,'amu*angstrom^2'), symmetry=1, barrier=(3.5344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153723,'amu*angstrom^2'), symmetry=1, barrier=(3.5344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153723,'amu*angstrom^2'), symmetry=1, barrier=(3.5344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153723,'amu*angstrom^2'), symmetry=1, barrier=(3.5344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153723,'amu*angstrom^2'), symmetry=1, barrier=(3.5344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153723,'amu*angstrom^2'), symmetry=1, barrier=(3.5344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153723,'amu*angstrom^2'), symmetry=1, barrier=(3.5344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153723,'amu*angstrom^2'), symmetry=1, barrier=(3.5344,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.07649,0.187513,-0.000234284,1.6232e-07,-4.5733e-11,-34856,55.9398], Tmin=(100,'K'), Tmax=(861.857,'K')), NASAPolynomial(coeffs=[20.9179,0.0715069,-3.23765e-05,6.13471e-09,-4.26656e-13,-39164.2,-60.9271], Tmin=(861.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-292.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CCsJOCs)"""),
)

species(
    label = 'CCCC(CCCO[CH][CH]O)OO(18059)',
    structure = SMILES('CCCC(CCCO[CH][CH]O)OO'),
    E0 = (-292.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,180,180,180,214.606,926.085,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153245,'amu*angstrom^2'), symmetry=1, barrier=(3.5234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.10223,0.189989,-0.000249355,1.85574e-07,-5.64541e-11,-34876.2,56.131], Tmin=(100,'K'), Tmax=(798.894,'K')), NASAPolynomial(coeffs=[18.7093,0.0757771,-3.49197e-05,6.63788e-09,-4.61049e-13,-38521.2,-48.8007], Tmin=(798.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-292.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CCsJOH)"""),
)

species(
    label = 'CCCC([CH]CCOC[CH]O)OO(18060)',
    structure = SMILES('CCCC([CH]CCOC[CH]O)OO'),
    E0 = (-272.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,180,180,180,180,976.009,1600,1821.46,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.09529,0.174741,-0.000196265,9.65272e-08,2.74212e-14,-32518.8,55.4994], Tmin=(100,'K'), Tmax=(585.159,'K')), NASAPolynomial(coeffs=[14.9492,0.0823395,-3.87307e-05,7.42713e-09,-5.18284e-13,-35160.4,-26.4126], Tmin=(585.159,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CCJCOOH)"""),
)

species(
    label = 'CCC[C](CCCOC[CH]O)OO(18061)',
    structure = SMILES('CCC[C](CCCOC[CH]O)OO'),
    E0 = (-285.87,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,360,370,350,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,3615,1310,387.5,850,1000,180,180,180,180,1012.91,1600,1786.14,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152239,'amu*angstrom^2'), symmetry=1, barrier=(3.50028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152239,'amu*angstrom^2'), symmetry=1, barrier=(3.50028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152239,'amu*angstrom^2'), symmetry=1, barrier=(3.50028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152239,'amu*angstrom^2'), symmetry=1, barrier=(3.50028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152239,'amu*angstrom^2'), symmetry=1, barrier=(3.50028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152239,'amu*angstrom^2'), symmetry=1, barrier=(3.50028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152239,'amu*angstrom^2'), symmetry=1, barrier=(3.50028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152239,'amu*angstrom^2'), symmetry=1, barrier=(3.50028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152239,'amu*angstrom^2'), symmetry=1, barrier=(3.50028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152239,'amu*angstrom^2'), symmetry=1, barrier=(3.50028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152239,'amu*angstrom^2'), symmetry=1, barrier=(3.50028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152239,'amu*angstrom^2'), symmetry=1, barrier=(3.50028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.87854,0.188352,-0.000266991,2.24845e-07,-7.72624e-11,-34113,56.1554], Tmin=(100,'K'), Tmax=(799.337,'K')), NASAPolynomial(coeffs=[13.0458,0.0856452,-4.04508e-05,7.71097e-09,-5.33486e-13,-36243.1,-18.1042], Tmin=(799.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-285.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CsJOOH) + radical(CCsJOH)"""),
)

species(
    label = 'CCCC(C[CH][CH]OCCO)OO(18062)',
    structure = SMILES('CCCC(C[CH][CH]OCCO)OO'),
    E0 = (-272.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,180,180,180,195.784,939.085,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15291,'amu*angstrom^2'), symmetry=1, barrier=(3.51569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15291,'amu*angstrom^2'), symmetry=1, barrier=(3.51569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15291,'amu*angstrom^2'), symmetry=1, barrier=(3.51569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15291,'amu*angstrom^2'), symmetry=1, barrier=(3.51569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15291,'amu*angstrom^2'), symmetry=1, barrier=(3.51569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15291,'amu*angstrom^2'), symmetry=1, barrier=(3.51569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15291,'amu*angstrom^2'), symmetry=1, barrier=(3.51569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15291,'amu*angstrom^2'), symmetry=1, barrier=(3.51569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15291,'amu*angstrom^2'), symmetry=1, barrier=(3.51569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15291,'amu*angstrom^2'), symmetry=1, barrier=(3.51569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15291,'amu*angstrom^2'), symmetry=1, barrier=(3.51569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15291,'amu*angstrom^2'), symmetry=1, barrier=(3.51569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.46447,0.173852,-0.000201206,1.31865e-07,-3.57282e-11,-32538.9,56.7982], Tmin=(100,'K'), Tmax=(888.371,'K')), NASAPolynomial(coeffs=[18.5032,0.0749391,-3.41933e-05,6.53128e-09,-4.57408e-13,-36442,-46.5827], Tmin=(888.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'CC[CH]C(CCCOC[CH]O)OO(18063)',
    structure = SMILES('CC[CH]C(CCCOC[CH]O)OO'),
    E0 = (-272.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,180,180,180,180,976.009,1600,1821.46,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152809,'amu*angstrom^2'), symmetry=1, barrier=(3.51338,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.09529,0.174741,-0.000196265,9.65272e-08,2.74212e-14,-32518.8,55.4994], Tmin=(100,'K'), Tmax=(585.159,'K')), NASAPolynomial(coeffs=[14.9492,0.0823395,-3.87307e-05,7.42713e-09,-5.18284e-13,-35160.4,-26.4126], Tmin=(585.159,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJCOOH) + radical(CCsJOH)"""),
)

species(
    label = 'CCCC(CCCOC[CH][O])OO(18064)',
    structure = SMILES('CCCC(CCCOC[CH][O])OO'),
    E0 = (-247.064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,180,180,180,180,729.762,1495.26,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154687,'amu*angstrom^2'), symmetry=1, barrier=(3.55655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154687,'amu*angstrom^2'), symmetry=1, barrier=(3.55655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154687,'amu*angstrom^2'), symmetry=1, barrier=(3.55655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154687,'amu*angstrom^2'), symmetry=1, barrier=(3.55655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154687,'amu*angstrom^2'), symmetry=1, barrier=(3.55655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154687,'amu*angstrom^2'), symmetry=1, barrier=(3.55655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154687,'amu*angstrom^2'), symmetry=1, barrier=(3.55655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154687,'amu*angstrom^2'), symmetry=1, barrier=(3.55655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154687,'amu*angstrom^2'), symmetry=1, barrier=(3.55655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154687,'amu*angstrom^2'), symmetry=1, barrier=(3.55655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154687,'amu*angstrom^2'), symmetry=1, barrier=(3.55655,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.63349,0.182853,-0.000253482,2.14405e-07,-7.48736e-11,-29454.4,55.843], Tmin=(100,'K'), Tmax=(784.801,'K')), NASAPolynomial(coeffs=[11.4985,0.0890921,-4.24778e-05,8.15279e-09,-5.671e-13,-31317.2,-10.229], Tmin=(784.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-247.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'CCCC(CCCOC[CH]O)O[O](18065)',
    structure = SMILES('CCCC(CCCOC[CH]O)O[O]'),
    E0 = (-320.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.68717,0.181393,-0.000243222,1.94157e-07,-6.40247e-11,-38313.9,55.9393], Tmin=(100,'K'), Tmax=(779.373,'K')), NASAPolynomial(coeffs=[14.5529,0.0811042,-3.73565e-05,7.07412e-09,-4.88984e-13,-40954.4,-26.2115], Tmin=(779.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-320.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOH)"""),
)

species(
    label = 'C[CH]CC(CCCOC[CH]O)OO(18066)',
    structure = SMILES('C[CH]CC(CCCOC[CH]O)OO'),
    E0 = (-278.323,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.7908,0.185729,-0.000260283,2.17181e-07,-7.40148e-11,-33207.8,58.2826], Tmin=(100,'K'), Tmax=(802.488,'K')), NASAPolynomial(coeffs=[13.1386,0.0844569,-3.94209e-05,7.4766e-09,-5.15855e-13,-35381.1,-16.2788], Tmin=(802.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-278.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(RCCJC)"""),
)

species(
    label = 'CCCC([CH]C[CH]OCCO)OO(18067)',
    structure = SMILES('CCCC([CH]C[CH]OCCO)OO'),
    E0 = (-272.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,180,180,180,191.064,949.093,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.7839,0.182944,-0.000232658,1.70287e-07,-5.13943e-11,-32468,57.7477], Tmin=(100,'K'), Tmax=(802.411,'K')), NASAPolynomial(coeffs=[17.3249,0.0777184,-3.59553e-05,6.86366e-09,-4.78535e-13,-35855.7,-39.4437], Tmin=(802.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]CCC(CCCOC[CH]O)OO(18068)',
    structure = SMILES('[CH2]CCC(CCCOC[CH]O)OO'),
    E0 = (-267.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.32758,0.177346,-0.000211099,1.31006e-07,-2.57611e-11,-31926.7,55.7601], Tmin=(100,'K'), Tmax=(625.805,'K')), NASAPolynomial(coeffs=[15.3638,0.0815553,-3.82548e-05,7.34081e-09,-5.13062e-13,-34729.8,-29.3586], Tmin=(625.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-267.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = 'CCC[C](CC[CH]OCCO)OO(18069)',
    structure = SMILES('CCC[C](CC[CH]OCCO)OO'),
    E0 = (-285.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,360,370,350,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,3615,1310,387.5,850,1000,180,180,180,180,970.149,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152962,'amu*angstrom^2'), symmetry=1, barrier=(3.5169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152962,'amu*angstrom^2'), symmetry=1, barrier=(3.5169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152962,'amu*angstrom^2'), symmetry=1, barrier=(3.5169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152962,'amu*angstrom^2'), symmetry=1, barrier=(3.5169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152962,'amu*angstrom^2'), symmetry=1, barrier=(3.5169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152962,'amu*angstrom^2'), symmetry=1, barrier=(3.5169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152962,'amu*angstrom^2'), symmetry=1, barrier=(3.5169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152962,'amu*angstrom^2'), symmetry=1, barrier=(3.5169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152962,'amu*angstrom^2'), symmetry=1, barrier=(3.5169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152962,'amu*angstrom^2'), symmetry=1, barrier=(3.5169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152962,'amu*angstrom^2'), symmetry=1, barrier=(3.5169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152962,'amu*angstrom^2'), symmetry=1, barrier=(3.5169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.87778,0.186147,-0.000252765,2.02596e-07,-6.69658e-11,-34091.7,56.0555], Tmin=(100,'K'), Tmax=(773.369,'K')), NASAPolynomial(coeffs=[15.1588,0.0815591,-3.80233e-05,7.23688e-09,-5.01618e-13,-36852.9,-29.7076], Tmin=(773.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-285.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(C2CsJOOH)"""),
)

species(
    label = 'CC[CH]C(CC[CH]OCCO)OO(18070)',
    structure = SMILES('CC[CH]C(CC[CH]OCCO)OO'),
    E0 = (-272.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,180,180,180,191.064,949.093,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152941,'amu*angstrom^2'), symmetry=1, barrier=(3.51643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.7839,0.182944,-0.000232658,1.70287e-07,-5.13943e-11,-32468,57.7477], Tmin=(100,'K'), Tmax=(802.411,'K')), NASAPolynomial(coeffs=[17.3249,0.0777184,-3.59553e-05,6.86366e-09,-4.78535e-13,-35855.7,-39.4437], Tmin=(802.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJCOOH) + radical(CCsJOCs)"""),
)

species(
    label = 'CCCC(CC[CH][O])OO(17462)',
    structure = SMILES('CCCC(CC[CH][O])OO'),
    E0 = (-70.4347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.184,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.590312,0.116872,-9.72251e-05,-8.20407e-09,5.32254e-11,-8321.76,39.241], Tmin=(100,'K'), Tmax=(521.994,'K')), NASAPolynomial(coeffs=[9.31624,0.0675691,-3.20149e-05,6.17245e-09,-4.32473e-13,-9718.52,-5.58476], Tmin=(521.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-70.4347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'CCCC(CCCOC=CO)OO(18071)',
    structure = SMILES('CCCC(CCCOC=CO)OO'),
    E0 = (-614.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.63184,0.167903,-0.000138396,4.04092e-08,1.8106e-12,-73536.9,53.8659], Tmin=(100,'K'), Tmax=(988.135,'K')), NASAPolynomial(coeffs=[37.548,0.041767,-1.46361e-05,2.59857e-09,-1.82023e-13,-84050.6,-160.144], Tmin=(988.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-614.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH)"""),
)

species(
    label = 'CCCC(CC=COCCO)OO(18072)',
    structure = SMILES('CCCC(CC=COCCO)OO'),
    E0 = (-576.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.91593,0.164862,-0.000156736,7.72425e-08,-1.52339e-11,-69021.6,54.3741], Tmin=(100,'K'), Tmax=(1221.99,'K')), NASAPolynomial(coeffs=[29.4081,0.0557803,-2.28371e-05,4.19231e-09,-2.88869e-13,-77165.8,-113.076], Tmin=(1221.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-576.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH)"""),
)

species(
    label = 'CCCC(CCCOCC=O)OO(18073)',
    structure = SMILES('CCCC(CCCOCC=O)OO'),
    E0 = (-559.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52345,0.152655,-0.000143585,7.65864e-08,-1.74771e-11,-67104.5,50.9802], Tmin=(100,'K'), Tmax=(1022.4,'K')), NASAPolynomial(coeffs=[16.0171,0.0801176,-3.7163e-05,7.19289e-09,-5.0879e-13,-70895.7,-38.878], Tmin=(1022.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-559.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH)"""),
)

species(
    label = 'CH2(S)(23)',
    structure = SMILES('[CH2]'),
    E0 = (419.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.36,2789.41,2993.36],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19195,-0.00230793,8.0509e-06,-6.60123e-09,1.95638e-12,50484.3,-0.754589], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.28556,0.00460255,-1.97412e-06,4.09548e-10,-3.34695e-14,50922.4,8.67684], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(419.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CCC(CC[CH]OC[CH]O)OO(18074)',
    structure = SMILES('CCC(CC[CH]OC[CH]O)OO'),
    E0 = (-268.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3615,1310,387.5,850,1000,180,180,180,325.428,1405.93,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149249,'amu*angstrom^2'), symmetry=1, barrier=(3.43153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149249,'amu*angstrom^2'), symmetry=1, barrier=(3.43153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149249,'amu*angstrom^2'), symmetry=1, barrier=(3.43153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149249,'amu*angstrom^2'), symmetry=1, barrier=(3.43153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149249,'amu*angstrom^2'), symmetry=1, barrier=(3.43153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149249,'amu*angstrom^2'), symmetry=1, barrier=(3.43153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149249,'amu*angstrom^2'), symmetry=1, barrier=(3.43153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149249,'amu*angstrom^2'), symmetry=1, barrier=(3.43153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149249,'amu*angstrom^2'), symmetry=1, barrier=(3.43153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149249,'amu*angstrom^2'), symmetry=1, barrier=(3.43153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149249,'amu*angstrom^2'), symmetry=1, barrier=(3.43153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (176.21,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.49743,0.175774,-0.000237847,1.81018e-07,-5.6017e-11,-32037.2,51.7093], Tmin=(100,'K'), Tmax=(787.334,'K')), NASAPolynomial(coeffs=[17.8737,0.0672017,-3.10037e-05,5.88017e-09,-4.07222e-13,-35402.5,-46.2848], Tmin=(787.334,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(636.057,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CCsJOCs)"""),
)

species(
    label = 'OH(5)',
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
    label = 'CCCC1CCC(OC[CH]O)O1(18075)',
    structure = SMILES('CCCC1CCC(OC[CH]O)O1'),
    E0 = (-479.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (173.229,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.5853,0.13049,-0.000103897,4.52096e-08,-8.04433e-12,-57374.4,45.4332], Tmin=(100,'K'), Tmax=(1336.24,'K')), NASAPolynomial(coeffs=[22.356,0.0558302,-2.00891e-05,3.39753e-09,-2.21769e-13,-64040,-82.124], Tmin=(1336.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-479.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(677.629,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C(CC(CCC)OO)OC[CH]O(18076)',
    structure = SMILES('[CH2]C(CC(CCC)OO)OC[CH]O'),
    E0 = (-272.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,180,180,180,334.244,1401.2,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149803,'amu*angstrom^2'), symmetry=1, barrier=(3.44426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149803,'amu*angstrom^2'), symmetry=1, barrier=(3.44426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149803,'amu*angstrom^2'), symmetry=1, barrier=(3.44426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149803,'amu*angstrom^2'), symmetry=1, barrier=(3.44426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149803,'amu*angstrom^2'), symmetry=1, barrier=(3.44426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149803,'amu*angstrom^2'), symmetry=1, barrier=(3.44426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149803,'amu*angstrom^2'), symmetry=1, barrier=(3.44426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149803,'amu*angstrom^2'), symmetry=1, barrier=(3.44426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149803,'amu*angstrom^2'), symmetry=1, barrier=(3.44426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149803,'amu*angstrom^2'), symmetry=1, barrier=(3.44426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149803,'amu*angstrom^2'), symmetry=1, barrier=(3.44426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149803,'amu*angstrom^2'), symmetry=1, barrier=(3.44426,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.67998,0.186419,-0.000228398,1.37915e-07,-2.27652e-11,-32551.2,54.6625], Tmin=(100,'K'), Tmax=(618.74,'K')), NASAPolynomial(coeffs=[17.5432,0.0782391,-3.64983e-05,6.96033e-09,-4.8373e-13,-35733.1,-42.0278], Tmin=(618.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)OC) + radical(CCsJOH)"""),
)

species(
    label = 'CCCC(CCC1OCC1O)OO(18036)',
    structure = SMILES('CCCC(CCC1OCC1O)OO'),
    E0 = (-525.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.37547,0.149766,-0.000130229,6.17785e-08,-1.19055e-11,-62876.5,51.8519], Tmin=(100,'K'), Tmax=(1244.32,'K')), NASAPolynomial(coeffs=[24.5732,0.0599227,-2.19245e-05,3.75277e-09,-2.4736e-13,-69832,-89.0935], Tmin=(1244.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-525.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(719.202,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane)"""),
)

species(
    label = 'CCCC([O])CCC(O)OC[CH]O(18077)',
    structure = SMILES('CCCC([O])CCC(O)OC[CH]O'),
    E0 = (-514.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.41793,0.174266,-0.000210132,1.46494e-07,-4.24189e-11,-61580.1,56.4434], Tmin=(100,'K'), Tmax=(832.416,'K')), NASAPolynomial(coeffs=[16.7191,0.077502,-3.57663e-05,6.84907e-09,-4.7947e-13,-64932.6,-37.0128], Tmin=(832.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-514.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C[CH]O(1827)',
    structure = SMILES('[O]C[CH]O'),
    E0 = (-10.3955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3025,407.5,1350,352.5,180,2261.05],'cm^-1')),
        HinderedRotor(inertia=(0.223718,'amu*angstrom^2'), symmetry=1, barrier=(5.14371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223266,'amu*angstrom^2'), symmetry=1, barrier=(5.13333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47466,0.0416325,-7.70032e-05,7.78856e-08,-2.90461e-11,-1203.08,16.406], Tmin=(100,'K'), Tmax=(877.571,'K')), NASAPolynomial(coeffs=[2.30008,0.0226909,-1.08908e-05,2.03334e-09,-1.36524e-13,-412.426,21.5557], Tmin=(877.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.3955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH]CCC(CCC)OO(17474)',
    structure = SMILES('[CH]CCC(CCC)OO'),
    E0 = (131.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.185,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.712022,0.105925,-9.11256e-05,4.29082e-08,-8.51487e-12,15996.1,37.3789], Tmin=(100,'K'), Tmax=(1171.5,'K')), NASAPolynomial(coeffs=[14.7379,0.0531728,-2.35826e-05,4.47195e-09,-3.12622e-13,12376.1,-39.6037], Tmin=(1171.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'HCOH(T)(285)',
    structure = SMILES('[CH]O'),
    E0 = (205.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,403.876,3308.82],'cm^-1')),
        HinderedRotor(inertia=(0.0103144,'amu*angstrom^2'), symmetry=1, barrier=(22.7121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75938,0.0029613,8.90411e-06,-1.35016e-08,5.39816e-12,24775.6,6.76286], Tmin=(100,'K'), Tmax=(940.429,'K')), NASAPolynomial(coeffs=[5.09112,0.00321239,-9.31686e-07,1.59615e-10,-1.15729e-14,24263.5,-0.971], Tmin=(940.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HCOH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]O[CH]CCC(CCC)OO(17553)',
    structure = SMILES('[CH2]O[CH]CCC(CCC)OO'),
    E0 = (-87.4643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (160.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.14421,0.140297,-0.000145454,8.42927e-08,-2.02529e-11,-10302.5,46.7062], Tmin=(100,'K'), Tmax=(993.983,'K')), NASAPolynomial(coeffs=[17.2695,0.0621715,-2.75568e-05,5.21852e-09,-3.64693e-13,-14161.9,-46.8366], Tmin=(993.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.4643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(CsJOCC) + radical(CCsJOCs)"""),
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
    E0 = (-292.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-131.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-177.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-201.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-220.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-185.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-163.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-164.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-175.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-126.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-132.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-223.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-245.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-243.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-196.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-231.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-217.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-163.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-175.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-138.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-163.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (73.0492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-228.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-283.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-283.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (151.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-237.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-114.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-284.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-227.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (121.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (118.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    products = ['CH2CHOH(42)', 'CCCC(CCC=O)OO(17459)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'CCCC(CC[CH]OCC=O)OO(18053)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2818 used for Od_CO-CsH;HJ
Exact match found for rate rule [Od_CO-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'CCCC(CC=COC[CH]O)OO(18054)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'CCCC(CC[CH]OC=CO)OO(18055)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.182e+10,'cm^3/(mol*s)'), n=0.859, Ea=(6.76971,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-OsH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]O(284)', 'CCCC(CCC=O)OO(17459)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(CCC)OO(331)', 'C=COC[CH]O(4708)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1071,'cm^3/(mol*s)'), n=2.72, Ea=(34.4762,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2780 used for Cds-HH_Cds-OsH;CsJ-CsHH
Exact match found for rate rule [Cds-HH_Cds-OsH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    products = ['CCCC(CC[CH]OCC[O])OO(18056)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4500,'s^-1'), n=2.62, Ea=(129.286,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 322 used for R2H_S;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CCCC(C[CH]COC[CH]O)OO(18057)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CCCC(CC[CH]O[CH]CO)OO(18058)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.253e-19,'s^-1'), n=8.985, Ea=(117.152,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_1H;Cs_H_out_H/NonDeO] + [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_1H] for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CCCC(CCCO[CH][CH]O)OO(18059)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.53164e+09,'s^-1'), n=1.09, Ea=(165.519,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CCCC([CH]CCOC[CH]O)OO(18060)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.73898e+09,'s^-1'), n=0.87, Ea=(139.955,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CCC[C](CCCOC[CH]O)OO(18061)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(304,'s^-1'), n=2.77, Ea=(62.3834,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_single;Cs_H_out_H/NonDeO] for rate rule [R4H_SSS;C_rad_out_NDMustO;Cs_H_out_H/NonDeO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CCCC(C[CH][CH]OCCO)OO(18062)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.0756983,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R5HJ_1;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CC[CH]C(CCCOC[CH]O)OO(18063)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.1016,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R5H_CCC;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CCCC(CCCOC[CH][O])OO(18064)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.19599e+09,'s^-1'), n=0.63, Ea=(50.8774,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_1;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    products = ['CCCC(CCCOC[CH]O)O[O](18065)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(46.1,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;C_rad_out_1H;XH_out] for rate rule [R6H_SSSSS;C_rad_out_H/NonDeO;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[CH]CC(CCCOC[CH]O)OO(18066)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(92.2,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R6H_SSSSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CCCC([CH]C[CH]OCCO)OO(18067)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [RnH;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO] for rate rule [R6HJ_2;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CCC(CCCOC[CH]O)OO(18068)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.58054e+06,'s^-1'), n=1.27667, Ea=(92.048,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;C_rad_out_2H;Cs_H_out_1H] for rate rule [R7H;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CCC[C](CC[CH]OCCO)OO(18069)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.31367e+08,'s^-1'), n=1.20512, Ea=(147.086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_NonDe;Cs_H_out_1H] for rate rule [R7HJ_3;C_rad_out_NDMustO;Cs_H_out_H/NonDeO]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC[CH]C(CC[CH]OCCO)OO(18070)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [RnH;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO] for rate rule [R8Hall;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]O(284)', 'CCCC(CC[CH][O])OO(17462)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    products = ['CCCC(CCCOC=CO)OO(18071)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    products = ['CCCC(CC=COCCO)OO(18072)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    products = ['CCCC(CCCOCC=O)OO(18073)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(S)(23)', 'CCC(CC[CH]OC[CH]O)OO(18074)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    products = ['OH(5)', 'CCCC1CCC(OC[CH]O)O1(18075)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.63e+10,'s^-1','*|/',1.41), n=0, Ea=(54.392,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4OO_SSS;C_sec_rad_intra;OOH] for rate rule [R4OO_SSS;C_rad/H/NonDeO_intra;OOH]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(CC(CCC)OO)OC[CH]O(18076)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    products = ['CCCC(CCC1OCC1O)OO(18036)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    products = ['CCCC([O])CCC(O)OC[CH]O(18077)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(65.2704,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4OOH_SSS;C_rad_out_1H] for rate rule [R4OOH_SSS;C_rad_out_H/NonDeO]
Euclidian distance = 1.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C[CH]O(1827)', '[CH]CCC(CCC)OO(17474)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['HCOH(T)(285)', '[CH2]O[CH]CCC(CCC)OO(17553)'],
    products = ['CCCC(CC[CH]OC[CH]O)OO(18034)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '3217',
    isomers = [
        'CCCC(CC[CH]OC[CH]O)OO(18034)',
    ],
    reactants = [
        ('CH2CHOH(42)', 'CCCC(CCC=O)OO(17459)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '3217',
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

