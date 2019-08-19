species(
    label = 'CCCC[C](CCOO)OC[CH]O(29159)',
    structure = SMILES('CCCC[C](CCOO)OC[CH]O'),
    E0 = (-289.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,360,370,350,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,3615,1310,387.5,850,1000,180,180,180,262.919,895.763,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154086,'amu*angstrom^2'), symmetry=1, barrier=(3.54274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154086,'amu*angstrom^2'), symmetry=1, barrier=(3.54274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154086,'amu*angstrom^2'), symmetry=1, barrier=(3.54274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154086,'amu*angstrom^2'), symmetry=1, barrier=(3.54274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154086,'amu*angstrom^2'), symmetry=1, barrier=(3.54274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154086,'amu*angstrom^2'), symmetry=1, barrier=(3.54274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154086,'amu*angstrom^2'), symmetry=1, barrier=(3.54274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154086,'amu*angstrom^2'), symmetry=1, barrier=(3.54274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154086,'amu*angstrom^2'), symmetry=1, barrier=(3.54274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154086,'amu*angstrom^2'), symmetry=1, barrier=(3.54274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154086,'amu*angstrom^2'), symmetry=1, barrier=(3.54274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154086,'amu*angstrom^2'), symmetry=1, barrier=(3.54274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.17089,0.197245,-0.000295247,2.55063e-07,-8.80277e-11,-34598.4,56.8544], Tmin=(100,'K'), Tmax=(825.501,'K')), NASAPolynomial(coeffs=[13.1183,0.0851889,-4.02442e-05,7.62548e-09,-5.23756e-13,-36489.2,-17.404], Tmin=(825.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CsJOCs) + radical(CCsJOH)"""),
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
    label = 'CCCCC(=O)CCOO(16766)',
    structure = SMILES('CCCCC(=O)CCOO'),
    E0 = (-419.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,375,552.5,462.5,1710,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4468.5,'J/mol'), sigma=(7.51027,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=697.97 K, Pc=23.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.021909,0.106478,-5.11605e-05,-1.11163e-07,1.3815e-10,-50327.2,35.0186], Tmin=(100,'K'), Tmax=(471.151,'K')), NASAPolynomial(coeffs=[8.18504,0.0687158,-3.25389e-05,6.25288e-09,-4.36324e-13,-51454.8,-2.15787], Tmin=(471.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-419.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs)"""),
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
    label = 'CCCC[C](CCOO)OCC=O(29337)',
    structure = SMILES('CCCC[C](CCOO)OCC=O'),
    E0 = (-377.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,2782.5,750,1395,475,1775,1000,360,370,350,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,336.996,1405.8,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149429,'amu*angstrom^2'), symmetry=1, barrier=(3.43566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149429,'amu*angstrom^2'), symmetry=1, barrier=(3.43566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149429,'amu*angstrom^2'), symmetry=1, barrier=(3.43566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149429,'amu*angstrom^2'), symmetry=1, barrier=(3.43566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149429,'amu*angstrom^2'), symmetry=1, barrier=(3.43566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149429,'amu*angstrom^2'), symmetry=1, barrier=(3.43566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149429,'amu*angstrom^2'), symmetry=1, barrier=(3.43566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149429,'amu*angstrom^2'), symmetry=1, barrier=(3.43566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149429,'amu*angstrom^2'), symmetry=1, barrier=(3.43566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149429,'amu*angstrom^2'), symmetry=1, barrier=(3.43566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149429,'amu*angstrom^2'), symmetry=1, barrier=(3.43566,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.229,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45244,0.159353,-0.000167768,7.70779e-08,3.2542e-12,-45130.7,50.7819], Tmin=(100,'K'), Tmax=(577.47,'K')), NASAPolynomial(coeffs=[12.1637,0.0841532,-4.00782e-05,7.76042e-09,-5.45798e-13,-47253,-15.4666], Tmin=(577.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-377.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(685.944,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(C2CsJOCs)"""),
)

species(
    label = 'CCCC=C(CCOO)OC[CH]O(29338)',
    structure = SMILES('CCCC=C(CCOO)OC[CH]O'),
    E0 = (-400.247,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,350,440,435,1725,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,180,180,180,230.764,1503.38,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.46997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.46997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.46997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.46997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.46997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.46997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.46997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.46997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.46997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.46997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.46997,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.229,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.43974,0.16802,-0.00018113,1.05367e-07,-2.49506e-11,-47874.3,53.004], Tmin=(100,'K'), Tmax=(1016.55,'K')), NASAPolynomial(coeffs=[22.7944,0.0647921,-2.88082e-05,5.4723e-09,-3.8342e-13,-53207.9,-73.9909], Tmin=(1016.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-400.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(685.944,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCsJOH)"""),
)

species(
    label = 'CCCCC(=CCOO)OC[CH]O(29339)',
    structure = SMILES('CCCCC(=CCOO)OC[CH]O'),
    E0 = (-394.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,350,440,435,1725,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,180,180,180,230.407,1505.29,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150977,'amu*angstrom^2'), symmetry=1, barrier=(3.47125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150977,'amu*angstrom^2'), symmetry=1, barrier=(3.47125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150977,'amu*angstrom^2'), symmetry=1, barrier=(3.47125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150977,'amu*angstrom^2'), symmetry=1, barrier=(3.47125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150977,'amu*angstrom^2'), symmetry=1, barrier=(3.47125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150977,'amu*angstrom^2'), symmetry=1, barrier=(3.47125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150977,'amu*angstrom^2'), symmetry=1, barrier=(3.47125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150977,'amu*angstrom^2'), symmetry=1, barrier=(3.47125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150977,'amu*angstrom^2'), symmetry=1, barrier=(3.47125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150977,'amu*angstrom^2'), symmetry=1, barrier=(3.47125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150977,'amu*angstrom^2'), symmetry=1, barrier=(3.47125,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.229,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.30479,0.166234,-0.000179527,1.05773e-07,-2.55014e-11,-47211.8,52.906], Tmin=(100,'K'), Tmax=(996.759,'K')), NASAPolynomial(coeffs=[21.473,0.0667995,-2.98892e-05,5.68916e-09,-3.98858e-13,-52151.3,-66.5517], Tmin=(996.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-394.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(685.944,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCsJOH)"""),
)

species(
    label = 'CCCC[C](CCOO)OC=CO(29340)',
    structure = SMILES('CCCC[C](CCOO)OC=CO'),
    E0 = (-413.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,360,370,350,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,180,180,180,185.297,1547.69,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151643,'amu*angstrom^2'), symmetry=1, barrier=(3.48658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151643,'amu*angstrom^2'), symmetry=1, barrier=(3.48658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151643,'amu*angstrom^2'), symmetry=1, barrier=(3.48658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151643,'amu*angstrom^2'), symmetry=1, barrier=(3.48658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151643,'amu*angstrom^2'), symmetry=1, barrier=(3.48658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151643,'amu*angstrom^2'), symmetry=1, barrier=(3.48658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151643,'amu*angstrom^2'), symmetry=1, barrier=(3.48658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151643,'amu*angstrom^2'), symmetry=1, barrier=(3.48658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151643,'amu*angstrom^2'), symmetry=1, barrier=(3.48658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151643,'amu*angstrom^2'), symmetry=1, barrier=(3.48658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151643,'amu*angstrom^2'), symmetry=1, barrier=(3.48658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.229,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.89637,0.17825,-0.000185178,9.61495e-08,-1.94392e-11,-49451.7,55.6244], Tmin=(100,'K'), Tmax=(1214.3,'K')), NASAPolynomial(coeffs=[36.9884,0.0402759,-1.47392e-05,2.57495e-09,-1.73795e-13,-59623.6,-154.577], Tmin=(1214.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-413.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(685.944,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C2CsJOC(O))"""),
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
    label = 'npropyl(83)',
    structure = SMILES('[CH2]CC'),
    E0 = (87.0621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0928812,'amu*angstrom^2'), symmetry=1, barrier=(2.13552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.092914,'amu*angstrom^2'), symmetry=1, barrier=(2.13628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0877,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.02815,0.0147023,2.4051e-05,-3.66738e-08,1.38611e-11,10512.1,12.4699], Tmin=(100,'K'), Tmax=(984.464,'K')), NASAPolynomial(coeffs=[6.16543,0.0184495,-6.79029e-06,1.23049e-09,-8.63866e-14,9095.06,-6.67607], Tmin=(984.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.0621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""npropyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C(CCOO)OC[CH]O(29341)',
    structure = SMILES('C=C(CCOO)OC[CH]O'),
    E0 = (-317.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,350,440,435,1725,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,200,800,1000,1200,1400,1600],'cm^-1')),
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
    molecularWeight = (147.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67308,0.127903,-0.00015188,9.45435e-08,-2.35033e-11,-38020.3,39.6952], Tmin=(100,'K'), Tmax=(979.236,'K')), NASAPolynomial(coeffs=[19.7196,0.0405188,-1.80263e-05,3.41694e-09,-2.38939e-13,-42210.1,-63.0638], Tmin=(979.236,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOH)"""),
)

species(
    label = 'CH2OOH(35)',
    structure = SMILES('[CH2]OO'),
    E0 = (52.1952,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(2.16183,'amu*angstrom^2'), symmetry=1, barrier=(49.7048,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0333,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3467.13,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.83127,-0.00351771,4.54551e-05,-5.66903e-08,2.21633e-11,6061.87,-0.579143], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.98746,0.00900484,-3.24367e-06,5.24325e-10,-3.13587e-14,5012.58,-10.2619], Tmin=(1000,'K'), Tmax=(2500,'K'))], Tmin=(200,'K'), Tmax=(2500,'K'), E0=(52.1952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), label="""CH2OOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=C(CCCC)OC[CH]O(29342)',
    structure = SMILES('C=C(CCCC)OC[CH]O'),
    E0 = (-277.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (143.203,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.87911,0.124698,-0.00011751,5.86883e-08,-1.18504e-11,-33104,41.4457], Tmin=(100,'K'), Tmax=(1188.15,'K')), NASAPolynomial(coeffs=[21.2537,0.0468191,-1.9189e-05,3.52004e-09,-2.42258e-13,-38601,-74.1446], Tmin=(1188.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-277.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOH)"""),
)

species(
    label = 'CCCC[C](CCOO)OCC[O](29343)',
    structure = SMILES('CCCC[C](CCOO)OCC[O]'),
    E0 = (-244.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.58137,0.186302,-0.000278529,2.51287e-07,-9.07547e-11,-29160.4,56.0473], Tmin=(100,'K'), Tmax=(820.35,'K')), NASAPolynomial(coeffs=[7.8481,0.0947341,-4.55673e-05,8.71354e-09,-6.01968e-13,-29829.7,10.52], Tmin=(820.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CsJOCs) + radical(CCOJ)"""),
)

species(
    label = 'CCC[CH]C(CCOO)OC[CH]O(29344)',
    structure = SMILES('CCC[CH]C(CCOO)OC[CH]O'),
    E0 = (-270.775,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.47855,0.176127,-0.000215271,1.5335e-07,-4.54749e-11,-32307.9,56.9514], Tmin=(100,'K'), Tmax=(813.088,'K')), NASAPolynomial(coeffs=[16.1756,0.0794321,-3.6874e-05,7.06854e-09,-4.94734e-13,-35503.8,-33.8003], Tmin=(813.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-270.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CCJCO)"""),
)

species(
    label = 'CCCCC([CH]COO)OC[CH]O(29345)',
    structure = SMILES('CCCCC([CH]COO)OC[CH]O'),
    E0 = (-270.262,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.09529,0.174741,-0.000196265,9.65272e-08,2.74212e-14,-32267.2,55.4994], Tmin=(100,'K'), Tmax=(585.159,'K')), NASAPolynomial(coeffs=[14.9492,0.0823395,-3.87307e-05,7.42713e-09,-5.18284e-13,-34908.8,-26.4126], Tmin=(585.159,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-270.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJCOOH) + radical(CCsJOH)"""),
)

species(
    label = 'CCCC[C](CCOO)O[CH]CO(29346)',
    structure = SMILES('CCCC[C](CCOO)O[CH]CO'),
    E0 = (-289.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,360,370,350,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,3615,1310,387.5,850,1000,180,180,180,312.495,841.664,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154633,'amu*angstrom^2'), symmetry=1, barrier=(3.55531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154633,'amu*angstrom^2'), symmetry=1, barrier=(3.55531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154633,'amu*angstrom^2'), symmetry=1, barrier=(3.55531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154633,'amu*angstrom^2'), symmetry=1, barrier=(3.55531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154633,'amu*angstrom^2'), symmetry=1, barrier=(3.55531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154633,'amu*angstrom^2'), symmetry=1, barrier=(3.55531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154633,'amu*angstrom^2'), symmetry=1, barrier=(3.55531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154633,'amu*angstrom^2'), symmetry=1, barrier=(3.55531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154633,'amu*angstrom^2'), symmetry=1, barrier=(3.55531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154633,'amu*angstrom^2'), symmetry=1, barrier=(3.55531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154633,'amu*angstrom^2'), symmetry=1, barrier=(3.55531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154633,'amu*angstrom^2'), symmetry=1, barrier=(3.55531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.21458,0.195626,-0.00028338,2.36371e-07,-7.94952e-11,-34575.2,56.9101], Tmin=(100,'K'), Tmax=(814.596,'K')), NASAPolynomial(coeffs=[15.3399,0.0809048,-3.76968e-05,7.12205e-09,-4.89389e-13,-37140.6,-29.6108], Tmin=(814.596,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CsJOCs) + radical(CCsJOCs)"""),
)

species(
    label = 'CCCCC(CCOO)O[CH][CH]O(29347)',
    structure = SMILES('CCCCC(CCOO)O[CH][CH]O'),
    E0 = (-290.221,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.10223,0.189989,-0.000249355,1.85574e-07,-5.64541e-11,-34624.6,56.131], Tmin=(100,'K'), Tmax=(798.894,'K')), NASAPolynomial(coeffs=[18.7093,0.0757771,-3.49197e-05,6.63788e-09,-4.61049e-13,-38269.5,-48.8007], Tmin=(798.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CCsJOH)"""),
)

species(
    label = 'CC[CH]CC(CCOO)OC[CH]O(29348)',
    structure = SMILES('CC[CH]CC(CCOO)OC[CH]O'),
    E0 = (-276.219,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.76645,0.186373,-0.000266145,2.26557e-07,-7.83361e-11,-32956.8,58.0222], Tmin=(100,'K'), Tmax=(809.509,'K')), NASAPolynomial(coeffs=[12.0939,0.0863222,-4.05793e-05,7.70828e-09,-5.31735e-13,-34814.2,-10.7554], Tmin=(809.509,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(RCCJCC)"""),
)

species(
    label = 'CCCCC(C[CH]OO)OC[CH]O(29349)',
    structure = SMILES('CCCCC(C[CH]OO)OC[CH]O'),
    E0 = (-282.095,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.09044,0.174988,-0.000194529,9.15043e-08,3.99803e-12,-33690.9,53.3345], Tmin=(100,'K'), Tmax=(580.732,'K')), NASAPolynomial(coeffs=[14.9095,0.0828875,-3.89834e-05,7.47191e-09,-5.21166e-13,-36319.1,-28.3509], Tmin=(580.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-282.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOOH) + radical(CCsJOH)"""),
)

species(
    label = 'C[CH]CCC(CCOO)OC[CH]O(29350)',
    structure = SMILES('C[CH]CCC(CCOO)OC[CH]O'),
    E0 = (-276.231,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.7908,0.185729,-0.000260283,2.17181e-07,-7.40148e-11,-32956.1,58.2826], Tmin=(100,'K'), Tmax=(802.488,'K')), NASAPolynomial(coeffs=[13.1386,0.0844569,-3.94209e-05,7.4766e-09,-5.15855e-13,-35129.5,-16.2788], Tmin=(802.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(RCCJC)"""),
)

species(
    label = 'CCC[CH][C](CCOO)OCCO(29351)',
    structure = SMILES('CCC[CH][C](CCOO)OCCO'),
    E0 = (-270.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,360,370,350,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,3615,1310,387.5,850,1000,180,180,180,253.241,903.027,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153924,'amu*angstrom^2'), symmetry=1, barrier=(3.53903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153924,'amu*angstrom^2'), symmetry=1, barrier=(3.53903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153924,'amu*angstrom^2'), symmetry=1, barrier=(3.53903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153924,'amu*angstrom^2'), symmetry=1, barrier=(3.53903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153924,'amu*angstrom^2'), symmetry=1, barrier=(3.53903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153924,'amu*angstrom^2'), symmetry=1, barrier=(3.53903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153924,'amu*angstrom^2'), symmetry=1, barrier=(3.53903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153924,'amu*angstrom^2'), symmetry=1, barrier=(3.53903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153924,'amu*angstrom^2'), symmetry=1, barrier=(3.53903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153924,'amu*angstrom^2'), symmetry=1, barrier=(3.53903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153924,'amu*angstrom^2'), symmetry=1, barrier=(3.53903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153924,'amu*angstrom^2'), symmetry=1, barrier=(3.53903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.6767,0.182892,-0.000253857,2.1113e-07,-7.20752e-11,-32254.9,58.0312], Tmin=(100,'K'), Tmax=(793.912,'K')), NASAPolynomial(coeffs=[12.9291,0.0843369,-3.95168e-05,7.51995e-09,-5.2029e-13,-34422.4,-15.2945], Tmin=(793.912,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-270.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'CCCC[C]([CH]COO)OCCO(29352)',
    structure = SMILES('CCCC[C]([CH]COO)OCCO'),
    E0 = (-269.856,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,360,370,350,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,3615,1310,387.5,850,1000,180,180,180,241.534,919.545,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15385,'amu*angstrom^2'), symmetry=1, barrier=(3.53731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15385,'amu*angstrom^2'), symmetry=1, barrier=(3.53731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15385,'amu*angstrom^2'), symmetry=1, barrier=(3.53731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15385,'amu*angstrom^2'), symmetry=1, barrier=(3.53731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15385,'amu*angstrom^2'), symmetry=1, barrier=(3.53731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15385,'amu*angstrom^2'), symmetry=1, barrier=(3.53731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15385,'amu*angstrom^2'), symmetry=1, barrier=(3.53731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15385,'amu*angstrom^2'), symmetry=1, barrier=(3.53731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15385,'amu*angstrom^2'), symmetry=1, barrier=(3.53731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15385,'amu*angstrom^2'), symmetry=1, barrier=(3.53731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15385,'amu*angstrom^2'), symmetry=1, barrier=(3.53731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15385,'amu*angstrom^2'), symmetry=1, barrier=(3.53731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.89693,0.190785,-0.000280918,2.4339e-07,-8.47954e-11,-32188.3,58.6265], Tmin=(100,'K'), Tmax=(818.593,'K')), NASAPolynomial(coeffs=[11.8169,0.0869792,-4.11885e-05,7.82894e-09,-5.39343e-13,-33855.6,-8.50844], Tmin=(818.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-269.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CsJOCs) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]CCCC(CCOO)OC[CH]O(29353)',
    structure = SMILES('[CH2]CCCC(CCOO)OC[CH]O'),
    E0 = (-265.431,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.32758,0.177346,-0.000211099,1.31006e-07,-2.57611e-11,-31675.1,55.7601], Tmin=(100,'K'), Tmax=(625.805,'K')), NASAPolynomial(coeffs=[15.3638,0.0815553,-3.82548e-05,7.34081e-09,-5.13062e-13,-34478.2,-29.3586], Tmin=(625.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-265.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = 'CCCCC(CCO[O])OC[CH]O(29354)',
    structure = SMILES('CCCCC(CCO[O])OC[CH]O'),
    E0 = (-318.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,492.5,1135,1000,180,180,180,451.77,703.698,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156381,'amu*angstrom^2'), symmetry=1, barrier=(3.5955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156381,'amu*angstrom^2'), symmetry=1, barrier=(3.5955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156381,'amu*angstrom^2'), symmetry=1, barrier=(3.5955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156381,'amu*angstrom^2'), symmetry=1, barrier=(3.5955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156381,'amu*angstrom^2'), symmetry=1, barrier=(3.5955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156381,'amu*angstrom^2'), symmetry=1, barrier=(3.5955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156381,'amu*angstrom^2'), symmetry=1, barrier=(3.5955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156381,'amu*angstrom^2'), symmetry=1, barrier=(3.5955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156381,'amu*angstrom^2'), symmetry=1, barrier=(3.5955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156381,'amu*angstrom^2'), symmetry=1, barrier=(3.5955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156381,'amu*angstrom^2'), symmetry=1, barrier=(3.5955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.68717,0.181393,-0.000243222,1.94157e-07,-6.40247e-11,-38062.3,55.9393], Tmin=(100,'K'), Tmax=(779.373,'K')), NASAPolynomial(coeffs=[14.5529,0.0811042,-3.73565e-05,7.07412e-09,-4.88984e-13,-40702.8,-26.2115], Tmin=(779.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-318.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(ROOJ)"""),
)

species(
    label = 'CCCCC(CCOO)OC[CH][O](29355)',
    structure = SMILES('CCCCC(CCOO)OC[CH][O]'),
    E0 = (-244.972,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.63349,0.182853,-0.000253482,2.14405e-07,-7.48736e-11,-29202.8,55.843], Tmin=(100,'K'), Tmax=(784.801,'K')), NASAPolynomial(coeffs=[11.4985,0.0890921,-4.24778e-05,8.15279e-09,-5.671e-13,-31065.6,-10.229], Tmin=(784.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'CC[CH]C[C](CCOO)OCCO(29356)',
    structure = SMILES('CC[CH]C[C](CCOO)OCCO'),
    E0 = (-275.813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,360,370,350,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,3615,1310,387.5,850,1000,180,180,180,180,1001.66,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153343,'amu*angstrom^2'), symmetry=1, barrier=(3.52566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153343,'amu*angstrom^2'), symmetry=1, barrier=(3.52566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153343,'amu*angstrom^2'), symmetry=1, barrier=(3.52566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153343,'amu*angstrom^2'), symmetry=1, barrier=(3.52566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153343,'amu*angstrom^2'), symmetry=1, barrier=(3.52566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153343,'amu*angstrom^2'), symmetry=1, barrier=(3.52566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153343,'amu*angstrom^2'), symmetry=1, barrier=(3.52566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153343,'amu*angstrom^2'), symmetry=1, barrier=(3.52566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153343,'amu*angstrom^2'), symmetry=1, barrier=(3.52566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153343,'amu*angstrom^2'), symmetry=1, barrier=(3.52566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153343,'amu*angstrom^2'), symmetry=1, barrier=(3.52566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153343,'amu*angstrom^2'), symmetry=1, barrier=(3.52566,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.68492,0.189448,-0.000289759,2.61433e-07,-9.33192e-11,-32915.6,58.1228], Tmin=(100,'K'), Tmax=(835.73,'K')), NASAPolynomial(coeffs=[8.33555,0.0921587,-4.37856e-05,8.29747e-09,-5.69015e-13,-33536.4,10.5944], Tmin=(835.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-275.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(C2CsJOCs)"""),
)

species(
    label = 'CCCC[C](C[CH]OO)OCCO(29357)',
    structure = SMILES('CCCC[C](C[CH]OO)OCCO'),
    E0 = (-281.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,360,370,350,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,3615,1310,387.5,850,1000,180,180,180,180,996.669,1600,1812.87,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152835,'amu*angstrom^2'), symmetry=1, barrier=(3.51399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152835,'amu*angstrom^2'), symmetry=1, barrier=(3.51399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152835,'amu*angstrom^2'), symmetry=1, barrier=(3.51399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152835,'amu*angstrom^2'), symmetry=1, barrier=(3.51399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152835,'amu*angstrom^2'), symmetry=1, barrier=(3.51399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152835,'amu*angstrom^2'), symmetry=1, barrier=(3.51399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152835,'amu*angstrom^2'), symmetry=1, barrier=(3.51399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152835,'amu*angstrom^2'), symmetry=1, barrier=(3.51399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152835,'amu*angstrom^2'), symmetry=1, barrier=(3.51399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152835,'amu*angstrom^2'), symmetry=1, barrier=(3.51399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152835,'amu*angstrom^2'), symmetry=1, barrier=(3.51399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152835,'amu*angstrom^2'), symmetry=1, barrier=(3.51399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.92848,0.191615,-0.000282242,2.44707e-07,-8.52889e-11,-33610.4,56.5841], Tmin=(100,'K'), Tmax=(819.298,'K')), NASAPolynomial(coeffs=[11.7628,0.0875521,-4.14557e-05,7.87719e-09,-5.42514e-13,-35260.2,-10.3667], Tmin=(819.298,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-281.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOOH) + radical(C2CsJOCs)"""),
)

species(
    label = 'C[CH]CC[C](CCOO)OCCO(29358)',
    structure = SMILES('C[CH]CC[C](CCOO)OCCO'),
    E0 = (-275.825,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,360,370,350,2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,3615,1310,387.5,850,1000,180,180,180,180,999.786,1600,1826.71,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153272,'amu*angstrom^2'), symmetry=1, barrier=(3.52403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153272,'amu*angstrom^2'), symmetry=1, barrier=(3.52403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153272,'amu*angstrom^2'), symmetry=1, barrier=(3.52403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153272,'amu*angstrom^2'), symmetry=1, barrier=(3.52403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153272,'amu*angstrom^2'), symmetry=1, barrier=(3.52403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153272,'amu*angstrom^2'), symmetry=1, barrier=(3.52403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153272,'amu*angstrom^2'), symmetry=1, barrier=(3.52403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153272,'amu*angstrom^2'), symmetry=1, barrier=(3.52403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153272,'amu*angstrom^2'), symmetry=1, barrier=(3.52403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153272,'amu*angstrom^2'), symmetry=1, barrier=(3.52403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153272,'amu*angstrom^2'), symmetry=1, barrier=(3.52403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153272,'amu*angstrom^2'), symmetry=1, barrier=(3.52403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.72165,0.188963,-0.000284528,2.52982e-07,-8.94419e-11,-32914.5,58.4267], Tmin=(100,'K'), Tmax=(833.347,'K')), NASAPolynomial(coeffs=[9.41411,0.0902323,-4.25906e-05,8.05684e-09,-5.52376e-13,-33864.9,4.88255], Tmin=(833.347,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-275.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CsJOCs) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CCC[C](CCOO)OCCO(29359)',
    structure = SMILES('[CH2]CCC[C](CCOO)OCCO'),
    E0 = (-265.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2757.69,2765.38,2773.08,2780.77,2788.46,2796.15,2803.85,2811.54,2819.23,2826.92,2834.62,2842.31,2850,1425,1429.17,1433.33,1437.5,1441.67,1445.83,1450,1225,1233.33,1241.67,1250,1258.33,1266.67,1275,1270,1281.67,1293.33,1305,1316.67,1328.33,1340,700,716.667,733.333,750,766.667,783.333,800,300,316.667,333.333,350,366.667,383.333,400,360,370,350,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,3615,1277.5,1000,180,180,180,206.918,953.285,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153417,'amu*angstrom^2'), symmetry=1, barrier=(3.52737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153417,'amu*angstrom^2'), symmetry=1, barrier=(3.52737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153417,'amu*angstrom^2'), symmetry=1, barrier=(3.52737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153417,'amu*angstrom^2'), symmetry=1, barrier=(3.52737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153417,'amu*angstrom^2'), symmetry=1, barrier=(3.52737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153417,'amu*angstrom^2'), symmetry=1, barrier=(3.52737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153417,'amu*angstrom^2'), symmetry=1, barrier=(3.52737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153417,'amu*angstrom^2'), symmetry=1, barrier=(3.52737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153417,'amu*angstrom^2'), symmetry=1, barrier=(3.52737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153417,'amu*angstrom^2'), symmetry=1, barrier=(3.52737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153417,'amu*angstrom^2'), symmetry=1, barrier=(3.52737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153417,'amu*angstrom^2'), symmetry=1, barrier=(3.52737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.86972,0.189302,-0.000274743,2.35478e-07,-8.15709e-11,-31607.3,58.0112], Tmin=(100,'K'), Tmax=(813.059,'K')), NASAPolynomial(coeffs=[12.2995,0.0860847,-4.06516e-05,7.72875e-09,-5.33006e-13,-33454.3,-11.8383], Tmin=(813.059,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-265.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(C2CsJOCs)"""),
)

species(
    label = 'CCCC[C](CCO[O])OCCO(29360)',
    structure = SMILES('CCCC[C](CCO[O])OCCO'),
    E0 = (-318.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.66311,0.185219,-0.000269833,2.33501e-07,-8.1196e-11,-38018.7,56.2414], Tmin=(100,'K'), Tmax=(826.37,'K')), NASAPolynomial(coeffs=[10.9413,0.0866746,-4.04024e-05,7.62411e-09,-5.22933e-13,-39481.5,-5.67712], Tmin=(826.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-318.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CsJOCs) + radical(ROOJ)"""),
)

species(
    label = 'CCCC[C]([O])CCOO(16758)',
    structure = SMILES('CCCC[C]([O])CCOO'),
    E0 = (-69.4491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,360,370,350,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38813,0.12719,-0.000146358,1.00264e-07,-2.90472e-11,-8166.49,40.7595], Tmin=(100,'K'), Tmax=(826.048,'K')), NASAPolynomial(coeffs=[11.9547,0.0625799,-2.90356e-05,5.57995e-09,-3.91711e-13,-10370.9,-21.0624], Tmin=(826.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.4491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'CCCCC(CCOO)OC=CO(29361)',
    structure = SMILES('CCCCC(CCOO)OC=CO'),
    E0 = (-614.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.63184,0.167903,-0.000138396,4.04092e-08,1.8106e-12,-73536.9,53.8659], Tmin=(100,'K'), Tmax=(988.135,'K')), NASAPolynomial(coeffs=[37.548,0.041767,-1.46361e-05,2.59857e-09,-1.82023e-13,-84050.6,-160.144], Tmin=(988.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-614.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH)"""),
)

species(
    label = 'CCCC=C(CCOO)OCCO(29362)',
    structure = SMILES('CCCC=C(CCOO)OCCO'),
    E0 = (-580.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.26875,0.160583,-0.000155049,8.03343e-08,-1.70464e-11,-69562.3,52.8305], Tmin=(100,'K'), Tmax=(1123.31,'K')), NASAPolynomial(coeffs=[23.2652,0.0660979,-2.88789e-05,5.45382e-09,-3.81139e-13,-75523.4,-78.2657], Tmin=(1123.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-580.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH)"""),
)

species(
    label = 'CCCCC(=CCOO)OCCO(29363)',
    structure = SMILES('CCCCC(=CCOO)OCCO'),
    E0 = (-574.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.10382,0.158455,-0.000152309,7.93498e-08,-1.70395e-11,-68901.1,52.6242], Tmin=(100,'K'), Tmax=(1106.58,'K')), NASAPolynomial(coeffs=[21.8491,0.0682563,-3.00426e-05,5.68945e-09,-3.98085e-13,-74423.6,-70.2862], Tmin=(1106.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-574.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH)"""),
)

species(
    label = 'CCCCC(CCOO)OCC=O(29364)',
    structure = SMILES('CCCCC(CCOO)OCC=O'),
    E0 = (-557.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52345,0.152655,-0.000143585,7.65864e-08,-1.74771e-11,-66852.9,50.9802], Tmin=(100,'K'), Tmax=(1022.4,'K')), NASAPolynomial(coeffs=[16.0171,0.0801176,-3.7163e-05,7.19289e-09,-5.0879e-13,-70644.1,-38.878], Tmin=(1022.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-557.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH)"""),
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
    label = 'CCC[C](CCOO)OC[CH]O(10466)',
    structure = SMILES('CCC[C](CCOO)OC[CH]O'),
    E0 = (-266.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,360,370,350,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,180,180,180,306.831,1448.58,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150343,'amu*angstrom^2'), symmetry=1, barrier=(3.45668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150343,'amu*angstrom^2'), symmetry=1, barrier=(3.45668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150343,'amu*angstrom^2'), symmetry=1, barrier=(3.45668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150343,'amu*angstrom^2'), symmetry=1, barrier=(3.45668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150343,'amu*angstrom^2'), symmetry=1, barrier=(3.45668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150343,'amu*angstrom^2'), symmetry=1, barrier=(3.45668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150343,'amu*angstrom^2'), symmetry=1, barrier=(3.45668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150343,'amu*angstrom^2'), symmetry=1, barrier=(3.45668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150343,'amu*angstrom^2'), symmetry=1, barrier=(3.45668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150343,'amu*angstrom^2'), symmetry=1, barrier=(3.45668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150343,'amu*angstrom^2'), symmetry=1, barrier=(3.45668,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (176.21,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.49897,0.182148,-0.000280174,2.45088e-07,-8.48617e-11,-31762.2,52.1977], Tmin=(100,'K'), Tmax=(834.407,'K')), NASAPolynomial(coeffs=[12.1523,0.0768503,-3.64712e-05,6.9027e-09,-4.72898e-13,-33320.4,-14.1634], Tmin=(834.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(636.057,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(C2CsJOCs)"""),
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
    label = 'CCCCC1(CCO1)OC[CH]O(29365)',
    structure = SMILES('CCCCC1(CCO1)OC[CH]O'),
    E0 = (-412.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (173.229,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.50148,0.153397,-0.000141318,6.68337e-08,-1.21415e-11,-49287.4,51.8848], Tmin=(100,'K'), Tmax=(1508.38,'K')), NASAPolynomial(coeffs=[33.4018,0.0365505,-8.87975e-06,1.12099e-09,-6.05171e-14,-58863.9,-140.398], Tmin=(1508.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-412.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(673.472,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C(CCC)(CCOO)OC[CH]O(29366)',
    structure = SMILES('[CH2]C(CCC)(CCOO)OC[CH]O'),
    E0 = (-269.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3615,1310,387.5,850,1000,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.645,0.203577,-0.00029277,2.36206e-07,-7.67177e-11,-32098.5,56.1353], Tmin=(100,'K'), Tmax=(807.888,'K')), NASAPolynomial(coeffs=[18.98,0.0756327,-3.48429e-05,6.55385e-09,-4.49533e-13,-35557.7,-50.5853], Tmin=(807.888,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-269.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-SQ) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C(CCCC)(COO)OC[CH]O(29367)',
    structure = SMILES('[CH2]C(CCCC)(COO)OC[CH]O'),
    E0 = (-272.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3615,1310,387.5,850,1000,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.645,0.203577,-0.00029277,2.36206e-07,-7.67177e-11,-32501.1,56.1353], Tmin=(100,'K'), Tmax=(807.888,'K')), NASAPolynomial(coeffs=[18.98,0.0756327,-3.48429e-05,6.55385e-09,-4.49533e-13,-35960.3,-50.5853], Tmin=(807.888,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(706.73,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-SQ) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)OC) + radical(CCsJOH)"""),
)

species(
    label = 'CCCCC1(CCOO)OCC1O(29161)',
    structure = SMILES('CCCCC1(CCOO)OCC1O'),
    E0 = (-530.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.51915,0.155474,-0.000144024,7.32337e-08,-1.50767e-11,-63515.7,50.478], Tmin=(100,'K'), Tmax=(1171.48,'K')), NASAPolynomial(coeffs=[24.411,0.0601054,-2.19103e-05,3.74014e-09,-2.46235e-13,-70059.6,-88.6891], Tmin=(1171.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-530.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(719.202,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane)"""),
)

species(
    label = 'CCCCC(O)(CC[O])OC[CH]O(29368)',
    structure = SMILES('CCCCC(O)(CC[O])OC[CH]O'),
    E0 = (-528.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.237,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.83518,0.175773,-0.000193979,1.15141e-07,-2.7633e-11,-63301.8,55.6485], Tmin=(100,'K'), Tmax=(1007.68,'K')), NASAPolynomial(coeffs=[24.2813,0.0641643,-2.78407e-05,5.22588e-09,-3.63698e-13,-68968.3,-80.2122], Tmin=(1007.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-528.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(710.887,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CCOJ)"""),
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
    label = 'CCCC[C]CCOO(28384)',
    structure = SMILES('CCCC[C]CCOO'),
    E0 = (144.248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.602271,0.106692,-9.65495e-05,4.98095e-08,-1.103e-11,17510.2,36.2426], Tmin=(100,'K'), Tmax=(1049.61,'K')), NASAPolynomial(coeffs=[12.3607,0.0572903,-2.59494e-05,4.96716e-09,-3.49172e-13,14789,-26.9238], Tmin=(1049.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
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
    label = '[CH2]O[C](CCCC)CCOO(28849)',
    structure = SMILES('[CH2]O[C](CCCC)CCOO'),
    E0 = (-89.3111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,360,370,350,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.65288,0.158133,-0.000217847,1.79592e-07,-6.09066e-11,-10513.4,46.7441], Tmin=(100,'K'), Tmax=(789.359,'K')), NASAPolynomial(coeffs=[12.0913,0.0723598,-3.38401e-05,6.43874e-09,-4.45662e-13,-12496.5,-18.718], Tmin=(789.359,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.3111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(C2CsJOCs) + radical(CsJOCC2)"""),
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
    E0 = (-289.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-129.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-184.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-178.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-195.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-242.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-213.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-191.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-160.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-116.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-115.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-172.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-129.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-131.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-120.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-166.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-243.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-242.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-206.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-254.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-172.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-167.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-164.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-167.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-138.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-170.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (74.0348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-226.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-281.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-281.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-281.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (153.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-214.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-112.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-115.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-281.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-179.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (133.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (116.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    products = ['CH2CHOH(42)', 'CCCCC(=O)CCOO(16766)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'CCCC[C](CCOO)OCC=O(29337)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2818 used for Od_CO-CsH;HJ
Exact match found for rate rule [Od_CO-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'CCCC=C(CCOO)OC[CH]O(29338)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(137.12,'m^3/(mol*s)'), n=1.63155, Ea=(4.2466,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds;HJ] for rate rule [Cds-CsH_Cds-OsCs;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'CCCCC(=CCOO)OC[CH]O(29339)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(137.12,'m^3/(mol*s)'), n=1.63155, Ea=(4.2466,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds;HJ] for rate rule [Cds-CsH_Cds-OsCs;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'CCCC[C](CCOO)OC=CO(29340)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.182e+10,'cm^3/(mol*s)'), n=0.859, Ea=(6.76971,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-OsH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]O(284)', 'CCCCC(=O)CCOO(16766)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(32300,'cm^3/(mol*s)'), n=2.98, Ea=(33.0536,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-CsCs;YJ] for rate rule [Od_CO-CsCs;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['npropyl(83)', 'C=C(CCOO)OC[CH]O(29341)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00238412,'m^3/(mol*s)'), n=2.47216, Ea=(17.7199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds;CsJ-CsHH] for rate rule [Cds-HH_Cds-OsCs;CsJ-CsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2OOH(35)', 'C=C(CCCC)OC[CH]O(29342)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(182.434,'m^3/(mol*s)'), n=0.88, Ea=(33.1163,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds;CsJ-OsHH] for rate rule [Cds-HH_Cds-OsCs;CsJ-OsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    products = ['CCCC[C](CCOO)OCC[O](29343)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4500,'s^-1'), n=2.62, Ea=(129.286,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 322 used for R2H_S;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CCC[CH]C(CCOO)OC[CH]O(29344)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CCCCC([CH]COO)OC[CH]O(29345)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CCCC[C](CCOO)O[CH]CO(29346)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.253e-19,'s^-1'), n=8.985, Ea=(117.152,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_1H;Cs_H_out_H/NonDeO] + [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_1H] for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CCCCC(CCOO)O[CH][CH]O(29347)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.97e+09,'s^-1'), n=1.01, Ea=(160.958,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CC[CH]CC(CCOO)OC[CH]O(29348)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.27e+09,'s^-1'), n=0.66, Ea=(144.766,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_NonDe] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CCCCC(C[CH]OO)OC[CH]O(29349)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.12723e+10,'s^-1'), n=0.640667, Ea=(161.697,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_1H;Cs_H_out_NonDe] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[CH]CCC(CCOO)OC[CH]O(29350)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.86e+10,'s^-1'), n=0.58, Ea=(109.579,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_single;Cs_H_out_NonDe] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_NDMustO]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CCC[CH][C](CCOO)OCCO(29351)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.0756983,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R5HJ_1;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CCCC[C]([CH]COO)OCCO(29352)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.0756983,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R5HJ_1;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CCCC(CCOO)OC[CH]O(29353)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.80589e+07,'s^-1'), n=1.02417, Ea=(59.2315,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_SSSS;Y_rad_out;Cs_H_out_NDMustO] + [R5H_SSSS;C_rad_out_2H;Cs_H_out_noH] + [R5H_CCC;Y_rad_out;Cs_H_out_NonDe] for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CCCCC(CCO[O])OC[CH]O(29354)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.9e+07,'s^-1'), n=1.1, Ea=(64.4336,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 276 used for R5H_SSSS_OCC;O_rad_out;Cs_H_out_NDMustO
Exact match found for rate rule [R5H_SSSS_OCC;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CCCCC(CCOO)OC[CH][O](29355)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.99186e+09,'s^-1'), n=0.575, Ea=(72.1112,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_Cs2] for rate rule [R5HJ_1;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC[CH]C[C](CCOO)OCCO(29356)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [RnH;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO] for rate rule [R6HJ_2;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CCCC[C](C[CH]OO)OCCO(29357)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.253e-19,'s^-1'), n=8.985, Ea=(117.152,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_1H;Cs_H_out_H/NonDeO] + [RnH;C_rad_out_H/NonDeO;Cs_H_out_1H] for rate rule [R6HJ_2;C_rad_out_H/NonDeO;Cs_H_out_H/NonDeO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C[CH]CC[C](CCOO)OCCO(29358)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [RnH;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO] for rate rule [R7HJ_3;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]CCC[C](CCOO)OCCO(29359)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4538.23,'s^-1'), n=2.54667, Ea=(126.817,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_2H;Cs_H_out_H/NonDeO] + [R8Hall;C_rad_out_2H;Cs_H_out_1H] for rate rule [R8Hall;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    products = ['CCCC[C](CCO[O])OCCO(29360)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.0900658,'s^-1'), n=3.79083, Ea=(119.293,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_1H;O_H_out] + [RnH;C_rad_out_H/NonDeO;XH_out] for rate rule [R8Hall;C_rad_out_H/NonDeO;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]O(284)', 'CCCC[C]([O])CCOO(16758)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    products = ['CCCCC(CCOO)OC=CO(29361)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    products = ['CCCC=C(CCOO)OCCO(29362)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    products = ['CCCCC(=CCOO)OCCO(29363)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    products = ['CCCCC(CCOO)OCC=O(29364)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(S)(23)', 'CCC[C](CCOO)OC[CH]O(10466)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    products = ['OH(5)', 'CCCCC1(CCO1)OC[CH]O(29365)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.31e+11,'s^-1','*|/',1.74), n=0, Ea=(75.8559,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3OO_SS;C_rad/NonDeC_intra;OOH] for rate rule [R3OO_SS;C_rad/NDMustO_intra;OOH]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(CCC)(CCOO)OC[CH]O(29366)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(CCCC)(COO)OC[CH]O(29367)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    products = ['CCCCC1(CCOO)OCC1O(29161)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_Cs2;Cpri_rad_out_H/NonDeO]
Euclidian distance = 3.60555127546
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    products = ['CCCCC(O)(CC[O])OC[CH]O(29368)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.79e+10,'s^-1'), n=0, Ea=(110.039,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3OOH_SS;C_rad_out_NonDe] for rate rule [R3OOH_SS;C_rad_out_NDMustO]
Euclidian distance = 1.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O]C[CH]O(1827)', 'CCCC[C]CCOO(28384)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['HCOH(T)(285)', '[CH2]O[C](CCCC)CCOO(28849)'],
    products = ['CCCC[C](CCOO)OC[CH]O(29159)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4669',
    isomers = [
        'CCCC[C](CCOO)OC[CH]O(29159)',
    ],
    reactants = [
        ('CH2CHOH(42)', 'CCCCC(=O)CCOO(16766)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4669',
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

