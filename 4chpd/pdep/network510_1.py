species(
    label = 'C=[C]OC1[CH]CCCCC1OO(2257)',
    structure = SMILES('C=[C]OC1[CH]CCCCC1OO'),
    E0 = (101.795,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,3615,1310,387.5,850,1000,2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.92042,0.0703018,0.000146123,-2.87632e-07,1.2955e-10,12512.2,47.4283], Tmin=(100,'K'), Tmax=(913.51,'K')), NASAPolynomial(coeffs=[51.2123,-0.00556244,1.32431e-05,-2.7746e-09,1.74867e-13,-3737.26,-239.907], Tmin=(913.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.795,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(C=CJO) + radical(CCJCO)"""),
)

species(
    label = 'CH2CO(43)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'OOC1C=CCCCC1(259)',
    structure = SMILES('OOC1C=CCCCC1'),
    E0 = (-124.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,3615,1310,387.5,850,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4468.79,'J/mol'), sigma=(7.53736,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=698.02 K, Pc=23.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.372801,0.0206158,0.000217263,-3.35136e-07,1.4146e-10,-14763.2,32.5035], Tmin=(100,'K'), Tmax=(913.852,'K')), NASAPolynomial(coeffs=[42.4349,-0.00964692,1.44126e-05,-2.93401e-09,1.83612e-13,-28874.9,-201.781], Tmin=(913.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene)"""),
)

species(
    label = 'H(19)',
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
    label = 'C=[C]OC1=CCCCCC1OO(3590)',
    structure = SMILES('C=[C]OC1=CCCCCC1OO'),
    E0 = (37.3499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,1685,370,3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32321,0.0623053,0.000141511,-2.66923e-07,1.18195e-10,4734.78,47.597], Tmin=(100,'K'), Tmax=(918.429,'K')), NASAPolynomial(coeffs=[45.3992,0.00116734,8.87219e-06,-1.88434e-09,1.12994e-13,-9851.18,-206.521], Tmin=(918.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.3499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC1C=CCCCC1OO(3591)',
    structure = SMILES('C=[C]OC1C=CCCCC1OO'),
    E0 = (14.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,3615,1310,387.5,850,1000,2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.64787,0.0610311,0.00016911,-3.12714e-07,1.38857e-10,1951.55,45.4161], Tmin=(100,'K'), Tmax=(914.207,'K')), NASAPolynomial(coeffs=[52.3901,-0.0103513,1.54172e-05,-3.1504e-09,1.98056e-13,-14826.2,-248.162], Tmin=(914.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene) + radical(C=CJO)"""),
)

species(
    label = 'C#COC1[CH]CCCCC1OO(3592)',
    structure = SMILES('C#COC1[CH]CCCCC1OO'),
    E0 = (71.3774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3615,1310,387.5,850,1000,2175,525,2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.13208,0.074204,0.000136832,-2.82637e-07,1.28999e-10,8862.17,42.5161], Tmin=(100,'K'), Tmax=(912.898,'K')), NASAPolynomial(coeffs=[53.6206,-0.0110249,1.55194e-05,-3.18549e-09,2.02402e-13,-7945.01,-257.679], Tmin=(912.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.3774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Ct-CtOs) + group(Ct-CtH) + ring(Cycloheptane) + radical(CCJCO)"""),
)

species(
    label = 'OOC1[CH][CH]CCCC1(365)',
    structure = SMILES('OOC1[CH][CH]CCCC1'),
    E0 = (133.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.113571,0.0375133,0.0001521,-2.5566e-07,1.10283e-10,16228.8,36.6232], Tmin=(100,'K'), Tmax=(914.695,'K')), NASAPolynomial(coeffs=[37.2426,4.92732e-05,8.71038e-06,-1.86634e-09,1.15238e-13,4211.41,-167.754], Tmin=(914.695,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(cycloheptane)"""),
)

species(
    label = 'C=[C]O[C]1CCCCCC1OO(3593)',
    structure = SMILES('C=[C]O[C]1CCCCCC1OO'),
    E0 = (102.109,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,3615,1310,387.5,850,1000,2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.19916,0.0838313,9.72559e-05,-2.29866e-07,1.07108e-10,12553,47.4296], Tmin=(100,'K'), Tmax=(915.508,'K')), NASAPolynomial(coeffs=[48.4572,-7.30184e-05,9.57201e-06,-2.05851e-09,1.27866e-13,-2481.27,-223.939], Tmin=(915.508,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(C2CsJOC(O)) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC1C[CH]CCCC1OO(3594)',
    structure = SMILES('C=[C]OC1C[CH]CCCC1OO'),
    E0 = (71.3218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93518,0.075938,0.000119723,-2.52227e-07,1.14811e-10,8842.62,47.3315], Tmin=(100,'K'), Tmax=(915.88,'K')), NASAPolynomial(coeffs=[47.9582,0.000619226,9.55729e-06,-2.05803e-09,1.26897e-13,-6276.91,-221.639], Tmin=(915.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.3218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(C=CJO) + radical(cycloheptane)"""),
)

species(
    label = '[CH]=COC1[CH]CCCCC1OO(3595)',
    structure = SMILES('[CH]=COC1[CH]CCCCC1OO'),
    E0 = (109.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,3120,650,792.5,1650,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.20946,0.0721195,0.000156231,-3.09017e-07,1.39926e-10,13411.2,45.7409], Tmin=(100,'K'), Tmax=(911.293,'K')), NASAPolynomial(coeffs=[55.6412,-0.0126031,1.71712e-05,-3.5362e-09,2.26379e-13,-4158.44,-266.53], Tmin=(911.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(Cds_P) + radical(CCJCO)"""),
)

species(
    label = 'C=CO[C]1[CH]CCCCC1OO(3596)',
    structure = SMILES('[CH2][CH]OC1=CCCCCC1OO'),
    E0 = (47.8014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52156,0.083136,0.000121061,-2.70379e-07,1.25465e-10,6040.33,47.2158], Tmin=(100,'K'), Tmax=(912.688,'K')), NASAPolynomial(coeffs=[55.1033,-0.0110445,1.55661e-05,-3.20103e-09,2.03931e-13,-11074.5,-261.662], Tmin=(912.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.8014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(611.114,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CJCO) + radical(CCsJOC(O))"""),
)

species(
    label = 'C=[C]OC1CCCCC[C]1OO(3597)',
    structure = SMILES('C=[C]OC1CCCCC[C]1OO'),
    E0 = (88.7921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.9684,0.0779255,0.000112812,-2.4381e-07,1.11519e-10,10943.8,46.0907], Tmin=(100,'K'), Tmax=(915.757,'K')), NASAPolynomial(coeffs=[47.3072,0.00205101,8.8237e-06,-1.92691e-09,1.18695e-13,-3924.5,-219.204], Tmin=(915.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.7921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(C=CJO) + radical(C2CsJOOH)"""),
)

species(
    label = 'C=[C]OC1CC[CH]CCC1OO(3598)',
    structure = SMILES('C=[C]OC1CC[CH]CCC1OO'),
    E0 = (71.3218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93518,0.075938,0.000119723,-2.52227e-07,1.14811e-10,8842.62,46.6384], Tmin=(100,'K'), Tmax=(915.88,'K')), NASAPolynomial(coeffs=[47.9582,0.000619226,9.55729e-06,-2.05803e-09,1.26897e-13,-6276.91,-222.332], Tmin=(915.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.3218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(C=CJO) + radical(cycloheptane)"""),
)

species(
    label = 'C=[C]OC1CCCC[CH]C1OO(3599)',
    structure = SMILES('C=[C]OC1CCCC[CH]C1OO'),
    E0 = (102.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07994,0.0774285,0.000121966,-2.59417e-07,1.18652e-10,12576.2,47.8092], Tmin=(100,'K'), Tmax=(914.22,'K')), NASAPolynomial(coeffs=[49.8785,-0.00252533,1.13361e-05,-2.40867e-09,1.51006e-13,-3083.05,-231.885], Tmin=(914.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(C=CJO) + radical(CCJCOOH)"""),
)

species(
    label = 'C=COC1[CH]CCCC[C]1OO(3600)',
    structure = SMILES('C=COC1[CH]CCCC[C]1OO'),
    E0 = (48.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99689,0.0710991,0.000147878,-2.91649e-07,1.31483e-10,6160.01,45.5614], Tmin=(100,'K'), Tmax=(913.24,'K')), NASAPolynomial(coeffs=[52.048,-0.00643261,1.37621e-05,-2.87561e-09,1.81594e-13,-10349.3,-246.613], Tmin=(913.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(C2CsJOOH) + radical(CCJCO)"""),
)

species(
    label = 'C=[C]OC1CCC[CH]CC1OO(3601)',
    structure = SMILES('C=[C]OC1CCC[CH]CC1OO'),
    E0 = (71.3218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93518,0.075938,0.000119723,-2.52227e-07,1.14811e-10,8842.62,47.3315], Tmin=(100,'K'), Tmax=(915.88,'K')), NASAPolynomial(coeffs=[47.9582,0.000619226,9.55729e-06,-2.05803e-09,1.26897e-13,-6276.91,-221.639], Tmin=(915.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.3218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(cycloheptane) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]OC1CCCCCC1OO(3602)',
    structure = SMILES('[CH]=[C]OC1CCCCCC1OO'),
    E0 = (148.989,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,3615,1310,387.5,850,1000,2750,2786.36,2822.73,2859.09,2895.45,2931.82,2968.18,3004.55,3040.91,3077.27,3113.64,3150,900,918.182,936.364,954.545,972.727,990.909,1009.09,1027.27,1045.45,1063.64,1081.82,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.1798,0.0789319,0.000121214,-2.61241e-07,1.19988e-10,18195,45.5729], Tmin=(100,'K'), Tmax=(913.33,'K')), NASAPolynomial(coeffs=[50.8946,-0.00410954,1.2227e-05,-2.58613e-09,1.63366e-13,2268.73,-239.782], Tmin=(913.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]OC1CCCCCC1O[O](3603)',
    structure = SMILES('C=[C]OC1CCCCCC1O[O]'),
    E0 = (53.8977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.82222,0.0715809,0.000134003,-2.70423e-07,1.22634e-10,6744.81,45.3386], Tmin=(100,'K'), Tmax=(911.275,'K')), NASAPolynomial(coeffs=[48.8536,-0.00257001,1.19698e-05,-2.57703e-09,1.64362e-13,-8648.19,-228.217], Tmin=(911.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.8977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(C=CJO) + radical(ROOJ)"""),
)

species(
    label = 'C=COC1[CH]CCC[CH]C1OO(3604)',
    structure = SMILES('C=COC1[CH]CCC[CH]C1OO'),
    E0 = (62.4658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10972,0.0706179,0.000156975,-3.07179e-07,1.38583e-10,7792.47,47.2845], Tmin=(100,'K'), Tmax=(912.054,'K')), NASAPolynomial(coeffs=[54.6245,-0.011018,1.62798e-05,-3.35862e-09,2.1401e-13,-9509.98,-259.323], Tmin=(912.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.4658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(CCJCO)"""),
)

species(
    label = 'C=COC1[CH][CH]CCCC1OO(3605)',
    structure = SMILES('C=COC1[CH][CH]CCCC1OO'),
    E0 = (31.4797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96444,0.0691215,0.000154752,-3.00012e-07,1.34751e-10,4058.83,46.8049], Tmin=(100,'K'), Tmax=(913.425,'K')), NASAPolynomial(coeffs=[52.701,-0.00786796,1.44978e-05,-3.00723e-09,1.89838e-13,-12702.5,-249.059], Tmin=(913.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.4797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(CCJCO) + radical(cycloheptane)"""),
)

species(
    label = 'C=COC1[CH]CCCCC1O[O](3606)',
    structure = SMILES('C=COC1[CH]CCCCC1O[O]'),
    E0 = (14.0556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.85179,0.0647671,0.000169028,-3.18215e-07,1.42583e-10,1961.03,44.8132], Tmin=(100,'K'), Tmax=(909.469,'K')), NASAPolynomial(coeffs=[53.6018,-0.0110664,1.69156e-05,-3.52747e-09,2.27406e-13,-15076.1,-255.668], Tmin=(909.469,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.0556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(CCJCO) + radical(ROOJ)"""),
)

species(
    label = 'C=COC1[CH]CC[CH]CC1OO(3607)',
    structure = SMILES('C=COC1[CH]CC[CH]CC1OO'),
    E0 = (31.4797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96444,0.0691215,0.000154752,-3.00012e-07,1.34751e-10,4058.83,46.8049], Tmin=(100,'K'), Tmax=(913.425,'K')), NASAPolynomial(coeffs=[52.701,-0.00786796,1.44978e-05,-3.00723e-09,1.89838e-13,-12702.5,-249.059], Tmin=(913.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.4797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(CCJCO) + radical(cycloheptane)"""),
)

species(
    label = 'C=COC1[CH]C[CH]CCC1OO(3608)',
    structure = SMILES('C=COC1[CH]C[CH]CCC1OO'),
    E0 = (31.4797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96444,0.0691215,0.000154752,-3.00012e-07,1.34751e-10,4058.83,46.8049], Tmin=(100,'K'), Tmax=(913.425,'K')), NASAPolynomial(coeffs=[52.701,-0.00786796,1.44978e-05,-3.00723e-09,1.89838e-13,-12702.5,-249.059], Tmin=(913.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.4797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(CCJCO) + radical(cycloheptane)"""),
)

species(
    label = '[CH2][C]=O(360)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97567e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02477e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2]CCC(OO)C(C=C)O[C]=C(3609)',
    structure = SMILES('[CH2]CCC(OO)C(C=C)O[C]=C'),
    E0 = (154.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1310,387.5,850,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.83886,0.139508,-0.000130768,6.21831e-08,-1.1748e-11,18827.7,53.5311], Tmin=(100,'K'), Tmax=(1277.59,'K')), NASAPolynomial(coeffs=[27.847,0.0434328,-1.79668e-05,3.32116e-09,-2.2977e-13,10987,-102.028], Tmin=(1277.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC=CCCCC[CH]OO(3610)',
    structure = SMILES('C=[C]OC=CCCCC[CH]OO'),
    E0 = (179.806,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,1685,370,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96049,0.137096,-0.000140786,8.06693e-08,-1.92543e-11,21835.4,50.6503], Tmin=(100,'K'), Tmax=(996.954,'K')), NASAPolynomial(coeffs=[16.6847,0.0622867,-2.82289e-05,5.40122e-09,-3.79639e-13,18117.7,-39.2452], Tmin=(996.954,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CCsJOOH)"""),
)

species(
    label = 'C=COC1=CCCCCC1OO(3611)',
    structure = SMILES('C=COC1=CCCCCC1OO'),
    E0 = (-202.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53937,0.0606166,0.000168187,-3.05517e-07,1.34381e-10,-24085.9,45.2445], Tmin=(100,'K'), Tmax=(917.124,'K')), NASAPolynomial(coeffs=[49.5647,-0.00300166,1.17441e-05,-2.44198e-09,1.49211e-13,-40157.9,-233.4], Tmin=(917.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-202.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene)"""),
)

species(
    label = 'C=COC1C=CCCCC1OO(3612)',
    structure = SMILES('C=COC1C=CCCCC1OO'),
    E0 = (-225.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.86475,0.0593508,0.000195758,-3.51273e-07,1.5503e-10,-26869.2,43.0661], Tmin=(100,'K'), Tmax=(913.485,'K')), NASAPolynomial(coeffs=[56.5597,-0.0145273,1.82931e-05,-3.70901e-09,2.34354e-13,-45134.7,-275.064], Tmin=(913.485,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-225.699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptene)"""),
)

species(
    label = 'OH(21)',
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
    label = 'C=[C]OC1C2CCCCC1O2(3613)',
    structure = SMILES('C=[C]OC1C2CCCCC1O2'),
    E0 = (-141.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.78288,0.130286,-0.000115672,5.11569e-08,-8.89642e-12,-16807.5,9.0553], Tmin=(100,'K'), Tmax=(1393.87,'K')), NASAPolynomial(coeffs=[30.1562,0.0357605,-1.39489e-05,2.50422e-09,-1.70214e-13,-25990.1,-160.795], Tmin=(1393.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-141.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + polycyclic(PolycyclicRing) + radical(C=CJO)"""),
)

species(
    label = 'C=C1OC2[CH]CCCCC2O1(3614)',
    structure = SMILES('C=C1OC2[CH]CCCCC2O1'),
    E0 = (-122.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.979389,0.0836268,8.91399e-06,-7.36664e-08,3.40779e-11,-14540.7,33.0149], Tmin=(100,'K'), Tmax=(999.005,'K')), NASAPolynomial(coeffs=[25.0839,0.0395551,-1.54306e-05,2.98504e-09,-2.20521e-13,-22756.5,-107.756], Tmin=(999.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + polycyclic(s2_5_7) + radical(CCJCO)"""),
)

species(
    label = 'C=[C]O[CH]C1CCCCC1OO(3615)',
    structure = SMILES('C=[C]O[CH]C1CCCCC1OO'),
    E0 = (69.2532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.47851,0.120155,-6.391e-05,-5.56382e-09,1.05812e-11,8581.6,42.7523], Tmin=(100,'K'), Tmax=(1067.14,'K')), NASAPolynomial(coeffs=[28.3821,0.0465884,-1.96913e-05,3.7884e-09,-2.72296e-13,-402.573,-119.372], Tmin=(1067.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.2532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(611.114,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(CCsJOC(O)) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1CCCC(OO)C1O[C]=C(3616)',
    structure = SMILES('[CH2]C1CCCC(OO)C1O[C]=C'),
    E0 = (67.7652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.06251,0.109939,-3.77817e-05,-3.09126e-08,1.94067e-11,8388.88,42.5251], Tmin=(100,'K'), Tmax=(1030.18,'K')), NASAPolynomial(coeffs=[26.5589,0.048099,-1.95109e-05,3.70965e-09,-2.66608e-13,-123.746,-109.102], Tmin=(1030.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.7652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(611.114,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(C=CJO) + radical(Isobutyl)"""),
)

species(
    label = 'C=C1OC2C(CCCCC12)OO(3617)',
    structure = SMILES('C=C1OC2C(CCCCC12)OO'),
    E0 = (-353.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.31996,0.141307,-0.000116538,4.13255e-08,-3.80741e-12,-42189.2,9.75809], Tmin=(100,'K'), Tmax=(1098.16,'K')), NASAPolynomial(coeffs=[31.6482,0.0410261,-1.65653e-05,3.09771e-09,-2.18621e-13,-51502.8,-169.655], Tmin=(1098.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-353.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + polycyclic(PolycyclicRing)"""),
)

species(
    label = 'C=[C]OC1C([O])CCCCC1O(3618)',
    structure = SMILES('C=[C]OC1C([O])CCCCC1O'),
    E0 = (-111.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96052,0.0716379,0.000141948,-2.85247e-07,1.29621e-10,-13166.8,45.7876], Tmin=(100,'K'), Tmax=(909.939,'K')), NASAPolynomial(coeffs=[51.644,-0.00758348,1.4691e-05,-3.09906e-09,1.99498e-13,-29397.8,-243.347], Tmin=(909.939,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(CC(C)OJ) + radical(C=CJO)"""),
)

species(
    label = 'H2CC(T)(358)',
    structure = SMILES('[C]=C'),
    E0 = (600.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94832,-0.00181667,2.0713e-05,-2.36697e-08,8.4549e-12,72206.8,5.31546], Tmin=(100,'K'), Tmax=(953.658,'K')), NASAPolynomial(coeffs=[4.2027,0.00473426,-1.57307e-06,2.8591e-10,-2.08647e-14,71811.9,2.28404], Tmin=(953.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]C1[CH]CCCCC1OO(1585)',
    structure = SMILES('[O]C1[CH]CCCCC1OO'),
    E0 = (21.3476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,3615,1310,387.5,850,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.47042,0.0435283,0.00016546,-2.85954e-07,1.24818e-10,2779.66,38.0372], Tmin=(100,'K'), Tmax=(912.732,'K')), NASAPolynomial(coeffs=[43.6489,-0.00670377,1.28079e-05,-2.66108e-09,1.6804e-13,-11235.6,-203.442], Tmin=(912.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.3476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cycloheptane) + radical(CC(C)OJ) + radical(CCJCO)"""),
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
    E0 = (101.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (260.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (231.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (298.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (102.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (272.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (260.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (214.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (243.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (250.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (250.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (148.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (209.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (148.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (294.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (128.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (178.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (301.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (162.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (171.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (301.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (293.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (181.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (226.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (201.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (141.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (177.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (156.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (261.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (347.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (110.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (213.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (621.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['CH2CO(43)', 'OOC1C=CCCCC1(259)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(19)', 'C=[C]OC1=CCCCCC1OO(3590)'],
    products = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', 'C=[C]OC1C=CCCCC1OO(3591)'],
    products = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'C#COC1[CH]CCCCC1OO(3592)'],
    products = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OOC1[CH][CH]CCCC1(365)', 'CH2CO(43)'],
    products = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(11.6997,'m^3/(mol*s)'), n=2.021, Ea=(29.883,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd;CJ]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C]O[C]1CCCCCC1OO(3593)'],
    products = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.03713e+09,'s^-1'), n=1.3215, Ea=(170.736,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_H/NonDeC] for rate rule [R2H_S;C_rad_out_NDMustO;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=[C]OC1C[CH]CCCC1OO(3594)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=COC1[CH]CCCCC1OO(3595)'],
    products = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=CO[C]1[CH]CCCCC1OO(3596)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=[C]OC1CCCCC[C]1OO(3597)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.05815e+09,'s^-1'), n=0.95, Ea=(148.741,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=[C]OC1CC[CH]CCC1OO(3598)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.84148e+07,'s^-1'), n=1.625, Ea=(148.323,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)] + [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=[C]OC1CCCC[CH]C1OO(3599)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.00228,'s^-1'), n=3.95, Ea=(46.7353,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=COC1[CH]CCCC[C]1OO(3600)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.19773e+08,'s^-1'), n=1.45, Ea=(108.157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS_OCs;Y_rad_out;Cs_H_out_OOH/Cs] for rate rule [R4H_SSS_OCs;Cd_rad_out_Cd;Cs_H_out_OOH/Cs]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=[C]OC1CCC[CH]CC1OO(3601)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00228,'s^-1'), n=3.95, Ea=(46.7353,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]OC1CCCCCC1OO(3602)'],
    products = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=[C]OC1CCCCCC1O[O](3603)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_H/NonDeC;XH_out] for rate rule [R5H_SSSS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=COC1[CH]CCC[CH]C1OO(3604)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.0994e+09,'s^-1'), n=0.735, Ea=(76.7137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS_OCC;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_SSSS_OCC;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=COC1[CH][CH]CCCC1OO(3605)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=COC1[CH]CCCCC1O[O](3606)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(564.492,'s^-1'), n=2.19647, Ea=(60.7245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;Y_rad_out;XH_out] for rate rule [R6H_SSSSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=COC1[CH]CC[CH]CC1OO(3607)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.24e+06,'s^-1'), n=1.22, Ea=(69.4544,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6H_SSSSS_OOCCC(Cs/Cs);Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6H_SSSSS_OOCCC(Cs/Cs);Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=COC1[CH]C[CH]CCC1OO(3608)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R6HJ_3;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['OOC1[CH][CH]CCCC1(365)', '[CH2][C]=O(360)'],
    products = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.35017e+06,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]CCC(OO)C(C=C)O[C]=C(3609)'],
    products = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.2, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 46 used for R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]OC=CCCCC[CH]OO(3610)'],
    products = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.1182e+07,'s^-1'), n=0.786667, Ea=(47.07,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_linear;doublebond_intra_pri;radadd_intra_csHNd] + [R7_linear;doublebond_intra_pri_HNd;radadd_intra_cs] for rate rule [R7_linear;doublebond_intra_pri_HNd_O;radadd_intra_csHO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=COC1=CCCCCC1OO(3611)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.47101e+09,'s^-1'), n=0.37, Ea=(99.6901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad_NDe] + [R3radExo;Y_rad;XH_Rrad_NDe] for rate rule [R3radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=COC1C=CCCCC1OO(3612)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['OH(21)', 'C=[C]OC1C2CCCCC1O2(3613)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1','*|/',1.74), n=0, Ea=(75.8559,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 9 used for R3OO_SS;C_rad/H/NonDeC_intra;OOH
Exact match found for rate rule [R3OO_SS;C_rad/H/NonDeC_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['OH(21)', 'C=C1OC2[CH]CCCCC2O1(3614)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OOH] for rate rule [R4OO;Cd_rad_out;OOH]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=[C]O[CH]C1CCCCC1OO(3615)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['[CH2]C1CCCC(OO)C1O[C]=C(3616)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=C1OC2C(CCCCC12)OO(3617)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    products = ['C=[C]OC1C([O])CCCCC1O(3618)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.88e+10,'s^-1'), n=0, Ea=(111.294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for R3OOH_SS;C_rad_out_H/NonDeC
Exact match found for rate rule [R3OOH_SS;C_rad_out_H/NonDeC]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H2CC(T)(358)', '[O]C1[CH]CCCCC1OO(1585)'],
    products = ['C=[C]OC1[CH]CCCCC1OO(2257)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

network(
    label = '510',
    isomers = [
        'C=[C]OC1[CH]CCCCC1OO(2257)',
    ],
    reactants = [
        ('CH2CO(43)', 'OOC1C=CCCCC1(259)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '510',
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

