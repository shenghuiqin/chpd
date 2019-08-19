species(
    label = '[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)',
    structure = SMILES('[CH2]C(C=C)C(O)([CH]C)C(=C)O'),
    E0 = (-23.0335,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,344.984,784.383,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.29486,0.132516,-0.000125856,6.01436e-08,-1.10024e-11,-2483.82,53.6144], Tmin=(100,'K'), Tmax=(1483.82,'K')), NASAPolynomial(coeffs=[31.0249,0.0273455,-6.7481e-06,8.82588e-10,-4.95903e-14,-11275.8,-120.808], Tmin=(1483.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.0335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCJCO)"""),
)

species(
    label = 'butadiene13(2459)',
    structure = SMILES('C=CC=C'),
    E0 = (96.4553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.30712,'amu*angstrom^2'), symmetry=1, barrier=(30.0532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80599,0.0102584,6.1726e-05,-9.01643e-08,3.59117e-11,11658.5,12.0621], Tmin=(100,'K'), Tmax=(946.047,'K')), NASAPolynomial(coeffs=[12.4694,0.0100554,-2.41207e-06,4.57077e-10,-3.93161e-14,8010.78,-43.6375], Tmin=(946.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""butadiene13""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C(O)C(O)=CC(5562)',
    structure = SMILES('C=C(O)C(O)=CC'),
    E0 = (-369.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.87493,'amu*angstrom^2'), symmetry=1, barrier=(20.1164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.877997,'amu*angstrom^2'), symmetry=1, barrier=(20.1869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.874595,'amu*angstrom^2'), symmetry=1, barrier=(20.1087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.875946,'amu*angstrom^2'), symmetry=1, barrier=(20.1397,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4302.09,'J/mol'), sigma=(6.83849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.98 K, Pc=30.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.884099,0.0877395,-9.34163e-05,4.76671e-08,-9.0713e-12,-44294.8,25.8217], Tmin=(100,'K'), Tmax=(1473.72,'K')), NASAPolynomial(coeffs=[23.5689,0.00887539,-4.29798e-07,-1.49469e-10,1.60597e-14,-50145.5,-97.0299], Tmin=(1473.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-369.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][CH]C=C(2458)',
    structure = SMILES('[CH2]C=C[CH2]'),
    E0 = (274.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.210311,'amu*angstrom^2'), symmetry=1, barrier=(25.2351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779031,'amu*angstrom^2'), symmetry=1, barrier=(93.4717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56318,0.0223429,1.87067e-05,-3.93099e-08,1.63982e-11,33100.5,13.4097], Tmin=(100,'K'), Tmax=(974.264,'K')), NASAPolynomial(coeffs=[9.82995,0.0151966,-5.22272e-06,9.67656e-10,-7.0786e-14,30607.7,-26.985], Tmin=(974.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)C1(O)C(C)C1([CH2])O(32620)',
    structure = SMILES('[CH2]C(C=C)C1(O)C(C)C1([CH2])O'),
    E0 = (46.3504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.13544,0.135742,-0.000133209,6.64952e-08,-1.28554e-11,5849.57,46.4308], Tmin=(100,'K'), Tmax=(1342.92,'K')), NASAPolynomial(coeffs=[30.0822,0.030594,-8.83083e-06,1.30868e-09,-7.95821e-14,-2512.54,-121.536], Tmin=(1342.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(46.3504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CJC(C)2O) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1CC1C(O)([CH]C)C(=C)O(32621)',
    structure = SMILES('[CH2]C1CC1C(O)([CH]C)C(=C)O'),
    E0 = (1.66957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.83897,0.10385,-2.80413e-05,-5.68662e-08,3.53074e-11,433.717,46.4453], Tmin=(100,'K'), Tmax=(927.39,'K')), NASAPolynomial(coeffs=[30.5496,0.027291,-6.33628e-06,9.44057e-10,-6.69941e-14,-8288.81,-122.008], Tmin=(927.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.66957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCJCO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1(O)CC(C=C)C1(O)[CH]C(32622)',
    structure = SMILES('[CH2]C1(O)CC(C=C)C1(O)[CH]C'),
    E0 = (45.7922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07533,0.115301,-7.19183e-05,1.02963e-09,1.06781e-11,5742.59,43.6277], Tmin=(100,'K'), Tmax=(982.239,'K')), NASAPolynomial(coeffs=[26.9872,0.0368096,-1.29244e-05,2.30437e-09,-1.61879e-13,-1889.53,-105.849], Tmin=(982.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.7922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CJC(C)2O) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1C([CH2])C(O)(C(=C)O)C1C(32623)',
    structure = SMILES('[CH2]C1C([CH2])C(O)(C(=C)O)C1C'),
    E0 = (-6.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6639,0.0985677,-1.03744e-05,-7.65972e-08,4.30263e-11,-567.042,43.4514], Tmin=(100,'K'), Tmax=(918.007,'K')), NASAPolynomial(coeffs=[30.2963,0.027449,-5.50787e-06,7.24726e-10,-5.01846e-14,-9306.18,-123.642], Tmin=(918.007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = 'C=CC(=C)C(O)([CH]C)C(=C)O(32624)',
    structure = SMILES('C=CC(=C)C(O)([CH]C)C(=C)O'),
    E0 = (-128.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,2980,3010,3040,3070,3100,1330,1380,1430,900,975,1050,1000,1025,1050,1600,1650,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,180,180,180,436.356,690.876,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15555,'amu*angstrom^2'), symmetry=1, barrier=(3.5764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15555,'amu*angstrom^2'), symmetry=1, barrier=(3.5764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15555,'amu*angstrom^2'), symmetry=1, barrier=(3.5764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15555,'amu*angstrom^2'), symmetry=1, barrier=(3.5764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15555,'amu*angstrom^2'), symmetry=1, barrier=(3.5764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15555,'amu*angstrom^2'), symmetry=1, barrier=(3.5764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15555,'amu*angstrom^2'), symmetry=1, barrier=(3.5764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.20573,0.11678,-7.57528e-05,-3.48887e-10,1.26302e-11,-15223.8,44.3274], Tmin=(100,'K'), Tmax=(970.21,'K')), NASAPolynomial(coeffs=[29.7823,0.0295193,-9.82674e-06,1.75275e-09,-1.25701e-13,-23530.8,-119.851], Tmin=(970.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-128.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C(C=C)C(O)(C=C)C(=C)O(32625)',
    structure = SMILES('[CH2]C(C=C)C(O)(C=C)C(=C)O'),
    E0 = (-101.359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,2980,3010,3040,3070,3100,1330,1380,1430,900,975,1050,1000,1025,1050,1600,1650,1700,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,311.254,815.024,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152917,'amu*angstrom^2'), symmetry=1, barrier=(3.51587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152917,'amu*angstrom^2'), symmetry=1, barrier=(3.51587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152917,'amu*angstrom^2'), symmetry=1, barrier=(3.51587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152917,'amu*angstrom^2'), symmetry=1, barrier=(3.51587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152917,'amu*angstrom^2'), symmetry=1, barrier=(3.51587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152917,'amu*angstrom^2'), symmetry=1, barrier=(3.51587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152917,'amu*angstrom^2'), symmetry=1, barrier=(3.51587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00307,0.115433,-8.15887e-05,1.20932e-08,6.91635e-12,-11959.8,46.173], Tmin=(100,'K'), Tmax=(977.439,'K')), NASAPolynomial(coeffs=[26.7016,0.0339239,-1.16876e-05,2.05592e-09,-1.43344e-13,-19288.9,-100.443], Tmin=(977.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = 'CH2COH(99)',
    structure = SMILES('C=[C]O'),
    E0 = (103.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.989114,'amu*angstrom^2'), symmetry=1, barrier=(22.7417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1624,0.0134245,5.56346e-06,-1.95511e-08,9.36369e-12,12455.2,10.1544], Tmin=(100,'K'), Tmax=(925.618,'K')), NASAPolynomial(coeffs=[8.19875,0.00453462,-8.93448e-07,1.26083e-10,-9.46513e-15,10971.3,-16.733], Tmin=(925.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2COH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(C=C)C(O)=CC(32626)',
    structure = SMILES('[CH2]C(C=C)C(O)=CC'),
    E0 = (9.67819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,316.765,317.123],'cm^-1')),
        HinderedRotor(inertia=(0.17317,'amu*angstrom^2'), symmetry=1, barrier=(12.4004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173206,'amu*angstrom^2'), symmetry=1, barrier=(12.4074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173745,'amu*angstrom^2'), symmetry=1, barrier=(12.4053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173204,'amu*angstrom^2'), symmetry=1, barrier=(12.4057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172899,'amu*angstrom^2'), symmetry=1, barrier=(12.4139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.476456,0.0901748,-8.36736e-05,4.1519e-08,-8.23352e-12,1332.17,32.1226], Tmin=(100,'K'), Tmax=(1222.82,'K')), NASAPolynomial(coeffs=[17.6645,0.0308338,-1.08818e-05,1.83393e-09,-1.20102e-13,-3104.46,-59.046], Tmin=(1222.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.67819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(O)=C(O)[CH]C(4609)',
    structure = SMILES('[CH2]C(O)=C(O)[CH]C'),
    E0 = (-122.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,322.653],'cm^-1')),
        HinderedRotor(inertia=(0.00160507,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208419,'amu*angstrom^2'), symmetry=1, barrier=(15.514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20977,'amu*angstrom^2'), symmetry=1, barrier=(15.5164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209392,'amu*angstrom^2'), symmetry=1, barrier=(15.5175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207973,'amu*angstrom^2'), symmetry=1, barrier=(15.5093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.141743,0.0688338,-2.99459e-05,-3.072e-08,2.30897e-11,-14621,29.0351], Tmin=(100,'K'), Tmax=(905.863,'K')), NASAPolynomial(coeffs=[24.2058,0.00628117,1.26066e-06,-4.23747e-10,2.91719e-14,-20774,-94.5791], Tmin=(905.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(CCJCO)"""),
)

species(
    label = 'C2H3(30)',
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
    label = 'C=CC(O)([CH]C)C(=C)O(32291)',
    structure = SMILES('C=CC(O)([CH]C)C(=C)O'),
    E0 = (-178.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,180,180,180,180,1600,1625.66,2860.38,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150819,'amu*angstrom^2'), symmetry=1, barrier=(3.46762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150819,'amu*angstrom^2'), symmetry=1, barrier=(3.46762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150819,'amu*angstrom^2'), symmetry=1, barrier=(3.46762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150819,'amu*angstrom^2'), symmetry=1, barrier=(3.46762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150819,'amu*angstrom^2'), symmetry=1, barrier=(3.46762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150819,'amu*angstrom^2'), symmetry=1, barrier=(3.46762,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.161,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.74119,0.0881993,-4.90513e-05,-8.40888e-09,1.18976e-11,-21318.4,37.6411], Tmin=(100,'K'), Tmax=(981.75,'K')), NASAPolynomial(coeffs=[23.1971,0.0264226,-9.29506e-06,1.69222e-09,-1.21472e-13,-27741.8,-86.1822], Tmin=(981.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-178.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO)"""),
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
    label = '[CH2]C(C=C)C(=CC)C(=C)O(32627)',
    structure = SMILES('[CH2]C(C=C)C(=CC)C(=C)O'),
    E0 = (60.5507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.199,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.92893,0.122141,-0.000123886,6.62913e-08,-1.4056e-11,7503.1,37.9504], Tmin=(100,'K'), Tmax=(1150.88,'K')), NASAPolynomial(coeffs=[22.4679,0.0373474,-1.33715e-05,2.27412e-09,-1.49956e-13,1887.51,-83.1786], Tmin=(1150.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.5507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = 'C=C[C](C)C(O)([CH]C)C(=C)O(32628)',
    structure = SMILES('[CH2]C=C(C)C(O)([CH]C)C(=C)O'),
    E0 = (-102.371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.30415,0.12211,-9.43383e-05,2.6951e-08,5.45538e-13,-12071,46.0955], Tmin=(100,'K'), Tmax=(1032.03,'K')), NASAPolynomial(coeffs=[26.9401,0.0378991,-1.42898e-05,2.59715e-09,-1.81658e-13,-19658.8,-103.43], Tmin=(1032.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC(O)(C(=C)O)C([CH2])C=C(20155)',
    structure = SMILES('[CH2]CC(O)(C(=C)O)C([CH2])C=C'),
    E0 = (-17.6892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,309.411,820.877,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.11048,0.134576,-0.00013193,6.5739e-08,-1.26611e-11,-1852.96,52.2328], Tmin=(100,'K'), Tmax=(1360.54,'K')), NASAPolynomial(coeffs=[29.9995,0.029864,-8.3614e-06,1.20958e-09,-7.23992e-14,-10180.5,-115.191], Tmin=(1360.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.6892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C=C)C([O])(CC)C(=C)O(19804)',
    structure = SMILES('[CH2]C(C=C)C([O])(CC)C(=C)O'),
    E0 = (6.16714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,891.493,1315.34,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4819.15,'J/mol'), sigma=(8.10355,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=752.74 K, Pc=20.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10588,0.119884,-9.2049e-05,2.58179e-08,1.27997e-12,974.314,48.0474], Tmin=(100,'K'), Tmax=(991.066,'K')), NASAPolynomial(coeffs=[24.8373,0.0397338,-1.4018e-05,2.44045e-09,-1.66611e-13,-5770.46,-88.7809], Tmin=(991.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.16714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2][C](C=C)C(O)(CC)C(=C)O(20150)',
    structure = SMILES('[CH2]C=C([CH2])C(O)(CC)C(=C)O'),
    E0 = (-150.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.46663,0.124453,-9.31103e-05,2.1403e-08,3.42576e-12,-17885.5,44.4623], Tmin=(100,'K'), Tmax=(1011.41,'K')), NASAPolynomial(coeffs=[28.4254,0.0365467,-1.356e-05,2.46658e-09,-1.73752e-13,-25887.1,-113.589], Tmin=(1011.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-150.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C(C)C(O)([CH]C)C(=C)O(26113)',
    structure = SMILES('C=[C]C(C)C(O)([CH]C)C(=C)O'),
    E0 = (9.726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,180,180,180,523.054,610.775,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157345,'amu*angstrom^2'), symmetry=1, barrier=(3.61767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.31615,0.123569,-9.79371e-05,2.65283e-08,2.40075e-12,1410.89,48.7761], Tmin=(100,'K'), Tmax=(974.711,'K')), NASAPolynomial(coeffs=[27.2311,0.0349815,-1.18811e-05,2.05402e-09,-1.41163e-13,-5900.93,-100.976], Tmin=(974.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCJCO)"""),
)

species(
    label = 'C=CC(C)C([O])([CH]C)C(=C)O(20153)',
    structure = SMILES('C=CC(C)C([O])([CH]C)C(=C)O'),
    E0 = (0.986866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.9345,0.11304,-7.10988e-05,4.06188e-09,8.6232e-12,347.804,48.1475], Tmin=(100,'K'), Tmax=(996.605,'K')), NASAPolynomial(coeffs=[25.6259,0.0387884,-1.40757e-05,2.53058e-09,-1.77273e-13,-6951.51,-93.7823], Tmin=(996.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.986866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C([C]=C)C(O)(CC)C(=C)O(20154)',
    structure = SMILES('[CH2]C([C]=C)C(O)(CC)C(=C)O'),
    E0 = (14.9063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,411.055,724.254,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.46676,0.130155,-0.000117921,4.69382e-08,-4.32561e-12,2036.52,48.6023], Tmin=(100,'K'), Tmax=(957.092,'K')), NASAPolynomial(coeffs=[26.3948,0.0360113,-1.18734e-05,1.9759e-09,-1.31513e-13,-4700.86,-95.7079], Tmin=(957.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.9063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C=C)C(O)(CC)C(=C)[O](20156)',
    structure = SMILES('[CH2]C(C=C)C(O)(CC)C(=C)[O]'),
    E0 = (-85.1307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.5603,0.127299,-0.000120168,5.93212e-08,-1.15225e-11,-9988.29,49.9702], Tmin=(100,'K'), Tmax=(1276.34,'K')), NASAPolynomial(coeffs=[26.1857,0.0360116,-1.14763e-05,1.81319e-09,-1.14218e-13,-17228.6,-95.3445], Tmin=(1276.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-85.1307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C([CH2])C=C(20157)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C([CH2])C=C'),
    E0 = (24.1606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1556.5,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.18286,0.136944,-0.000136814,6.93583e-08,-1.35791e-11,3182.36,51.1069], Tmin=(100,'K'), Tmax=(1335.15,'K')), NASAPolynomial(coeffs=[30.2989,0.0295609,-8.2254e-06,1.18296e-09,-7.05231e-14,-5127.73,-117.74], Tmin=(1335.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.1606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(C)C(O)([CH]C)C(=C)O(32629)',
    structure = SMILES('[CH]=CC(C)C(O)([CH]C)C(=C)O'),
    E0 = (18.9804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,180,180,180,180,1600,1646.61,2842.08,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.3426,0.122319,-8.91799e-05,1.37959e-08,7.82678e-12,2526.67,48.8008], Tmin=(100,'K'), Tmax=(960.63,'K')), NASAPolynomial(coeffs=[28.5894,0.0327701,-1.06395e-05,1.8228e-09,-1.26342e-13,-5227.19,-108.612], Tmin=(960.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.9804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]C(O)(C(=C)O)C(C)C=C(32630)',
    structure = SMILES('[CH2][CH]C(O)(C(=C)O)C(C)C=C'),
    E0 = (-22.8694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,455.8,671.802,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155957,'amu*angstrom^2'), symmetry=1, barrier=(3.58576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155957,'amu*angstrom^2'), symmetry=1, barrier=(3.58576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155957,'amu*angstrom^2'), symmetry=1, barrier=(3.58576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155957,'amu*angstrom^2'), symmetry=1, barrier=(3.58576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155957,'amu*angstrom^2'), symmetry=1, barrier=(3.58576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155957,'amu*angstrom^2'), symmetry=1, barrier=(3.58576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155957,'amu*angstrom^2'), symmetry=1, barrier=(3.58576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155957,'amu*angstrom^2'), symmetry=1, barrier=(3.58576,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.23042,0.119516,-8.29256e-05,8.6045e-09,9.33142e-12,-2510.45,49.7815], Tmin=(100,'K'), Tmax=(962.49,'K')), NASAPolynomial(coeffs=[28.1243,0.0333399,-1.09229e-05,1.88313e-09,-1.30944e-13,-10205.3,-105.121], Tmin=(962.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.8694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = 'C=CC(C)C(O)([CH]C)C(=C)[O](32631)',
    structure = SMILES('C=CC(C)C(O)([CH]C)C(=C)[O]'),
    E0 = (-90.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93585,0.115218,-8.14416e-05,1.53544e-08,4.908e-12,-10634.6,48.4385], Tmin=(100,'K'), Tmax=(985.056,'K')), NASAPolynomial(coeffs=[24.7985,0.0386454,-1.35488e-05,2.37084e-09,-1.63133e-13,-17453.4,-88.0139], Tmin=(985.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)([CH]C)C(C)C=C(32632)',
    structure = SMILES('[CH]=C(O)C(O)([CH]C)C(C)C=C'),
    E0 = (18.9804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,180,180,180,180,1600,1646.61,2842.08,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152762,'amu*angstrom^2'), symmetry=1, barrier=(3.5123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.3426,0.122319,-8.91799e-05,1.37959e-08,7.82678e-12,2526.67,48.8008], Tmin=(100,'K'), Tmax=(960.63,'K')), NASAPolynomial(coeffs=[28.5894,0.0327701,-1.06395e-05,1.8228e-09,-1.26342e-13,-5227.19,-108.612], Tmin=(960.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.9804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])C(O)(CC)C(=C)O(20162)',
    structure = SMILES('[CH]=CC([CH2])C(O)(CC)C(=C)O'),
    E0 = (24.1606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1556.5,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.18286,0.136944,-0.000136814,6.93583e-08,-1.35791e-11,3182.36,51.1069], Tmin=(100,'K'), Tmax=(1335.15,'K')), NASAPolynomial(coeffs=[30.2988,0.0295609,-8.22542e-06,1.18296e-09,-7.05235e-14,-5127.72,-117.739], Tmin=(1335.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.1606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C=C)C1(O)[C](O)CC1C(32633)',
    structure = SMILES('[CH2]C(C=C)C1(O)[C](O)CC1C'),
    E0 = (24.5263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.82835,0.113062,-7.62797e-05,1.24942e-08,5.25028e-12,3173.05,44.2579], Tmin=(100,'K'), Tmax=(991.904,'K')), NASAPolynomial(coeffs=[23.5989,0.0414554,-1.47696e-05,2.59198e-09,-1.77898e-13,-3392.89,-85.8775], Tmin=(991.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.5263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(O)C(O)([CH]C)C1[CH]CC1(25394)',
    structure = SMILES('C=C(O)C(O)([CH]C)C1[CH]CC1'),
    E0 = (-10.9522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47563,0.0969741,-2.03252e-05,-5.25099e-08,2.97847e-11,-1098.92,46.6924], Tmin=(100,'K'), Tmax=(962.122,'K')), NASAPolynomial(coeffs=[26.8602,0.0352799,-1.16211e-05,2.07536e-09,-1.49595e-13,-9148.48,-102.414], Tmin=(962.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.9522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCJCO) + radical(cyclobutane)"""),
)

species(
    label = 'C=CC1CC[C](O)C1(O)[CH]C(32634)',
    structure = SMILES('C=CC1CC[C](O)C1(O)[CH]C'),
    E0 = (-53.3791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32435,0.0975654,-3.31334e-05,-3.15252e-08,2.05753e-11,-6210.87,43.6359], Tmin=(100,'K'), Tmax=(976.363,'K')), NASAPolynomial(coeffs=[23.5515,0.0404479,-1.4202e-05,2.53837e-09,-1.78739e-13,-13203.6,-86.7148], Tmin=(976.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.3791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C2CsJOH) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1[CH]CC(C)C1(O)C(=C)O(32635)',
    structure = SMILES('[CH2]C1[CH]CC(C)C1(O)C(=C)O'),
    E0 = (-89.8763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19536,0.0874673,1.12041e-05,-9.06742e-08,4.55048e-11,-10597.9,43.7078], Tmin=(100,'K'), Tmax=(933.006,'K')), NASAPolynomial(coeffs=[27.7742,0.0313514,-8.03546e-06,1.28457e-09,-9.25178e-14,-18967,-109.925], Tmin=(933.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.8763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(Isobutyl) + radical(Cs_S)"""),
)

species(
    label = 'C=CC(=C)C(O)(CC)C(=C)O(20165)',
    structure = SMILES('C=CC(=C)C(O)(CC)C(=C)O'),
    E0 = (-328.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.41239,0.122132,-8.48259e-05,9.65445e-09,8.60192e-12,-39259.9,42.5728], Tmin=(100,'K'), Tmax=(983.692,'K')), NASAPolynomial(coeffs=[29.3514,0.0335909,-1.17542e-05,2.11114e-09,-1.49685e-13,-47474.4,-120.137], Tmin=(983.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(C)C(O)(C=C)C(=C)O(32636)',
    structure = SMILES('C=CC(C)C(O)(C=C)C(=C)O'),
    E0 = (-306.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.04395,0.114003,-6.99125e-05,5.6691e-10,1.01527e-11,-36622.1,44.5389], Tmin=(100,'K'), Tmax=(1000.72,'K')), NASAPolynomial(coeffs=[27.1033,0.036977,-1.36312e-05,2.49478e-09,-1.77194e-13,-44432.5,-105.978], Tmin=(1000.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-306.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(C=C)C(C)[C](O)C(=C)O(32187)',
    structure = SMILES('[CH2]C(O)=C(O)C(C)C([CH2])C=C'),
    E0 = (-100.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.86515,0.139886,-0.000137382,6.66597e-08,-1.22112e-11,-11829.7,51.6786], Tmin=(100,'K'), Tmax=(1520.55,'K')), NASAPolynomial(coeffs=[33.576,0.0231354,-4.20037e-06,3.71392e-10,-1.42037e-14,-21105.5,-137.704], Tmin=(1520.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C(C=C)[C](O)C(C)C(=C)O(32637)',
    structure = SMILES('[CH2]C(C=C)[C](O)C(C)C(=C)O'),
    E0 = (-28.9034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.63686,0.133744,-0.000134273,7.0058e-08,-1.43875e-11,-3226.93,48.7265], Tmin=(100,'K'), Tmax=(1191.5,'K')), NASAPolynomial(coeffs=[26.076,0.0373516,-1.29232e-05,2.16025e-09,-1.41249e-13,-10069.2,-94.8273], Tmin=(1191.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.9034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C2CsJOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=CCC(O)([CH]C)C(=C)O(20177)',
    structure = SMILES('[CH2]C=CCC(O)([CH]C)C(=C)O'),
    E0 = (-80.8901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,379.414,746.693,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25836,0.119608,-8.11805e-05,6.45143e-09,1.00411e-11,-9487.27,48.1827], Tmin=(100,'K'), Tmax=(965.735,'K')), NASAPolynomial(coeffs=[28.3598,0.0336602,-1.11661e-05,1.94135e-09,-1.35617e-13,-17307,-108.332], Tmin=(965.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.8901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(Allyl_P)"""),
)

species(
    label = 'C=CC[CH]C(O)([CH]C)C(=C)O(32638)',
    structure = SMILES('C=CC[CH]C(O)([CH]C)C(=C)O'),
    E0 = (-20.2418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1077.04,1130.01,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.16119,'amu*angstrom^2'), symmetry=1, barrier=(3.70608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16119,'amu*angstrom^2'), symmetry=1, barrier=(3.70608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16119,'amu*angstrom^2'), symmetry=1, barrier=(3.70608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16119,'amu*angstrom^2'), symmetry=1, barrier=(3.70608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16119,'amu*angstrom^2'), symmetry=1, barrier=(3.70608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16119,'amu*angstrom^2'), symmetry=1, barrier=(3.70608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16119,'amu*angstrom^2'), symmetry=1, barrier=(3.70608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16119,'amu*angstrom^2'), symmetry=1, barrier=(3.70608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12717,0.116868,-7.72231e-05,3.87145e-09,1.06898e-11,-2197.82,50.366], Tmin=(100,'K'), Tmax=(965.048,'K')), NASAPolynomial(coeffs=[27.7827,0.033618,-1.112e-05,1.93185e-09,-1.34932e-13,-9866.97,-102.692], Tmin=(965.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.2418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(CCJCO)"""),
)

species(
    label = 'C=CC1CC(C)C1(O)C(=C)O(32193)',
    structure = SMILES('C=CC1CC(C)C1(O)C(=C)O'),
    E0 = (-278.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.4649,0.0904109,1.26798e-05,-9.65895e-08,4.79578e-11,-33317.8,40.6531], Tmin=(100,'K'), Tmax=(943.857,'K')), NASAPolynomial(coeffs=[30.6343,0.0287942,-7.66367e-06,1.31352e-09,-9.95047e-14,-42692,-129.912], Tmin=(943.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-278.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C(C=C)C(O)([CH]C)C(C)=O(32639)',
    structure = SMILES('[CH2]C(C=C)C(O)([CH]C)C(C)=O'),
    E0 = (-32.3095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.92647,0.123347,-0.000116316,5.88072e-08,-1.19753e-11,-3666.65,47.7829], Tmin=(100,'K'), Tmax=(1184.4,'K')), NASAPolynomial(coeffs=[21.2355,0.045123,-1.72473e-05,3.0439e-09,-2.04869e-13,-9153.23,-67.8801], Tmin=(1184.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.3095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C=CCCC(O)=C(O)[CH]C(20034)',
    structure = SMILES('[CH2]C=CCCC(O)=C(O)[CH]C'),
    E0 = (-108.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.29519,0.117368,-6.79165e-05,-1.18977e-08,1.75545e-11,-12829.4,47.0514], Tmin=(100,'K'), Tmax=(958.062,'K')), NASAPolynomial(coeffs=[30.4474,0.0303824,-9.56779e-06,1.6662e-09,-1.19117e-13,-21385,-121.418], Tmin=(958.062,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCJCO) + radical(Allyl_P)"""),
)

species(
    label = 'CHCH3(T)(95)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438699,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82363,-0.000909515,3.2138e-05,-3.7348e-08,1.3309e-11,41371.4,7.10948], Tmin=(100,'K'), Tmax=(960.812,'K')), NASAPolynomial(coeffs=[4.30487,0.00943069,-3.27559e-06,5.95121e-10,-4.27307e-14,40709.1,1.84202], Tmin=(960.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(C=C)[C](O)C(=C)O(32640)',
    structure = SMILES('[CH2]C(O)=C(O)C([CH2])C=C'),
    E0 = (-45.9637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95109,0.114701,-0.000128306,7.0741e-08,-1.47573e-11,-5299.99,38.156], Tmin=(100,'K'), Tmax=(1302.41,'K')), NASAPolynomial(coeffs=[26.5428,0.016472,-2.83065e-06,1.95552e-10,-3.10703e-15,-11813.2,-103.35], Tmin=(1302.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.9637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(O)CJ)"""),
)

species(
    label = 'CH2(19)',
    structure = SMILES('[CH2]'),
    E0 = (381.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1032.72,2936.3,3459],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8328,0.000224446,4.68033e-06,-6.04743e-09,2.59009e-12,45920.8,1.40666], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16229,0.00281798,-7.56235e-07,5.05446e-11,5.65236e-15,46099.1,4.77656], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(381.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=C[CH]C(O)([CH]C)C(=C)O(32641)',
    structure = SMILES('C=C[CH]C(O)([CH]C)C(=C)O'),
    E0 = (-126.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
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
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31481,0.0918562,-1.7451e-05,-6.12016e-08,3.52131e-11,-14966.5,43.5086], Tmin=(100,'K'), Tmax=(940.955,'K')), NASAPolynomial(coeffs=[29.9439,0.0215879,-5.24657e-06,8.6838e-10,-6.67358e-14,-23620.9,-120.123], Tmin=(940.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-126.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(CCJCO)"""),
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
    E0 = (-23.0335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (47.3978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (14.2041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (45.7922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (31.7769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (83.2095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (113.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-23.0335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (120.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-4.44454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (121.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (104.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (100.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (130.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (82.5285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (121.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (161.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (60.6465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (59.2148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (61.0649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (68.4692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (63.2889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (31.3134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (19.6433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (52.0206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (62.3834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (151.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (101.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (101.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (92.0265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (35.1241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (40.3667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (1.93979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (134.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (134.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (136.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (225.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-14.7491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (181.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (124.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (297.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (255.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['butadiene13(2459)', 'C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C(C=C)C1(O)C(C)C1([CH2])O(32620)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.39513e+10,'s^-1'), n=0.560608, Ea=(70.4313,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S_D;doublebond_intra_2H;radadd_intra_csHNd] + [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_cs] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_csHNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C1CC1C(O)([CH]C)C(=C)O(32621)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C1(O)CC(C=C)C1(O)[CH]C(32622)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(68.8257,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 349 used for R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 68.1 to 68.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C1C([CH2])C(O)(C(=C)O)C1C(32623)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.73e+06,'s^-1'), n=1.31, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 344 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_csHNd
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=CC(=C)C(O)([CH]C)C(=C)O(32624)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.88727,'m^3/(mol*s)'), n=1.9687, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 18 used for Cds-CdCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -8.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2]C(C=C)C(O)(C=C)C(=C)O(32625)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.014e+08,'cm^3/(mol*s)'), n=1.733, Ea=(3.17984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2823 used for Cds-HH_Cds-Cs\O2s/H;HJ
Exact match found for rate rule [Cds-HH_Cds-Cs\O2s/H;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][CH]C=C(2458)', 'C=C(O)C(O)=CC(5562)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00784329,'m^3/(mol*s)'), n=2.41519, Ea=(72.1429,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-CdOs_Cds;CJ]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 65.4 to 72.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2COH(99)', '[CH2]C(C=C)C(O)=CC(32626)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00712612,'m^3/(mol*s)'), n=2.40979, Ea=(7.81798,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CdsJ] for rate rule [Cds-OsCs_Cds;CdsJ-O2s]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['butadiene13(2459)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.0136943,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C2H3(30)', 'C=CC(O)([CH]C)C(=C)O(32291)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6870,'cm^3/(mol*s)'), n=2.41, Ea=(13.7235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CdsJ-H] for rate rule [Cds-Cs\O2s/H_Cds-HH;CdsJ-H]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OH(5)', '[CH2]C(C=C)C(=CC)C(=C)O(32627)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.567844,'m^3/(mol*s)'), n=2.025, Ea=(15.9149,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdCs_Cds-CsH;YJ] for rate rule [Cds-CdCs_Cds-CsH;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[C](C)C(O)([CH]C)C(=C)O(32628)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CC(O)(C(=C)O)C([CH2])C=C(20155)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2][C](C=C)C(O)(CC)C(=C)O(20150)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(588307,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]C(C)C(O)([CH]C)C(=C)O(26113)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['C=CC(C)C([O])([CH]C)C(=C)O(20153)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([C]=C)C(O)(CC)C(=C)O(20154)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C(C=C)C(O)(CC)C(=C)[O](20156)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;O_H_out] for rate rule [R4H_SS(Cd)S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(O)C(O)(CC)C([CH2])C=C(20157)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=CC(C)C(O)([CH]C)C(=C)O(32629)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]C(O)(C(=C)O)C(C)C=C(32630)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['C=CC(C)C(O)([CH]C)C(=C)[O](32631)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);C_rad_out_2H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(O)C(O)([CH]C)C(C)C=C(32632)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=CC([CH2])C(O)(CC)C(=C)O(20162)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(315594,'s^-1'), n=1.73223, Ea=(38.2227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;Cd_rad_out_singleH;Cs_H_out] + [R5H_RSSR;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]C=C(2458)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C(C=C)C1(O)[C](O)CC1C(32633)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.71e+08,'s^-1'), n=0.99, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra_secNd_2H;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['C=C(O)C(O)([CH]C)C1[CH]CC1(25394)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['C=CC1CC[C](O)C1(O)[CH]C(32634)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.71e+11,'s^-1'), n=0.2, Ea=(115.06,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C1[CH]CC(C)C1(O)C(=C)O(32635)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['C=CC(=C)C(O)(CC)C(=C)O(20165)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['C=CC(C)C(O)(C=C)C(=C)O(32636)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C(C=C)C(C)[C](O)C(=C)O(32187)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C(C=C)[C](O)C(C)C(=C)O(32637)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C=CCC(O)([CH]C)C(=C)O(20177)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=CC[CH]C(O)([CH]C)C(=C)O(32638)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['C=CC1CC(C)C1(O)C(=C)O(32193)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(C)=O(32639)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C=CCCC(O)=C(O)[CH]C(20034)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CHCH3(T)(95)', '[CH2]C(C=C)[C](O)C(=C)O(32640)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction42',
    reactants = ['CH2(19)', 'C=C[CH]C(O)([CH]C)C(=C)O(32641)'],
    products = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '5575',
    isomers = [
        '[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)',
    ],
    reactants = [
        ('butadiene13(2459)', 'C=C(O)C(O)=CC(5562)'),
        ('[CH2][CH]C=C(2458)', 'C=C(O)C(O)=CC(5562)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '5575',
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

