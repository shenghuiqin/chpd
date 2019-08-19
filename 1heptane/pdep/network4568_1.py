species(
    label = '[CH]=C(C=C)C=[C]C=C(26596)',
    structure = SMILES('[CH]=C(C=C)C=[C]C=C'),
    E0 = (662.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.37975,'amu*angstrom^2'), symmetry=1, barrier=(31.7232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38561,'amu*angstrom^2'), symmetry=1, barrier=(31.8578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37601,'amu*angstrom^2'), symmetry=1, barrier=(31.6371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371162,0.0693608,-3.99084e-05,-2.74722e-09,7.74147e-12,79839.3,26.692], Tmin=(100,'K'), Tmax=(961.871,'K')), NASAPolynomial(coeffs=[16.7573,0.0252327,-8.54203e-06,1.46909e-09,-1.00661e-13,75576.1,-57.4993], Tmin=(961.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(662.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJC=C)"""),
)

species(
    label = 'CH2CHCCH(26391)',
    structure = SMILES('C#CC=C'),
    E0 = (274.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,3010,987.5,1337.5,450,1655,2175,525],'cm^-1')),
        HinderedRotor(inertia=(1.46338,'amu*angstrom^2'), symmetry=1, barrier=(33.6459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.87083,0.0182042,1.06711e-05,-2.72492e-08,1.19478e-11,33023.8,11.2934], Tmin=(100,'K'), Tmax=(955.249,'K')), NASAPolynomial(coeffs=[8.52653,0.0108962,-3.56564e-06,6.31243e-10,-4.51891e-14,31196.2,-19.6435], Tmin=(955.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""CH2CHCCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C1C=C1C=[C]C=C(28031)',
    structure = SMILES('[CH2]C1C=C1C=[C]C=C'),
    E0 = (747.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.586472,0.0880494,-9.04612e-05,4.77122e-08,-9.67689e-12,90036.1,27.551], Tmin=(100,'K'), Tmax=(1322.89,'K')), NASAPolynomial(coeffs=[19.8864,0.01929,-4.72239e-06,5.86758e-10,-3.07438e-14,85219.3,-74.6799], Tmin=(1322.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(747.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Isobutyl) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C1C=C(C=C)C1[CH2](28032)',
    structure = SMILES('[CH]=C1C=C(C=C)C1[CH2]'),
    E0 = (688.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.220199,0.0672244,-1.87323e-05,-3.73533e-08,2.35297e-11,82963.1,24.589], Tmin=(100,'K'), Tmax=(920.754,'K')), NASAPolynomial(coeffs=[21.5221,0.0162097,-3.27531e-06,4.37543e-10,-3.06656e-14,77280,-85.9808], Tmin=(920.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(688.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[C]=CC(C=C)=C1(28033)',
    structure = SMILES('[CH2]C1[C]=CC(C=C)=C1'),
    E0 = (635.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18838,0.0506356,-2.03277e-06,-2.97609e-08,1.4152e-11,76545.3,23.4251], Tmin=(100,'K'), Tmax=(1000.06,'K')), NASAPolynomial(coeffs=[12.1738,0.0319869,-1.19945e-05,2.16712e-09,-1.50999e-13,73083.5,-35.8966], Tmin=(1000.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(635.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentadiene) + radical(1,3-cyclopentadiene-vinyl-1) + radical(Isobutyl)"""),
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
    label = '[CH]=C(C#CC=C)C=C(28034)',
    structure = SMILES('[CH]=C(C#CC=C)C=C'),
    E0 = (626.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2100,2250,500,550,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.38383,'amu*angstrom^2'), symmetry=1, barrier=(31.8171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38442,'amu*angstrom^2'), symmetry=1, barrier=(31.8306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38269,'amu*angstrom^2'), symmetry=1, barrier=(31.7908,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.570824,0.0587501,-7.31768e-06,-4.03404e-08,2.15142e-11,75543.9,27.7802], Tmin=(100,'K'), Tmax=(969.523,'K')), NASAPolynomial(coeffs=[19.9984,0.0183224,-6.2313e-06,1.17492e-09,-8.85641e-14,69909.7,-74.9741], Tmin=(969.523,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)Ct) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-Ct(Cds-Cds)) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)C=C=C=C(28035)',
    structure = SMILES('[CH]=C(C=C)C=C=C=C'),
    E0 = (693.047,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,563.333,586.667,610,1970,2140,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.42967,'amu*angstrom^2'), symmetry=1, barrier=(32.871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42584,'amu*angstrom^2'), symmetry=1, barrier=(32.783,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.308611,0.0712082,-5.98117e-05,2.57301e-08,-4.40771e-12,83495.6,26.4146], Tmin=(100,'K'), Tmax=(1404.09,'K')), NASAPolynomial(coeffs=[17.0356,0.0235564,-8.90519e-06,1.55969e-09,-1.04165e-13,78798.3,-59.9605], Tmin=(1404.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(693.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Cds_P)"""),
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
    label = 'C#CC=[C]C=C(27472)',
    structure = SMILES('C#CC=[C]C=C'),
    E0 = (521.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.62249,'amu*angstrom^2'), symmetry=1, barrier=(37.3042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61979,'amu*angstrom^2'), symmetry=1, barrier=(37.2422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45098,0.0475986,-2.51853e-05,-8.16809e-09,8.69319e-12,62818.4,18.2704], Tmin=(100,'K'), Tmax=(937.719,'K')), NASAPolynomial(coeffs=[14.4661,0.0131635,-3.8272e-06,6.24214e-10,-4.31862e-14,59450.6,-48.6255], Tmin=(937.719,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(521.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=[C]C=C(4699)',
    structure = SMILES('[CH]=C=C[CH2]'),
    E0 = (451.584,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3000,3100,440,815,1455,1000,180,1024.85,1025.53,1026.61],'cm^-1')),
        HinderedRotor(inertia=(0.00938781,'amu*angstrom^2'), symmetry=1, barrier=(7.01846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76805,0.020302,8.75519e-06,-2.87666e-08,1.37354e-11,54363.7,13.5565], Tmin=(100,'K'), Tmax=(915.031,'K')), NASAPolynomial(coeffs=[9.46747,0.00887314,-1.78262e-06,2.38534e-10,-1.6263e-14,52390.1,-22.2544], Tmin=(915.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C=CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C([C]=CC=C)C=C(28036)',
    structure = SMILES('[CH]C(=C=CC=C)C=C'),
    E0 = (619.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0781588,0.0767216,-4.08012e-05,-7.77164e-09,1.02254e-11,74650.9,27.8369], Tmin=(100,'K'), Tmax=(973.353,'K')), NASAPolynomial(coeffs=[19.0312,0.0273638,-9.69385e-06,1.71369e-09,-1.19418e-13,69549,-70.9373], Tmin=(973.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(619.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C=C)C=C[C]=C(28037)',
    structure = SMILES('[CH]C(C=C)=CC=C=C'),
    E0 = (619.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0781588,0.0767216,-4.08012e-05,-7.77164e-09,1.02254e-11,74650.9,27.8369], Tmin=(100,'K'), Tmax=(973.353,'K')), NASAPolynomial(coeffs=[19.0312,0.0273638,-9.69385e-06,1.71369e-09,-1.19418e-13,69549,-70.9373], Tmin=(973.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(619.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]C(=C)C=[C]C=C(28038)',
    structure = SMILES('[CH2]C(=C=C)C=[C]C=C'),
    E0 = (599.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.145655,0.0807195,-6.12846e-05,1.1753e-08,4.44743e-12,72223,26.8477], Tmin=(100,'K'), Tmax=(946.999,'K')), NASAPolynomial(coeffs=[19.7998,0.021376,-6.7332e-06,1.11901e-09,-7.60755e-14,67328.7,-74.1871], Tmin=(946.999,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(599.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[C]=[C]C(=C)C=C(28039)',
    structure = SMILES('[CH2]C=C([CH2])C#CC=C'),
    E0 = (523.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736313,0.0603737,-2.12413e-05,-1.58434e-08,1.05237e-11,63101.4,27.4416], Tmin=(100,'K'), Tmax=(999.488,'K')), NASAPolynomial(coeffs=[14.9048,0.0284886,-1.06347e-05,1.92508e-09,-1.34711e-13,59029.5,-47.107], Tmin=(999.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(523.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-Ct(Cds-Cds)) + radical(Allyl_P) + radical(CTCC=CCJ)"""),
)

species(
    label = '[CH]=CC=CC(=[CH])C=C(28040)',
    structure = SMILES('[CH]=CC=CC(=[CH])C=C'),
    E0 = (710.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3115,3125,620,680,785,800,1600,1700,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.13046,'amu*angstrom^2'), symmetry=1, barrier=(25.9915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13174,'amu*angstrom^2'), symmetry=1, barrier=(26.021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12607,'amu*angstrom^2'), symmetry=1, barrier=(25.8905,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.379576,0.066288,-2.72767e-05,-1.85171e-08,1.37212e-11,85626.8,27.4684], Tmin=(100,'K'), Tmax=(965.393,'K')), NASAPolynomial(coeffs=[18.5222,0.0223396,-7.50526e-06,1.3318e-09,-9.47259e-14,80668.8,-66.9561], Tmin=(965.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(710.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([C]=C)C=CC=C(28041)',
    structure = SMILES('[CH]C(=C=C)C=CC=C'),
    E0 = (619.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0781588,0.0767216,-4.08012e-05,-7.77164e-09,1.02254e-11,74650.9,27.8369], Tmin=(100,'K'), Tmax=(973.353,'K')), NASAPolynomial(coeffs=[19.0312,0.0273638,-9.69385e-06,1.71369e-09,-1.19418e-13,69549,-70.9373], Tmin=(973.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(619.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC(=C)C=[C]C=C(28042)',
    structure = SMILES('[CH]=CC(=C)C=[C]C=C'),
    E0 = (662.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.37975,'amu*angstrom^2'), symmetry=1, barrier=(31.7232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38561,'amu*angstrom^2'), symmetry=1, barrier=(31.8578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37601,'amu*angstrom^2'), symmetry=1, barrier=(31.6371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371162,0.0693608,-3.99084e-05,-2.74722e-09,7.74147e-12,79839.3,26.692], Tmin=(100,'K'), Tmax=(961.871,'K')), NASAPolynomial(coeffs=[16.7573,0.0252327,-8.54203e-06,1.46909e-09,-1.00661e-13,75576.1,-57.4993], Tmin=(961.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(662.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CC(=[CH])C=CC=C(28043)',
    structure = SMILES('[CH]=CC(=[CH])C=CC=C'),
    E0 = (710.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3115,3125,620,680,785,800,1600,1700,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.13046,'amu*angstrom^2'), symmetry=1, barrier=(25.9915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13174,'amu*angstrom^2'), symmetry=1, barrier=(26.021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12607,'amu*angstrom^2'), symmetry=1, barrier=(25.8905,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.379576,0.066288,-2.72767e-05,-1.85171e-08,1.37212e-11,85626.8,27.4684], Tmin=(100,'K'), Tmax=(965.393,'K')), NASAPolynomial(coeffs=[18.5222,0.0223396,-7.50526e-06,1.3318e-09,-9.47259e-14,80668.8,-66.9561], Tmin=(965.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(710.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][C]=CC(=C)C=C(28044)',
    structure = SMILES('[CH2]C#CC=C([CH2])C=C'),
    E0 = (520.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.425434,0.0672574,-3.47703e-05,-6.82322e-09,8.60793e-12,62744.2,27.0824], Tmin=(100,'K'), Tmax=(983.131,'K')), NASAPolynomial(coeffs=[16.8429,0.0254472,-9.10161e-06,1.62188e-09,-1.13253e-13,58308.6,-57.9841], Tmin=(983.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(520.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Propargyl) + radical(CTCC=CCJ)"""),
)

species(
    label = '[CH]=C[C]=CC(=C)C=C(28045)',
    structure = SMILES('[CH]C=C=CC(=C)C=C'),
    E0 = (640.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.489027,0.0663625,-2.86695e-05,-6.41474e-09,6.20934e-12,77113.6,27.4115], Tmin=(100,'K'), Tmax=(1060.43,'K')), NASAPolynomial(coeffs=[13.988,0.0357909,-1.42073e-05,2.58775e-09,-1.78901e-13,73106.7,-43.8995], Tmin=(1060.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(640.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C=C)C=C1[CH]C1(28046)',
    structure = SMILES('[CH]=C(C=C)C=C1[CH]C1'),
    E0 = (697.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02472,0.0479685,2.16474e-05,-6.46805e-08,2.87353e-11,84048.5,23.8245], Tmin=(100,'K'), Tmax=(967.63,'K')), NASAPolynomial(coeffs=[16.8577,0.0250208,-8.66681e-06,1.59933e-09,-1.17197e-13,78994.6,-62.3214], Tmin=(967.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(697.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(Cds_P)"""),
)

species(
    label = 'C=C[C]=CC1[CH]CC=1(28047)',
    structure = SMILES('C=C[C]=CC1[CH]CC=1'),
    E0 = (609.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.753829,0.0509664,2.62305e-05,-8.00016e-08,3.71916e-11,73499.9,22.8142], Tmin=(100,'K'), Tmax=(937.071,'K')), NASAPolynomial(coeffs=[20.476,0.0184759,-4.51274e-06,7.42881e-10,-5.66918e-14,67533.9,-83.1629], Tmin=(937.071,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(609.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CJC=C) + radical(cyclobutene-allyl)"""),
)

species(
    label = '[CH]=C1[CH]CC(C=C)=C1(28048)',
    structure = SMILES('[CH]C1=CCC(C=C)=C1'),
    E0 = (505.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.70217,0.0504982,3.83219e-05,-9.11671e-08,3.99702e-11,60940.1,24.1235], Tmin=(100,'K'), Tmax=(953.212,'K')), NASAPolynomial(coeffs=[19.7443,0.0252798,-8.05364e-06,1.45688e-09,-1.08352e-13,54825.3,-79.8639], Tmin=(953.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(505.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentadiene) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC1C=[C][CH]CC=1(28049)',
    structure = SMILES('[CH2]C=C1[CH]CC=C=C1'),
    E0 = (436.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75442,0.0218868,0.000107088,-1.51629e-07,5.80456e-11,52579.6,18.4459], Tmin=(100,'K'), Tmax=(974.105,'K')), NASAPolynomial(coeffs=[16.8297,0.0292657,-1.0962e-05,2.17839e-09,-1.676e-13,46355.5,-70.7605], Tmin=(974.105,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(Cyclohexane) + radical(C=CC=CCJ) + radical(Allyl_S)"""),
)

species(
    label = 'C=C=C=CC(=C)C=C(28050)',
    structure = SMILES('C=C=C=CC(=C)C=C'),
    E0 = (445.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.533444,0.0660312,-3.79807e-05,2.53388e-09,3.6017e-12,53769.1,25.1965], Tmin=(100,'K'), Tmax=(1058.68,'K')), NASAPolynomial(coeffs=[15.1956,0.0289278,-1.13313e-05,2.07523e-09,-1.44572e-13,49639.4,-51.218], Tmin=(1058.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC1=CC(C=C)=C1(28051)',
    structure = SMILES('C=CC1=CC(C=C)=C1'),
    E0 = (565.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45976,0.0407616,2.84501e-05,-6.13606e-08,2.49507e-11,68091,20.5706], Tmin=(100,'K'), Tmax=(996.557,'K')), NASAPolynomial(coeffs=[12.7845,0.0315263,-1.21666e-05,2.28148e-09,-1.63969e-13,64035.3,-43.0492], Tmin=(996.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(cyclobutadiene_13)"""),
)

species(
    label = '[CH2]C=[C]C1C=C1C=C(28018)',
    structure = SMILES('[CH2]C=[C]C1C=C1C=C'),
    E0 = (754.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.153884,0.0753571,-5.95682e-05,1.96286e-08,-9.82161e-13,90877.4,26.6776], Tmin=(100,'K'), Tmax=(1034.67,'K')), NASAPolynomial(coeffs=[17.0605,0.0253762,-9.40549e-06,1.67364e-09,-1.15025e-13,86555.6,-59.4416], Tmin=(1034.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(754.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C1C=C=CCC1[CH2](28052)',
    structure = SMILES('[CH]=C1C=C=CCC1[CH2]'),
    E0 = (633.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20193,0.0380142,6.25694e-05,-1.10239e-07,4.50518e-11,76323.7,21.6618], Tmin=(100,'K'), Tmax=(965.943,'K')), NASAPolynomial(coeffs=[18.3152,0.0255846,-8.87497e-06,1.70037e-09,-1.29397e-13,70291.4,-74.4187], Tmin=(965.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclohexane) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)C=C=[C]C(28053)',
    structure = SMILES('[CH]C(C=C)=CC#CC'),
    E0 = (613.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.290503,0.0688901,-2.46667e-05,-2.02787e-08,1.38757e-11,73919,27.9688], Tmin=(100,'K'), Tmax=(963.786,'K')), NASAPolynomial(coeffs=[16.7277,0.030262,-1.06017e-05,1.84884e-09,-1.27448e-13,69376.3,-57.8546], Tmin=(963.786,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([C]=C=CC)C=C(28054)',
    structure = SMILES('[CH]=C(C#C[CH]C)C=C'),
    E0 = (652.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.472264,0.0655935,-2.64102e-05,-2.06233e-08,1.5515e-11,78560.5,29.0043], Tmin=(100,'K'), Tmax=(925.502,'K')), NASAPolynomial(coeffs=[17.6219,0.0223564,-6.3877e-06,1.00903e-09,-6.78004e-14,74063.4,-59.5505], Tmin=(925.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(652.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Ct) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Cds_P) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=C([C]=C)C=C=CC(28055)',
    structure = SMILES('[CH]C(=C=C)C=C=CC'),
    E0 = (672.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,540,563.333,586.667,610,1970,2140,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.99972,'amu*angstrom^2'), symmetry=1, barrier=(45.9775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.9997,'amu*angstrom^2'), symmetry=1, barrier=(45.9771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00375,'amu*angstrom^2'), symmetry=1, barrier=(46.0702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0479731,0.0820058,-6.92863e-05,3.11371e-08,-5.72793e-12,80986.9,27.4159], Tmin=(100,'K'), Tmax=(1286.53,'K')), NASAPolynomial(coeffs=[15.2653,0.0346934,-1.4124e-05,2.55281e-09,-1.73449e-13,77071.3,-49.8332], Tmin=(1286.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC(=[CH])C=C=CC(28056)',
    structure = SMILES('[CH]=CC(=[CH])C=C=CC'),
    E0 = (763.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,3115,3125,620,680,785,800,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680],'cm^-1')),
        HinderedRotor(inertia=(1.03559,'amu*angstrom^2'), symmetry=1, barrier=(23.8101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02922,'amu*angstrom^2'), symmetry=1, barrier=(23.6639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03355,'amu*angstrom^2'), symmetry=1, barrier=(23.7633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.38342,0.0729708,-6.04511e-05,2.61958e-08,-4.59586e-12,91968.1,27.489], Tmin=(100,'K'), Tmax=(1354.69,'K')), NASAPolynomial(coeffs=[15.3964,0.0286421,-1.13676e-05,2.04096e-09,-1.38233e-13,87900.5,-49.4972], Tmin=(1354.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(763.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C1[CH]C(C=C)=C1(28057)',
    structure = SMILES('[CH2]C=C1[CH]C(C=C)=C1'),
    E0 = (453.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12331,0.032393,9.49875e-05,-1.62268e-07,6.94201e-11,54698.1,21.9172], Tmin=(100,'K'), Tmax=(921.531,'K')), NASAPolynomial(coeffs=[24.554,0.0104315,9.3705e-07,-3.28625e-10,1.41321e-14,46993.7,-107.579], Tmin=(921.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.686,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CCJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=C(C=C)C1[C]=CC1(28058)',
    structure = SMILES('[CH]=C(C=C)C1[C]=CC1'),
    E0 = (765.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.467061,0.0655391,-3.04743e-05,-1.10072e-08,9.8941e-12,92261.5,25.9267], Tmin=(100,'K'), Tmax=(991.847,'K')), NASAPolynomial(coeffs=[17.1881,0.0247314,-9.02691e-06,1.64258e-09,-1.1637e-13,87634.8,-61.2079], Tmin=(991.847,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(765.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1[CH]CCC=C=C1(28059)',
    structure = SMILES('[CH]C1C=C=CCCC=1'),
    E0 = (628.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06553,0.047489,2.97765e-05,-6.80472e-08,2.80264e-11,75753.2,20.1183], Tmin=(100,'K'), Tmax=(992.505,'K')), NASAPolynomial(coeffs=[14.154,0.0355073,-1.37281e-05,2.56046e-09,-1.83332e-13,71147.2,-53.0429], Tmin=(992.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(628.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#C[CH]CC([CH2])C#C(26589)',
    structure = SMILES('C#C[CH]CC([CH2])C#C'),
    E0 = (690.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.413,0.0775862,-6.82273e-05,2.00101e-08,4.75277e-12,83180.4,29.4531], Tmin=(100,'K'), Tmax=(770.958,'K')), NASAPolynomial(coeffs=[14.2582,0.0259108,-6.9066e-06,8.99538e-10,-4.79033e-14,80446.5,-37.6258], Tmin=(770.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(690.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Isobutyl) + radical(Sec_Propargyl)"""),
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
    label = '[CH]=C=CC(=[CH])C=C(28060)',
    structure = SMILES('[CH]C(C=C)=CC#C'),
    E0 = (655.573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.02047,'amu*angstrom^2'), symmetry=1, barrier=(46.4547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02183,'amu*angstrom^2'), symmetry=1, barrier=(46.4859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02095,'amu*angstrom^2'), symmetry=1, barrier=(46.4656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.1225,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.690423,0.0609185,-2.60113e-05,-1.49293e-08,1.15088e-11,78977,23.1129], Tmin=(100,'K'), Tmax=(972.069,'K')), NASAPolynomial(coeffs=[16.7842,0.0219118,-7.82026e-06,1.39954e-09,-9.87604e-14,74562.2,-60.6886], Tmin=(972.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
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
    E0 = (662.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (749.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (716.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (675.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (852.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (921.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (827.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (752.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (874.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (877.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (884.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (854.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (933.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (831.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1066.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (743.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (971.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (853.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (903.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (800.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (800.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (700.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (675.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (687.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (670.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (754.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (850.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (870.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (900.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (814.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (895.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (800.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (781.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (689.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (739.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (1037.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['CH2CHCCH(26391)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH2]C1C=C1C=[C]C=C(28031)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(86.3586,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 14 used for R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH
Exact match found for rate rule [R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH]=C1C=C(C=C)C1[CH2](28032)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra_2H_pri;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra_2H_pri;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH2]C1[C]=CC(C=C)=C1(28033)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.4226e+09,'s^-1'), n=0.3735, Ea=(13.1942,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH]=C(C#CC=C)C=C(28034)'],
    products = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2371.94,'m^3/(mol*s)'), n=1.49517, Ea=(13.5032,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-De_Ct-De;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]=C(C=C)C=C=C=C(28035)'],
    products = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H3(30)', 'C#CC=[C]C=C(27472)'],
    products = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0420833,'m^3/(mol*s)'), n=2.41, Ea=(19.4765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-De_Ct-H;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]C=C(4699)', 'CH2CHCCH(26391)'],
    products = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0345031,'m^3/(mol*s)'), n=2.33733, Ea=(26.7584,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-Cd_Ct-H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH]=C([C]=CC=C)C=C(28036)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.0945e+07,'s^-1'), n=1.951, Ea=(212.263,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_single] for rate rule [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH]=C(C=C)C=C[C]=C(28037)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['C=[C]C(=C)C=[C]C=C(28038)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['C=C[C]=[C]C(=C)C=C(28039)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC=CC(=[CH])C=C(28040)'],
    products = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH]=C([C]=C)C=CC=C(28041)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.66043e+07,'s^-1'), n=1.66125, Ea=(169.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_RSR;Cd_rad_out_singleDe_Cd;XH_out] + [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cd_H_out_doubleC]
Euclidian distance = 2.82842712475
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC(=C)C=[C]C=C(28042)'],
    products = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.27529e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_singleH] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=CC(=[CH])C=CC=C(28043)'],
    products = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSD;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 3.60555127546
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['C=[C][C]=CC(=C)C=C(28044)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH]=C[C]=CC(=C)C=C(28045)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R6HJ_3;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]C=C(4699)', '[CH]=[C]C=C(4699)'],
    products = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH]=C(C=C)C=C1[CH]C1(28046)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['C=C[C]=CC1[CH]CC=1(28047)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R4_D_D;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH]=C1[CH]CC(C=C)=C1(28048)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.80657e+08,'s^-1'), n=0.835, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra_pri_2H;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra_pri_2H;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['C=CC1C=[C][CH]CC=1(28049)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.16959e+10,'s^-1'), n=0.31, Ea=(12.393,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['C=C=C=CC(=C)C=C(28050)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['C=CC1=CC(C=C)=C1(28051)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleH_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH2]C=[C]C1C=C1C=C(28018)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(92.3417,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH]=C1C=C=CCC1[CH2](28052)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(188.28,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 40 used for R7_SMMS;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R7_SMMS;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH]=C(C=C)C=C=[C]C(28053)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH]=C([C]=C=CC)C=C(28054)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.0029655,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C([C]=C)C=C=CC(28055)'],
    products = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6H;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=CC(=[CH])C=C=CC(28056)'],
    products = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R7H;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH2]C=C1[CH]C(C=C)=C1(28057)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH]=C(C=C)C1[C]=CC1(28058)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(118.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['[CH]=C1[CH]CCC=C=C1(28059)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.2, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 46 used for R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    products = ['C#C[CH]CC([CH2])C#C(26589)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.213e+11,'s^-1'), n=0.07, Ea=(76.6885,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for 1_4_5_hexatriene
Exact match found for rate rule [1_4_5_hexatriene]
Euclidian distance = 0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(19)', '[CH]=C=CC(=[CH])C=C(28060)'],
    products = ['[CH]=C(C=C)C=[C]C=C(26596)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4568',
    isomers = [
        '[CH]=C(C=C)C=[C]C=C(26596)',
    ],
    reactants = [
        ('CH2CHCCH(26391)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4568',
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

