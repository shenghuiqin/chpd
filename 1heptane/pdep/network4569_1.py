species(
    label = '[CH]=C(C=C)C(=[CH])C=C(26597)',
    structure = SMILES('[CH]=C(C=C)C(=[CH])C=C'),
    E0 = (729.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3115,3125,620,680,785,800,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,337.663],'cm^-1')),
        HinderedRotor(inertia=(0.355212,'amu*angstrom^2'), symmetry=1, barrier=(28.7257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279382,'amu*angstrom^2'), symmetry=1, barrier=(22.4526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.998459,'amu*angstrom^2'), symmetry=1, barrier=(80.0405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.794749,0.0598441,-2.34502e-05,-1.05305e-08,7.82432e-12,87906.3,25.6889], Tmin=(100,'K'), Tmax=(1032.31,'K')), NASAPolynomial(coeffs=[14.0191,0.0301847,-1.17143e-05,2.14326e-09,-1.4969e-13,84026,-44.1016], Tmin=(1032.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
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
    label = '[CH]=C(C=C)C1=CC1[CH2](27997)',
    structure = SMILES('[CH]=C(C=C)C1=CC1[CH2]'),
    E0 = (814.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.31549,0.0729094,-5.44487e-05,1.47765e-08,1.02819e-12,98082.3,25.5241], Tmin=(100,'K'), Tmax=(977.162,'K')), NASAPolynomial(coeffs=[15.7247,0.0266606,-9.28671e-06,1.58908e-09,-1.06939e-13,94267.4,-52.5712], Tmin=(977.162,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(814.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C(C=C)=CC1[CH2](27998)',
    structure = SMILES('[CH]=C1C(C=C)=CC1[CH2]'),
    E0 = (709.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.823671,0.0565134,-5.78462e-06,-3.63419e-08,1.92961e-11,85424.2,24.029], Tmin=(100,'K'), Tmax=(951.01,'K')), NASAPolynomial(coeffs=[15.9327,0.0255341,-8.2933e-06,1.42859e-09,-9.97137e-14,81077.6,-55.8478], Tmin=(951.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(709.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(Isobutyl)"""),
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
    label = '[CH]=C(C#C)C=C(27493)',
    structure = SMILES('[CH]=C(C#C)C=C'),
    E0 = (570.504,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(1.4333,'amu*angstrom^2'), symmetry=1, barrier=(32.9545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42794,'amu*angstrom^2'), symmetry=1, barrier=(32.8313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44695,0.045192,-1.51416e-05,-2.05694e-08,1.32788e-11,68717.7,19.1078], Tmin=(100,'K'), Tmax=(951.09,'K')), NASAPolynomial(coeffs=[15.9808,0.0108076,-3.08614e-06,5.41859e-10,-4.07907e-14,64743.6,-56.6388], Tmin=(951.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(570.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)Ct) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)C(=C)[C]=C(27999)',
    structure = SMILES('[CH]=C(C=C)C([CH2])=C=C'),
    E0 = (666.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.268788,0.0713297,-4.53682e-05,4.82858e-09,4.08826e-12,80290.4,26.5694], Tmin=(100,'K'), Tmax=(1017.64,'K')), NASAPolynomial(coeffs=[17.0173,0.0263965,-9.9421e-06,1.80135e-09,-1.25752e-13,75799.4,-59.8423], Tmin=(1017.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(666.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C([C]=C)C(=C)C=C(28000)',
    structure = SMILES('[CH]C(=C=C)C(=C)C=C'),
    E0 = (638.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.312772,0.0705457,-3.7816e-05,1.14349e-09,4.01283e-12,76931.5,26.8385], Tmin=(100,'K'), Tmax=(1063.15,'K')), NASAPolynomial(coeffs=[14.7408,0.0348595,-1.37064e-05,2.47959e-09,-1.70657e-13,72812.6,-48.5952], Tmin=(1063.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(638.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC(=C)C(=[CH])C=C(28001)',
    structure = SMILES('[CH]=CC(=C)C(=[CH])C=C'),
    E0 = (729.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.794749,0.0598441,-2.34502e-05,-1.05305e-08,7.82432e-12,87906.3,26.3821], Tmin=(100,'K'), Tmax=(1032.31,'K')), NASAPolynomial(coeffs=[14.0191,0.0301847,-1.17143e-05,2.14326e-09,-1.4969e-13,84026,-43.4085], Tmin=(1032.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(=[CH])C(=C)C=C(28002)',
    structure = SMILES('[CH]=CC(=[CH])C(=C)C=C'),
    E0 = (729.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.794749,0.0598441,-2.34502e-05,-1.05305e-08,7.82432e-12,87906.3,26.3821], Tmin=(100,'K'), Tmax=(1032.31,'K')), NASAPolynomial(coeffs=[14.0191,0.0301847,-1.17143e-05,2.14326e-09,-1.4969e-13,84026,-43.4085], Tmin=(1032.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)C1[CH]CC=1(28003)',
    structure = SMILES('[CH]=C(C=C)C1[CH]CC=1'),
    E0 = (677.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20719,0.0411352,4.35827e-05,-8.85805e-08,3.7424e-11,81565.6,22.3953], Tmin=(100,'K'), Tmax=(961.64,'K')), NASAPolynomial(coeffs=[17.4079,0.023974,-7.9937e-06,1.48892e-09,-1.11616e-13,76127.4,-67.2049], Tmin=(961.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(677.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(cyclobutene-allyl)"""),
)

species(
    label = '[CH]=C1[CH]CC=C1C=C(28004)',
    structure = SMILES('[CH]C1=CCC=C1C=C'),
    E0 = (526.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31475,0.0396466,5.19404e-05,-9.13325e-08,3.63977e-11,63400.8,23.5328], Tmin=(100,'K'), Tmax=(979.86,'K')), NASAPolynomial(coeffs=[14.2509,0.0344484,-1.29849e-05,2.42797e-09,-1.75779e-13,58580.1,-50.2765], Tmin=(979.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(526.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentadiene) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC1=CC=C1C=C(28005)',
    structure = SMILES('C=CC1=CC=C1C=C'),
    E0 = (565.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45976,0.0407616,2.84501e-05,-6.13606e-08,2.49507e-11,68091,21.2638], Tmin=(100,'K'), Tmax=(996.557,'K')), NASAPolynomial(coeffs=[12.7845,0.0315263,-1.21666e-05,2.28148e-09,-1.63969e-13,64035.3,-42.3561], Tmin=(996.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(cyclobutadiene_13)"""),
)

species(
    label = 'C#C[CH]CC[CH]C#C(26588)',
    structure = SMILES('C#C[CH]CC[CH]C#C'),
    E0 = (640.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00747775,0.0843511,-9.48732e-05,5.92416e-08,-1.44726e-11,77214,28.5753], Tmin=(100,'K'), Tmax=(1110.58,'K')), NASAPolynomial(coeffs=[13.9285,0.026959,-7.63376e-06,1.03634e-09,-5.62696e-14,74562.5,-38.1205], Tmin=(1110.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(640.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Sec_Propargyl)"""),
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
    E0 = (729.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (816.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (783.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (752.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (876.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (952.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (774.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1133.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (783.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (903.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (867.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (767.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (738.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (876.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    products = ['CH2CHCCH(26391)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    products = ['[CH]=C(C=C)C1=CC1[CH2](27997)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.08e+10,'s^-1'), n=0.69, Ea=(86.3586,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 14 used for R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH
Exact match found for rate rule [R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    products = ['[CH]=C1C(C=C)=CC1[CH2](27998)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.082e+11,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 41 used for R5_DS_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH
Exact match found for rate rule [R5_DS_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=[C]C=C(4699)', 'CH2CHCCH(26391)'],
    products = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0345031,'m^3/(mol*s)'), n=2.33733, Ea=(26.7584,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-Cd_Ct-H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H3(30)', '[CH]=C(C#C)C=C(27493)'],
    products = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0420833,'m^3/(mol*s)'), n=2.41, Ea=(19.4765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-De_Ct-H;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    products = ['[CH]=C(C=C)C(=C)[C]=C(27999)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(383,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    products = ['[CH]=C([C]=C)C(=C)C=C(28000)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    products = ['[CH]=CC(=C)C(=[CH])C=C(28001)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.55058e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_singleH] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    products = ['[CH]=CC(=[CH])C(=C)C=C(28002)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(769414,'s^-1'), n=1.8337, Ea=(53.4313,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;Cd_rad_out_singleH;XH_out] + [R5H_RSSR;Y_rad_out;Cd_H_out_singleH] for rate rule [R5H_DSSD;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]C=C(4699)', '[CH]=[C]C=C(4699)'],
    products = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.73038e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    products = ['[CH]=C(C=C)C1[CH]CC=1(28003)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.906e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R4_D_D;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    products = ['[CH]=C1[CH]CC=C1C=C(28004)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.58e+09,'s^-1'), n=0.62, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for R5_DS_D;doublebond_intra_pri_2H;radadd_intra_cdsingleH
Exact match found for rate rule [R5_DS_D;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    products = ['C=CC1=CC=C1C=C(28005)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleH_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(C=C)C(=[CH])C=C(26597)'],
    products = ['C#C[CH]CC[CH]C#C(26588)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

network(
    label = '4569',
    isomers = [
        '[CH]=C(C=C)C(=[CH])C=C(26597)',
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
    label = '4569',
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

