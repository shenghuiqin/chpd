species(
    label = 'C#C[CH]CC=[C]C=C(26590)',
    structure = SMILES('C#C[CH]CC=[C]C=C'),
    E0 = (626.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(3.10668,'amu*angstrom^2'), symmetry=1, barrier=(71.4287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843074,'amu*angstrom^2'), symmetry=1, barrier=(19.3839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.09628,'amu*angstrom^2'), symmetry=1, barrier=(71.1895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.09656,'amu*angstrom^2'), symmetry=1, barrier=(71.196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.26904,0.0773468,-6.6166e-05,2.27803e-08,3.35569e-13,75523.2,28.3503], Tmin=(100,'K'), Tmax=(870.628,'K')), NASAPolynomial(coeffs=[15.2039,0.0261562,-7.9932e-06,1.22521e-09,-7.63792e-14,72262.2,-45.4257], Tmin=(870.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJC=C) + radical(Sec_Propargyl)"""),
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
    label = 'C#CC1CC=[C]C1[CH2](28289)',
    structure = SMILES('C#CC1CC=[C]C1[CH2]'),
    E0 = (690.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21478,0.045401,2.42989e-05,-6.80867e-08,3.11673e-11,83106.8,26.0564], Tmin=(100,'K'), Tmax=(928.901,'K')), NASAPolynomial(coeffs=[15.5615,0.0241074,-6.69272e-06,1.07631e-09,-7.48538e-14,78694.8,-51.5018], Tmin=(928.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(690.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopentene) + radical(Isobutyl) + radical(cyclopentene-vinyl)"""),
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
    label = 'C#CC=CC=[C]C=C(28333)',
    structure = SMILES('C#CC=CC=[C]C=C'),
    E0 = (573.247,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,2175,525,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.46742,'amu*angstrom^2'), symmetry=1, barrier=(33.7389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46535,'amu*angstrom^2'), symmetry=1, barrier=(33.6912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45773,'amu*angstrom^2'), symmetry=1, barrier=(33.5161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.152956,0.0721106,-4.31333e-05,-7.1722e-09,1.14211e-11,69095.6,25.5647], Tmin=(100,'K'), Tmax=(939.149,'K')), NASAPolynomial(coeffs=[20.3463,0.0171111,-4.81361e-06,7.82054e-10,-5.47634e-14,63935.2,-77.869], Tmin=(939.149,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJC=C)"""),
)

species(
    label = 'C#C[CH]CC#CC=C(28334)',
    structure = SMILES('C#C[CH]CC#CC=C'),
    E0 = (598.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2175,2250,500,525,550,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,331.693,331.693,331.693],'cm^-1')),
        HinderedRotor(inertia=(0.85511,'amu*angstrom^2'), symmetry=1, barrier=(66.7609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855112,'amu*angstrom^2'), symmetry=1, barrier=(66.7609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855112,'amu*angstrom^2'), symmetry=1, barrier=(66.7609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855112,'amu*angstrom^2'), symmetry=1, barrier=(66.7609,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0732787,0.0712476,-6.28478e-05,2.95181e-08,-5.38353e-12,72145.5,30.1919], Tmin=(100,'K'), Tmax=(1505.47,'K')), NASAPolynomial(coeffs=[16.0853,0.0210939,-5.29381e-06,6.73712e-10,-3.60052e-14,68186.8,-50.743], Tmin=(1505.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(598.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH]CC=C=C=C(28335)',
    structure = SMILES('C#C[CH]CC=C=C=C'),
    E0 = (657.171,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,2175,525,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.19834,'amu*angstrom^2'), symmetry=1, barrier=(27.5523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72975,'amu*angstrom^2'), symmetry=1, barrier=(39.7704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73795,'amu*angstrom^2'), symmetry=1, barrier=(39.9589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.379595,0.0774022,-8.10641e-05,4.67561e-08,-1.08389e-11,79171.7,27.4362], Tmin=(100,'K'), Tmax=(1049.45,'K')), NASAPolynomial(coeffs=[13.3874,0.0278213,-1.01953e-05,1.73535e-09,-1.13771e-13,76441.5,-35.9466], Tmin=(1049.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(657.171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
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
    label = 'C#CC[CH]C=[C]C=C(28269)',
    structure = SMILES('C#CCC=C[C]=C[CH2]'),
    E0 = (590.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483784,0.069306,-4.6557e-05,8.99298e-09,2.52479e-12,71205.9,28.2061], Tmin=(100,'K'), Tmax=(977.511,'K')), NASAPolynomial(coeffs=[14.4963,0.0287684,-1.01348e-05,1.73708e-09,-1.16687e-13,67663.6,-43.1838], Tmin=(977.511,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C#C[CH]C[C]=CC=C(28336)',
    structure = SMILES('C#C[CH]C[C]=CC=C'),
    E0 = (665.627,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,1685,370,253.854,253.856,253.88],'cm^-1')),
        HinderedRotor(inertia=(1.89122,'amu*angstrom^2'), symmetry=1, barrier=(86.4843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359554,'amu*angstrom^2'), symmetry=1, barrier=(16.4414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.8911,'amu*angstrom^2'), symmetry=1, barrier=(86.4827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359516,'amu*angstrom^2'), symmetry=1, barrier=(16.4416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.049721,0.0797839,-7.76351e-05,4.03541e-08,-8.26393e-12,80210.2,30.3648], Tmin=(100,'K'), Tmax=(1240.46,'K')), NASAPolynomial(coeffs=[16.6893,0.0236515,-7.15141e-06,1.07275e-09,-6.48833e-14,76223.2,-53.3297], Tmin=(1240.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(665.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH]CC=C[C]=C(28337)',
    structure = SMILES('C#C[CH]CC=C[C]=C'),
    E0 = (626.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(3.10668,'amu*angstrom^2'), symmetry=1, barrier=(71.4287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843074,'amu*angstrom^2'), symmetry=1, barrier=(19.3839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.09628,'amu*angstrom^2'), symmetry=1, barrier=(71.1895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.09656,'amu*angstrom^2'), symmetry=1, barrier=(71.196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.26904,0.0773468,-6.6166e-05,2.27803e-08,3.35569e-13,75523.2,28.3503], Tmin=(100,'K'), Tmax=(870.628,'K')), NASAPolynomial(coeffs=[15.2039,0.0261562,-7.9932e-06,1.22521e-09,-7.63792e-14,72262.2,-45.4257], Tmin=(870.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=CJC=C)"""),
)

species(
    label = 'C#CCC[C]=[C]C=C(28338)',
    structure = SMILES('C#CCCC#C[CH][CH2]'),
    E0 = (682.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2100,2175,2250,500,525,550,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773949,0.0688122,-5.05599e-05,8.00011e-09,6.96533e-12,82227.2,31.2005], Tmin=(100,'K'), Tmax=(788.26,'K')), NASAPolynomial(coeffs=[12.2709,0.028501,-8.15957e-06,1.15687e-09,-6.7032e-14,79854.5,-25.0832], Tmin=(788.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(682.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJ) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH][CH]C=CC=C(28339)',
    structure = SMILES('[CH]=C=CC=CC=C[CH2]'),
    E0 = (522.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.399106,0.0646047,-1.62831e-05,-3.5574e-08,2.18595e-11,63018.9,28.4058], Tmin=(100,'K'), Tmax=(922.36,'K')), NASAPolynomial(coeffs=[19.4806,0.0197339,-4.91448e-06,7.34704e-10,-5.01675e-14,57887.6,-70.8442], Tmin=(922.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC=CC[CH]C#C(28340)',
    structure = SMILES('[CH]=CC=CC[CH]C#C'),
    E0 = (674.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.82026,'amu*angstrom^2'), symmetry=1, barrier=(18.8594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.40303,'amu*angstrom^2'), symmetry=1, barrier=(78.2424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.821233,'amu*angstrom^2'), symmetry=1, barrier=(18.8818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.41097,'amu*angstrom^2'), symmetry=1, barrier=(78.425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.296245,0.0810858,-7.76214e-05,3.87364e-08,-7.49531e-12,81335.6,31.1821], Tmin=(100,'K'), Tmax=(1380.06,'K')), NASAPolynomial(coeffs=[18.4563,0.0207306,-5.49661e-06,7.43436e-10,-4.19032e-14,76731.3,-63.2578], Tmin=(1380.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(674.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_P)"""),
)

species(
    label = '[C]#CCCC=[C]C=C(28341)',
    structure = SMILES('[C]#CCCC=[C]C=C'),
    E0 = (817.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2175,525,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,216.443,216.443,216.443,216.444],'cm^-1')),
        HinderedRotor(inertia=(0.423448,'amu*angstrom^2'), symmetry=1, barrier=(14.0771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.42344,'amu*angstrom^2'), symmetry=1, barrier=(14.0771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11583,'amu*angstrom^2'), symmetry=1, barrier=(70.3397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11585,'amu*angstrom^2'), symmetry=1, barrier=(70.3397,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.305711,0.0798694,-7.94352e-05,4.05669e-08,-6.9374e-12,98482.8,29.1245], Tmin=(100,'K'), Tmax=(843.183,'K')), NASAPolynomial(coeffs=[13.2743,0.0296324,-1.01407e-05,1.65178e-09,-1.05475e-13,95894.6,-33.6081], Tmin=(843.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(817.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJC=C) + radical(Acetyl)"""),
)

species(
    label = 'C#CCCC=[C][C]=C(28342)',
    structure = SMILES('C#CCC[CH]C#C[CH2]'),
    E0 = (632.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2100,2175,2250,500,525,550,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.590081,0.0724789,-6.2378e-05,2.54638e-08,-2.13479e-12,76205,29.5311], Tmin=(100,'K'), Tmax=(855.099,'K')), NASAPolynomial(coeffs=[12.2281,0.0303792,-1.01754e-05,1.64229e-09,-1.04662e-13,73763.5,-27.432], Tmin=(855.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Propargyl)"""),
)

species(
    label = '[C]#C[CH]CC=CC=C(28343)',
    structure = SMILES('[C]#C[CH]CC=CC=C'),
    E0 = (764.929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2175,525,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,331.939,331.949,332.12,332.882],'cm^-1')),
        HinderedRotor(inertia=(0.923106,'amu*angstrom^2'), symmetry=1, barrier=(72.9423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154093,'amu*angstrom^2'), symmetry=1, barrier=(12.1124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151466,'amu*angstrom^2'), symmetry=1, barrier=(12.111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925502,'amu*angstrom^2'), symmetry=1, barrier=(72.8805,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0379229,0.0792813,-7.75321e-05,4.09764e-08,-8.4987e-12,92153.2,30.5887], Tmin=(100,'K'), Tmax=(1268.62,'K')), NASAPolynomial(coeffs=[16.0973,0.0238813,-6.67772e-06,9.30446e-10,-5.3003e-14,88423.5,-49.6585], Tmin=(1268.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(764.929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=C[C]=CCCC#C(28344)',
    structure = SMILES('[CH]C=C=CCCC#C'),
    E0 = (705.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.373454,0.0764555,-6.16203e-05,2.76053e-08,-5.18102e-12,84930.1,29.6746], Tmin=(100,'K'), Tmax=(1246.53,'K')), NASAPolynomial(coeffs=[12.4495,0.0377042,-1.49891e-05,2.66597e-09,-1.7923e-13,81919.4,-31.2466], Tmin=(1246.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(705.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#C[CH]CC=C1[CH]C1(28345)',
    structure = SMILES('C#C[CH]CC=C1[CH]C1'),
    E0 = (661.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852348,0.0568657,-8.22391e-06,-3.37645e-08,1.86734e-11,79735.4,25.7297], Tmin=(100,'K'), Tmax=(937.048,'K')), NASAPolynomial(coeffs=[15.3751,0.0258095,-8.03419e-06,1.3346e-09,-9.11196e-14,75655.4,-50.6375], Tmin=(937.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(661.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Methylene_cyclopropane) + radical(Sec_Propargyl) + radical(Allyl_S)"""),
)

species(
    label = 'C=C[C]=CCC1[C]=C1(28346)',
    structure = SMILES('C=C[C]=CCC1[C]=C1'),
    E0 = (800.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279456,0.0790062,-7.61636e-05,3.97583e-08,-8.40296e-12,96388.9,26.6715], Tmin=(100,'K'), Tmax=(1139.52,'K')), NASAPolynomial(coeffs=[14.2067,0.0301176,-1.18087e-05,2.10765e-09,-1.4267e-13,93214.9,-42.3383], Tmin=(1139.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(800.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(C=CJC=C) + radical(cyclopropenyl-vinyl)"""),
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
    label = 'C#CC1C[CH][C]=CC1(28322)',
    structure = SMILES('C#CC1C[CH][C]=CC1'),
    E0 = (592.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48659,0.0366411,4.91591e-05,-8.92689e-08,3.66889e-11,71353.6,20.7371], Tmin=(100,'K'), Tmax=(959.453,'K')), NASAPolynomial(coeffs=[14.7287,0.0277497,-9.3492e-06,1.69777e-09,-1.23619e-13,66680.8,-53.7094], Tmin=(959.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(592.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclohexene) + radical(cyclohexene-allyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CC=CC=CC=C(28347)',
    structure = SMILES('C#CC=CC=CC=C'),
    E0 = (374.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.228344,0.0656007,-1.41703e-05,-3.98594e-08,2.30981e-11,45163.6,25.6971], Tmin=(100,'K'), Tmax=(949.617,'K')), NASAPolynomial(coeffs=[21.8794,0.0170357,-4.80164e-06,8.41278e-10,-6.3504e-14,39129.3,-87.7587], Tmin=(949.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CCCC=C=C=C(28348)',
    structure = SMILES('C#CCCC=C=C=C'),
    E0 = (510.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.436106,0.0758907,-7.00638e-05,3.54171e-08,-7.31922e-12,61584.8,27.3958], Tmin=(100,'K'), Tmax=(1157.43,'K')), NASAPolynomial(coeffs=[13.2995,0.0314353,-1.24504e-05,2.23227e-09,-1.51386e-13,58607.1,-36.5434], Tmin=(1157.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC([CH2])C=[C]C=C(26593)',
    structure = SMILES('C#CC([CH2])C=[C]C=C'),
    E0 = (652.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.14611,'amu*angstrom^2'), symmetry=1, barrier=(26.3512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.99228,'amu*angstrom^2'), symmetry=1, barrier=(45.8063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13279,'amu*angstrom^2'), symmetry=1, barrier=(26.045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14769,'amu*angstrom^2'), symmetry=1, barrier=(26.3876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.26,'J/mol'), sigma=(6.13635,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=566.26 K, Pc=35.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.106035,0.0879548,-9.62203e-05,5.49266e-08,-1.18564e-11,78637.3,28.4965], Tmin=(100,'K'), Tmax=(897.011,'K')), NASAPolynomial(coeffs=[15.6045,0.0269246,-9.25904e-06,1.51484e-09,-9.70545e-14,75455.7,-47.6144], Tmin=(897.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(652.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJC=C) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC1CC=C1C=C(28249)',
    structure = SMILES('C#CC1CC=C1C=C'),
    E0 = (407.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.476308,0.0607852,-5.91824e-06,-4.41525e-08,2.36362e-11,49204.3,22.3493], Tmin=(100,'K'), Tmax=(952.417,'K')), NASAPolynomial(coeffs=[19.9484,0.0199387,-6.05426e-06,1.06787e-09,-7.85786e-14,43638.7,-80.3891], Tmin=(952.417,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene)"""),
)

species(
    label = 'C3H2(81)',
    structure = SMILES('[CH]C#C'),
    E0 = (530.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.621997,'amu*angstrom^2'), symmetry=1, barrier=(14.3009,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6874,0.0353211,-7.31245e-05,7.03801e-08,-2.45313e-11,63833.3,8.93418], Tmin=(100,'K'), Tmax=(908.454,'K')), NASAPolynomial(coeffs=[4.78002,0.0100611,-4.92178e-06,8.868e-10,-5.66859e-14,64115.2,2.68365], Tmin=(908.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCCH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C=[C]C=C(26613)',
    structure = SMILES('[CH2]C=[C]C=C'),
    E0 = (374.947,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(2.07984,'amu*angstrom^2'), symmetry=1, barrier=(47.8195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07746,'amu*angstrom^2'), symmetry=1, barrier=(47.7649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22509,0.0302876,1.05277e-05,-3.68274e-08,1.72727e-11,45167.7,17.3269], Tmin=(100,'K'), Tmax=(923.516,'K')), NASAPolynomial(coeffs=[10.3879,0.0174192,-5.09531e-06,8.16549e-10,-5.51179e-14,42701,-26.5965], Tmin=(923.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C#CC1CC1[C]=C[CH2](28271)',
    structure = SMILES('C#CC1CC1[C]=C[CH2]'),
    E0 = (704.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.862538,0.0527609,9.61622e-06,-5.6315e-08,2.7414e-11,84900.2,27.2849], Tmin=(100,'K'), Tmax=(940.505,'K')), NASAPolynomial(coeffs=[17.6436,0.0221038,-6.42227e-06,1.08082e-09,-7.72796e-14,79943,-62.2176], Tmin=(940.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(704.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C1[CH]CC=C=CC1(28309)',
    structure = SMILES('[CH]C1=CCC=C=CC1'),
    E0 = (647.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34081,0.0427995,3.41334e-05,-6.59687e-08,2.56101e-11,77956.5,21.8992], Tmin=(100,'K'), Tmax=(1017.93,'K')), NASAPolynomial(coeffs=[11.7699,0.039302,-1.59487e-05,3.00638e-09,-2.14357e-13,73891.3,-38.1398], Tmin=(1017.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#C[CH]CC=C=[C]C(28349)',
    structure = SMILES('C#C[CH]C[CH]C#CC'),
    E0 = (622.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.192348,0.0781314,-8.19394e-05,4.9019e-08,-1.14396e-11,74995.3,29.7233], Tmin=(100,'K'), Tmax=(1198.97,'K')), NASAPolynomial(coeffs=[12.5743,0.0275013,-6.9357e-06,8.30145e-10,-3.9558e-14,72696.1,-29.4655], Tmin=(1198.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(622.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH]C[C]=C=CC(28350)',
    structure = SMILES('C#C[CH]CC#C[CH]C'),
    E0 = (623.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.182905,0.0759148,-7.56173e-05,4.28581e-08,-9.45286e-12,75152.6,30.6502], Tmin=(100,'K'), Tmax=(1275.33,'K')), NASAPolynomial(coeffs=[12.8684,0.0265293,-6.24264e-06,6.91811e-10,-3.0281e-14,72697.5,-30.575], Tmin=(1275.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(623.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH][CH]C=C=CC(28351)',
    structure = SMILES('[CH]=C=C[CH]C=C=CC'),
    E0 = (590.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.651761,0.063364,-2.87251e-05,-1.00002e-08,9.27771e-12,71100.4,28.6625], Tmin=(100,'K'), Tmax=(970.24,'K')), NASAPolynomial(coeffs=[14.9521,0.0277264,-9.67919e-06,1.68363e-09,-1.15432e-13,67227.9,-45.5523], Tmin=(970.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CCJC=C) + radical(C=C=CJ)"""),
)

species(
    label = '[C]#C[CH]CC=C=CC(28352)',
    structure = SMILES('[C]#C[CH]CC=C=CC'),
    E0 = (817.71,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.23147,'amu*angstrom^2'), symmetry=1, barrier=(5.32195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.2951,'amu*angstrom^2'), symmetry=1, barrier=(75.7609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231743,'amu*angstrom^2'), symmetry=1, barrier=(5.32822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.2937,'amu*angstrom^2'), symmetry=1, barrier=(75.7286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.673852,0.0770577,-7.64949e-05,3.56943e-08,-2.14609e-12,98464.2,28.1036], Tmin=(100,'K'), Tmax=(687.659,'K')), NASAPolynomial(coeffs=[10.2044,0.03483,-1.31974e-05,2.26378e-09,-1.48071e-13,96841.2,-16.5777], Tmin=(687.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(817.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC1C[CH]C1=C[CH2](28353)',
    structure = SMILES('C#CC1C[CH]C1=C[CH2]'),
    E0 = (570.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.4648,-0.00303949,0.000116657,-1.17352e-07,2.7137e-11,68223.7,-20.1432], Tmin=(100,'K'), Tmax=(1735.14,'K')), NASAPolynomial(coeffs=[80.4891,0.0242857,-7.01367e-05,1.71094e-08,-1.26909e-12,15509.4,-478.437], Tmin=(1735.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(570.192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(methylenecyclobutane) + radical(Allyl_P) + radical(Allyl_S)"""),
)

species(
    label = 'C#C[CH]CC1[C]=CC1(28321)',
    structure = SMILES('C#C[CH]CC1[C]=CC1'),
    E0 = (733.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.900158,0.056137,-8.17096e-06,-3.31069e-08,1.83468e-11,88397.2,27.0946], Tmin=(100,'K'), Tmax=(935.394,'K')), NASAPolynomial(coeffs=[15.0647,0.0257977,-7.99928e-06,1.32326e-09,-9.00366e-14,84424.7,-47.3649], Tmin=(935.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene) + radical(Sec_Propargyl) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[C]1[CH]CC=C=CCC=1(28354)',
    structure = SMILES('[C]1[CH]CC=C=CCC=1'),
    E0 = (634.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38129,0.011254,0.000114903,-1.49394e-07,5.53618e-11,76368.6,19.9135], Tmin=(100,'K'), Tmax=(975.891,'K')), NASAPolynomial(coeffs=[12.4877,0.0318019,-1.19356e-05,2.32571e-09,-1.75223e-13,71445.1,-43.717], Tmin=(975.891,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(634.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(Cyclooctane) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = 'C#CC=CC=C=CC(28355)',
    structure = SMILES('C#CC=CC=C=CC'),
    E0 = (427.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.192719,0.0728072,-4.95246e-05,8.30931e-09,2.98949e-12,51506.7,25.8564], Tmin=(100,'K'), Tmax=(1028.55,'K')), NASAPolynomial(coeffs=[17.7855,0.02488,-9.51235e-06,1.74405e-09,-1.22657e-13,46803.8,-64.7828], Tmin=(1028.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cdd-CdsCds) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
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
    label = 'C#C[CH]C[CH]C#C(28326)',
    structure = SMILES('C#C[CH]C[CH]C#C'),
    E0 = (664.547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,3000,3050,390,425,1340,1360,335,370,302.995,1819.84],'cm^-1')),
        HinderedRotor(inertia=(1.04658,'amu*angstrom^2'), symmetry=1, barrier=(69.1286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0415,'amu*angstrom^2'), symmetry=1, barrier=(69.0912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04252,'amu*angstrom^2'), symmetry=1, barrier=(69.179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05174,'amu*angstrom^2'), symmetry=1, barrier=(69.1917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.1225,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.659248,0.0693599,-8.04423e-05,5.06247e-08,-1.21856e-11,80050.4,23.9347], Tmin=(100,'K'), Tmax=(1175.95,'K')), NASAPolynomial(coeffs=[12.4762,0.0194204,-4.31177e-06,4.18405e-10,-1.40061e-14,77945,-32.1256], Tmin=(1175.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(664.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH2]C1[C]=CCC=C=C1(28279)',
    structure = SMILES('[CH2]C1[C]=CCC=C=C1'),
    E0 = (724.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15142,0.054436,-2.16257e-05,-2.82243e-09,3.0322e-12,87304.6,21.8676], Tmin=(100,'K'), Tmax=(1170.98,'K')), NASAPolynomial(coeffs=[10.5941,0.0353934,-1.4158e-05,2.56229e-09,-1.74734e-13,84187.3,-29.0462], Tmin=(1170.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(724.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = 'C=C=[C]CC=[C]C=C(28356)',
    structure = SMILES('[CH2]C#CCC=[C]C=C'),
    E0 = (623.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.320682,0.0721497,-5.9833e-05,2.63775e-08,-4.68988e-12,75084.2,29.9107], Tmin=(100,'K'), Tmax=(1347.77,'K')), NASAPolynomial(coeffs=[15.4278,0.0273137,-9.93246e-06,1.6944e-09,-1.11337e-13,71012,-47.4808], Tmin=(1347.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(623.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=CJC=C) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]CC=CC=C(28357)',
    structure = SMILES('[CH]=C=[C]CC=CC=C'),
    E0 = (673.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.932455,'amu*angstrom^2'), symmetry=1, barrier=(21.439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.931213,'amu*angstrom^2'), symmetry=1, barrier=(21.4104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.932753,'amu*angstrom^2'), symmetry=1, barrier=(21.4458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.260258,0.0740038,-5.8591e-05,1.9076e-08,-5.87858e-13,81135.5,30.8409], Tmin=(100,'K'), Tmax=(998.072,'K')), NASAPolynomial(coeffs=[16.2347,0.0255088,-9.04257e-06,1.56636e-09,-1.06111e-13,77173.4,-50.0699], Tmin=(998.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(673.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=CC=[C]C=C(28010)',
    structure = SMILES('[CH2]C=[C]C=CC=C=C'),
    E0 = (567.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.303211,0.0702234,-3.77927e-05,-7.90254e-09,1.04294e-11,68373.3,27.7194], Tmin=(100,'K'), Tmax=(942.661,'K')), NASAPolynomial(coeffs=[17.5212,0.0242367,-7.69896e-06,1.28316e-09,-8.71982e-14,63924.2,-60.7115], Tmin=(942.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C=CC[C]=[C]C=C(28358)',
    structure = SMILES('[CH2][CH]C#CCC=C=C'),
    E0 = (685.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2100,2250,500,550,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,240.987,242.034,242.306],'cm^-1')),
        HinderedRotor(inertia=(0.00287239,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.32107,'amu*angstrom^2'), symmetry=1, barrier=(96.8961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.34206,'amu*angstrom^2'), symmetry=1, barrier=(96.9025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372562,'amu*angstrom^2'), symmetry=1, barrier=(96.8103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579388,0.0691322,-5.76898e-05,2.69857e-08,-5.18069e-12,82573.9,32.0042], Tmin=(100,'K'), Tmax=(1244.08,'K')), NASAPolynomial(coeffs=[12.7134,0.0301184,-1.06503e-05,1.77849e-09,-1.15221e-13,79554.8,-29.1853], Tmin=(1244.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(685.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-(Cds-Cds)CtHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Cdd-CdsCds) + radical(Sec_Propargyl) + radical(RCCJ)"""),
)

species(
    label = 'C=[C][C]=CCC=C=C(28359)',
    structure = SMILES('[CH2]C#C[CH]CC=C=C'),
    E0 = (630.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.509211,0.0711871,-6.07181e-05,2.87241e-08,-5.56222e-12,76000.3,30.588], Tmin=(100,'K'), Tmax=(1234.07,'K')), NASAPolynomial(coeffs=[13.126,0.0302933,-1.10132e-05,1.87323e-09,-1.22862e-13,72886.2,-32.9343], Tmin=(1234.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(630.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C[C]=CCC=C=C(28360)',
    structure = SMILES('[CH]C=C=CCC=C=C'),
    E0 = (704.566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,563.333,586.667,610,1970,2140,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.08614,'amu*angstrom^2'), symmetry=1, barrier=(47.9646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08785,'amu*angstrom^2'), symmetry=1, barrier=(48.0038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0869,'amu*angstrom^2'), symmetry=1, barrier=(47.9818,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.433846,0.0732197,-5.21612e-05,1.92065e-08,-2.92412e-12,84872.1,30.2614], Tmin=(100,'K'), Tmax=(1505.72,'K')), NASAPolynomial(coeffs=[14.5319,0.0357676,-1.48513e-05,2.68731e-09,-1.81376e-13,80626.5,-43.5236], Tmin=(1505.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(704.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C=C=CCC=C=C(28361)',
    structure = SMILES('C=C=C=CCC=C=C'),
    E0 = (510.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.556096,0.0720136,-5.86563e-05,2.49021e-08,-4.32436e-12,61524.2,27.765], Tmin=(100,'K'), Tmax=(1351.67,'K')), NASAPolynomial(coeffs=[14.433,0.0309477,-1.3084e-05,2.42517e-09,-1.67116e-13,57772.7,-43.3649], Tmin=(1351.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]C=[C]C1C=C=CC1(28030)',
    structure = SMILES('[CH2]C=[C]C1C=C=CC1'),
    E0 = (860.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736605,0.0622291,-3.48862e-05,5.20316e-09,1.23395e-12,103668,24.5842], Tmin=(100,'K'), Tmax=(1189.8,'K')), NASAPolynomial(coeffs=[13.8758,0.031457,-1.29855e-05,2.39794e-09,-1.65638e-13,99593.2,-45.0741], Tmin=(1189.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(860.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=[C]CC=C=CC(28362)',
    structure = SMILES('[CH]=C=[C]CC=C=CC'),
    E0 = (726.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,563.333,586.667,610,1970,2140,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.7993,'amu*angstrom^2'), symmetry=1, barrier=(18.3775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800369,'amu*angstrom^2'), symmetry=1, barrier=(18.4021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.801238,'amu*angstrom^2'), symmetry=1, barrier=(18.422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.626943,0.0763578,-7.65235e-05,4.41094e-08,-1.06265e-11,87461.3,29.5648], Tmin=(100,'K'), Tmax=(989.446,'K')), NASAPolynomial(coeffs=[10.625,0.0359386,-1.52474e-05,2.82261e-09,-1.94615e-13,85482.9,-18.5639], Tmin=(989.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(726.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C(C=C)C([CH2])C#C(26594)',
    structure = SMILES('[CH]=C(C=C)C([CH2])C#C'),
    E0 = (699.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2175,525,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.05601,'amu*angstrom^2'), symmetry=1, barrier=(24.2798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05942,'amu*angstrom^2'), symmetry=1, barrier=(24.3582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05457,'amu*angstrom^2'), symmetry=1, barrier=(24.2466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05536,'amu*angstrom^2'), symmetry=1, barrier=(24.2649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3563.49,'J/mol'), sigma=(6.08916,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=556.61 K, Pc=35.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.319136,0.089632,-9.48754e-05,4.97278e-08,-9.46921e-12,84244.6,28.8601], Tmin=(100,'K'), Tmax=(948.213,'K')), NASAPolynomial(coeffs=[18.251,0.0228745,-7.58846e-06,1.2375e-09,-8.01996e-14,80202.3,-62.4874], Tmin=(948.213,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(699.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_P)"""),
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
    E0 = (626.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (642.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (734.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (793.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (822.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (885.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (728.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (753.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (862.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (841.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (807.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (819.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (897.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (964.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (756.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (953.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (749.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (903.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (764.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (852.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (649.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (689.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (719.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (666.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (812.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (634.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (905.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (734.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (696.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (834.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (864.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (702.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (945.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (751.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (745.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (705.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (644.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (1046.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (752.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (758.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (1067.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (671.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (797.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (935.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (895.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (651.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (860.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (835.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (703.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['CH2CHCCH(26391)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['[CH]=C1[CH]CC=C1C=C(28004)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.42e+11,'s^-1'), n=0.258, Ea=(15.8866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;triplebond_intra_H;radadd_intra_cdsingle] for rate rule [R6;triplebond_intra_H;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#CC1CC=[C]C1[CH2](28289)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;doublebond_intra_2H_pri;radadd_intra_csHDe] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_csHCt]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#CC=CC=[C]C=C(28333)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(105.863,'m^3/(mol*s)'), n=1.66278, Ea=(8.08939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#C[CH]CC#CC=C(28334)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.66e+08,'cm^3/(mol*s)'), n=1.64, Ea=(12.2591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2703 used for Ct-Cs_Ct-Cd;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C#C[CH]CC=C=C=C(28335)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]C=C(4699)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0132398,'m^3/(mol*s)'), n=2.333, Ea=(2.89605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CtH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#CC[CH]C=[C]C=C(28269)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.26084e+07,'s^-1'), n=1.66833, Ea=(126.789,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_H/OneDe;Cs_H_out_H/OneDe] + [R2H_S;C_rad_out_H/Ct;Cs_H_out_1H] for rate rule [R2H_S;C_rad_out_H/Ct;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#C[CH]C[C]=CC=C(28336)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#C[CH]CC=C[C]=C(28337)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#CCC[C]=[C]C=C(28338)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#C[CH][CH]C=CC=C(28339)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.56742e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC=CC[CH]C#C(28340)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[C]#CCCC=[C]C=C(28341)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#CCCC=[C][C]=C(28342)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.04875e+10,'s^-1'), n=0.94, Ea=(123.637,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;Cd_rad_out_Cd;Cs_H_out_H/Ct]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[C]#C[CH]CC=CC=C(28343)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.67387e+08,'s^-1'), n=1.29644, Ea=(188.354,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Y_rad_out;Cd_H_out_singleDe] for rate rule [R6HJ_2;Ct_rad_out;Cd_H_out_singleDe]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C[C]=CCCC#C(28344)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_H/Ct]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]C=C(4699)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#C[CH]CC=C1[CH]C1(28345)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C=C[C]=CCC1[C]=C1(28346)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_T;triplebond_intra_H;radadd_intra_cs] for rate rule [R3_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C=CC1C=[C][CH]CC=1(28049)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.926e+10,'s^-1'), n=0.198, Ea=(22.8237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;triplebond_intra_H;radadd_intra_cdsingle] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#CC1C[CH][C]=CC1(28322)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.0689e+08,'s^-1'), n=0.637531, Ea=(63.1312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra_csHDe] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_csHCt]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#CC=CC=CC=C(28347)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#CCCC=C=C=C(28348)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#CC1CC=C1C=C(28249)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriDe_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C3H2(81)', '[CH2]C=[C]C=C(26613)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.26928e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#CC1CC1[C]=C[CH2](28271)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHCt]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['[CH]=C1[CH]CC=C=CC1(28309)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.73116e+09,'s^-1'), n=0.571, Ea=(69.5506,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;triplebond_intra_H;radadd_intra_cs2H] + [R8;multiplebond_intra;radadd_intra_cs2H] for rate rule [R8;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#C[CH]CC=C=[C]C(28349)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#C[CH]C[C]=C=CC(28350)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.0029655,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#C[CH][CH]C=C=CC(28351)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.42353e+06,'s^-1'), n=1.41298, Ea=(75.8094,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;C_rad_out_2H;XH_out] for rate rule [R5H_SMMS;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[C]#C[CH]CC=C=CC(28352)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.87683e+06,'s^-1'), n=1.52702, Ea=(127.963,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;Cs_H_out_2H] + [R8Hall;Y_rad_out;Cs_H_out] for rate rule [R8Hall;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#CC1C[CH]C1=C[CH2](28353)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.78932e+07,'s^-1'), n=1.19089, Ea=(125.213,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_HH_D;doublebond_intra;radadd_intra_cs] for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_csHCt]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#C[CH]CC1[C]=CC1(28321)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(118.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['[C]1[CH]CC=C=CCC=1(28354)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.52039e+09,'s^-1'), n=0.455643, Ea=(78.7028,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;multiplebond_intra;radadd_intra_cs2H] + [R6plus;triplebond_intra_H;radadd_intra] for rate rule [R8_linear;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C#CC=CC=C=CC(28355)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CH2(19)', 'C#C[CH]C[CH]C#C(28326)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['[CH2]C1[C]=CCC=C=C1(28279)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(125.52,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7plus;doublebond_intra_2H_pri;radadd_intra_cdsingleH] for rate rule [R8;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C=C=[C]CC=[C]C=C(28356)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=C=[C]CC=CC=C(28357)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(9.01194e+11,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_double;Cd_H_out_singleDe]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C=[C]C=CC=[C]C=C(28010)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C=CC[C]=[C]C=C(28358)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(33613.8,'s^-1'), n=2.10442, Ea=(111.806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cd_H_out_singleH] for rate rule [R5H;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C=[C][C]=CCC=C=C(28359)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R7HJ_5;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=C[C]=CCC=C=C(28360)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R8Hall;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['C=C=C=CCC=C=C(28361)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['[CH2]C=[C]C1C=C=CC1(28030)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7.4226e+09,'s^-1'), n=0.3735, Ea=(234.124,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 231.6 to 234.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=C=[C]CC=C=CC(28362)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.3367e+07,'s^-1'), n=1.3423, Ea=(109.61,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_double;Cs_H_out_2H] + [R6H;Y_rad_out;Cs_H_out_2H] for rate rule [R6H;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(3.213e+11,'s^-1'), n=0.07, Ea=(76.6885,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for 1_4_5_hexatriene
Exact match found for rate rule [1_4_5_hexatriene]
Euclidian distance = 0
family: 6_membered_central_C-C_shift"""),
)

network(
    label = '4562',
    isomers = [
        'C#C[CH]CC=[C]C=C(26590)',
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
    label = '4562',
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

