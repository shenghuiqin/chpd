species(
    label = 'C#C[CH]CC([O])C=O(29456)',
    structure = SMILES('C#C[CH]CC([O])C=O'),
    E0 = (239.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,414.812,415.46,417.033],'cm^-1')),
        HinderedRotor(inertia=(0.000964841,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0025394,'amu*angstrom^2'), symmetry=1, barrier=(13.55,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.624137,'amu*angstrom^2'), symmetry=1, barrier=(76.3969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.628933,'amu*angstrom^2'), symmetry=1, barrier=(76.3932,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723544,0.0675002,-6.3529e-05,3.07116e-08,-5.42362e-12,28901.1,31.6329], Tmin=(100,'K'), Tmax=(975.151,'K')), NASAPolynomial(coeffs=[13.4857,0.0232898,-8.04314e-06,1.33747e-09,-8.72547e-14,26025.2,-31.6001], Tmin=(975.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=OCOJ)"""),
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
    label = 'C#C[CH]CC1OC1[O](29631)',
    structure = SMILES('C#C[CH]CC1OC1[O]'),
    E0 = (290.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493761,0.0724831,-7.92215e-05,4.87422e-08,-1.16e-11,35105.3,27.5474], Tmin=(100,'K'), Tmax=(1188.42,'K')), NASAPolynomial(coeffs=[11.9676,0.0240959,-5.81845e-06,6.48906e-10,-2.79398e-14,33068,-26.8854], Tmin=(1188.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Ethylene_oxide) + radical(Sec_Propargyl) + radical(CCOJ)"""),
)

species(
    label = 'C#CC1CC([O])C1[O](29632)',
    structure = SMILES('C#CC1CC([O])C1[O]'),
    E0 = (366.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05681,0.0486461,9.53648e-06,-5.45564e-08,2.65947e-11,44243,25.2073], Tmin=(100,'K'), Tmax=(943.287,'K')), NASAPolynomial(coeffs=[17.8754,0.0175306,-4.91546e-06,8.40781e-10,-6.22088e-14,39281.4,-64.4317], Tmin=(943.287,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C1[CH]CC(C=O)O1(29558)',
    structure = SMILES('[CH]C1=CCC(C=O)O1'),
    E0 = (83.7355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51481,0.0311804,7.26192e-05,-1.24168e-07,5.23954e-11,10182.2,24.0663], Tmin=(100,'K'), Tmax=(927.652,'K')), NASAPolynomial(coeffs=[18.4459,0.018287,-3.73381e-06,5.5886e-10,-4.39768e-14,4454.51,-70.2858], Tmin=(927.652,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.7355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(2,3-Dihydrofuran) + radical(AllylJ2_triplet)"""),
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
    label = 'C#CC=CC([O])C=O(29633)',
    structure = SMILES('C#CC=CC([O])C=O'),
    E0 = (201.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2782.5,750,1395,475,1775,1000,750,770,3400,2100,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,406.092,406.104],'cm^-1')),
        HinderedRotor(inertia=(0.167508,'amu*angstrom^2'), symmetry=1, barrier=(19.5956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167435,'amu*angstrom^2'), symmetry=1, barrier=(19.5956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.852283,'amu*angstrom^2'), symmetry=1, barrier=(19.5957,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.7068,0.0625552,-5.37809e-05,2.28866e-08,-3.83142e-12,24416.2,30.3993], Tmin=(100,'K'), Tmax=(1443.98,'K')), NASAPolynomial(coeffs=[16.9959,0.0174323,-6.90723e-06,1.24561e-09,-8.46397e-14,19712,-54.1707], Tmin=(1443.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-OdCsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=OCOJ)"""),
)

species(
    label = 'C#C[CH]CC(=O)C=O(29634)',
    structure = SMILES('[CH]=C=CCC(=O)C=O'),
    E0 = (93.9248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,204.528],'cm^-1')),
        HinderedRotor(inertia=(0.869594,'amu*angstrom^2'), symmetry=1, barrier=(25.6055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.864851,'amu*angstrom^2'), symmetry=1, barrier=(25.6244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.864148,'amu*angstrom^2'), symmetry=1, barrier=(25.6271,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29649,0.0599189,-4.687e-05,1.7838e-08,-2.7709e-12,11392.9,25.2842], Tmin=(100,'K'), Tmax=(1473.18,'K')), NASAPolynomial(coeffs=[13.4852,0.0268241,-1.31727e-05,2.58876e-09,-1.83099e-13,7801.67,-38.2413], Tmin=(1473.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(93.9248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ)"""),
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
    label = 'HCO(16)',
    structure = SMILES('[CH]=O'),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C#C[CH]CC=O(26397)',
    structure = SMILES('[CH]=C=CCC=O'),
    E0 = (201.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(1.07579,'amu*angstrom^2'), symmetry=1, barrier=(24.7345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07601,'amu*angstrom^2'), symmetry=1, barrier=(24.7395,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92723,0.037918,-1.30986e-05,-7.10902e-09,4.24726e-12,24366,20.5999], Tmin=(100,'K'), Tmax=(1147.22,'K')), NASAPolynomial(coeffs=[11.2197,0.0202251,-9.19453e-06,1.79682e-09,-1.28617e-13,21266,-29.7254], Tmin=(1147.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ)"""),
)

species(
    label = 'C#CC[CH]C([O])C=O(29635)',
    structure = SMILES('C#CC[CH]C([O])C=O'),
    E0 = (292.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,180,855.535,3311.6],'cm^-1')),
        HinderedRotor(inertia=(0.0350033,'amu*angstrom^2'), symmetry=1, barrier=(18.2141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0980158,'amu*angstrom^2'), symmetry=1, barrier=(2.25358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792727,'amu*angstrom^2'), symmetry=1, barrier=(18.2263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157492,'amu*angstrom^2'), symmetry=1, barrier=(81.9043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.851198,0.0622344,-4.90523e-05,1.66783e-08,-1.02986e-12,35356.3,33.8334], Tmin=(100,'K'), Tmax=(1021.57,'K')), NASAPolynomial(coeffs=[14.0662,0.0224313,-8.14135e-06,1.42213e-09,-9.64419e-14,32033.2,-33.2529], Tmin=(1021.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ) + radical(CCJCO)"""),
)

species(
    label = 'C#C[CH]C[C](O)C=O(29636)',
    structure = SMILES('C#C[CH]CC(O)=C[O]'),
    E0 = (94.1249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29549,0.0892995,-9.71148e-05,4.9784e-08,-9.30224e-12,11534.7,33.0609], Tmin=(100,'K'), Tmax=(1573.23,'K')), NASAPolynomial(coeffs=[23.5609,0.005437,2.54615e-06,-7.96906e-10,6.22226e-14,6270.91,-89.993], Tmin=(1573.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.1249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=COJ)"""),
)

species(
    label = 'C#CCC[C]([O])C=O(29637)',
    structure = SMILES('C#CCCC([O])=C[O]'),
    E0 = (85.7195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.606859,0.0806681,-8.25689e-05,4.10992e-08,-7.62371e-12,10493.3,32.2537], Tmin=(100,'K'), Tmax=(1534.77,'K')), NASAPolynomial(coeffs=[21.2657,0.0102976,-7.30351e-07,-1.23295e-10,1.52821e-14,5353.49,-77.5104], Tmin=(1534.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.7195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C#C[CH][CH]C(O)C=O(29638)',
    structure = SMILES('[CH]=C=C[CH]C(O)C=O'),
    E0 = (118.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.320358,0.0716916,-5.9263e-05,1.90948e-08,-4.42292e-13,14451.7,31.0605], Tmin=(100,'K'), Tmax=(1014.78,'K')), NASAPolynomial(coeffs=[17.7845,0.0196958,-7.30259e-06,1.31573e-09,-9.18885e-14,10040,-57.7233], Tmin=(1014.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(C=C=CJ)"""),
)

species(
    label = 'C#C[CH]CC(O)[C]=O(29639)',
    structure = SMILES('C#C[CH]CC(O)[C]=O'),
    E0 = (155.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.247302,0.0753351,-8.16954e-05,4.59647e-08,-1.00357e-11,18844.9,33.3779], Tmin=(100,'K'), Tmax=(1183.31,'K')), NASAPolynomial(coeffs=[16.4485,0.0177003,-4.99849e-06,7.05238e-10,-4.07182e-14,15211.6,-46.6615], Tmin=(1183.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCCJ=O) + radical(Sec_Propargyl)"""),
)

species(
    label = '[C]#CCCC([O])C=O(29640)',
    structure = SMILES('[C]#CCCC([O])C=O'),
    E0 = (430.217,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2175,525,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,388.425,388.425,388.425,2369.94],'cm^-1')),
        HinderedRotor(inertia=(0.0833626,'amu*angstrom^2'), symmetry=1, barrier=(8.92461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0833607,'amu*angstrom^2'), symmetry=1, barrier=(8.92468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00111734,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.665914,'amu*angstrom^2'), symmetry=1, barrier=(71.2942,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.849017,0.0689562,-7.29874e-05,4.34601e-08,-1.05152e-11,51856.9,32.0899], Tmin=(100,'K'), Tmax=(1001.48,'K')), NASAPolynomial(coeffs=[11.3046,0.0271955,-1.04389e-05,1.82279e-09,-1.21228e-13,49762.7,-18.3679], Tmin=(1001.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(C=OCOJ)"""),
)

species(
    label = 'C#CCCC([O])[C]=O(29641)',
    structure = SMILES('C#CCCC([O])[C]=O'),
    E0 = (253.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,1855,455,950,1380,1390,370,380,2900,435,2175,525,304.832,304.839,304.846],'cm^-1')),
        HinderedRotor(inertia=(0.260036,'amu*angstrom^2'), symmetry=1, barrier=(17.1455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00181441,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00151007,'amu*angstrom^2'), symmetry=1, barrier=(17.1455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31204,'amu*angstrom^2'), symmetry=1, barrier=(86.5147,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.617405,0.0700178,-7.13306e-05,3.85916e-08,-8.28947e-12,30558.5,32.9945], Tmin=(100,'K'), Tmax=(1135.73,'K')), NASAPolynomial(coeffs=[14.264,0.0219551,-7.85239e-06,1.33028e-09,-8.74008e-14,27458.7,-34.5792], Tmin=(1135.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[C]#C[CH]CC(O)C=O(29642)',
    structure = SMILES('[C]#C[CH]CC(O)C=O'),
    E0 = (332.694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.653292,0.0721287,-7.53845e-05,3.97139e-08,-7.1105e-12,40135.9,31.8532], Tmin=(100,'K'), Tmax=(848.935,'K')), NASAPolynomial(coeffs=[13.2232,0.0234072,-7.85972e-06,1.26357e-09,-8.00622e-14,37623.1,-28.9604], Tmin=(848.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl)"""),
)

species(
    label = '[O]C(C=O)CC1[C]=C1(29643)',
    structure = SMILES('[O]C(C=O)CC1[C]=C1'),
    E0 = (412.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02511,0.0656137,-6.05405e-05,2.9879e-08,-6.06981e-12,49754.3,28.9166], Tmin=(100,'K'), Tmax=(1165.22,'K')), NASAPolynomial(coeffs=[11.9104,0.0282462,-1.24371e-05,2.35731e-09,-1.64986e-13,47217.6,-25.2634], Tmin=(1165.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(C=OCOJ) + radical(cyclopropenyl-vinyl)"""),
)

species(
    label = 'C#C[CH]CC1[CH]OO1(29644)',
    structure = SMILES('C#C[CH]CC1[CH]OO1'),
    E0 = (504.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.786547,0.0727395,-8.23401e-05,5.49282e-08,-1.51197e-11,60739.2,23.6301], Tmin=(100,'K'), Tmax=(879.14,'K')), NASAPolynomial(coeffs=[9.72526,0.032068,-1.29437e-05,2.3021e-09,-1.5399e-13,59167.6,-18.3423], Tmin=(879.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(12dioxetane) + radical(Sec_Propargyl) + radical(CCsJOO)"""),
)

species(
    label = 'C#CC1CC([O])[CH]O1(29526)',
    structure = SMILES('C#CC1CC([O])[CH]O1'),
    E0 = (272.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.994898,0.0481663,2.30713e-05,-8.42637e-08,4.35053e-11,32956.9,22.3729], Tmin=(100,'K'), Tmax=(870.723,'K')), NASAPolynomial(coeffs=[20.9214,0.00818008,3.14389e-06,-1.00774e-09,7.74596e-14,27532.5,-82.2247], Tmin=(870.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Tetrahydrofuran) + radical(CCsJOCs) + radical(CC(C)OJ)"""),
)

species(
    label = '[O][CH]C1CC=C=CO1(29493)',
    structure = SMILES('[O][CH]C1CC=C=CO1'),
    E0 = (233.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59086,0.129459,-0.000192323,1.44991e-07,-4.28841e-11,28265.1,14.4783], Tmin=(100,'K'), Tmax=(832.77,'K')), NASAPolynomial(coeffs=[18.9415,0.0308354,-1.46794e-05,2.77882e-09,-1.91144e-13,24845.4,-80.821], Tmin=(832.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'C#CC=CC(O)C=O(29645)',
    structure = SMILES('C#CC=CC(O)C=O'),
    E0 = (-41.7796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.649932,0.0634781,-4.37697e-05,7.81997e-09,2.2117e-12,-4895.53,29.8121], Tmin=(100,'K'), Tmax=(1058.44,'K')), NASAPolynomial(coeffs=[16.8032,0.0204505,-8.32628e-06,1.57868e-09,-1.12886e-13,-9324.24,-53.8031], Tmin=(1058.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.7796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-OdCsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CCCC(=O)C=O(29646)',
    structure = SMILES('C#CCCC(=O)C=O'),
    E0 = (-64.9899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04985,0.0716323,-7.73791e-05,3.9848e-08,-2.94239e-12,-7716.52,25.2688], Tmin=(100,'K'), Tmax=(609.63,'K')), NASAPolynomial(coeffs=[8.25308,0.0355665,-1.61899e-05,3.06241e-09,-2.12372e-13,-8802.86,-7.62425], Tmin=(609.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.9899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + group(Ct-CtCs) + group(Ct-CtH)"""),
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
    label = 'C#C[CH]CC[O](26478)',
    structure = SMILES('C#C[CH]CC[O]'),
    E0 = (337.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2175,525,750,770,3400,2100,3025,407.5,1350,352.5,180,1836.18,1840.93],'cm^-1')),
        HinderedRotor(inertia=(2.82432,'amu*angstrom^2'), symmetry=1, barrier=(64.9367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.82373,'amu*angstrom^2'), symmetry=1, barrier=(64.9231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226791,'amu*angstrom^2'), symmetry=1, barrier=(5.21437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3643.1,'J/mol'), sigma=(6.20007,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.04 K, Pc=34.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63858,0.0561205,-7.21531e-05,5.87321e-08,-1.92566e-11,40719.7,22.0443], Tmin=(100,'K'), Tmax=(895.404,'K')), NASAPolynomial(coeffs=[5.42687,0.0297023,-1.19906e-05,2.09572e-09,-1.36966e-13,40421.9,6.312], Tmin=(895.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CCOJ)"""),
)

species(
    label = 'C#CC([CH2])C([O])C=O(29457)',
    structure = SMILES('C#CC([CH2])C([O])C=O'),
    E0 = (289.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2782.5,750,1395,475,1775,1000,750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,482.678,487.449],'cm^-1')),
        HinderedRotor(inertia=(0.54176,'amu*angstrom^2'), symmetry=1, barrier=(12.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00035764,'amu*angstrom^2'), symmetry=1, barrier=(1.71898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.481135,'amu*angstrom^2'), symmetry=1, barrier=(76.4926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.47078,'amu*angstrom^2'), symmetry=1, barrier=(76.5333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.527693,0.0684877,-6.66024e-05,3.45847e-08,-7.08117e-12,34894,33.9991], Tmin=(100,'K'), Tmax=(1234.58,'K')), NASAPolynomial(coeffs=[14.8807,0.0203933,-6.235e-06,9.42737e-10,-5.73343e-14,31471.3,-37.7793], Tmin=(1234.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC1CC(C=O)O1(29464)',
    structure = SMILES('C#CC1CC(C=O)O1'),
    E0 = (20.2394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49021,0.0346665,5.47604e-05,-1.10589e-07,5.06067e-11,2544.02,25.3224], Tmin=(100,'K'), Tmax=(890.104,'K')), NASAPolynomial(coeffs=[19.1387,0.0108414,1.40774e-06,-5.97554e-10,4.4514e-14,-2795.77,-70.1137], Tmin=(890.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.2394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane)"""),
)

species(
    label = 'C#C[CH]CO[CH]C=O(29460)',
    structure = SMILES('[CH]=C=CCOC=C[O]'),
    E0 = (151.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.27601,'amu*angstrom^2'), symmetry=1, barrier=(29.3379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27732,'amu*angstrom^2'), symmetry=1, barrier=(29.3681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27675,'amu*angstrom^2'), symmetry=1, barrier=(29.355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4082.95,'J/mol'), sigma=(6.5735,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=637.75 K, Pc=32.62 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0113538,0.0692281,-2.33584e-05,-4.2182e-08,2.8141e-11,18415.5,28.8598], Tmin=(100,'K'), Tmax=(908.122,'K')), NASAPolynomial(coeffs=[26.2347,0.00364152,2.52032e-06,-6.48694e-10,4.3257e-14,11594.4,-106.459], Tmin=(908.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=COJ)"""),
)

species(
    label = 'C#C[CH]CC([O])=CO(29647)',
    structure = SMILES('C#C[CH]CC([O])=CO'),
    E0 = (90.4671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2175,525,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,213.749,215.12,215.334],'cm^-1')),
        HinderedRotor(inertia=(0.627434,'amu*angstrom^2'), symmetry=1, barrier=(20.6405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.51233,'amu*angstrom^2'), symmetry=1, barrier=(82.1364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.631651,'amu*angstrom^2'), symmetry=1, barrier=(20.6618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.631741,'amu*angstrom^2'), symmetry=1, barrier=(20.6502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18092,0.0906059,-0.000103004,5.53427e-08,-1.08302e-11,11087.3,32.7274], Tmin=(100,'K'), Tmax=(1499.81,'K')), NASAPolynomial(coeffs=[23.1692,0.00511131,3.05658e-06,-9.38229e-10,7.41324e-14,6094.83,-86.9115], Tmin=(1499.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.4671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C([O])C=O(2583)',
    structure = SMILES('[CH2]C([O])C=O'),
    E0 = (82.1174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.347137,'amu*angstrom^2'), symmetry=1, barrier=(18.8168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0159131,'amu*angstrom^2'), symmetry=1, barrier=(18.8156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37294,0.0300135,-1.22984e-05,-6.17239e-09,4.64176e-12,9940.08,21.4594], Tmin=(100,'K'), Tmax=(1011.88,'K')), NASAPolynomial(coeffs=[9.84794,0.0127268,-4.85036e-06,8.96671e-10,-6.36509e-14,7799.56,-17.7935], Tmin=(1011.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.1174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CJCO)"""),
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
    label = 'C#C[CH]C[CH]C=O(27376)',
    structure = SMILES('C#C[CH]CC=C[O]'),
    E0 = (308.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,288.093,288.102,288.117],'cm^-1')),
        HinderedRotor(inertia=(0.583639,'amu*angstrom^2'), symmetry=1, barrier=(34.3834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.474064,'amu*angstrom^2'), symmetry=1, barrier=(27.9184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1292,'amu*angstrom^2'), symmetry=1, barrier=(66.4958,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0647921,0.0677002,-6.37812e-05,3.0016e-08,-5.2724e-12,37291.1,28.5977], Tmin=(100,'K'), Tmax=(1651.61,'K')), NASAPolynomial(coeffs=[17.2967,0.0123796,-1.48396e-06,4.11297e-12,6.9518e-15,33366.6,-58.3917], Tmin=(1651.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=COJ)"""),
)

species(
    label = 'C#CC1CC([CH][O])O1(29513)',
    structure = SMILES('C#CC1CC([CH][O])O1'),
    E0 = (346.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.331361,0.069596,-6.68626e-05,3.5165e-08,-7.14138e-12,41831,28.0954], Tmin=(100,'K'), Tmax=(1380.58,'K')), NASAPolynomial(coeffs=[13.9822,0.0211929,-4.6547e-06,4.81187e-10,-1.97156e-14,38905.4,-39.1091], Tmin=(1380.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.626,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O]C1C=C=CCC1[O](29581)',
    structure = SMILES('[O]C1C=C=CCC1[O]'),
    E0 = (284.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.760904,0.110313,-0.000143117,9.68426e-08,-2.62052e-11,34357.9,14.6815], Tmin=(100,'K'), Tmax=(900.895,'K')), NASAPolynomial(coeffs=[16.3498,0.0343404,-1.66206e-05,3.23366e-09,-2.2835e-13,31274.9,-66.082], Tmin=(900.895,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C=[C]CC([O])C=O(29648)',
    structure = SMILES('[CH2]C#CCC([O])C=O'),
    E0 = (231.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.8952,0.0607262,-4.92132e-05,2.10284e-08,-3.62715e-12,27910.5,33.4181], Tmin=(100,'K'), Tmax=(1384.02,'K')), NASAPolynomial(coeffs=[13.8411,0.0233108,-8.66257e-06,1.4957e-09,-9.88974e-14,24327,-33.2459], Tmin=(1384.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=OCOJ) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]CC(O)C=O(29649)',
    structure = SMILES('[CH]=C=[C]CC(O)C=O'),
    E0 = (239.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.642294,0.0707161,-7.05825e-05,3.70035e-08,-7.75419e-12,28978.6,33.1804], Tmin=(100,'K'), Tmax=(1155.74,'K')), NASAPolynomial(coeffs=[14.271,0.023547,-9.3626e-06,1.68965e-09,-1.15324e-13,25828.4,-34.5427], Tmin=(1155.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C]C=CC([O])C=O(29567)',
    structure = SMILES('C=C=C[CH]C([O])C=O'),
    E0 = (208.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.410395,0.0691934,-5.91236e-05,2.50845e-08,-4.20575e-12,25183.5,30.9039], Tmin=(100,'K'), Tmax=(1435.83,'K')), NASAPolynomial(coeffs=[17.8413,0.0206334,-8.39302e-06,1.52972e-09,-1.04478e-13,20178,-59.4954], Tmin=(1435.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C=CC[C]([O])C=O(29650)',
    structure = SMILES('C=C=CCC([O])=C[O]'),
    E0 = (85.2447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.358143,0.0669658,-3.74789e-05,-1.20911e-08,1.31299e-11,10395.7,29.5825], Tmin=(100,'K'), Tmax=(940.196,'K')), NASAPolynomial(coeffs=[20.7508,0.0128403,-3.19031e-06,5.12883e-10,-3.78769e-14,5118.77,-75.213], Tmin=(940.196,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.2447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C=CCC([O])[C]=O(29651)',
    structure = SMILES('C=C=CCC([O])[C]=O'),
    E0 = (251.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1855,455,950,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,308.665,308.668,308.67],'cm^-1')),
        HinderedRotor(inertia=(0.00176934,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151833,'amu*angstrom^2'), symmetry=1, barrier=(10.2654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151836,'amu*angstrom^2'), symmetry=1, barrier=(10.2654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.75248,0.0660439,-5.96059e-05,2.77117e-08,-5.15116e-12,30344.5,33.2854], Tmin=(100,'K'), Tmax=(1293.36,'K')), NASAPolynomial(coeffs=[14.8497,0.0224444,-9.0398e-06,1.64685e-09,-1.1289e-13,26698,-38.352], Tmin=(1293.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1[CH]OC=C=CC1(29500)',
    structure = SMILES('[O]C1[CH]OC=C=CC1'),
    E0 = (303.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.435897,0.0592969,-3.67295e-06,-5.14755e-08,2.7584e-11,36591.8,19.828], Tmin=(100,'K'), Tmax=(947.519,'K')), NASAPolynomial(coeffs=[23.0611,0.0111599,-2.46882e-06,4.47164e-10,-3.87801e-14,30177.5,-99.3284], Tmin=(947.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CCsJOC(O)) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C=CCC(=O)C=O(29652)',
    structure = SMILES('C=C=CCC(=O)C=O'),
    E0 = (-60.5519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18836,0.0597749,-4.06391e-05,1.23926e-08,-1.48255e-12,-7180.59,25.0629], Tmin=(100,'K'), Tmax=(1912.83,'K')), NASAPolynomial(coeffs=[18.9934,0.0225416,-1.14412e-05,2.21628e-09,-1.52537e-13,-13992.1,-72.384], Tmin=(1912.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.5519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    E0 = (239.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (341.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (366.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (365.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (416.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (338.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (258.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (265.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (328.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (418.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (393.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (440.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (314.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (314.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (577.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (355.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (407.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (432.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (465.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (505.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (395.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (270.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (317.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (328.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (523.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (448.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (247.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (465.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (239.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (612.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (551.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (361.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (382.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (371.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (384.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (283.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (272.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (383.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (330.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (264.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['OCHCHO(48)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C#C[CH]CC1OC1[O](29631)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C#CC1CC([O])C1[O](29632)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(127.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHCt]
Euclidian distance = 3.0
family: Intra_R_Add_Exocyclic
Ea raised from 125.5 to 127.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['[CH]=C1[CH]CC(C=O)O1(29558)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.03297e+11,'s^-1'), n=0.45637, Ea=(125.756,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;triplebond_intra_H;radadd_intra] for rate rule [R6;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#CC=CC([O])C=O(29633)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(2.76144,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2560 used for Cds-CsH_Cds-CtH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CtH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C#C[CH]CC(=O)C=O(29634)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(25.1243,'m^3/(mol*s)'), n=1.86, Ea=(32.426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-DeNd_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C=C[O](2537)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0264797,'m^3/(mol*s)'), n=2.333, Ea=(2.89605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CtH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OCHCHO(48)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(210.531,'m^3/(mol*s)'), n=1.35673, Ea=(38.9443,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-DeH_O;YJ] for rate rule [CO-DeH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['HCO(16)', 'C#C[CH]CC=O(26397)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;CO_pri_rad] for rate rule [CO-CsH_O;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#CC[CH]C([O])C=O(29635)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C#C[CH]C[C](O)C=O(29636)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C#CCC[C]([O])C=O(29637)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.07201e+08,'s^-1'), n=1.25033, Ea=(201.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/OneDe;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_H/Ct;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C#C[CH][CH]C(O)C=O(29638)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C#C[CH]CC(O)[C]=O(29639)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[C]#CCCC([O])C=O(29640)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#CCCC([O])[C]=O(29641)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(363473,'s^-1'), n=1.92229, Ea=(102.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;CO_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[C]#C[CH]CC(O)C=O(29642)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(19101.7,'s^-1'), n=1.79759, Ea=(74.4436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;XH_out] for rate rule [R6HJ_2;Ct_rad_out;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C=C[O](2537)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['[O]C(C=O)CC1[C]=C1(29643)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_T;triplebond_intra_H;radadd_intra_cs] for rate rule [R3_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C#C[CH]CC1[CH]OO1(29644)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(265.786,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C#CC1CC([O])[CH]O1(29526)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.8435e+15,'s^-1'), n=-1.17677, Ea=(156.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHDe] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_csHCt]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['[O][CH]C1CC=C=CO1(29493)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.23539e+09,'s^-1'), n=0.445806, Ea=(31.0324,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;triplebond_intra_H;radadd_intra] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C#CC=CC(O)C=O(29645)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C#CCCC(=O)C=O(29646)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CO(12)', 'C#C[CH]CC[O](26478)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C#CC1CC(C=O)O1(29464)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#C[CH]CC([O])=CO(29647)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(148.816,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol
Ea raised from 146.8 to 148.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C3H2(81)', '[CH2]C([O])C=O(2583)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['O(4)', 'C#C[CH]C[CH]C=O(27376)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C#CC1CC([CH][O])O1(29513)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['[O]C1C=C=CCC1[O](29581)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.73685e+10,'s^-1'), n=0.191667, Ea=(142.744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R7_MMSR;carbonylbond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C=C=[C]CC([O])C=O(29648)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['[CH]=C=[C]CC(O)C=O(29649)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C=[C]C=CC([O])C=O(29567)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C=C=CC[C]([O])C=O(29650)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C=CCC([O])[C]=O(29651)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(101399,'s^-1'), n=2.02226, Ea=(132.65,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;Cd_H_out_singleH] + [R6H;Y_rad_out;XH_out] for rate rule [R6H;CO_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['[O]C1[CH]OC=C=CC1(29500)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.41957e+11,'s^-1'), n=0.2055, Ea=(90.797,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_linear;multiplebond_intra;radadd_intra_cdsingleH] + [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C#C[CH]CC([O])C=O(29456)'],
    products = ['C=C=CCC(=O)C=O(29652)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

network(
    label = '4674',
    isomers = [
        'C#C[CH]CC([O])C=O(29456)',
    ],
    reactants = [
        ('OCHCHO(48)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4674',
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

