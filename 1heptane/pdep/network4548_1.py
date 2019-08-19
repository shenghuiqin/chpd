species(
    label = 'C#C[CH]CC(=C)[CH]C(26576)',
    structure = SMILES('C#C[CH]CC([CH2])=CC'),
    E0 = (452.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,280.958,4000],'cm^-1')),
        HinderedRotor(inertia=(0.31491,'amu*angstrom^2'), symmetry=1, barrier=(17.4021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.69638,'amu*angstrom^2'), symmetry=1, barrier=(93.7549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.69715,'amu*angstrom^2'), symmetry=1, barrier=(93.852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.69066,'amu*angstrom^2'), symmetry=1, barrier=(93.8213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.68589,'amu*angstrom^2'), symmetry=1, barrier=(93.8855,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0773253,0.0821768,-7.17924e-05,3.45969e-08,-6.76576e-12,54568.1,30.1996], Tmin=(100,'K'), Tmax=(1230.35,'K')), NASAPolynomial(coeffs=[15.1429,0.032694,-1.14643e-05,1.90791e-09,-1.23502e-13,50822.9,-46.3841], Tmin=(1230.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(Sec_Propargyl)"""),
)

species(
    label = 'CH3CHCCH2(18175)',
    structure = SMILES('C=C=CC'),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.759584,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74635,0.0218189,8.22353e-06,-2.14768e-08,8.55624e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.6,'K')), NASAPolynomial(coeffs=[6.82078,0.0192338,-7.45622e-06,1.36536e-09,-9.53195e-14,16028,-10.4333], Tmin=(1025.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C#CC1CC1([CH2])[CH]C(28161)',
    structure = SMILES('C#CC1CC1([CH2])[CH]C'),
    E0 = (586.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579915,0.0597992,-9.59294e-07,-4.24954e-08,2.08555e-11,70661,29.7444], Tmin=(100,'K'), Tmax=(979.017,'K')), NASAPolynomial(coeffs=[16.7225,0.0300351,-1.08047e-05,1.96693e-09,-1.40153e-13,65765.9,-56.6494], Tmin=(979.017,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(Neopentyl) + radical(Cs_S)"""),
)

species(
    label = '[CH]=C1[CH]CC(=CC)C1(27981)',
    structure = SMILES('[CH]C1=CCC(=CC)C1'),
    E0 = (410.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08528,0.0482722,3.06403e-05,-6.3622e-08,2.47034e-11,49547.5,28.1394], Tmin=(100,'K'), Tmax=(1027.01,'K')), NASAPolynomial(coeffs=[11.7778,0.0440546,-1.78651e-05,3.34989e-09,-2.37376e-13,45377.4,-33.3403], Tmin=(1027.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(4-Methylenecyclopentene) + radical(AllylJ2_triplet)"""),
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
    label = 'C#CC=CC([CH2])=CC(28162)',
    structure = SMILES('C#CC=CC([CH2])=CC'),
    E0 = (400.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2175,525,750,770,3400,2100,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.2085,'amu*angstrom^2'), symmetry=1, barrier=(27.7858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20844,'amu*angstrom^2'), symmetry=1, barrier=(27.7843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20924,'amu*angstrom^2'), symmetry=1, barrier=(27.8029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20905,'amu*angstrom^2'), symmetry=1, barrier=(27.7985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00542329,0.0745271,-3.99925e-05,-7.22065e-09,9.62764e-12,48308.1,26.5494], Tmin=(100,'K'), Tmax=(989.334,'K')), NASAPolynomial(coeffs=[19.4153,0.0251069,-9.11723e-06,1.66001e-09,-1.17998e-13,43045.5,-74.0704], Tmin=(989.334,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C][CH]C(18176)',
    structure = SMILES('[CH2][C]=CC'),
    E0 = (361.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.352622,'amu*angstrom^2'), symmetry=1, barrier=(8.10748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828631,'amu*angstrom^2'), symmetry=1, barrier=(19.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12013e-13,43485.7,14.8304], Tmin=(100,'K'), Tmax=(2065.83,'K')), NASAPolynomial(coeffs=[10.7464,0.014324,-5.20136e-06,8.69079e-10,-5.48385e-14,40045.6,-31.3799], Tmin=(2065.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Allyl_P)"""),
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
    label = 'C#CC[CH]C([CH2])=CC(28163)',
    structure = SMILES('C#CC[CH]C([CH2])=CC'),
    E0 = (447.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0596728,0.0780646,-6.02488e-05,2.46605e-08,-4.08945e-12,53957.4,29.2035], Tmin=(100,'K'), Tmax=(1434.13,'K')), NASAPolynomial(coeffs=[16.7108,0.0312897,-1.13261e-05,1.91866e-09,-1.25094e-13,49147.2,-57.7514], Tmin=(1434.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'C#C[CH][CH]C(C)=CC(28164)',
    structure = SMILES('[CH]=C=C[CH]C(C)=CC'),
    E0 = (410.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.422372,0.0680663,-2.93979e-05,-1.01148e-08,9.07127e-12,49503.9,29.3827], Tmin=(100,'K'), Tmax=(983.316,'K')), NASAPolynomial(coeffs=[14.6619,0.0333835,-1.19448e-05,2.08938e-09,-1.42731e-13,45579.9,-44.7886], Tmin=(983.316,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'C#C[CH]CC(C)=[C]C(28165)',
    structure = SMILES('C#C[CH]CC(C)=[C]C'),
    E0 = (538.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,231.072,2589.07],'cm^-1')),
        HinderedRotor(inertia=(0.303747,'amu*angstrom^2'), symmetry=1, barrier=(11.5091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303748,'amu*angstrom^2'), symmetry=1, barrier=(11.5091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303748,'amu*angstrom^2'), symmetry=1, barrier=(11.5091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.90699,'amu*angstrom^2'), symmetry=1, barrier=(72.2559,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.90697,'amu*angstrom^2'), symmetry=1, barrier=(72.256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.220367,0.085302,-9.10043e-05,5.87218e-08,-1.58321e-11,64933.8,29.3309], Tmin=(100,'K'), Tmax=(892.676,'K')), NASAPolynomial(coeffs=[10.1527,0.0407971,-1.62227e-05,2.87502e-09,-1.92217e-13,63160.4,-17.4596], Tmin=(892.676,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Sec_Propargyl)"""),
)

species(
    label = '[C]#CCCC([CH2])=CC(28166)',
    structure = SMILES('[C]#CCCC([CH2])=CC'),
    E0 = (643.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.209518,0.0817651,-7.49063e-05,3.9419e-08,-8.61659e-12,77516.9,30.0758], Tmin=(100,'K'), Tmax=(1091.08,'K')), NASAPolynomial(coeffs=[12.1771,0.03789,-1.45862e-05,2.56162e-09,-1.71249e-13,74905.4,-28.7036], Tmin=(1091.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(643.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(Acetyl)"""),
)

species(
    label = 'C#CCCC([CH2])=[C]C(28167)',
    structure = SMILES('C#CCCC([CH2])=[C]C'),
    E0 = (544.064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,250.223,2831.6],'cm^-1')),
        HinderedRotor(inertia=(0.364178,'amu*angstrom^2'), symmetry=1, barrier=(16.1837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374863,'amu*angstrom^2'), symmetry=1, barrier=(16.2652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76993,'amu*angstrom^2'), symmetry=1, barrier=(77.332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361866,'amu*angstrom^2'), symmetry=1, barrier=(16.2059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78494,'amu*angstrom^2'), symmetry=1, barrier=(77.319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.224699,0.0819724,-7.40829e-05,3.77435e-08,-7.99362e-12,65572.6,29.7535], Tmin=(100,'K'), Tmax=(1120.63,'K')), NASAPolynomial(coeffs=[12.5367,0.0380259,-1.52591e-05,2.74902e-09,-1.86749e-13,62813.1,-31.0468], Tmin=(1120.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(544.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C#C[CH]C[C](C)C=C(27577)',
    structure = SMILES('C#C[CH]C[C](C)C=C'),
    E0 = (452.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,360,370,350,2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.455224,0.069971,-3.84041e-05,-2.18379e-10,5.97447e-12,54553.4,27.9418], Tmin=(100,'K'), Tmax=(947.544,'K')), NASAPolynomial(coeffs=[13.3576,0.0343147,-1.17365e-05,1.96979e-09,-1.30516e-13,51263.9,-38.0656], Tmin=(947.544,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_T) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CCC[C]([CH2])C=C(27579)',
    structure = SMILES('C#CCCC([CH2])=C[CH2]'),
    E0 = (457.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0322804,'amu*angstrom^2'), symmetry=1, barrier=(14.8879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42833,'amu*angstrom^2'), symmetry=1, barrier=(32.8401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.65373,'amu*angstrom^2'), symmetry=1, barrier=(84.0065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0322825,'amu*angstrom^2'), symmetry=1, barrier=(14.8878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182152,'amu*angstrom^2'), symmetry=1, barrier=(83.9957,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0021803,0.0779156,-5.14823e-05,9.20494e-09,2.89897e-12,55203.8,30.3564], Tmin=(100,'K'), Tmax=(1003.52,'K')), NASAPolynomial(coeffs=[16.5062,0.031646,-1.14907e-05,2.01525e-09,-1.37362e-13,50908.8,-54.2198], Tmin=(1003.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[C]#C[CH]CC(C)=CC(28168)',
    structure = SMILES('[C]#C[CH]CC(C)=CC'),
    E0 = (638.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.484886,0.0812841,-7.55535e-05,3.36882e-08,-1.71803e-12,76866.1,28.6792], Tmin=(100,'K'), Tmax=(691.981,'K')), NASAPolynomial(coeffs=[9.68665,0.0408813,-1.56935e-05,2.72453e-09,-1.79969e-13,75286.4,-14.5382], Tmin=(691.981,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(638.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH]C[C]1CC1C(28169)',
    structure = SMILES('C#C[CH]C[C]1CC1C'),
    E0 = (530.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.607002,0.0670053,-3.69593e-05,3.31475e-09,3.31764e-12,63946.1,29.0427], Tmin=(100,'K'), Tmax=(998.557,'K')), NASAPolynomial(coeffs=[12.0294,0.0361399,-1.2962e-05,2.2267e-09,-1.48673e-13,60922.5,-29.7644], Tmin=(998.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(Tertalkyl) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH2]C(=CC)CC1[C]=C1(28170)',
    structure = SMILES('[CH2]C(=CC)CC1[C]=C1'),
    E0 = (625.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.339159,0.0789872,-6.44584e-05,2.83788e-08,-5.20538e-12,75416.1,27.0672], Tmin=(100,'K'), Tmax=(1270.19,'K')), NASAPolynomial(coeffs=[13.3556,0.0379966,-1.60518e-05,2.97234e-09,-2.04862e-13,72109.5,-38.8429], Tmin=(1270.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(625.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(Allyl_P) + radical(cyclopropenyl-vinyl)"""),
)

species(
    label = 'C#CC1C[C]([CH2])C1C(28171)',
    structure = SMILES('C#CC1C[C]([CH2])C1C'),
    E0 = (575.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.959326,0.0558984,-3.0174e-06,-3.24788e-08,1.61278e-11,69388.2,28.7633], Tmin=(100,'K'), Tmax=(962.861,'K')), NASAPolynomial(coeffs=[11.9939,0.0360465,-1.25779e-05,2.17312e-09,-1.47735e-13,66058.6,-30.3104], Tmin=(962.861,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane) + radical(Tertalkyl) + radical(Isobutyl)"""),
)

species(
    label = 'CC=C1C[CH][C]=CC1(27962)',
    structure = SMILES('CC=C1C[CH][C]=CC1'),
    E0 = (406.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6148,0.0311042,7.6155e-05,-1.09156e-07,3.98783e-11,48977.6,21.298], Tmin=(100,'K'), Tmax=(1021.95,'K')), NASAPolynomial(coeffs=[12.711,0.0413672,-1.77207e-05,3.49601e-09,-2.56978e-13,43905.8,-46.1942], Tmin=(1021.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexane) + radical(cyclohexene-allyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CC=CC(C)=CC(28172)',
    structure = SMILES('C#CC=CC(C)=CC'),
    E0 = (248.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0434438,0.077577,-5.03683e-05,8.30971e-09,2.78593e-12,30086.7,26.4457], Tmin=(100,'K'), Tmax=(1046.4,'K')), NASAPolynomial(coeffs=[17.591,0.0303814,-1.1691e-05,2.12973e-09,-1.48322e-13,25289.5,-64.7181], Tmin=(1046.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
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
    label = 'C#C[CH]CC([CH2])=C(28173)',
    structure = SMILES('C#C[CH]CC([CH2])=C'),
    E0 = (488.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,826.233,831.223],'cm^-1')),
        HinderedRotor(inertia=(3.29369,'amu*angstrom^2'), symmetry=1, barrier=(75.7283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.28752,'amu*angstrom^2'), symmetry=1, barrier=(75.5865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.545813,'amu*angstrom^2'), symmetry=1, barrier=(12.5493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242193,'amu*angstrom^2'), symmetry=1, barrier=(75.6123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.312219,0.0695724,-6.13802e-05,2.90648e-08,-5.43237e-12,58890.5,26.3061], Tmin=(100,'K'), Tmax=(1397.91,'K')), NASAPolynomial(coeffs=[15.5156,0.0222604,-6.5261e-06,9.55762e-10,-5.68192e-14,55012,-50.8031], Tmin=(1397.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Allyl_P)"""),
)

species(
    label = 'C#CC([CH2])C(=C)[CH]C(26577)',
    structure = SMILES('C#CC([CH2])C([CH2])=CC'),
    E0 = (478.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,750,770,3400,2100,273.559],'cm^-1')),
        HinderedRotor(inertia=(0.20589,'amu*angstrom^2'), symmetry=1, barrier=(10.9289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.2534,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25263,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205847,'amu*angstrom^2'), symmetry=1, barrier=(10.929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25256,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.119174,0.0887879,-8.75458e-05,4.76961e-08,-1.05787e-11,57667.9,29.155], Tmin=(100,'K'), Tmax=(1086.04,'K')), NASAPolynomial(coeffs=[14.3672,0.0354333,-1.38548e-05,2.46114e-09,-1.65913e-13,54521.3,-41.9294], Tmin=(1086.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = 'C#C[CH]CC[C]=CC(26584)',
    structure = SMILES('C#C[CH]CC[C]=CC'),
    E0 = (554.049,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3593.44,'J/mol'), sigma=(6.30673,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.29 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.346599,0.0809579,-7.75396e-05,4.41698e-08,-1.05714e-11,66767.8,30.364], Tmin=(100,'K'), Tmax=(996.474,'K')), NASAPolynomial(coeffs=[10.5847,0.0398595,-1.5672e-05,2.77757e-09,-1.86406e-13,64727.5,-18.9922], Tmin=(996.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(554.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC1CC(=CC)C1(27996)',
    structure = SMILES('C#CC1CC(=CC)C1'),
    E0 = (301.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.4345,-0.00487263,0.000126057,-1.23641e-07,2.84961e-11,35914.1,-16.3388], Tmin=(100,'K'), Tmax=(1723.76,'K')), NASAPolynomial(coeffs=[76.2092,0.0332752,-7.31529e-05,1.761e-08,-1.30181e-12,-15105.2,-451.692], Tmin=(1723.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(methylenecyclobutane)"""),
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
    label = '[CH2]C([CH2])=CC(24219)',
    structure = SMILES('[CH2]C([CH2])=CC'),
    E0 = (234.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.0177712,'amu*angstrom^2'), symmetry=1, barrier=(20.2255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.87837,'amu*angstrom^2'), symmetry=1, barrier=(20.1954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.61389,'amu*angstrom^2'), symmetry=1, barrier=(106.082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80244,0.0390483,-7.97741e-07,-2.45997e-08,1.16044e-11,28236,18.1778], Tmin=(100,'K'), Tmax=(1004.18,'K')), NASAPolynomial(coeffs=[10.9852,0.023482,-8.93213e-06,1.6381e-09,-1.15442e-13,25332.4,-31.4366], Tmin=(1004.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
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
    label = '[CH]=C=CC[C]=CC(27971)',
    structure = SMILES('C#C[CH]C[C]=CC'),
    E0 = (577.83,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,180,2463.27],'cm^-1')),
        HinderedRotor(inertia=(0.581389,'amu*angstrom^2'), symmetry=1, barrier=(13.3673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.581494,'amu*angstrom^2'), symmetry=1, barrier=(13.3697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.14472,'amu*angstrom^2'), symmetry=1, barrier=(72.3033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.14547,'amu*angstrom^2'), symmetry=1, barrier=(72.3205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.966803,0.0665056,-6.49262e-05,3.77707e-08,-9.15604e-12,69606.2,25.8907], Tmin=(100,'K'), Tmax=(989.891,'K')), NASAPolynomial(coeffs=[9.62785,0.0315073,-1.18923e-05,2.0534e-09,-1.35446e-13,67891.6,-15.8058], Tmin=(989.891,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=C1[CH]CC(=C)C1C(28070)',
    structure = SMILES('[CH]C1=CCC(=C)C1C'),
    E0 = (413.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803804,0.0569721,5.78996e-06,-3.79323e-08,1.57157e-11,49850.6,25.5452], Tmin=(100,'K'), Tmax=(1058.1,'K')), NASAPolynomial(coeffs=[11.8521,0.0445446,-1.81843e-05,3.37815e-09,-2.36425e-13,45870.2,-36.1413], Tmin=(1058.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(4-Methylenecyclopentene) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#C[CH]CC(=C)C=C(27575)',
    structure = SMILES('C#C[CH]CC(=C)C=C'),
    E0 = (426.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,252.368,253.651,254.43],'cm^-1')),
        HinderedRotor(inertia=(0.379159,'amu*angstrom^2'), symmetry=1, barrier=(17.2126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.60486,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.63777,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.65163,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.560284,0.08356,-7.60547e-05,3.61026e-08,-6.63641e-12,51440.9,30.5272], Tmin=(100,'K'), Tmax=(1473.29,'K')), NASAPolynomial(coeffs=[19.2502,0.0221348,-5.73796e-06,7.64566e-10,-4.27438e-14,46432.7,-69.9099], Tmin=(1473.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH]CC(=C)C[CH2](28174)',
    structure = SMILES('C#C[CH]CC(=C)C[CH2]'),
    E0 = (519.559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0123319,0.0813589,-7.15197e-05,3.48173e-08,-6.88735e-12,62638.8,32.5783], Tmin=(100,'K'), Tmax=(1216.05,'K')), NASAPolynomial(coeffs=[14.7761,0.0327148,-1.15171e-05,1.9225e-09,-1.24702e-13,59042.1,-41.6599], Tmin=(1216.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(519.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJ) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH][CH]C(=C)CC(28175)',
    structure = SMILES('[CH]=C=CC=C([CH2])CC'),
    E0 = (409.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.180364,0.0722758,-3.25195e-05,-1.33011e-08,1.18889e-11,49376.3,30.2233], Tmin=(100,'K'), Tmax=(952.705,'K')), NASAPolynomial(coeffs=[16.8723,0.0300952,-1.00378e-05,1.70765e-09,-1.16233e-13,44929.6,-56.1424], Tmin=(952.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=C(CC)C[CH]C#C(28176)',
    structure = SMILES('[CH]=C(CC)C[CH]C#C'),
    E0 = (561.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,422.647,422.656],'cm^-1')),
        HinderedRotor(inertia=(0.0740902,'amu*angstrom^2'), symmetry=1, barrier=(9.3921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0740967,'amu*angstrom^2'), symmetry=1, barrier=(9.39313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0741095,'amu*angstrom^2'), symmetry=1, barrier=(9.39238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.611023,'amu*angstrom^2'), symmetry=1, barrier=(77.4098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.610557,'amu*angstrom^2'), symmetry=1, barrier=(77.4098,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.100004,0.0838844,-7.68632e-05,3.89252e-08,-7.97579e-12,67674.9,31.5089], Tmin=(100,'K'), Tmax=(1178.76,'K')), NASAPolynomial(coeffs=[15.0601,0.0324398,-1.13983e-05,1.90017e-09,-1.23196e-13,64100.9,-44.1229], Tmin=(1178.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]C)CCC#C(28177)',
    structure = SMILES('[CH]C(=CC)CCC#C'),
    E0 = (525.407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.177346,0.0807729,-6.09558e-05,2.57263e-08,-4.60703e-12,63332.1,30.2746], Tmin=(100,'K'), Tmax=(1287.94,'K')), NASAPolynomial(coeffs=[12.2918,0.0431487,-1.71368e-05,3.04464e-09,-2.0433e-13,60211.6,-31.2359], Tmin=(1287.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]#C[CH]CC(=C)CC(28178)',
    structure = SMILES('[C]#C[CH]CC(=C)CC'),
    E0 = (651.457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0943572,0.082762,-7.88588e-05,4.35047e-08,-9.84281e-12,78495.4,31.1503], Tmin=(100,'K'), Tmax=(1066.27,'K')), NASAPolynomial(coeffs=[12.6444,0.0356837,-1.26331e-05,2.1e-09,-1.35375e-13,75819,-30.2023], Tmin=(1066.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Acetyl)"""),
)

species(
    label = 'C#CC1C[C]([CH]C)C1(28179)',
    structure = SMILES('C#CC1C[C]([CH]C)C1'),
    E0 = (574.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00134,0.0567975,-1.59065e-05,-1.0203e-08,5.58415e-12,69204.2,30.3673], Tmin=(100,'K'), Tmax=(1130.37,'K')), NASAPolynomial(coeffs=[10.3932,0.0402243,-1.60232e-05,2.90539e-09,-1.98865e-13,66016.5,-20.8022], Tmin=(1130.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(574.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane) + radical(Tertalkyl) + radical(Cs_S)"""),
)

species(
    label = 'C=C1C[CH][C]=CC1C(28085)',
    structure = SMILES('C=C1C[CH][C]=CC1C'),
    E0 = (408.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34643,0.039656,5.17881e-05,-8.40352e-08,3.11047e-11,49280.2,19.3496], Tmin=(100,'K'), Tmax=(1038.46,'K')), NASAPolynomial(coeffs=[12.6846,0.0420227,-1.8133e-05,3.54585e-09,-2.57793e-13,44442.9,-47.7309], Tmin=(1038.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(Cds_S) + radical(cyclohexene-allyl)"""),
)

species(
    label = 'C#CC=CC(=C)CC(28180)',
    structure = SMILES('C#CC=CC(=C)CC'),
    E0 = (262.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00630962,0.0732859,-3.03189e-05,-1.7942e-08,1.34066e-11,31697.3,27.3719], Tmin=(100,'K'), Tmax=(985.233,'K')), NASAPolynomial(coeffs=[19.3941,0.0271938,-9.81031e-06,1.78757e-09,-1.27387e-13,26293.7,-73.9097], Tmin=(985.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CCCC(=C)C=C(27589)',
    structure = SMILES('C#CCCC(=C)C=C'),
    E0 = (280.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0606653,0.0755317,-4.29141e-05,-3.01587e-09,8.32733e-12,33829.3,28.4525], Tmin=(100,'K'), Tmax=(960.435,'K')), NASAPolynomial(coeffs=[17.4693,0.0286318,-9.653e-06,1.65262e-09,-1.12719e-13,29304.4,-60.9788], Tmin=(960.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#C[CH]CC(C)[C]=C(26580)',
    structure = SMILES('C#C[CH]CC(C)[C]=C'),
    E0 = (558.323,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3534.47,'J/mol'), sigma=(6.26728,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=552.08 K, Pc=32.58 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0763857,0.0807808,-7.15244e-05,3.5528e-08,-7.21357e-12,67296.6,31.5209], Tmin=(100,'K'), Tmax=(1182.29,'K')), NASAPolynomial(coeffs=[13.9108,0.0339756,-1.21418e-05,2.04361e-09,-1.33186e-13,64025.4,-37.5386], Tmin=(1182.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CC1CC(=C)C1C(28160)',
    structure = SMILES('C#CC1CC(=C)C1C'),
    E0 = (305.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[9.92201,-0.00250806,0.000124471,-1.23939e-07,2.88394e-11,36454.4,-15.6752], Tmin=(100,'K'), Tmax=(1708.66,'K')), NASAPolynomial(coeffs=[74.672,0.0341723,-7.30004e-05,1.75912e-08,-1.30315e-12,-13154.4,-443.162], Tmin=(1708.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(methylenecyclobutane)"""),
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
    label = 'C#C[CH]C[C]=C(27523)',
    structure = SMILES('C#C[CH]C[C]=C'),
    E0 = (613.855,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,1685,370,276.5,276.957],'cm^-1')),
        HinderedRotor(inertia=(0.00217879,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44516,'amu*angstrom^2'), symmetry=1, barrier=(78.3056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44431,'amu*angstrom^2'), symmetry=1, barrier=(78.2948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61257,0.0509025,-4.40529e-05,1.85449e-08,-1.84925e-12,73917.3,21.7687], Tmin=(100,'K'), Tmax=(854.459,'K')), NASAPolynomial(coeffs=[9.60387,0.0217389,-7.3325e-06,1.18963e-09,-7.60246e-14,72250.6,-17.2893], Tmin=(854.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH2]C1([CH]C)C=C=CC1(28181)',
    structure = SMILES('[CH2]C1([CH]C)C=C=CC1'),
    E0 = (747.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.791425,0.0568324,-2.46446e-06,-3.08311e-08,1.36261e-11,89993.3,28.1857], Tmin=(100,'K'), Tmax=(1075.01,'K')), NASAPolynomial(coeffs=[14.2247,0.0365383,-1.55743e-05,2.98977e-09,-2.13653e-13,85389.6,-45.5727], Tmin=(1075.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(747.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene) + radical(Cs_S) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C(=CC)C[C]=C=C(28182)',
    structure = SMILES('[CH2]C#CCC([CH2])=CC'),
    E0 = (448.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.265845,0.0733324,-5.16755e-05,1.87586e-08,-2.77888e-12,54116.7,30.7277], Tmin=(100,'K'), Tmax=(1567.66,'K')), NASAPolynomial(coeffs=[16.1296,0.0328546,-1.29445e-05,2.28771e-09,-1.52191e-13,49142.9,-52.9377], Tmin=(1567.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(448.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Allyl_P) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]CC(C)=CC(28183)',
    structure = SMILES('[CH]=C=[C]CC(C)=CC'),
    E0 = (546.566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1685,370,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.54953,'amu*angstrom^2'), symmetry=1, barrier=(12.6348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549472,'amu*angstrom^2'), symmetry=1, barrier=(12.6335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549243,'amu*angstrom^2'), symmetry=1, barrier=(12.6282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549284,'amu*angstrom^2'), symmetry=1, barrier=(12.6291,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.450065,0.0804381,-7.50082e-05,4.1138e-08,-9.6044e-12,65862.6,30.0969], Tmin=(100,'K'), Tmax=(1009.52,'K')), NASAPolynomial(coeffs=[10.1936,0.041833,-1.76485e-05,3.26009e-09,-2.24529e-13,63895.3,-17.0025], Tmin=(1009.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2]C([CH]C=C=C)=CC(28184)',
    structure = SMILES('[CH2]C([CH]C)=CC=C=C'),
    E0 = (395.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.35884,0.0661228,-1.58305e-05,-2.7435e-08,1.55174e-11,47764.3,27.8502], Tmin=(100,'K'), Tmax=(987.888,'K')), NASAPolynomial(coeffs=[16.6374,0.0317915,-1.16554e-05,2.10819e-09,-1.48297e-13,43007,-58.286], Tmin=(987.888,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC=CCJ) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C(=[C]C)CC=C=C(28185)',
    structure = SMILES('[CH2]C(=[C]C)CC=C=C'),
    E0 = (543.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1685,370,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.791856,'amu*angstrom^2'), symmetry=1, barrier=(18.2063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792101,'amu*angstrom^2'), symmetry=1, barrier=(18.212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.335738,'amu*angstrom^2'), symmetry=1, barrier=(7.71927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0771764,'amu*angstrom^2'), symmetry=1, barrier=(18.2094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.379724,0.0777298,-6.1599e-05,2.60785e-08,-4.59768e-12,65510.3,29.9935], Tmin=(100,'K'), Tmax=(1315.27,'K')), NASAPolynomial(coeffs=[13.5281,0.0377428,-1.59958e-05,2.96373e-09,-2.04131e-13,62051.6,-37.043], Tmin=(1315.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C](C=C)CC=C=C(27596)',
    structure = SMILES('[CH2]C=C([CH2])CC=C=C'),
    E0 = (457.246,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,366.309,1042.98],'cm^-1')),
        HinderedRotor(inertia=(0.95395,'amu*angstrom^2'), symmetry=1, barrier=(21.9332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16646,'amu*angstrom^2'), symmetry=1, barrier=(26.8192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.952068,'amu*angstrom^2'), symmetry=1, barrier=(21.8899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.961552,'amu*angstrom^2'), symmetry=1, barrier=(22.108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0624414,0.0747112,-4.23516e-05,1.64994e-09,4.59371e-12,55146,30.943], Tmin=(100,'K'), Tmax=(1056.32,'K')), NASAPolynomial(coeffs=[17.0253,0.032143,-1.2669e-05,2.33302e-09,-1.63221e-13,50353.6,-57.5437], Tmin=(1056.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C=[C]CC(=C)CC(28186)',
    structure = SMILES('[CH]=C=[C]CC(=C)CC'),
    E0 = (559.946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,292.687,292.747],'cm^-1')),
        HinderedRotor(inertia=(0.328062,'amu*angstrom^2'), symmetry=1, barrier=(19.9643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187263,'amu*angstrom^2'), symmetry=1, barrier=(11.3965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187314,'amu*angstrom^2'), symmetry=1, barrier=(11.3964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187392,'amu*angstrom^2'), symmetry=1, barrier=(11.3972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.205152,0.0797035,-6.76564e-05,3.14917e-08,-6.04278e-12,67485.7,32.0733], Tmin=(100,'K'), Tmax=(1233.91,'K')), NASAPolynomial(coeffs=[13.8209,0.0355648,-1.39994e-05,2.50151e-09,-1.69133e-13,64125.6,-36.4769], Tmin=(1233.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([CH]C)CC=C=C(28187)',
    structure = SMILES('[CH]C(=CC)CC=C=C'),
    E0 = (524.932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.216171,0.0777674,-5.21836e-05,1.80521e-08,-2.59346e-12,63275,30.9403], Tmin=(100,'K'), Tmax=(1574.8,'K')), NASAPolynomial(coeffs=[14.9443,0.0403575,-1.65504e-05,2.96725e-09,-1.98721e-13,58636.3,-46.8028], Tmin=(1574.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C=CCC(=C)C=C(27600)',
    structure = SMILES('C=C=CCC(=C)C=C'),
    E0 = (279.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.119854,0.0723755,-3.41348e-05,-9.82397e-09,9.56411e-12,33771.4,29.0407], Tmin=(100,'K'), Tmax=(1004.98,'K')), NASAPolynomial(coeffs=[17.8141,0.0294074,-1.09848e-05,2.00542e-09,-1.41408e-13,28828.3,-63.3106], Tmin=(1004.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    E0 = (452.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (586.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (521.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (617.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (638.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (629.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (559.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (655.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (700.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (790.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (588.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (602.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (546.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (713.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (812.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (683.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (678.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (577.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (511.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (530.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (908.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (638.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (724.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (460.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (764.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (959.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (536.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (638.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (633.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (661.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (706.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (569.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (712.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (577.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (478.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (515.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (491.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (652.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (460.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (957.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (747.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (584.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (628.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (496.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (725.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (605.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (649.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (716.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (477.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['CH3CHCCH2(18175)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C#CC1CC1([CH2])[CH]C(28161)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(133.938,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_csHCt]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic
Ea raised from 131.2 to 133.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['[CH]=C1[CH]CC(=CC)C1(27981)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(69.2285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;triplebond_intra_H;radadd_intra_cs2H] for rate rule [R6;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#CC=CC([CH2])=CC(28162)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.56e+08,'cm^3/(mol*s)'), n=1.64, Ea=(4.97896,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2586 used for Cds-CdH_Cds-CtH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CtH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C][CH]C(18176)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0132398,'m^3/(mol*s)'), n=2.333, Ea=(2.89605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CtH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]C=C(4699)', 'CH3CHCCH2(18175)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C#CC[CH]C([CH2])=CC(28163)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.03429e+08,'s^-1'), n=1.36289, Ea=(107.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/OneDe;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_H/Ct;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C#C[CH][CH]C(C)=CC(28164)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.09894e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#C[CH]CC(C)=[C]C(28165)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[C]#CCCC([CH2])=CC(28166)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#CCCC([CH2])=[C]C(28167)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_H/OneDe] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_H/Ct]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#C[CH]C[C](C)C=C(27577)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#CCC[C]([CH2])C=C(27579)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(252000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;C_rad_out_2H;Cs_H_out_H/OneDe] for rate rule [R5H_SMSS;C_rad_out_2H;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[C]#C[CH]CC(C)=CC(28168)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.31883e+06,'s^-1'), n=1.02765, Ea=(75.0925,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R6HJ_2;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]C=C(4699)', 'C=[C][CH]C(18176)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C#C[CH]C[C]1CC1C(28169)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['[CH2]C(=CC)CC1[C]=C1(28170)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_T;triplebond_intra_H;radadd_intra_cs] for rate rule [R3_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C#CC1C[C]([CH2])C1C(28171)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.78932e+07,'s^-1'), n=1.19089, Ea=(125.213,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_HH_D;doublebond_intra;radadd_intra_cs] for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_csHCt]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['CC=C1C[CH][C]=CC1(27962)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.33254e+09,'s^-1'), n=0.487896, Ea=(59.5573,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6_linear;multiplebond_intra;radadd_intra_cs2H] + [R6_linear;triplebond_intra_H;radadd_intra] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C#CC=CC(C)=CC(28172)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(S)(23)', 'C#C[CH]CC([CH2])=C(28173)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C#CC1CC(=CC)C1(27996)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C3H2(81)', '[CH2]C([CH2])=CC(24219)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.26928e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(19)', '[CH]=C=CC[C]=CC(27971)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['[CH]=C1[CH]CC(=C)C1C(28070)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.98513e+09,'s^-1'), n=0.768, Ea=(84.4457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;triplebond_intra_H;radadd_intra_cs] for rate rule [R6;triplebond_intra_H;radadd_intra_csHNd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', 'C#C[CH]CC(=C)C=C(27575)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#C[CH]CC(=C)C[CH2](28174)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C#C[CH][CH]C(=C)CC(28175)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.23617e+10,'s^-1'), n=1.04667, Ea=(209.2,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C(CC)C[CH]C#C(28176)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C([CH]C)CCC#C(28177)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[C]#C[CH]CC(=C)CC(28178)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(8.93443e+06,'s^-1'), n=0.836584, Ea=(61.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6HJ_2;Ct_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C#CC1C[C]([CH]C)C1(28179)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.78932e+07,'s^-1'), n=1.19089, Ea=(125.213,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_HH_D;doublebond_intra;radadd_intra_cs] for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_csHCt]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C=C1C[CH][C]=CC1C(28085)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(487000,'s^-1'), n=1.17, Ea=(26.3592,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHCs] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C#CC=CC(=C)CC(28180)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C#CCCC(=C)C=C(27589)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C#C[CH]CC(C)[C]=C(26580)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C#CC1CC(=C)C1C(28160)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] + [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CHCH3(T)(95)', 'C#C[CH]C[C]=C(27523)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction41',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['[CH2]C1([CH]C)C=C=CC1(28181)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(7.4226e+09,'s^-1'), n=0.3735, Ea=(294.758,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra;radadd_intra_cdsingleH] for rate rule [R6;doublebond_intra_HNd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 291.4 to 294.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction42',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['[CH2]C(=CC)C[C]=C=C(28182)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=C=[C]CC(C)=CC(28183)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(32400,'s^-1'), n=2.04, Ea=(82.4248,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SS(Cd)S;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_SS(Cd)S;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['[CH2]C([CH]C=C=C)=CC(28184)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C(=[C]C)CC=C=C(28185)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.74e+08,'s^-1'), n=1.713, Ea=(181.895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;Cd_H_out_singleH] for rate rule [R6H;Cd_rad_out_Cs;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2][C](C=C)CC=C=C(27596)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(14131.6,'s^-1'), n=2.22233, Ea=(148.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_2H;Cd_H_out_singleH] + [R7H;C_rad_out_2H;XH_out] for rate rule [R7H;C_rad_out_2H;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=C=[C]CC(=C)CC(28186)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(199605,'s^-1'), n=1.94475, Ea=(89.3818,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;Cs_H_out_H/NonDeC] + [R4H_SS(Cd)S;Y_rad_out;Cs_H_out_1H] for rate rule [R4H_SS(Cd)S;Cd_rad_out_double;Cs_H_out_H/NonDeC]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=C([CH]C)CC=C=C(28187)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R6H;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C#C[CH]CC(=C)[CH]C(26576)'],
    products = ['C=C=CCC(=C)C=C(27600)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

network(
    label = '4548',
    isomers = [
        'C#C[CH]CC(=C)[CH]C(26576)',
    ],
    reactants = [
        ('CH3CHCCH2(18175)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4548',
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

