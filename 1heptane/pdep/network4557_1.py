species(
    label = 'C#CC([CH2])C[C]=CC(26585)',
    structure = SMILES('C#CC([CH2])C[C]=CC'),
    E0 = (603.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,750,770,3400,2100,218.855,219.126],'cm^-1')),
        HinderedRotor(inertia=(0.00350253,'amu*angstrom^2'), symmetry=1, barrier=(0.119786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161609,'amu*angstrom^2'), symmetry=1, barrier=(5.50317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.9251,'amu*angstrom^2'), symmetry=1, barrier=(101.848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170685,'amu*angstrom^2'), symmetry=1, barrier=(5.50478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.01599,'amu*angstrom^2'), symmetry=1, barrier=(102.07,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.301188,0.0802179,-7.47905e-05,4.08001e-08,-9.26938e-12,72754.1,32.1875], Tmin=(100,'K'), Tmax=(1051.89,'K')), NASAPolynomial(coeffs=[11.4468,0.0378349,-1.43523e-05,2.4957e-09,-1.65688e-13,70409.3,-22.1473], Tmin=(1051.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_S)"""),
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
    label = '[CH]=C1CC1C[C]=CC(27984)',
    structure = SMILES('[CH]=C1CC1C[C]=CC'),
    E0 = (684.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.398251,0.0681712,-3.34164e-05,-2.6329e-09,5.13843e-12,82444.5,30.0589], Tmin=(100,'K'), Tmax=(1068.97,'K')), NASAPolynomial(coeffs=[14.9979,0.0336819,-1.32833e-05,2.43745e-09,-1.69667e-13,78172.4,-46.7325], Tmin=(1068.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(684.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C(=CC)CC1[CH2](27985)',
    structure = SMILES('[CH]=C1C(=CC)CC1[CH2]'),
    E0 = (570.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.643732,0.0601405,-5.19823e-06,-3.60039e-08,1.83076e-11,68755.5,28.0056], Tmin=(100,'K'), Tmax=(973.561,'K')), NASAPolynomial(coeffs=[15.2841,0.03205,-1.1316e-05,2.01158e-09,-1.40562e-13,64385.4,-50.0367], Tmin=(973.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(570.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(Isobutyl) + radical(Cds_P)"""),
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
    label = 'C#CC(=C)C[C]=CC(27986)',
    structure = SMILES('C#CC(=C)C[C]=CC'),
    E0 = (514.694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2175,525,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,882.17],'cm^-1')),
        HinderedRotor(inertia=(0.520476,'amu*angstrom^2'), symmetry=1, barrier=(11.9668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694339,'amu*angstrom^2'), symmetry=1, barrier=(15.9642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694309,'amu*angstrom^2'), symmetry=1, barrier=(15.9635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694275,'amu*angstrom^2'), symmetry=1, barrier=(15.9627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.661453,0.072114,-5.80992e-05,2.5423e-08,-4.65261e-12,62024.5,28.3827], Tmin=(100,'K'), Tmax=(1269.79,'K')), NASAPolynomial(coeffs=[12.2653,0.035561,-1.49198e-05,2.75311e-09,-1.89356e-13,59077.6,-30.3708], Tmin=(1269.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(514.694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Cds_S)"""),
)

species(
    label = 'C#CC([CH2])C=C=CC(27987)',
    structure = SMILES('C#CC([CH2])C=C=CC'),
    E0 = (506.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.755999,'amu*angstrom^2'), symmetry=1, barrier=(17.3819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75499,'amu*angstrom^2'), symmetry=1, barrier=(17.3587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76393,'amu*angstrom^2'), symmetry=1, barrier=(40.5563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75689,'amu*angstrom^2'), symmetry=1, barrier=(17.4024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.192023,0.0855445,-9.15098e-05,5.57727e-08,-1.39794e-11,61037.4,27.8663], Tmin=(100,'K'), Tmax=(961.439,'K')), NASAPolynomial(coeffs=[11.9784,0.0365084,-1.50059e-05,2.72467e-09,-1.8552e-13,58771,-28.5325], Tmin=(961.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(506.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC([CH2])CC#CC(27988)',
    structure = SMILES('C#CC([CH2])CC#CC'),
    E0 = (525.902,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2175,2250,500,525,550,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180,3893.25],'cm^-1')),
        HinderedRotor(inertia=(3.8157,'amu*angstrom^2'), symmetry=1, barrier=(87.7304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200788,'amu*angstrom^2'), symmetry=1, barrier=(87.7394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.81524,'amu*angstrom^2'), symmetry=1, barrier=(87.7198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.81503,'amu*angstrom^2'), symmetry=1, barrier=(87.715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.81573,'amu*angstrom^2'), symmetry=1, barrier=(87.7312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.1075,0.0769229,-7.14609e-05,3.80914e-08,-8.0606e-12,63399,31.8569], Tmin=(100,'K'), Tmax=(1278.16,'K')), NASAPolynomial(coeffs=[13.3751,0.029069,-7.8693e-06,1.04671e-09,-5.67094e-14,60524.7,-33.3843], Tmin=(1278.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.902,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl)"""),
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
    label = 'C2H(33)',
    structure = SMILES('[C]#C'),
    E0 = (557.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89868,0.0132988,-2.80733e-05,2.89485e-08,-1.07502e-11,67061.6,6.18548], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.6627,0.00382492,-1.36633e-06,2.13455e-10,-1.23217e-14,67168.4,3.92206], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(557.301,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""C2H""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=CC[C]=CC(24188)',
    structure = SMILES('C=CC[C]=CC'),
    E0 = (290.566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,294.596,295.516],'cm^-1')),
        HinderedRotor(inertia=(0.158995,'amu*angstrom^2'), symmetry=1, barrier=(9.8665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163218,'amu*angstrom^2'), symmetry=1, barrier=(9.87073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156473,'amu*angstrom^2'), symmetry=1, barrier=(9.84833,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42465,0.0495641,-2.49826e-05,3.78863e-09,4.90012e-13,35045.5,24.3622], Tmin=(100,'K'), Tmax=(1316.94,'K')), NASAPolynomial(coeffs=[10.9212,0.0293798,-1.18562e-05,2.13691e-09,-1.44306e-13,31793.3,-26.9186], Tmin=(1316.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = 'CH3(17)',
    structure = SMILES('[CH3]'),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C#CCC([CH2])C#C(27989)',
    structure = SMILES('C#CCC([CH2])C#C'),
    E0 = (568.084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,756.667,763.333,770,3350,3450,2000,2200,2100,2250,500,550,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,251.497],'cm^-1')),
        HinderedRotor(inertia=(1.46184,'amu*angstrom^2'), symmetry=1, barrier=(67.9331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183916,'amu*angstrom^2'), symmetry=1, barrier=(8.84531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144444,'amu*angstrom^2'), symmetry=1, barrier=(67.9519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49853,'amu*angstrom^2'), symmetry=1, barrier=(67.9276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.59172,0.0679603,-6.93556e-05,3.89913e-08,-8.54012e-12,68453.4,26.6984], Tmin=(100,'K'), Tmax=(1248.07,'K')), NASAPolynomial(coeffs=[13.2009,0.0211117,-5.31429e-06,6.50826e-10,-3.2446e-14,65807.3,-34.9194], Tmin=(1248.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(568.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Isobutyl)"""),
)

species(
    label = 'C#C[C](C)C[C]=CC(27990)',
    structure = SMILES('C#C[C](C)C[C]=CC'),
    E0 = (533.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,360,370,350,2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.225815,'amu*angstrom^2'), symmetry=1, barrier=(5.19193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04843,'amu*angstrom^2'), symmetry=1, barrier=(24.1055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225873,'amu*angstrom^2'), symmetry=1, barrier=(5.19328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225909,'amu*angstrom^2'), symmetry=1, barrier=(5.19409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.17152,'amu*angstrom^2'), symmetry=1, barrier=(72.9194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.292636,0.0834632,-8.64102e-05,5.43793e-08,-1.43693e-11,64337.2,29.8571], Tmin=(100,'K'), Tmax=(908.122,'K')), NASAPolynomial(coeffs=[9.93613,0.0409886,-1.62558e-05,2.88045e-09,-1.9272e-13,62585.6,-15.7382], Tmin=(908.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(533.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Tert_Propargyl)"""),
)

species(
    label = 'C#CC([CH2])[CH]C=CC(27703)',
    structure = SMILES('C#CC([CH2])C=C[CH]C'),
    E0 = (484.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0680949,0.0803632,-6.75223e-05,3.08321e-08,-5.75547e-12,58389,29.25], Tmin=(100,'K'), Tmax=(1274.72,'K')), NASAPolynomial(coeffs=[14.8678,0.0339226,-1.28744e-05,2.25183e-09,-1.50285e-13,54615.9,-45.7426], Tmin=(1274.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Allyl_S)"""),
)

species(
    label = 'C#CC([CH2])CC=[C]C(27700)',
    structure = SMILES('C#CC([CH2])CC=[C]C'),
    E0 = (603.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,750,770,3400,2100,218.855,219.126],'cm^-1')),
        HinderedRotor(inertia=(0.00350253,'amu*angstrom^2'), symmetry=1, barrier=(0.119786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161609,'amu*angstrom^2'), symmetry=1, barrier=(5.50317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.9251,'amu*angstrom^2'), symmetry=1, barrier=(101.848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170685,'amu*angstrom^2'), symmetry=1, barrier=(5.50478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.01599,'amu*angstrom^2'), symmetry=1, barrier=(102.07,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.301188,0.0802179,-7.47905e-05,4.08001e-08,-9.26938e-12,72754.1,32.1875], Tmin=(100,'K'), Tmax=(1051.89,'K')), NASAPolynomial(coeffs=[11.4468,0.0378349,-1.43523e-05,2.4957e-09,-1.65688e-13,70409.3,-22.1473], Tmin=(1051.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC(C)C=[C][CH]C(26416)',
    structure = SMILES('C#CC(C)C=[C][CH]C'),
    E0 = (517.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.347132,0.0794049,-6.63539e-05,3.02794e-08,-5.77196e-12,62314.6,26.9399], Tmin=(100,'K'), Tmax=(1225.55,'K')), NASAPolynomial(coeffs=[12.9023,0.0384269,-1.61992e-05,2.99662e-09,-2.06539e-13,59237.2,-36.1851], Tmin=(1225.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(517.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'C#C[C]([CH2])CC=CC(27706)',
    structure = SMILES('[CH]=C=C([CH2])CC=CC'),
    E0 = (460.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0714832,0.0798727,-6.43495e-05,2.75111e-08,-4.75967e-12,55506.7,31.7737], Tmin=(100,'K'), Tmax=(1378.23,'K')), NASAPolynomial(coeffs=[16.6259,0.0314123,-1.16075e-05,1.99916e-09,-1.32004e-13,50904.2,-54.1381], Tmin=(1378.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Allyl_P)"""),
)

species(
    label = 'C#CC([CH2])C[CH]C=C(26565)',
    structure = SMILES('C#CC([CH2])CC=C[CH2]'),
    E0 = (517.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,383.441,383.446],'cm^-1')),
        HinderedRotor(inertia=(0.00114658,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842492,'amu*angstrom^2'), symmetry=1, barrier=(87.9005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842493,'amu*angstrom^2'), symmetry=1, barrier=(87.9005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0861232,'amu*angstrom^2'), symmetry=1, barrier=(8.98551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842475,'amu*angstrom^2'), symmetry=1, barrier=(87.9004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3619.53,'J/mol'), sigma=(6.34117,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=565.36 K, Pc=32.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.218991,0.074501,-4.63476e-05,4.60551e-09,4.93605e-12,62379.2,32.2873], Tmin=(100,'K'), Tmax=(952.65,'K')), NASAPolynomial(coeffs=[15.1165,0.0319601,-1.08728e-05,1.82972e-09,-1.21895e-13,58632.8,-43.6277], Tmin=(952.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(517.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[C]#CC(C)C[C]=CC(27991)',
    structure = SMILES('[C]#CC(C)C[C]=CC'),
    E0 = (735.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.419932,0.0806603,-7.9616e-05,4.80541e-08,-1.23069e-11,88630.5,30.7274], Tmin=(100,'K'), Tmax=(930.746,'K')), NASAPolynomial(coeffs=[9.50466,0.0416174,-1.66937e-05,2.9845e-09,-2.01053e-13,86939.4,-12.4491], Tmin=(930.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(735.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Acetyl)"""),
)

species(
    label = '[C]#CC([CH2])CC=CC(27707)',
    structure = SMILES('[C]#CC([CH2])CC=CC'),
    E0 = (703.098,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283818,0.0800483,-7.57904e-05,4.2746e-08,-1.00176e-11,84698.4,32.5167], Tmin=(100,'K'), Tmax=(1024.87,'K')), NASAPolynomial(coeffs=[11.1717,0.0375534,-1.35946e-05,2.28809e-09,-1.48502e-13,82466.7,-20.2784], Tmin=(1024.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(703.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Acetyl)"""),
)

species(
    label = 'C#CC(C)C[C]=[C]C(27992)',
    structure = SMILES('C#CC(C)C[C]=[C]C'),
    E0 = (636.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,185.118,2290.2],'cm^-1')),
        HinderedRotor(inertia=(0.563046,'amu*angstrom^2'), symmetry=1, barrier=(13.6251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.498806,'amu*angstrom^2'), symmetry=1, barrier=(13.6562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.550213,'amu*angstrom^2'), symmetry=1, barrier=(13.638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.527872,'amu*angstrom^2'), symmetry=1, barrier=(13.6107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.5401,'amu*angstrom^2'), symmetry=1, barrier=(70.8941,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.448284,0.0806892,-7.80791e-05,4.53609e-08,-1.12229e-11,76685.7,30.3597], Tmin=(100,'K'), Tmax=(957.745,'K')), NASAPolynomial(coeffs=[9.7043,0.0420297,-1.7528e-05,3.21038e-09,-2.19768e-13,74912.8,-13.8951], Tmin=(957.745,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(C)C[C]=C[CH2](27702)',
    structure = SMILES('C#CC(C)C[C]=C[CH2]'),
    E0 = (550.213,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0209563,'amu*angstrom^2'), symmetry=1, barrier=(14.0974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.612973,'amu*angstrom^2'), symmetry=1, barrier=(14.0935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.612763,'amu*angstrom^2'), symmetry=1, barrier=(14.0886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.69703,'amu*angstrom^2'), symmetry=1, barrier=(85.0021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126326,'amu*angstrom^2'), symmetry=1, barrier=(84.9442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0100579,0.0792003,-6.44596e-05,2.82527e-08,-5.03521e-12,66326.1,31.7335], Tmin=(100,'K'), Tmax=(1337.34,'K')), NASAPolynomial(coeffs=[15.6875,0.0323091,-1.18652e-05,2.03439e-09,-1.34018e-13,62132.9,-48.4582], Tmin=(1337.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(550.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'CC=[C]CC1[C]=CC1(27970)',
    structure = SMILES('CC=[C]CC1[C]=CC1'),
    E0 = (647.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.810406,0.0580837,-9.95914e-06,-2.28331e-08,1.12264e-11,77970,29.7746], Tmin=(100,'K'), Tmax=(1055.65,'K')), NASAPolynomial(coeffs=[13.3104,0.0359905,-1.44742e-05,2.69493e-09,-1.89506e-13,73922.8,-37.8768], Tmin=(1055.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Cds_S) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH2]C1[C]=CC(=CC)C1(27993)',
    structure = SMILES('[CH2]C1[C]=CC(=CC)C1'),
    E0 = (491.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.961058,0.0518598,1.47115e-05,-5.37401e-08,2.38393e-11,59268.8,27.5104], Tmin=(100,'K'), Tmax=(974.944,'K')), NASAPolynomial(coeffs=[14.1969,0.0333638,-1.19238e-05,2.14542e-09,-1.51251e-13,54986.2,-44.7369], Tmin=(974.944,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(3-Methylenecyclopentene) + radical(Isobutyl) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C#CC(=C)CC=CC(27709)',
    structure = SMILES('C#CC(=C)CC=CC'),
    E0 = (276.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.164423,0.0734581,-5.27267e-05,1.95054e-08,-2.92228e-12,33444.7,29.7132], Tmin=(100,'K'), Tmax=(1567.06,'K')), NASAPolynomial(coeffs=[16.9824,0.03053,-1.16362e-05,2.02471e-09,-1.3354e-13,28173.7,-58.9784], Tmin=(1567.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC(C)C=C=CC(26421)',
    structure = SMILES('C#CC(C)C=C=CC'),
    E0 = (301.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.300839,0.0824127,-7.41365e-05,3.71326e-08,-7.79039e-12,36368.4,25.6906], Tmin=(100,'K'), Tmax=(1121.32,'K')), NASAPolynomial(coeffs=[12.2681,0.0397227,-1.70294e-05,3.18005e-09,-2.20589e-13,33684.6,-33.4146], Tmin=(1121.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH)"""),
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
    label = 'C#CC([CH2])C[C]=C(27994)',
    structure = SMILES('C#CC([CH2])C[C]=C'),
    E0 = (639.822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,288.058,288.059],'cm^-1')),
        HinderedRotor(inertia=(0.00203163,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0791242,'amu*angstrom^2'), symmetry=1, barrier=(4.6589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.365,'amu*angstrom^2'), symmetry=1, barrier=(80.36,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36479,'amu*angstrom^2'), symmetry=1, barrier=(80.36,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.766294,0.066825,-6.204e-05,3.27229e-08,-7.02112e-12,77072.9,28.7085], Tmin=(100,'K'), Tmax=(1126.75,'K')), NASAPolynomial(coeffs=[11.9273,0.0272029,-9.29226e-06,1.51336e-09,-9.64085e-14,74557.8,-26.4688], Tmin=(1126.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_S)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.346599,0.0809579,-7.75396e-05,4.41698e-08,-1.05714e-11,66767.8,30.364], Tmin=(100,'K'), Tmax=(996.474,'K')), NASAPolynomial(coeffs=[10.5847,0.0398595,-1.5672e-05,2.77757e-09,-1.86406e-13,64727.5,-18.9922], Tmin=(996.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(554.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC[CH]C[C]=CC(27995)',
    structure = SMILES('C#CC[CH]C[C]=CC'),
    E0 = (602.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.568446,0.0810105,-9.46252e-05,7.51152e-08,-2.54608e-11,72557.9,32.2917], Tmin=(100,'K'), Tmax=(818.037,'K')), NASAPolynomial(coeffs=[5.72931,0.0482894,-2.08999e-05,3.84578e-09,-2.6158e-13,71964,9.96092], Tmin=(818.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJCC) + radical(Cds_S)"""),
)

species(
    label = 'C#CC([CH2])C(=C)[CH]C(26577)',
    structure = SMILES('C#CC([CH2])C([CH2])=CC'),
    E0 = (478.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.119174,0.0887879,-8.75458e-05,4.76961e-08,-1.05787e-11,57667.9,29.155], Tmin=(100,'K'), Tmax=(1086.04,'K')), NASAPolynomial(coeffs=[14.3672,0.0354333,-1.38548e-05,2.46114e-09,-1.65913e-13,54521.3,-41.9294], Tmin=(1086.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Allyl_P)"""),
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
    label = '[C]=CC(24199)',
    structure = SMILES('[C]=CC'),
    E0 = (564.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.40488,'amu*angstrom^2'), symmetry=1, barrier=(9.30899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28451,0.0144722,-2.9151e-06,-2.04635e-09,7.91551e-13,67868.8,10.0416], Tmin=(100,'K'), Tmax=(1455.86,'K')), NASAPolynomial(coeffs=[5.67821,0.0121697,-4.94658e-06,9.00486e-10,-6.07653e-14,66718.9,-3.96135], Tmin=(1455.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'C#CC([CH2])[CH2](26629)',
    structure = SMILES('C#CC([CH2])[CH2]'),
    E0 = (526.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00245333,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00268778,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19638,'amu*angstrom^2'), symmetry=1, barrier=(57.3712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59865,0.0445704,-4.03451e-05,2.09026e-08,-4.1805e-12,63457.8,21.3678], Tmin=(100,'K'), Tmax=(1443.27,'K')), NASAPolynomial(coeffs=[9.29133,0.0160372,-3.19385e-06,2.79115e-10,-8.33799e-15,61988.6,-15.9649], Tmin=(1443.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(526.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    E0 = (603.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (684.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (684.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (742.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (726.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (756.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (652.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (867.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (613.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (734.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (737.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (761.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (837.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (721.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (745.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (765.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (831.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (747.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (798.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (656.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (812.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (734.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (662.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (692.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (682.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1059.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (763.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (763.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (773.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (612.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (959.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (1090.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['CH3CHCCH2(18175)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['[CH]=C1CC1C[C]=CC(27984)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(80.5299,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for R4_S_T;triplebond_intra_H;radadd_intra_cs2H
Exact match found for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 78.7 to 80.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['[CH]=C1C(=CC)CC1[CH2](27985)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.98674e+07,'s^-1'), n=1.31443, Ea=(80.4773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#CC(=C)C[C]=CC(27986)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(15.6063,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2632 used for Cds-CtCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CtCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#CC([CH2])C=C=CC(27987)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.23e+08,'cm^3/(mol*s)'), n=1.64, Ea=(8.49352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2566 used for Cds-CsH_Ca;HJ
Exact match found for rate rule [Cds-CsH_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C#CC([CH2])CC#CC(27988)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.02e+09,'cm^3/(mol*s)'), n=1.64, Ea=(18.4933,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2702 used for Ct-Cs_Ct-Cs;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C][CH]C(18176)', 'CH2CHCCH(26391)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00294841,'m^3/(mol*s)'), n=2.48333, Ea=(17.6885,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CtH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2H(33)', 'C=CC[C]=CC(24188)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CJ] for rate rule [Cds-CsH_Cds-HH;CtJ_Ct]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]C=C(4699)', 'CH3CHCCH2(18175)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00429749,'m^3/(mol*s)'), n=2.45395, Ea=(16.6104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH3(17)', 'C#CCC([CH2])C#C(27989)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(178000,'cm^3/(mol*s)'), n=2.41, Ea=(30.1248,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2271 used for Ct-H_Ct-Cs;CsJ-HHH
Exact match found for rate rule [Ct-H_Ct-Cs;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#C[C](C)C[C]=CC(27990)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#CC([CH2])[CH]C=CC(27703)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#CC([CH2])CC=[C]C(27700)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#CC(C)C=[C][CH]C(26416)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#C[C]([CH2])CC=CC(27706)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#CC([CH2])C[CH]C=C(26565)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[C]#CC(C)C[C]=CC(27991)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.39293e+07,'s^-1'), n=1.32074, Ea=(96.0416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_TSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[C]#CC([CH2])CC=CC(27707)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_TSSS;Ct_rad_out;Cd_H_out_doubleC]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#CC(C)C[C]=[C]C(27992)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#CC(C)C[C]=C[CH2](27702)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5436.63,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R6HJ_3;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C]C=C(4699)', 'C=[C][CH]C(18176)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['CC=[C]CC1[C]=CC1(27970)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(130.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['[CH2]C1[C]=CC(=CC)C1(27993)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 13 used for R5_SS_T;triplebond_intra_H;radadd_intra_cddouble
Exact match found for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#CC(=C)CC=CC(27709)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#CC(C)C=C=CC(26421)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(S)(23)', 'C#CC([CH2])C[C]=C(27994)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#CC[CH]C[C]=CC(27995)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#CC([CH2])C(=C)[CH]C(26577)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#CC1CC(=CC)C1(27996)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH2(19)', '[CH]=C=CC[C]=CC(27971)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['[C]=CC(24199)', 'C#CC([CH2])[CH2](26629)'],
    products = ['C#CC([CH2])C[C]=CC(26585)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4557',
    isomers = [
        'C#CC([CH2])C[C]=CC(26585)',
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
    label = '4557',
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

