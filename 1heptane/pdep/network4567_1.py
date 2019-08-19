species(
    label = 'C=C[C]=CC=[C]C=C(26595)',
    structure = SMILES('C=C[C]=CC=[C]C=C'),
    E0 = (595.456,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44439,'amu*angstrom^2'), symmetry=1, barrier=(33.2095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43535,'amu*angstrom^2'), symmetry=1, barrier=(33.0015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4485,'amu*angstrom^2'), symmetry=1, barrier=(33.3039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.79274,0.0874901,-8.59447e-05,4.2632e-08,-8.05122e-12,71804.7,28.9726], Tmin=(100,'K'), Tmax=(1460.09,'K')), NASAPolynomial(coeffs=[21.1639,0.0173709,-3.66883e-06,3.8996e-10,-1.78827e-14,66455.4,-81.6279], Tmin=(1460.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(595.456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
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
    label = '[CH2]C1[C]=CC=C1C=C(28006)',
    structure = SMILES('[CH2]C1[C]=CC=C1C=C'),
    E0 = (614.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.595677,0.061251,-1.4823e-05,-3.06625e-08,1.81721e-11,74083.7,23.9447], Tmin=(100,'K'), Tmax=(946.974,'K')), NASAPolynomial(coeffs=[17.569,0.022984,-7.1584e-06,1.21844e-09,-8.54266e-14,69370.2,-64.9312], Tmin=(946.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentadiene) + radical(Isobutyl) + radical(1,3-cyclopentadiene-vinyl-1)"""),
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
    label = 'C=C[C]=CC#CC=C(28007)',
    structure = SMILES('C=C[C]=CC#CC=C'),
    E0 = (577.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2100,2250,500,550,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.5507,'amu*angstrom^2'), symmetry=1, barrier=(35.6535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5618,'amu*angstrom^2'), symmetry=1, barrier=(35.9089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55291,'amu*angstrom^2'), symmetry=1, barrier=(35.7044,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.568973,0.0612327,-1.76608e-05,-2.74969e-08,1.67133e-11,69644.9,26.9635], Tmin=(100,'K'), Tmax=(965.388,'K')), NASAPolynomial(coeffs=[18.4863,0.0206725,-6.96839e-06,1.25623e-09,-9.08664e-14,64616.1,-66.9743], Tmin=(965.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-Ct(Cds-Cds)) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C=C=CC=[C]C=C(28008)',
    structure = SMILES('C=C=C=CC=[C]C=C'),
    E0 = (625.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,563.333,586.667,610,1970,2140,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44949,'amu*angstrom^2'), symmetry=1, barrier=(33.3266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44945,'amu*angstrom^2'), symmetry=1, barrier=(33.3258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.24134,0.082304,-8.22074e-05,4.18004e-08,-8.27208e-12,75433.9,27.1719], Tmin=(100,'K'), Tmax=(1266.82,'K')), NASAPolynomial(coeffs=[19.3685,0.0191687,-6.0098e-06,9.42872e-10,-5.94069e-14,70563.2,-71.6866], Tmin=(1266.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(625.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CJC=C)"""),
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
    label = 'C=C[C]=C[C]=CC=C(28009)',
    structure = SMILES('[CH2]C=[C]C=C=CC=C'),
    E0 = (567.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.303211,0.0702234,-3.77927e-05,-7.90254e-09,1.04294e-11,68373.3,27.7194], Tmin=(100,'K'), Tmax=(942.661,'K')), NASAPolynomial(coeffs=[17.5212,0.0242367,-7.69896e-06,1.28316e-09,-8.71982e-14,63924.2,-60.7115], Tmin=(942.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
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
    label = 'C=C[C]=[C]C=CC=C(28011)',
    structure = SMILES('C=CC#C[CH][CH]C=C'),
    E0 = (573.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.904837,0.0554434,-5.64721e-06,-3.4347e-08,1.80795e-11,69108.4,28.2521], Tmin=(100,'K'), Tmax=(954.583,'K')), NASAPolynomial(coeffs=[14.9718,0.0270556,-9.05617e-06,1.56778e-09,-1.08778e-13,65030.6,-46.2514], Tmin=(954.583,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Sec_Propargyl) + radical(Allyl_S)"""),
)

species(
    label = '[CH]=CC=CC=[C]C=C(28012)',
    structure = SMILES('[CH]=CC=CC=[C]C=C'),
    E0 = (643.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3120,650,792.5,1650,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24028,'amu*angstrom^2'), symmetry=1, barrier=(28.5165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24055,'amu*angstrom^2'), symmetry=1, barrier=(28.5226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23963,'amu*angstrom^2'), symmetry=1, barrier=(28.5016,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0553591,0.0758877,-4.37495e-05,-1.1187e-08,1.40738e-11,77560.3,27.8221], Tmin=(100,'K'), Tmax=(924.348,'K')), NASAPolynomial(coeffs=[21.5508,0.0169116,-4.0658e-06,5.9576e-10,-4.0643e-14,72091.2,-82.6929], Tmin=(924.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(643.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][C]=CC=CC=C(28013)',
    structure = SMILES('[CH2]C#CC=CC=C[CH2]'),
    E0 = (518.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.758841,0.0558626,3.33153e-07,-4.35112e-08,2.14756e-11,62520,28.999], Tmin=(100,'K'), Tmax=(971.941,'K')), NASAPolynomial(coeffs=[17.1262,0.0251686,-8.88251e-06,1.62288e-09,-1.1704e-13,57606.5,-58.4073], Tmin=(971.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(518.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=CC=CCJ) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C[C]=CC=CC=C(28014)',
    structure = SMILES('[CH]C=C=CC=CC=C'),
    E0 = (620.931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,540,610,2055,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.96296,'amu*angstrom^2'), symmetry=1, barrier=(45.1322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.96256,'amu*angstrom^2'), symmetry=1, barrier=(45.1232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.96279,'amu*angstrom^2'), symmetry=1, barrier=(45.1284,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0885534,0.0726542,-3.20749e-05,-1.47637e-08,1.21712e-11,74833.5,28.4439], Tmin=(100,'K'), Tmax=(977.71,'K')), NASAPolynomial(coeffs=[18.3075,0.0282452,-1.01657e-05,1.81497e-09,-1.2709e-13,69830.9,-66.4057], Tmin=(977.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(620.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C[C]=CC=C1[CH]C1(28015)',
    structure = SMILES('C=C[C]=CC=C1[CH]C1'),
    E0 = (630.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.575743,0.0577491,4.46449e-06,-5.63042e-08,2.85792e-11,75982.6,24.2277], Tmin=(100,'K'), Tmax=(937.859,'K')), NASAPolynomial(coeffs=[19.8974,0.0195709,-5.21346e-06,8.59786e-10,-6.28098e-14,70413.3,-78.1185], Tmin=(937.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(630.596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=CJC=C) + radical(Allyl_S)"""),
)

species(
    label = 'C=CC1=CC=[C][CH]C1(28016)',
    structure = SMILES('[CH2]C=C1[CH]C=C=CC1'),
    E0 = (448.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88206,0.0163432,0.000125511,-1.72866e-07,6.6101e-11,54077.3,20.254], Tmin=(100,'K'), Tmax=(967.879,'K')), NASAPolynomial(coeffs=[17.8092,0.0270486,-9.68272e-06,1.94644e-09,-1.53303e-13,47409.6,-74.583], Tmin=(967.879,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(448.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(Cyclohexane) + radical(C=CCJC=C) + radical(Allyl_P)"""),
)

species(
    label = 'C=C=C=CC=CC=C(28017)',
    structure = SMILES('C=C=C=CC=CC=C'),
    E0 = (426.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.13238,0.0723242,-4.13613e-05,-5.89461e-09,9.62027e-12,51488.9,26.2313], Tmin=(100,'K'), Tmax=(972.117,'K')), NASAPolynomial(coeffs=[19.548,0.0213291,-7.26035e-06,1.29572e-09,-9.22147e-14,46348.8,-73.9114], Tmin=(972.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    label = '[CH2]C1[C]=CC=C=CC1(28019)',
    structure = SMILES('[CH2]C1[C]=CC=C=CC1'),
    E0 = (706.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21388,0.0464091,1.818e-05,-5.61421e-08,2.48128e-11,85136.7,21.4365], Tmin=(100,'K'), Tmax=(967.016,'K')), NASAPolynomial(coeffs=[14.255,0.0285038,-9.94751e-06,1.78807e-09,-1.27105e-13,80929.5,-49.7544], Tmin=(967.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(706.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = 'C=C[C]=CC=C=[C]C(28020)',
    structure = SMILES('[CH2]C=[C]C=CC#CC'),
    E0 = (561.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.671796,0.0623863,-2.16031e-05,-2.05429e-08,1.41688e-11,67641.4,27.852], Tmin=(100,'K'), Tmax=(932.756,'K')), NASAPolynomial(coeffs=[15.245,0.0270899,-8.58159e-06,1.41248e-09,-9.47509e-14,63739.5,-47.7836], Tmin=(932.756,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C[C]=C[C]=C=CC(28021)',
    structure = SMILES('C=C[C]=CC#C[CH]C'),
    E0 = (603.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.296914,0.076968,-6.71966e-05,3.09118e-08,-5.50925e-12,72695.2,30.9518], Tmin=(100,'K'), Tmax=(1542.26,'K')), NASAPolynomial(coeffs=[17.6212,0.022047,-5.56323e-06,7.17704e-10,-3.9006e-14,68173.1,-59.9982], Tmin=(1542.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=CJC=C) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C=C[C]=[C]C=C=CC(28022)',
    structure = SMILES('C=CC#CC=[C][CH]C'),
    E0 = (647.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,355.914,357.099,357.498],'cm^-1')),
        HinderedRotor(inertia=(0.30608,'amu*angstrom^2'), symmetry=1, barrier=(27.4733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304033,'amu*angstrom^2'), symmetry=1, barrier=(27.5256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.305144,'amu*angstrom^2'), symmetry=1, barrier=(27.5201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.306601,'amu*angstrom^2'), symmetry=1, barrier=(27.4806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.767801,0.0575123,-1.10229e-05,-2.60388e-08,1.34911e-11,77997.6,28.106], Tmin=(100,'K'), Tmax=(1028.38,'K')), NASAPolynomial(coeffs=[15.8567,0.0283397,-1.15254e-05,2.19742e-09,-1.582e-13,73333.3,-52.7002], Tmin=(1028.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-Ct(Cds-Cds)) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = 'C=[C][C]=CC=C=CC(28023)',
    structure = SMILES('[CH2]C#CC=C[C]=CC'),
    E0 = (599.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.137338,0.0742517,-6.15164e-05,2.64667e-08,-4.55169e-12,72272.3,29.2787], Tmin=(100,'K'), Tmax=(1397.71,'K')), NASAPolynomial(coeffs=[17.1036,0.0256973,-9.40878e-06,1.61299e-09,-1.06273e-13,67529.4,-58.2549], Tmin=(1397.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(599.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=CJC=C) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C[C]=CC=C=CC(28024)',
    structure = SMILES('[CH]C=C=CC=C=CC'),
    E0 = (673.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.99253,'amu*angstrom^2'), symmetry=1, barrier=(45.8121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98892,'amu*angstrom^2'), symmetry=1, barrier=(45.7293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98492,'amu*angstrom^2'), symmetry=1, barrier=(45.6372,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.169722,0.0784434,-6.21841e-05,2.60143e-08,-4.46868e-12,81171.4,28.1856], Tmin=(100,'K'), Tmax=(1365.9,'K')), NASAPolynomial(coeffs=[15.1269,0.0346419,-1.40824e-05,2.53694e-09,-1.71655e-13,77085.4,-48.638], Tmin=(1365.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(673.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C=C1[CH]C=C1C=C(28025)',
    structure = SMILES('[CH2]C=C1[CH]C=C1C=C'),
    E0 = (474.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74334,0.021483,0.000108644,-1.62182e-07,6.55789e-11,57158.5,21.9913], Tmin=(100,'K'), Tmax=(933.985,'K')), NASAPolynomial(coeffs=[18.8896,0.0198844,-4.15538e-06,6.80071e-10,-5.63846e-14,50822.5,-76.3313], Tmin=(933.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CCJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C[C]=CC1[C]=CC1(28026)',
    structure = SMILES('C=C[C]=CC1[C]=CC1'),
    E0 = (719.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.613599,0.0646829,-3.48589e-05,-1.6223e-09,5.61018e-12,86657,25.7994], Tmin=(100,'K'), Tmax=(1007.12,'K')), NASAPolynomial(coeffs=[14.7615,0.0284016,-1.04759e-05,1.86714e-09,-1.2882e-13,82797.6,-47.5694], Tmin=(1007.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(719.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(C=CJC=C)"""),
)

species(
    label = '[C]1[CH]CCC=C=CC=1(28027)',
    structure = SMILES('[C]1=CCCC=[C]C=C1'),
    E0 = (559.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06344,0.0459548,3.40192e-05,-8.23559e-08,3.68929e-11,67472.2,17.1279], Tmin=(100,'K'), Tmax=(932.49,'K')), NASAPolynomial(coeffs=[17.2916,0.0237205,-6.42635e-06,1.04594e-09,-7.47385e-14,62385.9,-71.0739], Tmin=(932.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cyclooctatriene) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C=C=CC=C=CC(28028)',
    structure = SMILES('C=C=C=CC=C=CC'),
    E0 = (479.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.266901,0.0775067,-6.9491e-05,3.25757e-08,-6.16292e-12,57824.6,25.7806], Tmin=(100,'K'), Tmax=(1264.01,'K')), NASAPolynomial(coeffs=[15.6183,0.0289265,-1.18407e-05,2.16965e-09,-1.49091e-13,53943.7,-51.8777], Tmin=(1264.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(479.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    label = '[CH]=C=CC=[C]C=C(28029)',
    structure = SMILES('C#CC=C[C]=C[CH2]'),
    E0 = (603.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,750,770,3400,2100,2175,525,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.8259,'amu*angstrom^2'), symmetry=1, barrier=(41.981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82498,'amu*angstrom^2'), symmetry=1, barrier=(41.96,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82933,'amu*angstrom^2'), symmetry=1, barrier=(42.0598,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.1225,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06706,0.0544769,-2.32027e-05,-1.48007e-08,1.16032e-11,72699.6,23.0123], Tmin=(100,'K'), Tmax=(939.868,'K')), NASAPolynomial(coeffs=[15.2958,0.018748,-5.80419e-06,9.64012e-10,-6.61249e-14,68928.4,-50.5841], Tmin=(939.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
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
    label = 'C#CC([CH2])C([CH2])C#C(26592)',
    structure = SMILES('C#CC([CH2])C([CH2])C#C'),
    E0 = (740.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,343.401],'cm^-1')),
        HinderedRotor(inertia=(0.0868022,'amu*angstrom^2'), symmetry=1, barrier=(7.26059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.831556,'amu*angstrom^2'), symmetry=1, barrier=(69.5361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.831076,'amu*angstrom^2'), symmetry=1, barrier=(69.5377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.831467,'amu*angstrom^2'), symmetry=1, barrier=(69.5362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0247503,'amu*angstrom^2'), symmetry=1, barrier=(69.536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.422844,0.0761775,-6.30966e-05,1.35801e-08,7.2947e-12,89164.4,30.3868], Tmin=(100,'K'), Tmax=(783.028,'K')), NASAPolynomial(coeffs=[14.8814,0.0243001,-5.82893e-06,6.75401e-10,-3.20101e-14,86226.2,-40.1344], Tmin=(783.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(740.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    E0 = (595.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (661.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (803.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (854.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (735.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (807.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (809.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (788.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (865.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (863.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (725.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (903.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (733.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (607.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (634.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (603.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (754.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (722.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (803.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (833.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (680.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (783.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (805.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (733.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (719.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (652.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (620.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (985.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (860.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (884.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['CH2CHCCH(26391)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['[CH2]C1[C]=CC=C1C=C(28006)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.82843e+12,'s^-1'), n=0.00925, Ea=(66.3007,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_2H_pri;radadd_intra_cdsingleDe] + [R6;doublebond_intra_2H_pri;radadd_intra_cdsingle] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C[C]=CC#CC=C(28007)'],
    products = ['C=C[C]=CC=[C]C=C(26595)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2371.94,'m^3/(mol*s)'), n=1.49517, Ea=(13.5032,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-De_Ct-De;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C=C=CC=[C]C=C(28008)'],
    products = ['C=C[C]=CC=[C]C=C(26595)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=[C]C=C(4699)', 'CH2CHCCH(26391)'],
    products = ['C=C[C]=CC=[C]C=C(26595)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0117197,'m^3/(mol*s)'), n=2.33233, Ea=(9.74423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=C[C]=C[C]=CC=C(28009)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.189e+07,'s^-1'), n=1.951, Ea=(212.263,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_single] for rate rule [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_singleDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=[C]C=CC=[C]C=C(28010)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.612e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=C[C]=[C]C=CC=C(28011)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.56742e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CC=CC=[C]C=C(28012)'],
    products = ['C=C[C]=CC=[C]C=C(26595)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=[C][C]=CC=CC=C(28013)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.8431e+08,'s^-1'), n=1.50489, Ea=(267.746,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_single;Cd_H_out_doubleC] + [RnH;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleDe_Cd;Cd_H_out_doubleC]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C[C]=CC=CC=C(28014)'],
    products = ['C=C[C]=CC=[C]C=C(26595)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.80899e+07,'s^-1'), n=1.452, Ea=(104.082,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleDe] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]C=C(4699)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC=[C]C=C(26595)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.73038e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=C[C]=CC=C1[CH]C1(28015)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.121e+13,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=CC1=CC=[C][CH]C1(28016)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.33917e+10,'s^-1'), n=0.31, Ea=(12.393,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra_cdsingle] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=C=C=CC=CC=C(28017)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=CC1=CC=C1C=C(28005)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleDe_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['[CH2]C=[C]C1C=C1C=C(28018)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.252e+14,'s^-1'), n=-0.355, Ea=(158.923,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 158.0 to 158.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['[CH2]C1[C]=CC=C=CC1(28019)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.57638e+09,'s^-1'), n=0.5525, Ea=(127.194,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra_2H_pri;radadd_intra_cs2H] for rate rule [R8;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=C[C]=CC=C=[C]C(28020)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.26e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=C[C]=C[C]=C=CC(28021)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.005931,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C[C]=[C]C=C=CC(28022)'],
    products = ['C=C[C]=CC=[C]C=C(26595)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out;Cs_H_out_2H] for rate rule [R5H_SMMS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=[C][C]=CC=C=CC(28023)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(8.23124e+07,'s^-1'), n=1.49763, Ea=(188.308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_2H;Cd_H_out_doubleC] for rate rule [R7HJ_5;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C[C]=CC=C=CC(28024)'],
    products = ['C=C[C]=CC=[C]C=C(26595)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R8Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['[CH2]C=C1[CH]C=C1C=C(28025)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=C[C]=CC1[C]=CC1(28026)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.07328e+10,'s^-1'), n=0.43543, Ea=(123.966,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 121.6 to 124.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['[C]1[CH]CCC=C=CC=1(28027)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.04798e+07,'s^-1'), n=0.841667, Ea=(57.3368,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_2H;radadd_intra_cs2H] for rate rule [R8_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['C=C=C=CC=C=CC(28028)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(19)', '[CH]=C=CC=[C]C=C(28029)'],
    products = ['C=C[C]=CC=[C]C=C(26595)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C[C]=CC=[C]C=C(26595)'],
    products = ['[CH2]C=[C]C1C=C=CC1(28030)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2e+11,'s^-1'), n=0.21, Ea=(265.449,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R6_SMM;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 262.6 to 265.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#CC([CH2])C([CH2])C#C(26592)'],
    products = ['C=C[C]=CC=[C]C=C(26595)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.30946e+10,'s^-1'), n=0.360276, Ea=(144.706,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for 1_5_hexadiyne
Exact match found for rate rule [1_5_hexadiyne]
Euclidian distance = 0
family: 6_membered_central_C-C_shift"""),
)

network(
    label = '4567',
    isomers = [
        'C=C[C]=CC=[C]C=C(26595)',
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
    label = '4567',
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

