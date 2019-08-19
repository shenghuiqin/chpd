species(
    label = 'C=C[C]=CC[CH]C=C(26566)',
    structure = SMILES('[CH2]C=CCC=[C]C=C'),
    E0 = (454.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,187.6,187.632,187.639],'cm^-1')),
        HinderedRotor(inertia=(0.803805,'amu*angstrom^2'), symmetry=1, barrier=(20.0799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.80385,'amu*angstrom^2'), symmetry=1, barrier=(20.0801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.803801,'amu*angstrom^2'), symmetry=1, barrier=(20.08,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43893,'amu*angstrom^2'), symmetry=1, barrier=(35.9474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.164635,0.0727761,-3.71724e-05,-5.30304e-09,7.85742e-12,54871,30.9083], Tmin=(100,'K'), Tmax=(997.04,'K')), NASAPolynomial(coeffs=[16.6052,0.0312354,-1.14101e-05,2.03301e-09,-1.40736e-13,50379,-54.4454], Tmin=(997.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Allyl_P)"""),
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
    label = '[CH2][CH]C1CC=C1C=C(27622)',
    structure = SMILES('[CH2][CH]C1CC=C1C=C'),
    E0 = (538.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.753105,0.0529293,2.09149e-05,-6.60601e-08,2.92009e-11,64951.8,30.2778], Tmin=(100,'K'), Tmax=(976.736,'K')), NASAPolynomial(coeffs=[17.5122,0.0287961,-1.03628e-05,1.93332e-09,-1.41263e-13,59555.3,-61.0471], Tmin=(976.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(RCCJ) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C1[C]=CCC=CC1(27623)',
    structure = SMILES('[CH2]C1[C]=CCC=CC1'),
    E0 = (501.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51728,-0.0040754,0.000252019,-3.6467e-07,1.52621e-10,60463.5,29.5066], Tmin=(100,'K'), Tmax=(901.705,'K')), NASAPolynomial(coeffs=[39.3389,-0.018215,1.99602e-05,-4.13921e-09,2.73406e-13,47396.7,-183.682], Tmin=(901.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(Isobutyl)"""),
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
    label = '[CH2]C=CCC#CC=C(27624)',
    structure = SMILES('C=CC#CC[CH]C=C'),
    E0 = (427.368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2100,2250,500,550,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,484.087,487.658,488.493,490.331],'cm^-1')),
        HinderedRotor(inertia=(0.230109,'amu*angstrom^2'), symmetry=1, barrier=(39.8124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115862,'amu*angstrom^2'), symmetry=1, barrier=(20.0031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116235,'amu*angstrom^2'), symmetry=1, barrier=(20.0006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.542739,'amu*angstrom^2'), symmetry=1, barrier=(89.5385,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.85218,0.0552235,7.88763e-07,-3.96427e-08,1.89438e-11,51526.3,28.6033], Tmin=(100,'K'), Tmax=(985.476,'K')), NASAPolynomial(coeffs=[14.9645,0.0305316,-1.12317e-05,2.04594e-09,-1.44844e-13,47162.4,-47.3029], Tmin=(985.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Allyl_S)"""),
)

species(
    label = 'C=C=CCC=[C]C=C(27625)',
    structure = SMILES('C=C=CCC=[C]C=C'),
    E0 = (480.096,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0876,'amu*angstrom^2'), symmetry=1, barrier=(25.0062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08811,'amu*angstrom^2'), symmetry=1, barrier=(25.0179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08736,'amu*angstrom^2'), symmetry=1, barrier=(25.0004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.152332,0.0755026,-5.67055e-05,1.87823e-08,-1.3481e-12,57888.5,29.7268], Tmin=(100,'K'), Tmax=(1068.26,'K')), NASAPolynomial(coeffs=[15.9033,0.0297892,-1.11423e-05,1.9713e-09,-1.34076e-13,53766.4,-50.8452], Tmin=(1068.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(480.096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=CCC=C=C=C(27626)',
    structure = SMILES('[CH2]C=CCC=C=C=C'),
    E0 = (485.381,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,540,563.333,586.667,610,1970,2140,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05167,'amu*angstrom^2'), symmetry=1, barrier=(24.1799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05193,'amu*angstrom^2'), symmetry=1, barrier=(24.186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05098,'amu*angstrom^2'), symmetry=1, barrier=(24.1642,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0835025,0.0747025,-5.67167e-05,2.18441e-08,-3.37243e-12,58528.3,30.7072], Tmin=(100,'K'), Tmax=(1536.1,'K')), NASAPolynomial(coeffs=[18.1543,0.027646,-1.07661e-05,1.90155e-09,-1.26776e-13,52976.6,-64.231], Tmin=(1536.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(485.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC[C]=CC=C(27627)',
    structure = SMILES('[CH2]C=CC[C]=CC=C'),
    E0 = (493.837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,273.828,273.828,273.828],'cm^-1')),
        HinderedRotor(inertia=(0.336441,'amu*angstrom^2'), symmetry=1, barrier=(17.9016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.336442,'amu*angstrom^2'), symmetry=1, barrier=(17.9016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.336441,'amu*angstrom^2'), symmetry=1, barrier=(17.9016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.336441,'amu*angstrom^2'), symmetry=1, barrier=(17.9016,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.209523,0.0708184,-3.27388e-05,-9.20246e-09,8.84488e-12,59542.3,31.6251], Tmin=(100,'K'), Tmax=(1015.44,'K')), NASAPolynomial(coeffs=[17.0403,0.0305112,-1.15928e-05,2.1221e-09,-1.49247e-13,54784.1,-56.4301], Tmin=(1015.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[C]=CCC=[C]C(27628)',
    structure = SMILES('C=C[C]=CCC=[C]C'),
    E0 = (541.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283233,0.0780811,-6.41802e-05,2.88981e-08,-5.38961e-12,65244.1,30.6758], Tmin=(100,'K'), Tmax=(1261.83,'K')), NASAPolynomial(coeffs=[13.5051,0.036168,-1.43563e-05,2.5746e-09,-1.74308e-13,61907.4,-36.1871], Tmin=(1261.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=CCC=C[C]=C(27629)',
    structure = SMILES('[CH2]C=CCC=C[C]=C'),
    E0 = (454.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.895016,'amu*angstrom^2'), symmetry=1, barrier=(20.5782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.89596,'amu*angstrom^2'), symmetry=1, barrier=(20.5999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895515,'amu*angstrom^2'), symmetry=1, barrier=(20.5897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.68982,'amu*angstrom^2'), symmetry=1, barrier=(84.8362,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.164636,0.0727761,-3.71724e-05,-5.30307e-09,7.85744e-12,54871,30.9083], Tmin=(100,'K'), Tmax=(997.04,'K')), NASAPolynomial(coeffs=[16.6052,0.0312354,-1.14101e-05,2.03301e-09,-1.40736e-13,50379,-54.4454], Tmin=(997.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=C[CH]C=CC=C(27630)',
    structure = SMILES('[CH2]C=CC=CC=C[CH2]'),
    E0 = (309.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.617249,0.0551877,2.34324e-05,-7.49889e-08,3.44961e-11,37394.1,28.8067], Tmin=(100,'K'), Tmax=(942.793,'K')), NASAPolynomial(coeffs=[18.6196,0.0268942,-8.05693e-06,1.37576e-09,-9.85656e-14,31862.5,-68.3175], Tmin=(942.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C[C]=CC[C]=CC(26586)',
    structure = SMILES('C=C[C]=CC[C]=CC'),
    E0 = (541.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.713501,'amu*angstrom^2'), symmetry=1, barrier=(16.4048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714239,'amu*angstrom^2'), symmetry=1, barrier=(16.4218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714825,'amu*angstrom^2'), symmetry=1, barrier=(16.4352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714532,'amu*angstrom^2'), symmetry=1, barrier=(16.4285,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283233,0.0780811,-6.41802e-05,2.88981e-08,-5.38961e-12,65244.1,30.6758], Tmin=(100,'K'), Tmax=(1261.83,'K')), NASAPolynomial(coeffs=[13.5051,0.036168,-1.43563e-05,2.5746e-09,-1.74308e-13,61907.4,-36.1871], Tmin=(1261.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC=CCC=C[CH2](27631)',
    structure = SMILES('[CH]=CC=CCC=C[CH2]'),
    E0 = (503.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.839344,'amu*angstrom^2'), symmetry=1, barrier=(19.2982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.839268,'amu*angstrom^2'), symmetry=1, barrier=(19.2964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.838324,'amu*angstrom^2'), symmetry=1, barrier=(19.2747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.838132,'amu*angstrom^2'), symmetry=1, barrier=(19.2703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.183585,0.0695769,-2.40899e-05,-2.16659e-08,1.40924e-11,60658.1,31.647], Tmin=(100,'K'), Tmax=(990.125,'K')), NASAPolynomial(coeffs=[18.3261,0.0284172,-1.04165e-05,1.90592e-09,-1.35646e-13,55490.3,-63.6543], Tmin=(990.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(503.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=[C]CC=CC=C(27632)',
    structure = SMILES('[CH2]C=[C]CC=CC=C'),
    E0 = (493.837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,273.009,274.056,276.281],'cm^-1')),
        HinderedRotor(inertia=(0.335716,'amu*angstrom^2'), symmetry=1, barrier=(17.8949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.337766,'amu*angstrom^2'), symmetry=1, barrier=(17.8956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343325,'amu*angstrom^2'), symmetry=1, barrier=(17.9158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333656,'amu*angstrom^2'), symmetry=1, barrier=(17.9001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.209523,0.0708184,-3.27388e-05,-9.20246e-09,8.84488e-12,59542.3,31.6251], Tmin=(100,'K'), Tmax=(1015.44,'K')), NASAPolynomial(coeffs=[17.0403,0.0305112,-1.15928e-05,2.1221e-09,-1.49247e-13,54784.1,-56.4301], Tmin=(1015.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[C]=C[CH]C=CC(27633)',
    structure = SMILES('[CH2]C=[C]C=CC=CC'),
    E0 = (390.693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242035,0.0708486,-2.97476e-05,-1.49244e-08,1.20942e-11,47135.5,28.8848], Tmin=(100,'K'), Tmax=(958.681,'K')), NASAPolynomial(coeffs=[16.4661,0.0308431,-1.04746e-05,1.79904e-09,-1.22891e-13,42752.4,-55.3383], Tmin=(958.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=CCC=CC=C(27634)',
    structure = SMILES('[CH2][C]=CCC=CC=C'),
    E0 = (493.837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,273.009,274.056,276.281],'cm^-1')),
        HinderedRotor(inertia=(0.335716,'amu*angstrom^2'), symmetry=1, barrier=(17.8949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.337766,'amu*angstrom^2'), symmetry=1, barrier=(17.8956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343325,'amu*angstrom^2'), symmetry=1, barrier=(17.9158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333656,'amu*angstrom^2'), symmetry=1, barrier=(17.9001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.209523,0.0708184,-3.27388e-05,-9.20246e-09,8.84488e-12,59542.3,31.6251], Tmin=(100,'K'), Tmax=(1015.44,'K')), NASAPolynomial(coeffs=[17.0403,0.0305112,-1.15928e-05,2.1221e-09,-1.49247e-13,54784.1,-56.4301], Tmin=(1015.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C[C]=[C]CC=CC(27635)',
    structure = SMILES('[CH2][CH]C#CCC=CC'),
    E0 = (508.888,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2100,2250,500,550,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493117,0.0699971,-5.01943e-05,2.02208e-08,-3.41467e-12,61337.1,33.2631], Tmin=(100,'K'), Tmax=(1378.71,'K')), NASAPolynomial(coeffs=[12.3229,0.0356754,-1.28527e-05,2.16428e-09,-1.4045e-13,58075.2,-27.6075], Tmin=(1378.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-(Cds-Cds)CtHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(RCCJ)"""),
)

species(
    label = 'C=[C][C]=CCC=CC(27636)',
    structure = SMILES('[CH2]C#C[CH]CC=CC'),
    E0 = (454.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.429812,0.0719829,-5.30289e-05,2.17637e-08,-3.73202e-12,54763.2,31.8213], Tmin=(100,'K'), Tmax=(1360.36,'K')), NASAPolynomial(coeffs=[12.7027,0.0358959,-1.32379e-05,2.26364e-09,-1.48434e-13,51424.1,-31.1655], Tmin=(1360.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Propargyl) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=C[C]=CCC=CC(27637)',
    structure = SMILES('[CH]C=C=CCC=CC'),
    E0 = (527.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.267132,0.0749054,-4.69852e-05,1.47664e-08,-1.90335e-12,63639.1,31.818], Tmin=(100,'K'), Tmax=(1746.79,'K')), NASAPolynomial(coeffs=[16.3192,0.0381473,-1.54201e-05,2.71944e-09,-1.79187e-13,58031.2,-54.5776], Tmin=(1746.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(527.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C[CH2](2461)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,229.711,230.18,230.787],'cm^-1')),
        HinderedRotor(inertia=(1.33306,'amu*angstrom^2'), symmetry=1, barrier=(50.5153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,44095,-3.44606], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
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
    label = 'C=C[C]=CCC1[CH]C1(27638)',
    structure = SMILES('C=C[C]=CCC1[CH]C1'),
    E0 = (566.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.666168,0.0592061,-2.90431e-06,-3.78868e-08,1.8715e-11,68225.1,30.3858], Tmin=(100,'K'), Tmax=(981.531,'K')), NASAPolynomial(coeffs=[15.3851,0.0320545,-1.15854e-05,2.08876e-09,-1.47017e-13,63754.1,-48.4065], Tmin=(981.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(566.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C=CJC=C) + radical(cyclopropane)"""),
)

species(
    label = '[CH2]C=CCC=C1[CH]C1(27639)',
    structure = SMILES('[CH2]C=CCC=C1[CH]C1'),
    E0 = (490.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.837321,0.0511543,2.52013e-05,-6.83116e-08,2.93136e-11,59079.4,27.9725], Tmin=(100,'K'), Tmax=(984.443,'K')), NASAPolynomial(coeffs=[16.6253,0.0311602,-1.16137e-05,2.18186e-09,-1.58815e-13,53831.3,-58.815], Tmin=(984.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1[CH]CC=C1C=C(27292)',
    structure = SMILES('[CH2]C1[CH]CC=C1C=C'),
    E0 = (432.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00111,0.0465303,3.85492e-05,-8.48524e-08,3.64854e-11,52147.6,28.2938], Tmin=(100,'K'), Tmax=(953.364,'K')), NASAPolynomial(coeffs=[16.7752,0.0286105,-9.19274e-06,1.63318e-09,-1.18228e-13,46946.6,-58.5571], Tmin=(953.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentene) + radical(cyclopentene-4) + radical(Isobutyl)"""),
)

species(
    label = '[C]1[CH]CCC=CCC=1(27640)',
    structure = SMILES('[C]1[CH]CC=CCCC=1'),
    E0 = (456.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52092,0.034587,6.35308e-05,-9.97166e-08,3.82778e-11,54952.9,20.7636], Tmin=(100,'K'), Tmax=(991.464,'K')), NASAPolynomial(coeffs=[13.2416,0.0366199,-1.41602e-05,2.69517e-09,-1.96481e-13,50204.8,-47.9058], Tmin=(991.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-cyclooctadiene) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'C=C=CCC=CC=C(27641)',
    structure = SMILES('C=C=CCC=CC=C'),
    E0 = (281.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.289511,0.0682706,-2.5263e-05,-1.70266e-08,1.16099e-11,33953.9,29.6373], Tmin=(100,'K'), Tmax=(1007.73,'K')), NASAPolynomial(coeffs=[17.0895,0.0302914,-1.14586e-05,2.1072e-09,-1.49125e-13,29110.4,-58.774], Tmin=(1007.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C=C=CCC=CC(27642)',
    structure = SMILES('C=C=C=CCC=CC'),
    E0 = (333.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449316,0.0731152,-5.19279e-05,1.90028e-08,-2.86436e-12,40288.2,29.0972], Tmin=(100,'K'), Tmax=(1522.24,'K')), NASAPolynomial(coeffs=[14.861,0.0352455,-1.46114e-05,2.65992e-09,-1.80324e-13,35900.6,-46.4862], Tmin=(1522.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC1=CCC=CC1(26573)',
    structure = SMILES('C=CC1=CCC=CC1'),
    E0 = (144.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03741,0.0417415,5.93487e-05,-1.08894e-07,4.5087e-11,17450,24.7165], Tmin=(100,'K'), Tmax=(960.56,'K')), NASAPolynomial(coeffs=[18.5174,0.0272195,-8.96643e-06,1.67155e-09,-1.25797e-13,11403.7,-72.9039], Tmin=(960.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,4-Cyclohexadiene)"""),
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
    label = '[CH]=CCC=[C]C=C(27643)',
    structure = SMILES('[CH]=CCC=[C]C=C'),
    E0 = (586.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.976949,'amu*angstrom^2'), symmetry=1, barrier=(22.462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.972333,'amu*angstrom^2'), symmetry=1, barrier=(22.3558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978992,'amu*angstrom^2'), symmetry=1, barrier=(22.5089,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.757375,0.0629921,-4.10763e-05,6.13658e-09,3.00445e-12,70677.4,27.0971], Tmin=(100,'K'), Tmax=(995.53,'K')), NASAPolynomial(coeffs=[14.4388,0.0248332,-8.91303e-06,1.56194e-09,-1.06757e-13,67120.2,-43.0309], Tmin=(995.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C1[C]=CCC1C=C(27315)',
    structure = SMILES('[CH2]C1[C]=CCC1C=C'),
    E0 = (525.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1824,0.0423047,4.72432e-05,-9.13014e-08,3.80935e-11,63272.2,28.3639], Tmin=(100,'K'), Tmax=(957.359,'K')), NASAPolynomial(coeffs=[15.8323,0.029936,-9.90199e-06,1.7809e-09,-1.29113e-13,58229,-53.3645], Tmin=(957.359,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopentene) + radical(cyclopentene-vinyl) + radical(Isobutyl)"""),
)

species(
    label = 'C=C[C]=CC=CC=C(27644)',
    structure = SMILES('C=C[C]=CC=CC=C'),
    E0 = (396.461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23226,'amu*angstrom^2'), symmetry=1, barrier=(28.3321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22615,'amu*angstrom^2'), symmetry=1, barrier=(28.1917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22442,'amu*angstrom^2'), symmetry=1, barrier=(28.1519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0129934,0.0724339,-2.73583e-05,-2.81833e-08,1.98046e-11,47840.8,27.1733], Tmin=(100,'K'), Tmax=(930.345,'K')), NASAPolynomial(coeffs=[21.3147,0.0197366,-5.09495e-06,7.93318e-10,-5.54065e-14,42194.2,-83.1021], Tmin=(930.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=C[CH]CC=C(27262)',
    structure = SMILES('[CH2]C=[C]C=CCC=C'),
    E0 = (421.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.453282,0.0662759,-2.21848e-05,-1.82024e-08,1.19665e-11,50838.5,31.1494], Tmin=(100,'K'), Tmax=(981.824,'K')), NASAPolynomial(coeffs=[14.9947,0.0333072,-1.19566e-05,2.10816e-09,-1.45196e-13,46716.8,-45.1869], Tmin=(981.824,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]CCC=[C]C=C(27645)',
    structure = SMILES('C=[C]CCC=[C]C=C'),
    E0 = (552.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1670,1700,300,440,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,257.3,257.3,257.304,257.306],'cm^-1')),
        HinderedRotor(inertia=(0.318575,'amu*angstrom^2'), symmetry=1, barrier=(14.9663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.318533,'amu*angstrom^2'), symmetry=1, barrier=(14.9663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.318571,'amu*angstrom^2'), symmetry=1, barrier=(14.9663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.54978,'amu*angstrom^2'), symmetry=1, barrier=(25.8292,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0960851,0.0795927,-6.55705e-05,2.90772e-08,-5.27044e-12,66573.1,31.7607], Tmin=(100,'K'), Tmax=(1308.6,'K')), NASAPolynomial(coeffs=[15.0997,0.0337312,-1.30014e-05,2.29595e-09,-1.54052e-13,62646.4,-44.6585], Tmin=(1308.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(552.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=[C]CCC=C(27646)',
    structure = SMILES('[CH2][CH]C#CCCC=C'),
    E0 = (516.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.609026,0.0678373,-3.70864e-05,8.00153e-11,5.80484e-12,62259.3,33.9325], Tmin=(100,'K'), Tmax=(925.481,'K')), NASAPolynomial(coeffs=[12.185,0.0350155,-1.1784e-05,1.94712e-09,-1.27409e-13,59379.6,-25.0004], Tmin=(925.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(RCCJ) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=CCCC=[C]C=C(27647)',
    structure = SMILES('[CH]=CCCC=[C]C=C'),
    E0 = (561.564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.81523,'amu*angstrom^2'), symmetry=1, barrier=(18.7438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.810812,'amu*angstrom^2'), symmetry=1, barrier=(18.6422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812241,'amu*angstrom^2'), symmetry=1, barrier=(18.675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.815087,'amu*angstrom^2'), symmetry=1, barrier=(18.7405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.199787,0.0814171,-6.71185e-05,2.91326e-08,-5.07793e-12,67700.8,32.7591], Tmin=(100,'K'), Tmax=(1376.78,'K')), NASAPolynomial(coeffs=[17.6138,0.0296629,-1.07327e-05,1.82954e-09,-1.2018e-13,62795.7,-58.8773], Tmin=(1376.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][C]=CCCC=C(27648)',
    structure = SMILES('[CH2]C#C[CH]CCC=C'),
    E0 = (466.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202359,0.0738064,-5.51891e-05,2.26058e-08,-3.80456e-12,56246.9,33.0817], Tmin=(100,'K'), Tmax=(1403.94,'K')), NASAPolynomial(coeffs=[14.4964,0.033081,-1.16773e-05,1.94401e-09,-1.25307e-13,52233.3,-40.7283], Tmin=(1403.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C[CH]CC=CC=C(27649)',
    structure = SMILES('[CH]C=CCC=CC=C'),
    E0 = (475.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,2950,3100,1380,975,1025,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.245066,0.0686045,-1.59733e-05,-2.59288e-08,1.41892e-11,57298.4,31.8523], Tmin=(100,'K'), Tmax=(1011.88,'K')), NASAPolynomial(coeffs=[15.8682,0.0371762,-1.43465e-05,2.62243e-09,-1.83671e-13,52583.9,-51.3777], Tmin=(1011.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C[C]=CCCC=C(27650)',
    structure = SMILES('[CH]C=C=CCCC=C'),
    E0 = (538.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0280778,0.0769032,-4.96115e-05,1.60727e-08,-2.1199e-12,64970.8,33.0993], Tmin=(100,'K'), Tmax=(1736.16,'K')), NASAPolynomial(coeffs=[17.8162,0.0359203,-1.4203e-05,2.4761e-09,-1.62024e-13,58794.3,-62.5311], Tmin=(1736.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC1=CC[CH][CH]C1(27621)',
    structure = SMILES('C=CC1=CC[CH][CH]C1'),
    E0 = (417.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43573,0.0430629,2.63301e-05,-5.44015e-08,2.10207e-11,50360.2,26.9744], Tmin=(100,'K'), Tmax=(1027.29,'K')), NASAPolynomial(coeffs=[10.1508,0.0406295,-1.61129e-05,2.99193e-09,-2.1076e-13,46907.4,-23.3953], Tmin=(1027.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(417.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclohexene) + radical(RCCJCC) + radical(RCCJCC)"""),
)

species(
    label = 'C=CC1C[CH][C]=CC1(27588)',
    structure = SMILES('C=CC1C[CH][C]=CC1'),
    E0 = (427.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46503,0.0333903,7.27858e-05,-1.13609e-07,4.42179e-11,51518.6,23.0074], Tmin=(100,'K'), Tmax=(979.574,'K')), NASAPolynomial(coeffs=[15.062,0.0334784,-1.25035e-05,2.38982e-09,-1.76867e-13,46186.7,-55.9282], Tmin=(979.574,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclohexene) + radical(cyclohexene-allyl) + radical(Cds_S)"""),
)

species(
    label = 'C=CC=CC=CC=C(27651)',
    structure = SMILES('C=CC=CC=CC=C'),
    E0 = (197.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0864509,0.0659507,1.49008e-06,-6.06862e-08,3.13845e-11,23908.9,26.6193], Tmin=(100,'K'), Tmax=(941.384,'K')), NASAPolynomial(coeffs=[22.8418,0.0196705,-5.08803e-06,8.53682e-10,-6.42382e-14,17391,-93.6509], Tmin=(941.384,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C=CCCC=C(27652)',
    structure = SMILES('C=C=C=CCCC=C'),
    E0 = (344.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.22411,0.0749939,-5.42882e-05,2.01008e-08,-3.0283e-12,41619.2,30.325], Tmin=(100,'K'), Tmax=(1544.87,'K')), NASAPolynomial(coeffs=[16.5634,0.0326881,-1.32112e-05,2.37456e-09,-1.59752e-13,36570.7,-55.6093], Tmin=(1544.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]C(C=C)C=[C]C=C(26570)',
    structure = SMILES('[CH2]C(C=C)C=[C]C=C'),
    E0 = (511.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.550407,'amu*angstrom^2'), symmetry=1, barrier=(12.6549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.548943,'amu*angstrom^2'), symmetry=1, barrier=(12.6213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.14273,'amu*angstrom^2'), symmetry=1, barrier=(72.2574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05201,'amu*angstrom^2'), symmetry=1, barrier=(24.1878,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3606.27,'J/mol'), sigma=(6.17911,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.29 K, Pc=34.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00957748,0.0846348,-7.96538e-05,4.20261e-08,-9.04714e-12,61610.1,30.815], Tmin=(100,'K'), Tmax=(1118.3,'K')), NASAPolynomial(coeffs=[13.9114,0.034841,-1.28637e-05,2.20926e-09,-1.45848e-13,58496.5,-37.9018], Tmin=(1118.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(511.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC1=CCC1C=C(27296)',
    structure = SMILES('C=CC1=CCC1C=C'),
    E0 = (266.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.583793,0.0572995,1.13794e-05,-5.81862e-08,2.69995e-11,32176.7,24.6308], Tmin=(100,'K'), Tmax=(971.368,'K')), NASAPolynomial(coeffs=[18.0582,0.0281974,-9.8596e-06,1.81029e-09,-1.31393e-13,26760,-69.5729], Tmin=(971.368,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH2]C=[C]C1CC=CC1(27653)',
    structure = SMILES('[CH2]C=[C]C1CC=CC1'),
    E0 = (447.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17342,0.0414409,5.10771e-05,-9.37615e-08,3.79951e-11,53949.6,27.3957], Tmin=(100,'K'), Tmax=(977.676,'K')), NASAPolynomial(coeffs=[16.1049,0.0312802,-1.14713e-05,2.1707e-09,-1.59939e-13,48596,-56.7507], Tmin=(977.676,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CCC=C=[C]C(27654)',
    structure = SMILES('[CH2]C=CC[CH]C#CC'),
    E0 = (449.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.534286,0.0681602,-3.29532e-05,-6.76743e-09,8.78063e-12,54171.2,30.6693], Tmin=(100,'K'), Tmax=(923.946,'K')), NASAPolynomial(coeffs=[13.2024,0.0340489,-1.12325e-05,1.84575e-09,-1.21087e-13,50945.4,-34.2335], Tmin=(923.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Allyl_P) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH2]C=CC[C]=C=CC(27655)',
    structure = SMILES('C=C[CH]CC#C[CH]C'),
    E0 = (452.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.764482,0.0619466,-1.79301e-05,-2.03206e-08,1.30635e-11,54542.4,29.788], Tmin=(100,'K'), Tmax=(930.085,'K')), NASAPolynomial(coeffs=[12.4952,0.0347208,-1.14766e-05,1.90078e-09,-1.25789e-13,51355.8,-31.356], Tmin=(930.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.453,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C=C[CH]C=C=CC(27656)',
    structure = SMILES('[CH2]C=CC=C[C]=CC'),
    E0 = (390.693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.241999,0.070849,-2.97492e-05,-1.49222e-08,1.20932e-11,47135.5,28.8849], Tmin=(100,'K'), Tmax=(958.69,'K')), NASAPolynomial(coeffs=[16.4662,0.0308428,-1.04745e-05,1.79901e-09,-1.22889e-13,42752.3,-55.339], Tmin=(958.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=[C]CC=C=CC(27657)',
    structure = SMILES('[CH2]C=[C]CC=C=CC'),
    E0 = (546.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.328597,'amu*angstrom^2'), symmetry=1, barrier=(7.55509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243112,'amu*angstrom^2'), symmetry=1, barrier=(17.4589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.759703,'amu*angstrom^2'), symmetry=1, barrier=(17.4671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75971,'amu*angstrom^2'), symmetry=1, barrier=(17.4672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.503535,0.0741618,-5.45351e-05,2.10444e-08,-3.38084e-12,65870.7,30.5981], Tmin=(100,'K'), Tmax=(1425.5,'K')), NASAPolynomial(coeffs=[13.6464,0.0372832,-1.57298e-05,2.89665e-09,-1.982e-13,62123.6,-37.4685], Tmin=(1425.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC=C=CC(27658)',
    structure = SMILES('[CH2][C]=CCC=C=CC'),
    E0 = (546.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.328597,'amu*angstrom^2'), symmetry=1, barrier=(7.55509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243112,'amu*angstrom^2'), symmetry=1, barrier=(17.4589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.759703,'amu*angstrom^2'), symmetry=1, barrier=(17.4671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75971,'amu*angstrom^2'), symmetry=1, barrier=(17.4672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.503535,0.0741618,-5.45351e-05,2.10444e-08,-3.38084e-12,65870.7,30.5981], Tmin=(100,'K'), Tmax=(1425.5,'K')), NASAPolynomial(coeffs=[13.6464,0.0372832,-1.57298e-05,2.89665e-09,-1.982e-13,62123.6,-37.4685], Tmin=(1425.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CCC1[C]=CC1(27659)',
    structure = SMILES('[CH2]C=CCC1[C]=CC1'),
    E0 = (560.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.884113,0.050552,2.46843e-05,-6.67919e-08,2.85974e-11,67588.4,29.314], Tmin=(100,'K'), Tmax=(985.455,'K')), NASAPolynomial(coeffs=[16.222,0.0313793,-1.17129e-05,2.19679e-09,-1.59518e-13,62473.5,-55.0723], Tmin=(985.455,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(560.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Allyl_P) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH2]C1[CH]CC=C=CC1(27599)',
    structure = SMILES('[CH2]C1[CH]CC=C=CC1'),
    E0 = (563.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2584,0.048066,1.30348e-05,-4.07246e-08,1.6134e-11,67848.8,25.1817], Tmin=(100,'K'), Tmax=(1052.97,'K')), NASAPolynomial(coeffs=[10.4149,0.0410669,-1.65745e-05,3.08104e-09,-2.161e-13,64380.2,-26.7798], Tmin=(1052.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(563.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = 'C=C=CCC=C=CC(27660)',
    structure = SMILES('C=C=CCC=C=CC'),
    E0 = (333.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449316,0.0731152,-5.19279e-05,1.90028e-08,-2.86436e-12,40288.2,29.0972], Tmin=(100,'K'), Tmax=(1522.24,'K')), NASAPolynomial(coeffs=[14.861,0.0352455,-1.46114e-05,2.65992e-09,-1.80324e-13,35900.6,-46.4862], Tmin=(1522.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C#CC([CH2])C([CH2])C=C(26569)',
    structure = SMILES('C#CC([CH2])C([CH2])C=C'),
    E0 = (575.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,750,770,3400,2100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,335.718,335.875],'cm^-1')),
        HinderedRotor(inertia=(3.70374e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.904772,'amu*angstrom^2'), symmetry=1, barrier=(72.8523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.906402,'amu*angstrom^2'), symmetry=1, barrier=(72.8519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128953,'amu*angstrom^2'), symmetry=1, barrier=(10.327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.911175,'amu*angstrom^2'), symmetry=1, barrier=(72.8883,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3560.51,'J/mol'), sigma=(6.30174,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=556.14 K, Pc=32.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.297611,0.0815261,-7.15621e-05,3.45681e-08,-6.60905e-12,69359.6,35.8367], Tmin=(100,'K'), Tmax=(1387.22,'K')), NASAPolynomial(coeffs=[16.3164,0.0281223,-7.87185e-06,1.1031e-09,-6.32268e-14,65279.1,-47.8477], Tmin=(1387.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C#C[CH]C[CH]C=C(27592)',
    structure = SMILES('C#C[CH]CC=C[CH2]'),
    E0 = (491.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,433.817,438.145],'cm^-1')),
        HinderedRotor(inertia=(0.0632098,'amu*angstrom^2'), symmetry=1, barrier=(8.75477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.643662,'amu*angstrom^2'), symmetry=1, barrier=(84.8187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608377,'amu*angstrom^2'), symmetry=1, barrier=(84.7284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658565,'amu*angstrom^2'), symmetry=1, barrier=(84.9366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924384,0.0603159,-3.48015e-05,-6.7064e-10,6.04832e-12,59229.7,25.848], Tmin=(100,'K'), Tmax=(934.069,'K')), NASAPolynomial(coeffs=[13.2645,0.0256867,-8.44297e-06,1.39437e-09,-9.22156e-14,56129.7,-37.0982], Tmin=(934.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=[C]C1CC1C=C(27265)',
    structure = SMILES('[CH2]C=[C]C1CC1C=C'),
    E0 = (539.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.832513,0.0496259,3.27554e-05,-7.9885e-08,3.45454e-11,65065.6,29.5846], Tmin=(100,'K'), Tmax=(966.979,'K')), NASAPolynomial(coeffs=[17.9491,0.0278758,-9.59988e-06,1.7781e-09,-1.30944e-13,59461.9,-64.2771], Tmin=(966.979,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(539.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C[CH]CC=C=CC(27661)',
    structure = SMILES('[CH]C=CCC=C=CC'),
    E0 = (527.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.267132,0.0749054,-4.69852e-05,1.47664e-08,-1.90335e-12,63639.1,31.818], Tmin=(100,'K'), Tmax=(1746.79,'K')), NASAPolynomial(coeffs=[16.3192,0.0381473,-1.54201e-05,2.71944e-09,-1.79187e-13,58031.2,-54.5776], Tmin=(1746.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(527.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC=CC=C=CC(27662)',
    structure = SMILES('C=CC=CC=C=CC'),
    E0 = (250.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0730811,0.0729175,-3.31534e-05,-1.32141e-08,1.1456e-11,30251,27.3905], Tmin=(100,'K'), Tmax=(986.052,'K')), NASAPolynomial(coeffs=[18.5254,0.0278831,-1.00069e-05,1.80491e-09,-1.27367e-13,25162.4,-68.7225], Tmin=(986.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    E0 = (454.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (538.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (582.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (651.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (703.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (713.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (591.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (690.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (662.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (669.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (602.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (703.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (725.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (888.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (572.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (567.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (570.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (643.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (660.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (751.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (726.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (680.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (592.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (493.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (512.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (494.216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (479.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (463.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (968.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (563.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (616.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (552.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (571.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (709.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (641.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (706.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (584.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (579.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (583.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (467.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (518.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (548.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (494.216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (670.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (462.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (534.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (662.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (693.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (514.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (689.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (627.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (573.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (563.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (479.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (775.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (873.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (563.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (660.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (472.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['butadiene13(2459)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2][CH]C1CC=C1C=C(27622)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(83.9398,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 81.3 to 83.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2]C1[C]=CCC=CC1(27623)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.78819e+09,'s^-1'), n=0.5525, Ea=(127.194,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra_2H_pri;radadd_intra_cs2H] for rate rule [R8;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C=CCC#CC=C(27624)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.66e+08,'cm^3/(mol*s)'), n=1.64, Ea=(12.2591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2703 used for Ct-Cs_Ct-Cd;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=C=CCC=[C]C=C(27625)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH2]C=CCC=C=C=C(27626)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH]C=C(2458)', 'CH2CHCCH(26391)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(45800,'cm^3/(mol*s)'), n=2.41, Ea=(42.8442,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2293 used for Ct-H_Ct-Cd;CsJ-CdHH
Exact match found for rate rule [Ct-H_Ct-Cd;CsJ-CdHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C=CC[C]=CC=C(27627)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=C[C]=CCC=[C]C(27628)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C=CCC=C[C]=C(27629)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2]C=C[CH]C=CC=C(27630)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[C]=CC[C]=CC(26586)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC=CCC=C[CH2](27631)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C=[C]CC=CC=C(27632)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.01194e+11,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_Cd;Cd_H_out_singleDe]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=C[C]=C[CH]C=CC(27633)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(256000,'s^-1'), n=2, Ea=(117.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C]=CCC=CC=C(27634)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(272058,'s^-1'), n=1.7475, Ea=(73.8225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cd_H_out_single] for rate rule [R5H_DSSD;Cd_rad_out;Cd_H_out_singleDe]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C[C]=[C]CC=CC(27635)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(770185,'s^-1'), n=1.81245, Ea=(61.132,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_RSMS;Cd_rad_out;Cs_H_out_2H] + [R5H_SSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_SSMS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=[C][C]=CCC=CC(27636)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.11562e+07,'s^-1'), n=1.49763, Ea=(188.308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_2H;Cd_H_out_doubleC] for rate rule [R7HJ_5;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C[C]=CCC=CC(27637)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R8Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C[CH2](2461)', '[CH2]C=[C]C=C(26613)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.33799e+08,'m^3/(mol*s)'), n=-0.455312, Ea=(0.377199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;C_rad/H2/Cd] + [Cd_rad;C_pri_rad] for rate rule [Cd_rad;C_rad/H2/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C=C(2458)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.13324e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H2/Cd;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=C[C]=CCC1[CH]C1(27638)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2]C=CCC=C1[CH]C1(27639)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2]C1[CH]CC=C1C=C(27292)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.80657e+08,'s^-1'), n=0.835, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra_pri;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra_pri;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[C]1[CH]CCC=CCC=1(27640)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.04798e+07,'s^-1'), n=0.841667, Ea=(57.3368,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_2H;radadd_intra_cs2H] for rate rule [R8_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=C=CCC=CC=C(27641)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=C=C=CCC=CC(27642)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=CC1=CCC=CC1(26573)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Cpri_rad_out_2H] for rate rule [R6;CdsingleDe_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(19)', '[CH]=CCC=[C]C=C(27643)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2]C1[C]=CCC1C=C(27315)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Rn;doublebond_intra_2H_pri;radadd_intra_csHCd] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_csHCd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', 'C=C[C]=CC=CC=C(27644)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(105.863,'m^3/(mol*s)'), n=1.66278, Ea=(8.08939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['butadiene13(2459)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.0534234,'m^3/(mol*s)'), n=2.459, Ea=(4.91982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=C[C]=C[CH]CC=C(27262)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.169e+11,'s^-1'), n=0.707, Ea=(116.068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C]CCC=[C]C=C(27645)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C[C]=[C]CCC=C(27646)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=CCCC=[C]C=C(27647)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=[C][C]=CCCC=C(27648)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.17074e+10,'s^-1'), n=0.96, Ea=(118.093,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_H/Cd] for rate rule [R5HJ_1;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C[CH]CC=CC=C(27649)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(7.80899e+07,'s^-1'), n=1.452, Ea=(104.082,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleDe] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=C[C]=CCCC=C(27650)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=CC1=CC[CH][CH]C1(27621)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.16959e+10,'s^-1'), n=0.31, Ea=(12.393,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra_cdsingle] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=CC1C[CH][C]=CC1(27588)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.0689e+08,'s^-1'), n=0.637531, Ea=(63.1312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=CC=CC=CC=C(27651)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=C=C=CCCC=C(27652)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C=C)C=[C]C=C(26570)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=CC1=CCC1C=C(27296)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriDe_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2]C=[C]C1CC=CC1(27653)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(79.3435,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_SMS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2]C=CCC=C=[C]C(27654)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2]C=CC[C]=C=CC(27655)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(0.0029655,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2]C=C[CH]C=C=CC(27656)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(62296.1,'s^-1'), n=1.86, Ea=(59.4128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;C_rad_out_2H;Cs_H_out_H/Cd] for rate rule [R5H_SMMS;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C=[C]CC=C=CC(27657)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6H;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2][C]=CCC=C=CC(27658)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(244756,'s^-1'), n=1.235, Ea=(80.8557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;Y_rad_out;Cs_H_out_2H] for rate rule [R7H;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2]C=CCC1[C]=CC1(27659)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(118.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2]C1[CH]CC=C=CC1(27599)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(163933,'s^-1'), n=1.17455, Ea=(108.229,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7_linear;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 104.6 to 108.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction54',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=C=CCC=C=CC(27660)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#CC([CH2])C([CH2])C=C(26569)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(2.214e+09,'s^-1'), n=0.749, Ea=(200.242,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 2 used for hex_1_ene_5_yne
Exact match found for rate rule [hex_1_ene_5_yne]
Euclidian distance = 0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction56',
    reactants = ['CH2(19)', 'C#C[CH]C[CH]C=C(27592)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction57',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['[CH2]C=[C]C1CC1C=C(27265)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH]=C[CH]CC=C=CC(27661)'],
    products = ['C=C[C]=CC[CH]C=C(26566)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R8Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['C=C[C]=CC[CH]C=C(26566)'],
    products = ['C=CC=CC=C=CC(27662)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

network(
    label = '4536',
    isomers = [
        'C=C[C]=CC[CH]C=C(26566)',
    ],
    reactants = [
        ('butadiene13(2459)', 'CH2CHCCH(26391)'),
        ('[CH2][CH]C=C(2458)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4536',
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

