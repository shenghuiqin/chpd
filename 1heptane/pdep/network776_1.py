species(
    label = 'C[CH]OO[CH]C=O(1192)',
    structure = SMILES('C[CH]OO[CH]C=O'),
    E0 = (35.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,308.229],'cm^-1')),
        HinderedRotor(inertia=(0.00341453,'amu*angstrom^2'), symmetry=1, barrier=(9.45894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.660503,'amu*angstrom^2'), symmetry=1, barrier=(44.7411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.662622,'amu*angstrom^2'), symmetry=1, barrier=(44.7368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.661317,'amu*angstrom^2'), symmetry=1, barrier=(44.7397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.662212,'amu*angstrom^2'), symmetry=1, barrier=(44.7375,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21585,0.0671781,-7.38746e-05,4.89805e-08,-1.41378e-11,4382.38,24.9422], Tmin=(100,'K'), Tmax=(817.093,'K')), NASAPolynomial(coeffs=[7.40276,0.036891,-1.82747e-05,3.61693e-09,-2.58386e-13,3371.31,-3.65632], Tmin=(817.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(CCsJOOC)"""),
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
    label = 'CH3CHO(37)',
    structure = SMILES('CC=O'),
    E0 = (-178.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1305.64,1305.66,1305.67,3976.84],'cm^-1')),
        HinderedRotor(inertia=(0.136163,'amu*angstrom^2'), symmetry=1, barrier=(3.13064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-178.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CC1OO[CH]C1[O](4459)',
    structure = SMILES('CC1OO[CH]C1[O]'),
    E0 = (83.4875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48856,0.0412275,1.89331e-05,-6.49233e-08,3.23355e-11,10144.9,20.9446], Tmin=(100,'K'), Tmax=(885.351,'K')), NASAPolynomial(coeffs=[16.0514,0.0143163,-1.35044e-06,-4.41421e-11,7.95149e-15,6042.36,-56.1457], Tmin=(885.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.4875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(CC(C)OJ) + radical(CCsJOOC)"""),
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
    label = 'C=COO[CH]C=O(4460)',
    structure = SMILES('C=COO[CH]C=O'),
    E0 = (-7.87623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.981489,'amu*angstrom^2'), symmetry=1, barrier=(22.5664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23905,'amu*angstrom^2'), symmetry=1, barrier=(28.4882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.53444,'amu*angstrom^2'), symmetry=1, barrier=(58.2717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0921527,'amu*angstrom^2'), symmetry=1, barrier=(28.4869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14262,0.0596679,-5.28045e-05,2.34668e-08,-4.1925e-12,-841.914,25.5714], Tmin=(100,'K'), Tmax=(1330.34,'K')), NASAPolynomial(coeffs=[13.8682,0.0214049,-9.66146e-06,1.84657e-09,-1.29547e-13,-4227.76,-39.4545], Tmin=(1330.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.87623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(OCJC=O)"""),
)

species(
    label = 'C[CH]OOC=C=O(4461)',
    structure = SMILES('C[CH]OOC=C=O'),
    E0 = (79.3517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2120,512.5,787.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.164091,'amu*angstrom^2'), symmetry=1, barrier=(11.6682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.506886,'amu*angstrom^2'), symmetry=1, barrier=(11.6543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47483,'amu*angstrom^2'), symmetry=1, barrier=(33.9092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.388524,'amu*angstrom^2'), symmetry=1, barrier=(67.0507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20384,0.0649858,-7.32633e-05,4.52659e-08,-1.15086e-11,9641.57,24.3551], Tmin=(100,'K'), Tmax=(944.426,'K')), NASAPolynomial(coeffs=[10.2208,0.0267941,-1.26027e-05,2.44443e-09,-1.72917e-13,7938.45,-18.6308], Tmin=(944.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.3517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + radical(CCsJOOC)"""),
)

species(
    label = '[CH2]COO[CH]C=O(4462)',
    structure = SMILES('[CH2]COO[CH]C=O'),
    E0 = (63.1688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15586,0.0679361,-7.2651e-05,4.41696e-08,-1.14184e-11,7695.32,26.8553], Tmin=(100,'K'), Tmax=(914.038,'K')), NASAPolynomial(coeffs=[8.96261,0.0337725,-1.65866e-05,3.27849e-09,-2.34312e-13,6268.18,-10.1061], Tmin=(914.038,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.1688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(CJCOOH)"""),
)

species(
    label = 'C[CH]OOC[C]=O(4463)',
    structure = SMILES('C[CH]OOC[C]=O'),
    E0 = (51.8331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,2750,2850,1437.5,1250,1305,750,350,251.124],'cm^-1')),
        HinderedRotor(inertia=(0.769382,'amu*angstrom^2'), symmetry=1, barrier=(34.5928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00514394,'amu*angstrom^2'), symmetry=1, barrier=(34.6105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.69955,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268738,'amu*angstrom^2'), symmetry=1, barrier=(11.7141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786156,'amu*angstrom^2'), symmetry=1, barrier=(34.5959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36149,0.0637316,-6.41439e-05,3.7367e-08,-9.46713e-12,6324.29,25.3403], Tmin=(100,'K'), Tmax=(919.801,'K')), NASAPolynomial(coeffs=[7.96188,0.0350288,-1.73369e-05,3.44239e-09,-2.46734e-13,5110.05,-5.95115], Tmin=(919.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.8331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(CCsJOOC) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2][CH]OOCC=O(4464)',
    structure = SMILES('[CH2][CH]OOCC=O'),
    E0 = (111.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.214307,'amu*angstrom^2'), symmetry=1, barrier=(16.262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72895,'amu*angstrom^2'), symmetry=1, barrier=(39.752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103232,'amu*angstrom^2'), symmetry=1, barrier=(2.37352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0116573,'amu*angstrom^2'), symmetry=1, barrier=(39.7562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.484351,'amu*angstrom^2'), symmetry=1, barrier=(39.7518,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21176,0.0665203,-6.94407e-05,4.08094e-08,-1.01918e-11,13461.3,27.1918], Tmin=(100,'K'), Tmax=(943.57,'K')), NASAPolynomial(coeffs=[9.10743,0.0330478,-1.62275e-05,3.21119e-09,-2.29746e-13,11971.3,-10.4414], Tmin=(943.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCOOH) + radical(CCsJOOC)"""),
)

species(
    label = 'CCOO[CH][C]=O(4465)',
    structure = SMILES('CCOO[CH][C]=O'),
    E0 = (3.8772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29761,0.0652428,-6.76935e-05,4.11782e-08,-1.08931e-11,558.682,25.0323], Tmin=(100,'K'), Tmax=(884.924,'K')), NASAPolynomial(coeffs=[7.8239,0.0357433,-1.76907e-05,3.50856e-09,-2.51215e-13,-596.386,-5.65545], Tmin=(884.924,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.8772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'C[CH][O](176)',
    structure = SMILES('C[CH][O]'),
    E0 = (157.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1642.51],'cm^-1')),
        HinderedRotor(inertia=(0.123965,'amu*angstrom^2'), symmetry=1, barrier=(2.85019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65562,0.0114444,2.34936e-06,-4.83164e-09,1.17966e-12,18963.9,10.3625], Tmin=(100,'K'), Tmax=(1718.65,'K')), NASAPolynomial(coeffs=[6.06294,0.0136322,-6.35953e-06,1.18407e-09,-7.90642e-14,16985.9,-5.90233], Tmin=(1718.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CCsJOH) + radical(CCOJ)"""),
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
    label = 'C[CH]OOC1[CH]O1(4466)',
    structure = SMILES('C[CH]OOC1[CH]O1'),
    E0 = (136.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.727301,0.0678351,-6.1343e-05,2.38368e-08,-1.26744e-12,16575.8,25.1134], Tmin=(100,'K'), Tmax=(881.048,'K')), NASAPolynomial(coeffs=[14.1466,0.0213196,-6.68085e-06,1.03751e-09,-6.51826e-14,13651.9,-41.1013], Tmin=(881.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCsJOOC) + radical(CCsJO)"""),
)

species(
    label = 'CC1O[CH][CH]OO1(4443)',
    structure = SMILES('CC1O[CH][CH]OO1'),
    E0 = (39.4351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14305,0.0452353,2.61745e-05,-8.99987e-08,4.72924e-11,4862.89,20.5244], Tmin=(100,'K'), Tmax=(853.028,'K')), NASAPolynomial(coeffs=[21.2185,0.00390968,5.97681e-06,-1.63571e-09,1.24922e-13,-483.529,-84.3992], Tmin=(853.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.4351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(124trioxane) + radical(CCsJOOC) + radical(CCsJOCs)"""),
)

species(
    label = 'C=COOCC=O(4467)',
    structure = SMILES('C=COOCC=O'),
    E0 = (-146.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.802186,0.060749,-4.95929e-05,1.98722e-08,-3.14525e-12,-17480.1,27.0517], Tmin=(100,'K'), Tmax=(1514.52,'K')), NASAPolynomial(coeffs=[16.8505,0.0183636,-7.61382e-06,1.39364e-09,-9.50122e-14,-22341.2,-57.0338], Tmin=(1514.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CCOOC=C=O(4468)',
    structure = SMILES('CCOOC=C=O'),
    E0 = (-107.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15415,0.0637473,-6.17147e-05,3.19383e-08,-6.79891e-12,-12778.3,24.2304], Tmin=(100,'K'), Tmax=(1115.52,'K')), NASAPolynomial(coeffs=[11.4276,0.0269091,-1.21796e-05,2.33461e-09,-1.64408e-13,-15070.3,-26.456], Tmin=(1115.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-107.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH)"""),
)

species(
    label = 'CC1OOC1C=O(1193)',
    structure = SMILES('CC1OOC1C=O'),
    E0 = (-168.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67368,0.0406657,-2.94025e-07,-2.52242e-08,1.12547e-11,-20205.3,21.0267], Tmin=(100,'K'), Tmax=(1053.68,'K')), NASAPolynomial(coeffs=[11.9615,0.0247103,-1.04643e-05,2.0165e-09,-1.45049e-13,-23655.7,-35.2289], Tmin=(1053.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(12dioxetane)"""),
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
    label = '[O]O[CH]C=O(3604)',
    structure = SMILES('[O]O[CH]C=O'),
    E0 = (38.8345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.304071,'amu*angstrom^2'), symmetry=1, barrier=(6.9912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17371,'amu*angstrom^2'), symmetry=1, barrier=(56.7408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0355,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48924,0.0359331,-5.02567e-05,4.15061e-08,-1.40856e-11,4722.53,16.7287], Tmin=(100,'K'), Tmax=(782.809,'K')), NASAPolynomial(coeffs=[5.94539,0.0159797,-7.6283e-06,1.46013e-09,-1.01349e-13,4251.69,1.34978], Tmin=(782.809,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.8345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(ROOJ)"""),
)

species(
    label = 'C[CH]O[O](97)',
    structure = SMILES('C[CH]O[O]'),
    E0 = (162.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.659845,'amu*angstrom^2'), symmetry=1, barrier=(15.1711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000844769,'amu*angstrom^2'), symmetry=1, barrier=(5.4037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88063,0.0212973,-8.38042e-06,-1.15683e-09,1.18565e-12,19629.9,16.0899], Tmin=(100,'K'), Tmax=(1195.85,'K')), NASAPolynomial(coeffs=[6.93453,0.0132692,-5.24916e-06,9.65724e-10,-6.67381e-14,18264.7,-5.84667], Tmin=(1195.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CCsJOOH) + radical(ROOJ)"""),
)

species(
    label = 'CHCHO(47)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (245.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06864,0.0187233,-1.21319e-05,-3.33727e-10,2.32882e-12,29739.4,14.7866], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.96288,0.00799899,-4.30606e-06,1.11076e-09,-1.11415e-13,28725.6,-5.17392], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(245.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CC1OOC1[CH][O](4469)',
    structure = SMILES('CC1OOC1[CH][O]'),
    E0 = (157.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24759,0.0668744,-9.11793e-05,8.06572e-08,-2.96644e-11,19050.1,21.1777], Tmin=(100,'K'), Tmax=(786.669,'K')), NASAPolynomial(coeffs=[4.98458,0.0383026,-1.84512e-05,3.55888e-09,-2.48275e-13,18758.2,5.9277], Tmin=(786.669,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxetane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'C[CH]OOC=[C]O(4470)',
    structure = SMILES('C[CH]OOC=[C]O'),
    E0 = (137.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,500,795,815,214.345],'cm^-1')),
        HinderedRotor(inertia=(0.57474,'amu*angstrom^2'), symmetry=1, barrier=(18.7385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574752,'amu*angstrom^2'), symmetry=1, barrier=(18.7385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.8912,'amu*angstrom^2'), symmetry=1, barrier=(94.2533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574735,'amu*angstrom^2'), symmetry=1, barrier=(18.7384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574768,'amu*angstrom^2'), symmetry=1, barrier=(18.7383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.323795,0.0747766,-8.08379e-05,4.3608e-08,-9.14869e-12,16689.4,29.8257], Tmin=(100,'K'), Tmax=(1171.95,'K')), NASAPolynomial(coeffs=[17.4486,0.0163283,-6.02952e-06,1.05342e-09,-7.10479e-14,12675.5,-55.5087], Tmin=(1171.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CCsJOOC)"""),
)

species(
    label = 'C[CH]OO[C]=CO(4471)',
    structure = SMILES('C[CH]OO[C]=CO'),
    E0 = (137.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,500,795,815,214.345],'cm^-1')),
        HinderedRotor(inertia=(0.57474,'amu*angstrom^2'), symmetry=1, barrier=(18.7385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574752,'amu*angstrom^2'), symmetry=1, barrier=(18.7385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.8912,'amu*angstrom^2'), symmetry=1, barrier=(94.2533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574735,'amu*angstrom^2'), symmetry=1, barrier=(18.7384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574768,'amu*angstrom^2'), symmetry=1, barrier=(18.7383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.323795,0.0747766,-8.08379e-05,4.3608e-08,-9.14869e-12,16689.4,29.8257], Tmin=(100,'K'), Tmax=(1171.95,'K')), NASAPolynomial(coeffs=[17.4486,0.0163283,-6.02952e-06,1.05342e-09,-7.10479e-14,12675.5,-55.5087], Tmin=(1171.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CCsJOOC)"""),
)

species(
    label = 'CCOO[C]=C[O](4472)',
    structure = SMILES('CCOO[C]=C[O]'),
    E0 = (92.6421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,222.25,325.84],'cm^-1')),
        HinderedRotor(inertia=(0.234808,'amu*angstrom^2'), symmetry=1, barrier=(16.833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23308,'amu*angstrom^2'), symmetry=1, barrier=(16.8345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232782,'amu*angstrom^2'), symmetry=1, barrier=(16.8362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67373,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.789982,0.0651822,-6.03244e-05,2.83049e-08,-5.28076e-12,11262.3,29.27], Tmin=(100,'K'), Tmax=(1292.33,'K')), NASAPolynomial(coeffs=[15.2174,0.0205263,-8.49193e-06,1.56609e-09,-1.08115e-13,7533.38,-44.0335], Tmin=(1292.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.6421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[CH2][CH]OOC=CO(4473)',
    structure = SMILES('[CH2][CH]OOC=CO'),
    E0 = (111.839,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.06057,'amu*angstrom^2'), symmetry=1, barrier=(24.3847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05935,'amu*angstrom^2'), symmetry=1, barrier=(24.3565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06033,'amu*angstrom^2'), symmetry=1, barrier=(24.3791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06021,'amu*angstrom^2'), symmetry=1, barrier=(24.3764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05699,'amu*angstrom^2'), symmetry=1, barrier=(24.3023,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.50007,0.0821297,-8.91529e-05,4.56364e-08,-8.74484e-12,13627.5,31.6542], Tmin=(100,'K'), Tmax=(1433.25,'K')), NASAPolynomial(coeffs=[23.1709,0.00682097,-6.59861e-07,-2.66872e-11,5.24423e-15,7791.93,-87.7515], Tmin=(1433.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CJCOOH) + radical(CCsJOOC)"""),
)

species(
    label = 'C=COOC=CO(4474)',
    structure = SMILES('C=COOC=CO'),
    E0 = (-145.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.400966,0.0610024,-1.61252e-05,-4.34154e-08,2.69775e-11,-17370.8,26.8066], Tmin=(100,'K'), Tmax=(919.782,'K')), NASAPolynomial(coeffs=[24.8874,0.00210334,2.31955e-06,-5.32558e-10,3.23674e-14,-23888.3,-100.222], Tmin=(919.782,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C1=COO1(3609)',
    structure = SMILES('C1=COO1'),
    E0 = (143.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31094,-0.00195743,7.27275e-05,-1.01582e-07,4.10226e-11,17349.4,10.0685], Tmin=(100,'K'), Tmax=(927.669,'K')), NASAPolynomial(coeffs=[13.7806,-0.00309387,3.40676e-06,-6.26613e-10,3.47414e-14,13513.3,-49.8617], Tmin=(927.669,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclobutene)"""),
)

species(
    label = 'CC1OC=COO1(4439)',
    structure = SMILES('CC1OC=COO1'),
    E0 = (-221.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42236,0.0301566,7.11353e-05,-1.33832e-07,5.97258e-11,-26577.2,19.1779], Tmin=(100,'K'), Tmax=(909.348,'K')), NASAPolynomial(coeffs=[24.0119,0.000488408,5.1058e-06,-1.13761e-09,7.31019e-14,-33567.3,-103.502], Tmin=(909.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(124trioxene)"""),
)

species(
    label = 'C[CH]OC([O])C=O(1191)',
    structure = SMILES('C[CH]OC([O])C=O'),
    E0 = (-130.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,363.391,363.391,363.393,4000],'cm^-1')),
        HinderedRotor(inertia=(0.107763,'amu*angstrom^2'), symmetry=1, barrier=(10.0983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00127658,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107763,'amu*angstrom^2'), symmetry=1, barrier=(10.0983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276613,'amu*angstrom^2'), symmetry=1, barrier=(25.9209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4012.74,'J/mol'), sigma=(6.59519,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=626.78 K, Pc=31.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11567,0.0676797,-8.0674e-05,5.45494e-08,-1.5257e-11,-15586.2,27.6479], Tmin=(100,'K'), Tmax=(861.407,'K')), NASAPolynomial(coeffs=[9.41464,0.0291425,-1.3567e-05,2.61286e-09,-1.837e-13,-17016,-11.1516], Tmin=(861.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(C=OCOJ)"""),
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
    label = '[CH]=COO[CH]C(3205)',
    structure = SMILES('[CH]=COO[CH]C'),
    E0 = (353.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(0.474872,'amu*angstrom^2'), symmetry=1, barrier=(24.9912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.427888,'amu*angstrom^2'), symmetry=1, barrier=(9.838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08667,'amu*angstrom^2'), symmetry=1, barrier=(24.9848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.78583,'amu*angstrom^2'), symmetry=1, barrier=(87.0438,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3992,0.0569617,-5.26171e-05,2.56104e-08,-5.10509e-12,42641.9,23.3828], Tmin=(100,'K'), Tmax=(1188.88,'K')), NASAPolynomial(coeffs=[11.3176,0.0235909,-1.05129e-05,2.00015e-09,-1.40222e-13,40283.6,-26.1834], Tmin=(1188.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCsJOOC) + radical(Cds_P)"""),
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
    E0 = (35.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (83.4875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (210.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (316.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (223.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (204.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (170.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (62.5292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (138.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (219.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (62.0062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (74.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (44.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (43.9313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (382.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (408.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (157.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (35.8107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (300.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (287.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (136.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (192.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (60.6202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (302.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (43.1782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (349.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (596.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['OCHCHO(48)', 'CH3CHO(37)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['CC1OO[CH]C1[O](4459)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(13013.2,'s^-1'), n=1.81618, Ea=(47.8405,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_csHNd] for rate rule [R6;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 44.1 to 47.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=COO[CH]C=O(4460)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.67e+12,'cm^3/(mol*s)'), n=0.1, Ea=(6.4601,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2816 used for Cds-HH_Cds-OsH;HJ
Exact match found for rate rule [Cds-HH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C[CH]OOC=C=O(4461)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.3e-15,'cm^3/(molecule*s)'), n=1.43, Ea=(25.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ck_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['[CH2]COO[CH]C=O(4462)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.41e+13,'s^-1','+|-',2), n=0, Ea=(188.28,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1500,'K'), comment="""From training reaction 346 used for R2H_S;C_rad_out_H/NonDeO;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C[CH]OOC[C]=O(4463)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.67823e+08,'s^-1'), n=1.48018, Ea=(152.494,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out] for rate rule [R2H_S;CO_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH]OOCC=O(4464)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(61977.9,'s^-1'), n=1.86063, Ea=(59.8007,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['CCOO[CH][C]=O(4465)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R5HJ_3;C_rad_out_H/NonDeC;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[CH][O](176)', '[O]C=C[O](2537)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['C[CH]OOC1[CH]O1(4466)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra_H;radadd_intra_csHO]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['CC1O[CH][CH]OO1(4443)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(487000,'s^-1'), n=1.17, Ea=(26.3592,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHCs] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['C=COOCC=O(4467)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['CCOOC=C=O(4468)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['CC1OOC1C=O(1193)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CHCH3(T)(95)', '[O]O[CH]C=O(3604)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[CH]O[O](97)', 'CHCHO(47)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['CC1OOC1[CH][O](4469)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.73e+06,'s^-1'), n=1.31, Ea=(121.971,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 121.1 to 122.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH3CHO(37)', '[O]C=C[O](2537)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(45607,'m^3/(mol*s)'), n=0.695, Ea=(233.422,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_rad/OneDe] + [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;O_rad/OneDe]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C[CH]OOC=[C]O(4470)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C[CH]OO[C]=CO(4471)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CCOO[C]=C[O](4472)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]OOC=CO(4473)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(60975.7,'s^-1'), n=1.58648, Ea=(80.9836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;C_rad_out_2H;XH_out] for rate rule [R7HJ_1;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['C=COOC=CO(4474)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['C[CH][O](176)', 'C1=COO1(3609)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(267.037,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO_intra] for rate rule [R3OO_SD;Y_rad_intra;OO_intra]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['CC1OC=COO1(4439)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_1H;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C[CH]OO[CH]C=O(1192)'],
    products = ['C[CH]OC([O])C=O(1191)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_R]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O(4)', '[CH]=COO[CH]C(3205)'],
    products = ['C[CH]OO[CH]C=O(1192)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '776',
    isomers = [
        'C[CH]OO[CH]C=O(1192)',
    ],
    reactants = [
        ('OCHCHO(48)', 'CH3CHO(37)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '776',
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

