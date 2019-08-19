species(
    label = 'S(754)(753)',
    structure = SMILES('[CH2]C=CC(C=O)C=C[O]'),
    E0 = (21.0141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,216.278,216.278,216.278],'cm^-1')),
        HinderedRotor(inertia=(0.437863,'amu*angstrom^2'), symmetry=1, barrier=(14.5344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.624656,'amu*angstrom^2'), symmetry=1, barrier=(20.7346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.624657,'amu*angstrom^2'), symmetry=1, barrier=(20.7346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0518089,'amu*angstrom^2'), symmetry=1, barrier=(20.7346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4437.15,'J/mol'), sigma=(7.04607,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=693.07 K, Pc=28.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.275399,0.0842065,-7.3822e-05,3.29303e-08,-5.83953e-12,2689.5,34.3085], Tmin=(100,'K'), Tmax=(1358.38,'K')), NASAPolynomial(coeffs=[19.1669,0.0269551,-1.06018e-05,1.90301e-09,-1.29181e-13,-2592.5,-65.4443], Tmin=(1358.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.0141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = 'HCO(14)(15)',
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
    label = 'C6H7O(506)(505)',
    structure = SMILES('[CH2]C=CC=CC=O'),
    E0 = (52.4633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.09775,'amu*angstrom^2'), symmetry=1, barrier=(25.2395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09742,'amu*angstrom^2'), symmetry=1, barrier=(25.2318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09757,'amu*angstrom^2'), symmetry=1, barrier=(25.2352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3745.2,'J/mol'), sigma=(6.02434,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=584.99 K, Pc=38.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19817,0.051949,-2.18837e-05,-5.89878e-09,4.90491e-12,6418.88,24.4996], Tmin=(100,'K'), Tmax=(1111.19,'K')), NASAPolynomial(coeffs=[13.2992,0.0262162,-1.12126e-05,2.13734e-09,-1.51482e-13,2628.95,-40.109], Tmin=(1111.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.4633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ)"""),
)

species(
    label = 'CHCHO(45)(45)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (245.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06864,0.0187233,-1.21319e-05,-3.33727e-10,2.32882e-12,29739.4,14.7866], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.96288,0.00799899,-4.30606e-06,1.11076e-09,-1.11415e-13,28725.6,-5.17392], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(245.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C5H6O(217)(216)',
    structure = SMILES('C=CC=CC=O'),
    E0 = (-29.5668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.9508,'amu*angstrom^2'), symmetry=1, barrier=(21.8608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.95303,'amu*angstrom^2'), symmetry=1, barrier=(21.912,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3519.98,'J/mol'), sigma=(5.63814,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.81 K, Pc=44.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58677,0.0449815,-2.33922e-05,-1.17435e-10,2.4605e-12,-3462.46,19.7432], Tmin=(100,'K'), Tmax=(1159.43,'K')), NASAPolynomial(coeffs=[12.6828,0.020301,-9.05764e-06,1.75766e-09,-1.25367e-13,-6949.62,-39.3724], Tmin=(1159.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.5668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C5H6O(219)(218)',
    structure = SMILES('[CH2]C=CC=C[O]'),
    E0 = (108.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.59878,'amu*angstrom^2'), symmetry=1, barrier=(36.759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59525,'amu*angstrom^2'), symmetry=1, barrier=(36.6779,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3794.25,'J/mol'), sigma=(6.17562,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.65 K, Pc=36.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.673,0.0325848,3.83615e-05,-8.44753e-08,3.82236e-11,13165.2,20.9814], Tmin=(100,'K'), Tmax=(922.022,'K')), NASAPolynomial(coeffs=[18.3994,0.00812125,-9.23824e-08,-9.0725e-11,1.79375e-15,8036.17,-69.4433], Tmin=(922.022,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'H(3)(3)',
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
    label = 'C=C=CC(C=O)C=C[O](3013)',
    structure = SMILES('C=C=CC(C=O)C=C[O]'),
    E0 = (46.1192,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.893878,'amu*angstrom^2'), symmetry=1, barrier=(20.552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.894485,'amu*angstrom^2'), symmetry=1, barrier=(20.566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.893483,'amu*angstrom^2'), symmetry=1, barrier=(20.5429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0978797,0.0825128,-7.85294e-05,3.87522e-08,-7.6859e-12,5690.11,31.7359], Tmin=(100,'K'), Tmax=(1210.91,'K')), NASAPolynomial(coeffs=[16.2796,0.0290597,-1.23151e-05,2.29798e-09,-1.59704e-13,1771.17,-49.4285], Tmin=(1210.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(46.1192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=CC(C=O)C=C=O(3014)',
    structure = SMILES('[CH2]C=CC(C=O)C=C=O'),
    E0 = (2.45453,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2120,512.5,787.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,231.176,1321.16],'cm^-1')),
        HinderedRotor(inertia=(0.232176,'amu*angstrom^2'), symmetry=1, barrier=(8.7102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.44809,'amu*angstrom^2'), symmetry=1, barrier=(16.8635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858997,'amu*angstrom^2'), symmetry=1, barrier=(32.5817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.44237,'amu*angstrom^2'), symmetry=1, barrier=(16.9093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.638001,0.0737311,-6.21896e-05,2.76754e-08,-5.10621e-12,416.249,31.7648], Tmin=(100,'K'), Tmax=(1263.5,'K')), NASAPolynomial(coeffs=[13.1376,0.0341602,-1.52124e-05,2.88884e-09,-2.01929e-13,-2742.44,-31.4621], Tmin=(1263.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.45453,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC(C=O)=CC=O(3042)',
    structure = SMILES('[CH2]C=CC(C=O)=CC=O'),
    E0 = (-61.4902,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(1.03471,'amu*angstrom^2'), symmetry=1, barrier=(23.79,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03481,'amu*angstrom^2'), symmetry=1, barrier=(23.7922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03441,'amu*angstrom^2'), symmetry=1, barrier=(23.7831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03468,'amu*angstrom^2'), symmetry=1, barrier=(23.7893,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.16087,0.0782262,-6.34218e-05,2.47446e-08,-3.84737e-12,-7252.6,29.4638], Tmin=(100,'K'), Tmax=(1514.44,'K')), NASAPolynomial(coeffs=[19.473,0.0272182,-1.29001e-05,2.50456e-09,-1.76028e-13,-13102,-71.7215], Tmin=(1514.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.4902,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ)"""),
)

species(
    label = 'S(525)(524)',
    structure = SMILES('O=CC=CC=O'),
    E0 = (-204.827,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.770505,'amu*angstrom^2'), symmetry=1, barrier=(17.7154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.772944,'amu*angstrom^2'), symmetry=1, barrier=(17.7715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3741.66,'J/mol'), sigma=(5.7541,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=584.44 K, Pc=44.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03135,0.0345588,-1.9969e-05,4.287e-09,-3.16505e-13,-24613.1,14.2636], Tmin=(100,'K'), Tmax=(2498.65,'K')), NASAPolynomial(coeffs=[27.6611,0.000737979,-3.0321e-06,6.66313e-10,-4.41138e-14,-38671.9,-130.618], Tmin=(2498.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.827,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = '[CH]=C[CH2](891)',
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
    label = '[O][CH]C=CC=O(2788)',
    structure = SMILES('[O]C=CC=C[O]'),
    E0 = (-40.7388,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.42008,'amu*angstrom^2'), symmetry=1, barrier=(32.6505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42346,0.0345738,3.51007e-05,-9.27285e-08,4.46226e-11,-4786.29,19.0904], Tmin=(100,'K'), Tmax=(909.432,'K')), NASAPolynomial(coeffs=[24.1228,-0.00681144,6.94748e-06,-1.41392e-09,9.1781e-14,-11332.3,-101.556], Tmin=(909.432,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.7388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'C=C=CC(C=O)C=CO(3032)',
    structure = SMILES('C=C=CC(C=O)C=CO'),
    E0 = (-95.3434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.602305,0.0928931,-9.39097e-05,4.81936e-08,-9.70181e-12,-11294.4,32.8383], Tmin=(100,'K'), Tmax=(1214.56,'K')), NASAPolynomial(coeffs=[20.2773,0.0241281,-8.98324e-06,1.57753e-09,-1.0647e-13,-16366.3,-71.9525], Tmin=(1214.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.3434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CC=CC(C=O)C=C=O(3033)',
    structure = SMILES('CC=CC(C=O)C=C=O'),
    E0 = (-149.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.850088,0.073643,-6.14265e-05,2.87361e-08,-5.83834e-12,-17816.2,30.7307], Tmin=(100,'K'), Tmax=(1121.25,'K')), NASAPolynomial(coeffs=[9.66271,0.0422047,-1.9369e-05,3.72999e-09,-2.62909e-13,-19792.5,-12.7938], Tmin=(1121.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'CC=CC(C=O)=CC=O(3051)',
    structure = SMILES('CC=CC(C=O)=CC=O'),
    E0 = (-179.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.360685,0.0816733,-6.8505e-05,2.83217e-08,-4.78516e-12,-21465.3,27.1756], Tmin=(100,'K'), Tmax=(1368.67,'K')), NASAPolynomial(coeffs=[16.2764,0.0351591,-1.75276e-05,3.49113e-09,-2.49637e-13,-25821.9,-54.6034], Tmin=(1368.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=C=CC(C=O)CC=O(3052)',
    structure = SMILES('C=C=CC(C=O)CC=O'),
    E0 = (-95.9976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1011,0.0756528,-4.61093e-05,-4.85241e-08,7.03809e-11,-11453.6,29.2728], Tmin=(100,'K'), Tmax=(486.624,'K')), NASAPolynomial(coeffs=[6.88878,0.0481185,-2.30079e-05,4.45411e-09,-3.12724e-13,-12254.1,3.08132], Tmin=(486.624,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.9976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC=C(C=O)C=CO(3061)',
    structure = SMILES('C=CC=C(C=O)C=CO'),
    E0 = (-177.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.93856,0.0864871,-3.76812e-05,-3.67559e-08,2.67366e-11,-21096.4,29.4533], Tmin=(100,'K'), Tmax=(937.108,'K')), NASAPolynomial(coeffs=[30.7281,0.00711695,5.06223e-08,-6.02856e-11,-3.68039e-15,-29481.4,-134.335], Tmin=(937.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CCC(C=O)C=C=O(2844)',
    structure = SMILES('C=CCC(C=O)C=C=O'),
    E0 = (-135.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,264.65,264.65,264.651],'cm^-1')),
        HinderedRotor(inertia=(0.00204861,'amu*angstrom^2'), symmetry=1, barrier=(23.26,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.467991,'amu*angstrom^2'), symmetry=1, barrier=(23.26,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.467989,'amu*angstrom^2'), symmetry=1, barrier=(23.26,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209396,'amu*angstrom^2'), symmetry=1, barrier=(10.4074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.512859,0.0739656,-5.75625e-05,2.3238e-08,-3.8681e-12,-16189.2,33.0393], Tmin=(100,'K'), Tmax=(1391.68,'K')), NASAPolynomial(coeffs=[14.3058,0.0343223,-1.48341e-05,2.7698e-09,-1.9126e-13,-20028.3,-38.0624], Tmin=(1391.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-135.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CCC(C=O)=CC=O(3064)',
    structure = SMILES('C=CCC(C=O)=CC=O'),
    E0 = (-165.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36301,0.0697783,-4.80254e-05,1.55599e-08,-2.11549e-12,-19826.9,27.9577], Tmin=(100,'K'), Tmax=(1545.65,'K')), NASAPolynomial(coeffs=[11.6476,0.043163,-2.21964e-05,4.41957e-09,-3.13624e-13,-23006.2,-26.1381], Tmin=(1545.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=CC=C(C=O)CC=O(3065)',
    structure = SMILES('C=CC=C(C=O)CC=O'),
    E0 = (-172.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0484139,0.0750053,-5.21291e-05,1.62694e-08,-1.96605e-12,-20555,31.175], Tmin=(100,'K'), Tmax=(1952.08,'K')), NASAPolynomial(coeffs=[26.674,0.0204494,-1.02096e-05,1.95391e-09,-1.32767e-13,-30950.5,-115.089], Tmin=(1952.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'O=CC1C=CCOC=C1(3035)',
    structure = SMILES('O=CC1C=CCOC=C1'),
    E0 = (-156.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82574,0.0103025,0.000222042,-3.42962e-07,1.46788e-10,-18639.5,27.254], Tmin=(100,'K'), Tmax=(903.912,'K')), NASAPolynomial(coeffs=[44.2633,-0.0256727,2.24603e-05,-4.53436e-09,2.97685e-13,-32875.3,-213.227], Tmin=(903.912,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-156.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(Cycloheptane)"""),
)

species(
    label = 'S(790)(789)',
    structure = SMILES('O=CC1C=CCC1C=O'),
    E0 = (-203.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4330.58,'J/mol'), sigma=(6.80044,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=676.43 K, Pc=31.24 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.866885,0.0581741,-1.79512e-05,-1.1148e-08,6.30878e-12,-24391.9,28.1439], Tmin=(100,'K'), Tmax=(1138.24,'K')), NASAPolynomial(coeffs=[12.987,0.0358288,-1.51864e-05,2.86058e-09,-2.00508e-13,-28462.7,-37.6597], Tmin=(1138.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + ring(Cyclopentene)"""),
)

species(
    label = 'S(651)(650)',
    structure = SMILES('C=CC1OC=CC1C=O'),
    E0 = (-174.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4159.18,'J/mol'), sigma=(6.6066,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=649.65 K, Pc=32.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.724483,0.0491554,4.09266e-05,-1.0001e-07,4.56656e-11,-20879.2,25.7074], Tmin=(100,'K'), Tmax=(925.525,'K')), NASAPolynomial(coeffs=[21.71,0.017362,-3.01104e-06,4.03535e-10,-3.23453e-14,-27286.6,-87.5408], Tmin=(925.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran)"""),
)

species(
    label = 'S(756)(755)',
    structure = SMILES('C=CC1C(C=O)C1C=O'),
    E0 = (-110.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4142.9,'J/mol'), sigma=(6.50707,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=647.11 K, Pc=34.12 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605121,0.0625765,-2.22475e-05,-1.36289e-08,8.88868e-12,-13140.4,30.71], Tmin=(100,'K'), Tmax=(1050.9,'K')), NASAPolynomial(coeffs=[15.1698,0.0315098,-1.26895e-05,2.37454e-09,-1.67898e-13,-17547.3,-46.682], Tmin=(1050.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclopropane)"""),
)

species(
    label = '[O][CH]C1CC=CC1C=O(3010)',
    structure = SMILES('[O][CH]C1CC=CC1C=O'),
    E0 = (121.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07395,0.0652969,-4.08181e-05,1.25015e-08,-1.5816e-12,14761.7,28.3182], Tmin=(100,'K'), Tmax=(1722.37,'K')), NASAPolynomial(coeffs=[13.1939,0.0371498,-1.63052e-05,3.01351e-09,-2.04439e-13,10586.7,-36.7434], Tmin=(1722.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopentene) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O]C=CC1C=CCC1[O](3011)',
    structure = SMILES('[O]C=CC1C=CCC1[O]'),
    E0 = (84.2593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.415789,0.0577845,1.71835e-05,-7.11295e-08,3.31984e-11,10282.1,27.9106], Tmin=(100,'K'), Tmax=(961.688,'K')), NASAPolynomial(coeffs=[21.4233,0.0218676,-7.06112e-06,1.32006e-09,-1.00391e-13,3861.94,-84.9899], Tmin=(961.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.2593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopentene) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC1OC=CC1[CH][O](2894)',
    structure = SMILES('C=CC1OC=CC1[CH][O]'),
    E0 = (153.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.464159,0.0582861,1.37669e-05,-7.07671e-08,3.50072e-11,18564.1,28.2839], Tmin=(100,'K'), Tmax=(925.996,'K')), NASAPolynomial(coeffs=[21.1023,0.0192555,-4.19537e-06,6.15017e-10,-4.48634e-14,12593.1,-81.2988], Tmin=(925.996,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C=CC1C=COC1[O](3012)',
    structure = SMILES('[CH2]C=CC1C=COC1[O]'),
    E0 = (84.1258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.54339,0.0594331,3.18583e-06,-5.30092e-08,2.67602e-11,10257.6,25.798], Tmin=(100,'K'), Tmax=(941.023,'K')), NASAPolynomial(coeffs=[18.5067,0.0245371,-7.2782e-06,1.22444e-09,-8.6673e-14,5041.12,-69.5266], Tmin=(941.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.1258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CCOJ) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=CC(C=O)C=C[O](3015)',
    structure = SMILES('C[C]=CC(C=O)C=C[O]'),
    E0 = (107.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.270622,0.084493,-8.35564e-05,4.54363e-08,-1.02401e-11,13044.3,32.5435], Tmin=(100,'K'), Tmax=(1054.62,'K')), NASAPolynomial(coeffs=[12.6627,0.0374922,-1.67072e-05,3.17869e-09,-2.22897e-13,10430.5,-27.9002], Tmin=(1054.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC(C=O)C=[C]O(3016)',
    structure = SMILES('[CH2]C=CC(C=O)C=[C]O'),
    E0 = (119.296,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,300.238,300.238],'cm^-1')),
        HinderedRotor(inertia=(0.206494,'amu*angstrom^2'), symmetry=1, barrier=(13.2089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206494,'amu*angstrom^2'), symmetry=1, barrier=(13.2089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206493,'amu*angstrom^2'), symmetry=1, barrier=(13.2089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.358991,'amu*angstrom^2'), symmetry=1, barrier=(22.9638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206494,'amu*angstrom^2'), symmetry=1, barrier=(13.2089,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.126078,0.0888906,-9.04715e-05,4.86347e-08,-1.05048e-11,14498.1,35.4863], Tmin=(100,'K'), Tmax=(1118.53,'K')), NASAPolynomial(coeffs=[16.0615,0.0310018,-1.284e-05,2.36477e-09,-1.63177e-13,10876.9,-44.4227], Tmin=(1118.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = 'CC=[C]C(C=O)C=C[O](3017)',
    structure = SMILES('CC=[C]C(C=O)C=C[O]'),
    E0 = (107.357,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,226.13,226.131,226.131],'cm^-1')),
        HinderedRotor(inertia=(0.361157,'amu*angstrom^2'), symmetry=1, barrier=(13.1051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361157,'amu*angstrom^2'), symmetry=1, barrier=(13.1051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361145,'amu*angstrom^2'), symmetry=1, barrier=(13.1051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361153,'amu*angstrom^2'), symmetry=1, barrier=(13.1051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.270622,0.084493,-8.35564e-05,4.54363e-08,-1.02401e-11,13044.3,32.5435], Tmin=(100,'K'), Tmax=(1054.62,'K')), NASAPolynomial(coeffs=[12.6627,0.0374922,-1.67072e-05,3.17869e-09,-2.22897e-13,10430.5,-27.9002], Tmin=(1054.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC([C]=CO)C=O(3018)',
    structure = SMILES('[CH2]C=CC([C]=CO)C=O'),
    E0 = (117.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.777378,'amu*angstrom^2'), symmetry=1, barrier=(17.8735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778043,'amu*angstrom^2'), symmetry=1, barrier=(17.8887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778291,'amu*angstrom^2'), symmetry=1, barrier=(17.8944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777974,'amu*angstrom^2'), symmetry=1, barrier=(17.8871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777129,'amu*angstrom^2'), symmetry=1, barrier=(17.8677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.599602,0.0944832,-9.81307e-05,5.1953e-08,-1.08074e-11,14290.5,34.5284], Tmin=(100,'K'), Tmax=(1175.24,'K')), NASAPolynomial(coeffs=[19.7996,0.0250532,-9.5147e-06,1.68464e-09,-1.1414e-13,9495.67,-67.1799], Tmin=(1175.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=C[C](C=O)C=C[O](3019)',
    structure = SMILES('CC=C[C](C=O)C=C[O]'),
    E0 = (-39.0601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.179177,0.0707381,-3.22626e-05,-1.02387e-08,9.15294e-12,-4548.71,30.3792], Tmin=(100,'K'), Tmax=(1026.84,'K')), NASAPolynomial(coeffs=[17.7884,0.0293675,-1.15998e-05,2.16734e-09,-1.54121e-13,-9600.39,-62.0306], Tmin=(1026.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.0601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CCJ(C=O)C=C)"""),
)

species(
    label = '[CH2]C=C[C](C=O)C=CO(3020)',
    structure = SMILES('C=C[CH]C(C=CO)=C[O]'),
    E0 = (-35.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.570257,0.0687487,2.95942e-05,-1.204e-07,6.10156e-11,-4068.25,31.459], Tmin=(100,'K'), Tmax=(908.224,'K')), NASAPolynomial(coeffs=[34.7697,-0.00235278,7.39451e-06,-1.61207e-09,1.05466e-13,-13974.4,-154.83], Tmin=(908.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.441,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC=CC([C]=C[O])C=O(3021)',
    structure = SMILES('CC=CC([C]=C[O])C=O'),
    E0 = (107.357,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,226.13,226.131,226.131],'cm^-1')),
        HinderedRotor(inertia=(0.361157,'amu*angstrom^2'), symmetry=1, barrier=(13.1051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361157,'amu*angstrom^2'), symmetry=1, barrier=(13.1051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361145,'amu*angstrom^2'), symmetry=1, barrier=(13.1051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361153,'amu*angstrom^2'), symmetry=1, barrier=(13.1051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.270622,0.084493,-8.35564e-05,4.54363e-08,-1.02401e-11,13044.3,32.5435], Tmin=(100,'K'), Tmax=(1054.62,'K')), NASAPolynomial(coeffs=[12.6627,0.0374922,-1.67072e-05,3.17869e-09,-2.22897e-13,10430.5,-27.9002], Tmin=(1054.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'CC=CC([C]=O)C=C[O](3022)',
    structure = SMILES('CC=CC([C]=O)C=C[O]'),
    E0 = (28.2035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1380,1390,370,380,2900,435,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,234.302,234.306,234.316],'cm^-1')),
        HinderedRotor(inertia=(0.364648,'amu*angstrom^2'), symmetry=1, barrier=(14.206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364644,'amu*angstrom^2'), symmetry=1, barrier=(14.206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36467,'amu*angstrom^2'), symmetry=1, barrier=(14.2061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364686,'amu*angstrom^2'), symmetry=1, barrier=(14.2061,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0533031,0.0822806,-7.3985e-05,3.46702e-08,-6.55833e-12,3537.9,33.9986], Tmin=(100,'K'), Tmax=(1263.2,'K')), NASAPolynomial(coeffs=[16.3384,0.0307125,-1.27495e-05,2.35228e-09,-1.6225e-13,-576.337,-48.3725], Tmin=(1263.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.2035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CC(C)CJ=O) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=[C]C(C=O)C=CO(3023)',
    structure = SMILES('[CH2]C=[C]C(C=O)C=CO'),
    E0 = (117.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.777378,'amu*angstrom^2'), symmetry=1, barrier=(17.8735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778043,'amu*angstrom^2'), symmetry=1, barrier=(17.8887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778291,'amu*angstrom^2'), symmetry=1, barrier=(17.8944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777974,'amu*angstrom^2'), symmetry=1, barrier=(17.8871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777129,'amu*angstrom^2'), symmetry=1, barrier=(17.8677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.599602,0.0944832,-9.81307e-05,5.1953e-08,-1.08074e-11,14290.5,34.5284], Tmin=(100,'K'), Tmax=(1175.24,'K')), NASAPolynomial(coeffs=[19.7996,0.0250532,-9.5147e-06,1.68464e-09,-1.1414e-13,9495.67,-67.1799], Tmin=(1175.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC([C]=O)C=CO(3024)',
    structure = SMILES('[CH2]C=CC([C]=O)C=CO'),
    E0 = (38.2401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.80953,'amu*angstrom^2'), symmetry=1, barrier=(18.6127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.809907,'amu*angstrom^2'), symmetry=1, barrier=(18.6213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.809706,'amu*angstrom^2'), symmetry=1, barrier=(18.6167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.809638,'amu*angstrom^2'), symmetry=1, barrier=(18.6152,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.809833,'amu*angstrom^2'), symmetry=1, barrier=(18.6197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.532337,0.0888511,-7.64481e-05,2.55318e-08,-5.64514e-13,4772.03,34.9687], Tmin=(100,'K'), Tmax=(986.667,'K')), NASAPolynomial(coeffs=[21.3722,0.0217959,-7.56742e-06,1.32965e-09,-9.23825e-14,-609.022,-75.7786], Tmin=(986.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.2401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'CC=CC(C=O)C=[C][O](3025)',
    structure = SMILES('CC=CC([CH][C]=O)C=O'),
    E0 = (87.2607,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,218.614,2741.03],'cm^-1')),
        HinderedRotor(inertia=(0.00351519,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285433,'amu*angstrom^2'), symmetry=1, barrier=(9.72356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286452,'amu*angstrom^2'), symmetry=1, barrier=(9.72268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286871,'amu*angstrom^2'), symmetry=1, barrier=(9.72298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286445,'amu*angstrom^2'), symmetry=1, barrier=(9.72287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07599,0.0754749,-5.27292e-05,-3.15805e-08,5.6035e-11,10589.1,33.6015], Tmin=(100,'K'), Tmax=(502.112,'K')), NASAPolynomial(coeffs=[7.30298,0.0453994,-2.12286e-05,4.06353e-09,-2.83249e-13,9717.56,5.39812], Tmin=(502.112,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.2607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]=CC(C=O)C=CO(3026)',
    structure = SMILES('[CH2][C]=CC(C=O)C=CO'),
    E0 = (117.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.777378,'amu*angstrom^2'), symmetry=1, barrier=(17.8735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778043,'amu*angstrom^2'), symmetry=1, barrier=(17.8887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778291,'amu*angstrom^2'), symmetry=1, barrier=(17.8944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777974,'amu*angstrom^2'), symmetry=1, barrier=(17.8871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777129,'amu*angstrom^2'), symmetry=1, barrier=(17.8677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.599602,0.0944832,-9.81307e-05,5.1953e-08,-1.08074e-11,14290.5,34.5284], Tmin=(100,'K'), Tmax=(1175.24,'K')), NASAPolynomial(coeffs=[19.7996,0.0250532,-9.5147e-06,1.68464e-09,-1.1414e-13,9495.67,-67.1799], Tmin=(1175.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[O]C=CC(C=O)C1[CH]C1(3027)',
    structure = SMILES('[O]C=CC(C=O)C1[CH]C1'),
    E0 = (134.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.368297,0.0660071,-2.25174e-05,-1.98837e-08,1.26449e-11,16329.2,33.964], Tmin=(100,'K'), Tmax=(1008.06,'K')), NASAPolynomial(coeffs=[17.5636,0.0279333,-1.07384e-05,2.00411e-09,-1.43437e-13,11330.2,-56.7317], Tmin=(1008.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(cyclopropane) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=CC(C=O)C1[CH]O1(3028)',
    structure = SMILES('[CH2]C=CC(C=O)C1[CH]O1'),
    E0 = (162.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.595589,0.0888792,-8.67998e-05,4.5206e-08,-9.1796e-12,19721.6,33.9719], Tmin=(100,'K'), Tmax=(1313.16,'K')), NASAPolynomial(coeffs=[18.2881,0.0250646,-6.71661e-06,8.9965e-10,-4.97255e-14,15304.7,-60.2094], Tmin=(1313.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(Allyl_P) + radical(CCsJO)"""),
)

species(
    label = '[O]C1[CH]C(C=O)C=CC1(3029)',
    structure = SMILES('[O]C1[CH]C(C=O)C=CC1'),
    E0 = (123.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.938135,0.0514305,1.34791e-05,-5.06341e-08,2.16662e-11,14923.5,28.9486], Tmin=(100,'K'), Tmax=(1016.47,'K')), NASAPolynomial(coeffs=[15.092,0.0326028,-1.31458e-05,2.51282e-09,-1.81755e-13,10141.3,-48.9361], Tmin=(1016.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclohexene) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C=CC1[CH]OCC=C1(3030)',
    structure = SMILES('[O]C=CC1[CH]OCC=C1'),
    E0 = (64.5817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.179281,0.0663368,-4.24126e-06,-5.44679e-08,2.99455e-11,7921.39,25.5282], Tmin=(100,'K'), Tmax=(917.734,'K')), NASAPolynomial(coeffs=[21.4385,0.0200196,-4.28221e-06,5.84652e-10,-4.01024e-14,2067.77,-85.8425], Tmin=(917.734,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.5817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(36dihydro2hpyran) + radical(C=COJ) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C1[CH]C(C=O)C=CO1(3031)',
    structure = SMILES('[CH2]C1[CH]C(C=O)C=CO1'),
    E0 = (106.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.365091,0.0634506,-4.95229e-06,-4.84224e-08,2.63285e-11,12910.3,27.472], Tmin=(100,'K'), Tmax=(929.341,'K')), NASAPolynomial(coeffs=[19.8402,0.0217658,-5.68477e-06,8.92895e-10,-6.25417e-14,7470.76,-74.8474], Tmin=(929.341,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(3,4-Dihydro-2H-pyran) + radical(CJC(C)OC) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C=CC1[CH]OOC=C1(2747)',
    structure = SMILES('[CH2]C=CC1[CH]OOC=C1'),
    E0 = (335.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.353545,0.0676095,-2.33438e-05,-2.11552e-08,1.41998e-11,40534.3,26.149], Tmin=(100,'K'), Tmax=(960.651,'K')), NASAPolynomial(coeffs=[16.6425,0.029436,-1.00368e-05,1.74026e-09,-1.20059e-13,36036.5,-58.9032], Tmin=(960.651,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(34dihydro12dioxin) + radical(Allyl_P) + radical(CCsJOOC)"""),
)

species(
    label = 'CO(10)(11)',
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
    label = '[CH2]C=CCC=C[O](3034)',
    structure = SMILES('[CH2]C=CCC=C[O]'),
    E0 = (136.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,288.83,288.868,289.433],'cm^-1')),
        HinderedRotor(inertia=(0.382724,'amu*angstrom^2'), symmetry=1, barrier=(22.4549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.376442,'amu*angstrom^2'), symmetry=1, barrier=(22.4702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.379495,'amu*angstrom^2'), symmetry=1, barrier=(22.4623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.949875,0.0501072,1.01307e-05,-5.57741e-08,2.67736e-11,16589.7,27.1272], Tmin=(100,'K'), Tmax=(955.034,'K')), NASAPolynomial(coeffs=[18.4853,0.0184131,-5.66331e-06,1.02492e-09,-7.696e-14,11336.3,-66.6323], Tmin=(955.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = 'S(650)(649)',
    structure = SMILES('C=C[CH]OC=CC=C[O]'),
    E0 = (19.4331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3025,407.5,1350,352.5,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4211.59,'J/mol'), sigma=(6.78269,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=657.84 K, Pc=30.63 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.483862,0.0752234,-9.10794e-06,-6.2386e-08,3.50091e-11,2520.25,32.7202], Tmin=(100,'K'), Tmax=(933.469,'K')), NASAPolynomial(coeffs=[28.4671,0.0106677,-9.87178e-07,1.00325e-10,-1.4161e-14,-5477.11,-118.844], Tmin=(933.469,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.4331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=CC=COC=C[O](3036)',
    structure = SMILES('C=CC=C[CH]OC=C[O]'),
    E0 = (19.4331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3025,407.5,1350,352.5,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.483862,0.0752234,-9.10794e-06,-6.2386e-08,3.50091e-11,2520.25,32.7202], Tmin=(100,'K'), Tmax=(933.469,'K')), NASAPolynomial(coeffs=[28.4671,0.0106677,-9.87178e-07,1.00325e-10,-1.4161e-14,-5477.11,-118.844], Tmin=(933.469,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.4331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH2]C=CC(C=C[O])=CO(3037)',
    structure = SMILES('C=C[CH]C(C=C[O])=CO'),
    E0 = (-35.441,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29615,'amu*angstrom^2'), symmetry=1, barrier=(29.801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29606,'amu*angstrom^2'), symmetry=1, barrier=(29.799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29614,'amu*angstrom^2'), symmetry=1, barrier=(29.8007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29641,'amu*angstrom^2'), symmetry=1, barrier=(29.807,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.570257,0.0687487,2.95942e-05,-1.204e-07,6.10156e-11,-4068.25,31.459], Tmin=(100,'K'), Tmax=(908.224,'K')), NASAPolynomial(coeffs=[34.7697,-0.00235278,7.39451e-06,-1.61207e-09,1.05466e-13,-13974.4,-154.83], Tmin=(908.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.441,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'CH2(17)(18)',
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
    label = '[CH]=CC(C=O)C=C[O](3038)',
    structure = SMILES('[CH]=CC(C=O)C=C[O]'),
    E0 = (152.637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.777964,'amu*angstrom^2'), symmetry=1, barrier=(17.8869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778742,'amu*angstrom^2'), symmetry=1, barrier=(17.9048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777458,'amu*angstrom^2'), symmetry=1, barrier=(17.8753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.53486,0.0719582,-6.96043e-05,3.45807e-08,-6.8438e-12,18486.3,29.7106], Tmin=(100,'K'), Tmax=(1221.03,'K')), NASAPolynomial(coeffs=[15.4802,0.0229988,-9.45945e-06,1.74255e-09,-1.20388e-13,14836.5,-45.3764], Tmin=(1221.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'O(4)(4)',
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
    label = '[CH]=CC(C=O)C=C[CH2](3039)',
    structure = SMILES('[CH]=CC(C=O)C=C[CH2]'),
    E0 = (335.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,293.72],'cm^-1')),
        HinderedRotor(inertia=(0.163306,'amu*angstrom^2'), symmetry=1, barrier=(9.99672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163274,'amu*angstrom^2'), symmetry=1, barrier=(9.9967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163233,'amu*angstrom^2'), symmetry=1, barrier=(9.99689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.546491,'amu*angstrom^2'), symmetry=1, barrier=(33.4717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.914813,0.0714194,-6.31986e-05,3.15858e-08,-6.77485e-12,40452.3,29.1719], Tmin=(100,'K'), Tmax=(1079.49,'K')), NASAPolynomial(coeffs=[9.79513,0.0385142,-1.7476e-05,3.34903e-09,-2.35559e-13,38535,-14.3498], Tmin=(1079.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]C1C(C=O)C1C=O(3040)',
    structure = SMILES('[CH2][CH]C1C(C=O)C1C=O'),
    E0 = (160.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.734603,0.0676307,-4.58586e-05,1.55365e-08,-2.1545e-12,19433.6,34.0811], Tmin=(100,'K'), Tmax=(1636.72,'K')), NASAPolynomial(coeffs=[14.9448,0.032902,-1.40308e-05,2.5724e-09,-1.74287e-13,14782.1,-41.4761], Tmin=(1636.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC1C([O])C1C=O(3041)',
    structure = SMILES('[CH2]C=CC1C([O])C1C=O'),
    E0 = (175.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.379845,0.0638823,-1.17787e-05,-3.3508e-08,1.80013e-11,21233.6,31.6301], Tmin=(100,'K'), Tmax=(991.639,'K')), NASAPolynomial(coeffs=[18.4856,0.0267937,-1.00486e-05,1.8825e-09,-1.36404e-13,15875.4,-64.4785], Tmin=(991.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C[C](C=O)CC=O(3043)',
    structure = SMILES('[CH2]C=CC(=C[O])CC=O'),
    E0 = (-28.2317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0833961,0.065652,-9.8627e-07,-5.27529e-08,2.61419e-11,-3236.22,31.7872], Tmin=(100,'K'), Tmax=(987.236,'K')), NASAPolynomial(coeffs=[22.5471,0.0226685,-8.65841e-06,1.71092e-09,-1.30109e-13,-10012.3,-88.1534], Tmin=(987.236,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.2317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=CC(C=O)C[C]=O(3044)',
    structure = SMILES('[CH2]C=CC(C=O)C[C]=O'),
    E0 = (38.8578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,252.343,252.343],'cm^-1')),
        HinderedRotor(inertia=(0.123948,'amu*angstrom^2'), symmetry=1, barrier=(5.60268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123996,'amu*angstrom^2'), symmetry=1, barrier=(5.60271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124009,'amu*angstrom^2'), symmetry=1, barrier=(5.60269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124002,'amu*angstrom^2'), symmetry=1, barrier=(5.60275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.497373,'amu*angstrom^2'), symmetry=1, barrier=(22.4714,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.280662,0.0880056,-0.000105607,7.71522e-08,-2.38956e-11,4801.78,34.0953], Tmin=(100,'K'), Tmax=(775.008,'K')), NASAPolynomial(coeffs=[8.80091,0.0440291,-2.0489e-05,3.9306e-09,-2.75139e-13,3481.17,-4.83802], Tmin=(775.008,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.8578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=[C]C(C=O)CC=O(3045)',
    structure = SMILES('[CH2]C=[C]C(C=O)CC=O'),
    E0 = (116.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,192.026,192.032],'cm^-1')),
        HinderedRotor(inertia=(0.272438,'amu*angstrom^2'), symmetry=1, barrier=(7.12924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272413,'amu*angstrom^2'), symmetry=1, barrier=(7.12926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54422,'amu*angstrom^2'), symmetry=1, barrier=(40.4136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000627906,'amu*angstrom^2'), symmetry=1, barrier=(7.12928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54438,'amu*angstrom^2'), symmetry=1, barrier=(40.4137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19266,0.0749876,-3.3265e-05,-9.21512e-08,1.12012e-10,14127.8,30.7049], Tmin=(100,'K'), Tmax=(463.219,'K')), NASAPolynomial(coeffs=[7.11101,0.0478216,-2.28197e-05,4.38844e-09,-3.05878e-13,13322.6,3.93434], Tmin=(463.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC([C]=O)CC=O(3046)',
    structure = SMILES('[CH2]C=CC([C]=O)CC=O'),
    E0 = (37.5859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,329.427,2066.64],'cm^-1')),
        HinderedRotor(inertia=(0.134966,'amu*angstrom^2'), symmetry=1, barrier=(10.3937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134966,'amu*angstrom^2'), symmetry=1, barrier=(10.3937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134967,'amu*angstrom^2'), symmetry=1, barrier=(10.3937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134966,'amu*angstrom^2'), symmetry=1, barrier=(10.3937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.56404,'amu*angstrom^2'), symmetry=1, barrier=(43.4366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.557167,0.0818763,-8.55154e-05,5.40012e-08,-1.47927e-11,4639.18,33.4483], Tmin=(100,'K'), Tmax=(860.223,'K')), NASAPolynomial(coeffs=[8.49926,0.0449466,-2.11211e-05,4.09701e-09,-2.89693e-13,3272.76,-3.67207], Tmin=(860.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.5859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2][C]=CC(C=O)CC=O(3047)',
    structure = SMILES('[CH2][C]=CC(C=O)CC=O'),
    E0 = (116.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,192.026,192.032],'cm^-1')),
        HinderedRotor(inertia=(0.272438,'amu*angstrom^2'), symmetry=1, barrier=(7.12924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272413,'amu*angstrom^2'), symmetry=1, barrier=(7.12926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54422,'amu*angstrom^2'), symmetry=1, barrier=(40.4136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000627906,'amu*angstrom^2'), symmetry=1, barrier=(7.12928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54438,'amu*angstrom^2'), symmetry=1, barrier=(40.4137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19266,0.0749876,-3.3265e-05,-9.21512e-08,1.12012e-10,14127.8,30.7049], Tmin=(100,'K'), Tmax=(463.219,'K')), NASAPolynomial(coeffs=[7.11101,0.0478216,-2.28197e-05,4.38844e-09,-3.05878e-13,13322.6,3.93434], Tmin=(463.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1[CH]C(C=O)C1C=O(3048)',
    structure = SMILES('[CH2]C1[CH]C(C=O)C1C=O'),
    E0 = (161.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.812216,0.0579572,-1.06263e-05,-2.47663e-08,1.28682e-11,19531.9,33.7481], Tmin=(100,'K'), Tmax=(1014.25,'K')), NASAPolynomial(coeffs=[13.9373,0.0327297,-1.25605e-05,2.29994e-09,-1.6136e-13,15504.7,-36.4872], Tmin=(1014.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(Isobutyl) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]C=CC1[CH]OC1C=O(3049)',
    structure = SMILES('[CH2]C=CC1[CH]OC1C=O'),
    E0 = (155.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.56363,0.054131,3.19621e-05,-9.49529e-08,4.56975e-11,18813.2,31.7529], Tmin=(100,'K'), Tmax=(903.388,'K')), NASAPolynomial(coeffs=[21.8584,0.0171052,-1.63947e-06,9.0865e-12,8.65567e-16,12629.1,-81.7511], Tmin=(903.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Oxetane) + radical(Allyl_P) + radical(CCsJOCs)"""),
)

species(
    label = 'O=CC1[CH][CH]OCC=C1(3050)',
    structure = SMILES('O=CC1[CH][CH]OCC=C1'),
    E0 = (148.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.757803,0.0228786,0.000164638,-2.68559e-07,1.16577e-10,18003.9,32.1786], Tmin=(100,'K'), Tmax=(904.474,'K')), NASAPolynomial(coeffs=[38.1085,-0.0149099,1.60356e-05,-3.30347e-09,2.16831e-13,6036.49,-173.073], Tmin=(904.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptane) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C=C[CH]C(C=O)C=O(3053)',
    structure = SMILES('[CH2][CH]C=CC(C=O)C=O'),
    E0 = (75.1916,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0179188,'amu*angstrom^2'), symmetry=1, barrier=(16.8218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00787879,'amu*angstrom^2'), symmetry=1, barrier=(7.39218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127256,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.669895,0.0671028,-4.33969e-05,1.36203e-08,-1.72419e-12,9168.03,34.5647], Tmin=(100,'K'), Tmax=(1796.53,'K')), NASAPolynomial(coeffs=[16.9843,0.0307781,-1.30674e-05,2.36534e-09,-1.5796e-13,3306.23,-53.7008], Tmin=(1796.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.1916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(Allyl_S) + radical(RCCJ)"""),
)

species(
    label = 'C=CC1C([O])C1[CH]C=O(2760)',
    structure = SMILES('C=CC1C([O])C1C=C[O]'),
    E0 = (177.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.32858,0.0542341,4.01074e-05,-1.05623e-07,4.83617e-11,21460.7,31.5755], Tmin=(100,'K'), Tmax=(939.776,'K')), NASAPolynomial(coeffs=[25.5935,0.013753,-2.30735e-06,3.89433e-10,-3.74003e-14,13750.9,-104.498], Tmin=(939.776,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'C=CC[C](C=O)C=C[O](3054)',
    structure = SMILES('C=CC[C](C=O)C=C[O]'),
    E0 = (2.59555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.300219,0.0653177,-1.29065e-05,-3.36545e-08,1.83097e-11,459.639,32.5193], Tmin=(100,'K'), Tmax=(989.687,'K')), NASAPolynomial(coeffs=[18.9535,0.0266551,-9.97421e-06,1.86769e-09,-1.35457e-13,-5031.26,-66.366], Tmin=(989.687,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.59555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(C)C=O) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]CC(C=O)C=C[O](3055)',
    structure = SMILES('C=[C]CC(C=O)C=C[O]'),
    E0 = (120.737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,357.798,357.806,358.05,358.381],'cm^-1')),
        HinderedRotor(inertia=(0.158523,'amu*angstrom^2'), symmetry=1, barrier=(14.4152,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158589,'amu*angstrom^2'), symmetry=1, barrier=(14.413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158449,'amu*angstrom^2'), symmetry=1, barrier=(14.4146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158549,'amu*angstrom^2'), symmetry=1, barrier=(14.416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00378177,0.0840878,-7.72703e-05,3.70533e-08,-7.16909e-12,14668.7,34.6268], Tmin=(100,'K'), Tmax=(1237.06,'K')), NASAPolynomial(coeffs=[16.3964,0.031058,-1.29685e-05,2.40001e-09,-1.65902e-13,10611.1,-47.9833], Tmin=(1237.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C=CCC([C]=C[O])C=O(3056)',
    structure = SMILES('C=CCC([C]=C[O])C=O'),
    E0 = (120.737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,357.97,357.985,357.988,358.006],'cm^-1')),
        HinderedRotor(inertia=(0.158514,'amu*angstrom^2'), symmetry=1, barrier=(14.4148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158518,'amu*angstrom^2'), symmetry=1, barrier=(14.4146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158513,'amu*angstrom^2'), symmetry=1, barrier=(14.4148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158489,'amu*angstrom^2'), symmetry=1, barrier=(14.4146,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00378152,0.0840878,-7.72703e-05,3.70532e-08,-7.16909e-12,14668.7,34.6268], Tmin=(100,'K'), Tmax=(1237.06,'K')), NASAPolynomial(coeffs=[16.3964,0.031058,-1.29685e-05,2.40001e-09,-1.65902e-13,10611.1,-47.9834], Tmin=(1237.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'C=CCC([C]=O)C=C[O](3057)',
    structure = SMILES('C=CCC([C]=O)C=C[O]'),
    E0 = (41.5836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,357.764,357.764,357.764,357.764],'cm^-1')),
        HinderedRotor(inertia=(0.169986,'amu*angstrom^2'), symmetry=1, barrier=(15.4395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169986,'amu*angstrom^2'), symmetry=1, barrier=(15.4395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169986,'amu*angstrom^2'), symmetry=1, barrier=(15.4395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169986,'amu*angstrom^2'), symmetry=1, barrier=(15.4395,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.1815,0.0812992,-6.53196e-05,2.29495e-08,-2.04786e-12,5161.01,35.949], Tmin=(100,'K'), Tmax=(1084.23,'K')), NASAPolynomial(coeffs=[18.7269,0.0265367,-1.03031e-05,1.87731e-09,-1.30357e-13,179.411,-60.8668], Tmin=(1084.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.5836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH]=CCC(C=O)C=C[O](3058)',
    structure = SMILES('[CH]=CCC(C=O)C=C[O]'),
    E0 = (129.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,298.806,298.953,298.969],'cm^-1')),
        HinderedRotor(inertia=(0.24838,'amu*angstrom^2'), symmetry=1, barrier=(15.7495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248654,'amu*angstrom^2'), symmetry=1, barrier=(15.7493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248444,'amu*angstrom^2'), symmetry=1, barrier=(15.7497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248322,'amu*angstrom^2'), symmetry=1, barrier=(15.749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.270055,0.0856176,-7.80044e-05,3.62986e-08,-6.7139e-12,15795,35.5149], Tmin=(100,'K'), Tmax=(1305.06,'K')), NASAPolynomial(coeffs=[18.8203,0.0271053,-1.07515e-05,1.9433e-09,-1.32681e-13,10812.2,-61.6675], Tmin=(1305.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC(C=O)C=[C][O](3059)',
    structure = SMILES('C=CCC([CH][C]=O)C=O'),
    E0 = (100.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.583247,0.0802365,-8.2429e-05,5.00634e-08,-1.30114e-11,12223,36.3349], Tmin=(100,'K'), Tmax=(909.367,'K')), NASAPolynomial(coeffs=[9.30428,0.0418756,-1.91526e-05,3.67477e-09,-2.58409e-13,10636.9,-4.91052], Tmin=(909.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C[CH]C(C=O)C=CO(3060)',
    structure = SMILES('[CH]C=CC(C=O)C=CO'),
    E0 = (98.7367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.690108,0.0937147,-8.6202e-05,4.11474e-08,-7.82454e-12,12052.1,35.2102], Tmin=(100,'K'), Tmax=(1272.31,'K')), NASAPolynomial(coeffs=[19.4622,0.0303584,-1.15081e-05,2.00943e-09,-1.34255e-13,6924.03,-66.8665], Tmin=(1272.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.7367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC1C([O])[CH]C1C=O(2781)',
    structure = SMILES('C=CC1C([O])[CH]C1C=O'),
    E0 = (231.575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808395,0.0520043,1.88893e-05,-6.24743e-08,2.75004e-11,27983.2,32.7604], Tmin=(100,'K'), Tmax=(986.324,'K')), NASAPolynomial(coeffs=[17.4412,0.0279424,-1.05077e-05,1.99898e-09,-1.46972e-13,22591.4,-57.9539], Tmin=(986.324,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = 'C=CC1O[CH]C1C=C[O](2914)',
    structure = SMILES('C=CC1O[CH]C1C=C[O]'),
    E0 = (156.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.229461,0.0571972,3.88677e-05,-1.14847e-07,5.63283e-11,18958.5,31.2771], Tmin=(100,'K'), Tmax=(894.374,'K')), NASAPolynomial(coeffs=[26.7676,0.00842785,3.39502e-06,-9.9574e-10,7.06634e-14,11415,-109.425], Tmin=(894.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(CCsJOCs) + radical(C=COJ)"""),
)

species(
    label = 'C=CC([CH]C=O)C=C[O](3062)',
    structure = SMILES('C=CC(C=C[O])C=C[O]'),
    E0 = (30.2956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02111,'amu*angstrom^2'), symmetry=1, barrier=(23.4774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01821,'amu*angstrom^2'), symmetry=1, barrier=(23.4107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02047,'amu*angstrom^2'), symmetry=1, barrier=(23.4626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.270555,0.0769583,-3.23369e-05,-2.61716e-08,1.92402e-11,3812.94,32.8179], Tmin=(100,'K'), Tmax=(948.157,'K')), NASAPolynomial(coeffs=[23.4749,0.0182827,-5.16431e-06,8.84808e-10,-6.53016e-14,-2555.36,-90.3133], Tmin=(948.157,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.2956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'C=CC([CH][O])C=CC=O(2786)',
    structure = SMILES('C=CC([CH][O])C=CC=O'),
    E0 = (177.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,210.331,218.61,224.227,241.853],'cm^-1')),
        HinderedRotor(inertia=(0.13923,'amu*angstrom^2'), symmetry=1, barrier=(5.79174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245669,'amu*angstrom^2'), symmetry=1, barrier=(5.64879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.410593,'amu*angstrom^2'), symmetry=1, barrier=(14.7571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131194,'amu*angstrom^2'), symmetry=1, barrier=(5.71066,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.13059,0.10302,-0.000162166,1.52588e-07,-5.72235e-11,21529.4,33.4214], Tmin=(100,'K'), Tmax=(801.784,'K')), NASAPolynomial(coeffs=[5.27823,0.0534804,-2.72893e-05,5.35434e-09,-3.7537e-13,21387.1,13.043], Tmin=(801.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=C[CH]C(C=O)CC=O(3063)',
    structure = SMILES('[CH]C=CC(C=O)CC=O'),
    E0 = (98.0825,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35819,0.070844,-6.95888e-06,-1.28581e-07,1.32455e-10,11878.4,30.487], Tmin=(100,'K'), Tmax=(449.643,'K')), NASAPolynomial(coeffs=[5.49977,0.0552847,-2.60563e-05,5.0069e-09,-3.50364e-13,11290.8,11.4241], Tmin=(449.643,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.0825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(AllylJ2_triplet)"""),
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
    label = 'Ar(8)',
    structure = SMILES('[Ar]'),
    E0 = (-6.19426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,-745,4.3663], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,-745,4.3663], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.19426,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: BurkeH2O2"""),
)

transitionState(
    label = 'TS1',
    E0 = (269.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (246.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (168.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (181.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (91.6899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (225.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (335.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (354.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (45.9873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (45.9873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (38.7961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (60.2391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (38.7961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (60.2391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (48.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (48.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (29.1031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (28.5453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (28.5453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (28.5453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (121.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (100.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (153.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (100.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (228.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (282.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (269.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (267.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (154.693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (134.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (140.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (117.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (150.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (132.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (168.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (202.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (246.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (234.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (123.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (120.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (106.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (335.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (361.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (333.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (333.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (108.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (534.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (578.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (160.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (175.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (164.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (195.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (258.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (162.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (161.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (161.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (155.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (148.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (197.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (177.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (184.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (278.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (262.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (166.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (275.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (144.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (250.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (231.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (156.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (199.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (300.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (142.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction427',
    reactants = ['H(3)(3)', 'C=C=CC(C=O)C=C[O](3013)'],
    products = ['S(754)(753)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction428',
    reactants = ['H(3)(3)', '[CH2]C=CC(C=O)C=C=O(3014)'],
    products = ['S(754)(753)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction459',
    reactants = ['H(3)(3)', '[CH2]C=CC(C=O)=CC=O(3042)'],
    products = ['S(754)(753)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(454.257,'m^3/(mol*s)'), n=1.51715, Ea=(17.7318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-TwoDe_Cds;HJ] for rate rule [Cds-CdCO_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction460',
    reactants = ['S(525)(524)', '[CH]=C[CH2](891)'],
    products = ['S(754)(753)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.53019,'m^3/(mol*s)'), n=1.97633, Ea=(9.23266,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-OneDeH_Cds;CJ] + [Cds-COH_Cds;YJ] for rate rule [Cds-COH_Cds;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction461',
    reactants = ['HCO(14)(15)', 'C6H7O(506)(505)'],
    products = ['S(754)(753)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00435403,'m^3/(mol*s)'), n=2.44862, Ea=(6.74833,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-OneDeH;CJ] for rate rule [Cds-CdH_Cds-COH;CO_pri_rad]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction475',
    reactants = ['CHCHO(45)(45)', 'C5H6O(217)(216)'],
    products = ['S(754)(753)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.7651,'m^3/(mol*s)'), n=1.97633, Ea=(9.23266,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-OneDeH_Cds;CJ] + [Cds-COH_Cds;YJ] for rate rule [Cds-COH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction441',
    reactants = ['[CH]=C[CH2](891)', '[O][CH]C=CC=O(2788)'],
    products = ['S(754)(753)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.07337e+06,'m^3/(mol*s)'), n=-0.0618178, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/TwoDe] for rate rule [Cd_rad;C_rad/H/TwoDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -8.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction442',
    reactants = ['CHCHO(45)(45)', 'C5H6O(219)(218)'],
    products = ['S(754)(753)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.53668e+06,'m^3/(mol*s)'), n=-0.0618178, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/TwoDe;Y_rad] for rate rule [C_rad/H/TwoDe;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -8.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction449',
    reactants = ['S(754)(753)'],
    products = ['C=C=CC(C=O)C=CO(3032)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction450',
    reactants = ['S(754)(753)'],
    products = ['CC=CC(C=O)C=C=O(3033)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction470',
    reactants = ['S(754)(753)'],
    products = ['CC=CC(C=O)=CC=O(3051)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction471',
    reactants = ['S(754)(753)'],
    products = ['C=C=CC(C=O)CC=O(3052)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction485',
    reactants = ['S(754)(753)'],
    products = ['C=CC=C(C=O)C=CO(3061)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction486',
    reactants = ['S(754)(753)'],
    products = ['C=CCC(C=O)C=C=O(2844)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction490',
    reactants = ['S(754)(753)'],
    products = ['C=CCC(C=O)=CC=O(3064)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction491',
    reactants = ['S(754)(753)'],
    products = ['C=CC=C(C=O)CC=O(3065)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction452',
    reactants = ['S(754)(753)'],
    products = ['O=CC1C=CCOC=C1(3035)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Ypri_rad_out] for rate rule [R7;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction473',
    reactants = ['S(754)(753)'],
    products = ['S(790)(789)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] for rate rule [R5_SSDS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction408',
    reactants = ['S(754)(753)'],
    products = ['S(651)(650)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R5_SSDS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction492',
    reactants = ['S(754)(753)'],
    products = ['S(756)(755)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] for rate rule [R3_SS;C_rad_out_H/OneDe;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction423',
    reactants = ['S(754)(753)'],
    products = ['[O][CH]C1CC=CC1C=O(3010)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(100.862,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_SMS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 98.8 to 100.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction424',
    reactants = ['S(754)(753)'],
    products = ['[O]C=CC1C=CCC1[O](3011)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(79.3435,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_SMS_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction425',
    reactants = ['S(754)(753)'],
    products = ['C=CC1OC=CC1[CH][O](2894)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(132.128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS_D;doublebond_intra;radadd_intra] for rate rule [R6_SMS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 129.1 to 132.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction426',
    reactants = ['S(754)(753)'],
    products = ['[CH2]C=CC1C=COC1[O](3012)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(79.3435,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra] for rate rule [R6_SMS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction429',
    reactants = ['S(754)(753)'],
    products = ['C[C]=CC(C=O)C=C[O](3015)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction430',
    reactants = ['[CH2]C=CC(C=O)C=[C]O(3016)'],
    products = ['S(754)(753)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction431',
    reactants = ['CC=[C]C(C=O)C=C[O](3017)'],
    products = ['S(754)(753)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction432',
    reactants = ['[CH2]C=CC([C]=CO)C=O(3018)'],
    products = ['S(754)(753)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction433',
    reactants = ['S(754)(753)'],
    products = ['CC=C[C](C=O)C=C[O](3019)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(184752,'s^-1'), n=1.905, Ea=(133.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_2H;Cs_H_out] for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction434',
    reactants = ['S(754)(753)'],
    products = ['[CH2]C=C[C](C=O)C=CO(3020)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.11e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4H_SDS;O_rad_out;Cs_H_out] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction435',
    reactants = ['CC=CC([C]=C[O])C=O(3021)'],
    products = ['S(754)(753)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSMS;Cd_rad_out;Cs_H_out_2H] for rate rule [R5H_SSMS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction436',
    reactants = ['CC=CC([C]=O)C=C[O](3022)'],
    products = ['S(754)(753)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.45388e+06,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_SSMS;CO_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction437',
    reactants = ['[CH2]C=[C]C(C=O)C=CO(3023)'],
    products = ['S(754)(753)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSMS;Cd_rad_out;XH_out] for rate rule [R5H_SSMS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction438',
    reactants = ['[CH2]C=CC([C]=O)C=CO(3024)'],
    products = ['S(754)(753)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(722272,'s^-1'), n=1.6737, Ea=(94.6126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;CO_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction439',
    reactants = ['CC=CC(C=O)C=[C][O](3025)'],
    products = ['S(754)(753)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(37753.8,'s^-1'), n=1.925, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_RSSMS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction440',
    reactants = ['[CH2][C]=CC(C=O)C=CO(3026)'],
    products = ['S(754)(753)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(16140,'s^-1'), n=1.92259, Ea=(84.7802,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction443',
    reactants = ['S(754)(753)'],
    products = ['[O]C=CC(C=O)C1[CH]C1(3027)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction444',
    reactants = ['S(754)(753)'],
    products = ['[CH2]C=CC(C=O)C1[CH]O1(3028)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(9.85157e+11,'s^-1'), n=0.224969, Ea=(213.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;doublebond_intra_pri_HNd_Cs;radadd_intra] + [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction445',
    reactants = ['S(754)(753)'],
    products = ['[O]C1[CH]C(C=O)C=CC1(3029)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(5.041e+08,'s^-1'), n=0.7, Ea=(102.034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_SDS_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 97.8 to 102.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction446',
    reactants = ['S(754)(753)'],
    products = ['[O]C=CC1[CH]OCC=C1(3030)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.39072e+10,'s^-1'), n=0.346137, Ea=(99.0339,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_SMS_CO;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction447',
    reactants = ['S(754)(753)'],
    products = ['[CH2]C1[CH]C(C=O)C=CO1(3031)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(5.041e+08,'s^-1'), n=0.7, Ea=(85.1135,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SDS_D;doublebond_intra_pri;radadd_intra] for rate rule [R6_SDS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 84.7 to 85.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction448',
    reactants = ['S(754)(753)'],
    products = ['[CH2]C=CC1[CH]OOC=C1(2747)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.39072e+10,'s^-1'), n=0.346137, Ea=(314.821,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra] for rate rule [R6_SMS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 313.0 to 314.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction451',
    reactants = ['CO(10)(11)', '[CH2]C=CCC=C[O](3034)'],
    products = ['S(754)(753)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.532e+06,'cm^3/(mol*s)'), n=2.07, Ea=(343.925,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO;C_sec] for rate rule [CO;C/H2/TwoDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction363',
    reactants = ['S(650)(649)'],
    products = ['S(754)(753)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction453',
    reactants = ['[CH2]C=CC=COC=C[O](3036)'],
    products = ['S(754)(753)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction454',
    reactants = ['[CH2]C=CC(C=C[O])=CO(3037)'],
    products = ['S(754)(753)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction455',
    reactants = ['CH2(17)(18)', '[CH]=CC(C=O)C=C[O](3038)'],
    products = ['S(754)(753)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction456',
    reactants = ['O(4)(4)', '[CH]=CC(C=O)C=C[CH2](3039)'],
    products = ['S(754)(753)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction457',
    reactants = ['S(754)(753)'],
    products = ['[CH2][CH]C1C(C=O)C1C=O(3040)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(139.565,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHDe]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 137.4 to 139.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction458',
    reactants = ['S(754)(753)'],
    products = ['[CH2]C=CC1C([O])C1C=O(3041)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(2.90568e+10,'s^-1'), n=0.237, Ea=(154.332,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHDe] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_csHDe]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 151.9 to 154.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction462',
    reactants = ['S(754)(753)'],
    products = ['[CH2]C=C[C](C=O)CC=O(3043)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(2.33935e+09,'s^-1'), n=1.12899, Ea=(143.262,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/OneDe;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction463',
    reactants = ['[CH2]C=CC(C=O)C[C]=O(3044)'],
    products = ['S(754)(753)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction464',
    reactants = ['[CH2]C=[C]C(C=O)CC=O(3045)'],
    products = ['S(754)(753)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction465',
    reactants = ['[CH2]C=CC([C]=O)CC=O(3046)'],
    products = ['S(754)(753)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;CO_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction466',
    reactants = ['[CH2][C]=CC(C=O)CC=O(3047)'],
    products = ['S(754)(753)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_DSS;Cd_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction467',
    reactants = ['S(754)(753)'],
    products = ['[CH2]C1[CH]C(C=O)C1C=O(3048)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(2.93362e+07,'s^-1'), n=1.19915, Ea=(140.341,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R4_Cs_RR_D;doublebond_intra_pri;radadd_intra_csHCO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 136.8 to 140.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction468',
    reactants = ['S(754)(753)'],
    products = ['[CH2]C=CC1[CH]OC1C=O(3049)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(6.89861e+07,'s^-1'), n=1.13751, Ea=(134.212,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_csHCO]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic
Ea raised from 130.9 to 134.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction469',
    reactants = ['S(754)(753)'],
    products = ['O=CC1[CH][CH]OCC=C1(3050)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(2.97949e+09,'s^-1'), n=0.649948, Ea=(127.33,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_linear;multiplebond_intra;radadd_intra_cs2H] + [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 120.3 to 127.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction472',
    reactants = ['[CH2]C=C[CH]C(C=O)C=O(3053)'],
    products = ['S(754)(753)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.7778e+12,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ-CdH;CO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction474',
    reactants = ['S(754)(753)'],
    products = ['C=CC1C([O])C1[CH]C=O(2760)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(156.119,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCd] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 152.8 to 156.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction476',
    reactants = ['S(754)(753)'],
    products = ['C=CC[C](C=O)C=C[O](3054)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(7.95861e+10,'s^-1'), n=0.595417, Ea=(163.653,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction477',
    reactants = ['C=[C]CC(C=O)C=C[O](3055)'],
    products = ['S(754)(753)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction478',
    reactants = ['C=CCC([C]=C[O])C=O(3056)'],
    products = ['S(754)(753)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction479',
    reactants = ['C=CCC([C]=O)C=C[O](3057)'],
    products = ['S(754)(753)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction480',
    reactants = ['[CH]=CCC(C=O)C=C[O](3058)'],
    products = ['S(754)(753)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction481',
    reactants = ['C=CCC(C=O)C=[C][O](3059)'],
    products = ['S(754)(753)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction482',
    reactants = ['[CH]=C[CH]C(C=O)C=CO(3060)'],
    products = ['S(754)(753)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction483',
    reactants = ['S(754)(753)'],
    products = ['C=CC1C([O])[CH]C1C=O(2781)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(2.93362e+07,'s^-1'), n=1.19915, Ea=(210.561,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R4_Cs_RR_D;doublebond_intra_pri;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 206.5 to 210.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction484',
    reactants = ['S(754)(753)'],
    products = ['C=CC1O[CH]C1C=C[O](2914)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(6.89861e+07,'s^-1'), n=1.13751, Ea=(135.285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_csHCd]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic
Ea raised from 132.6 to 135.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction487',
    reactants = ['C=CC([CH]C=O)C=C[O](3062)'],
    products = ['S(754)(753)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(7.33094e+10,'s^-1'), n=0.732, Ea=(169.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ-OneDeH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction488',
    reactants = ['C=CC([CH][O])C=CC=O(2786)'],
    products = ['S(754)(753)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ-CdH;CO]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction489',
    reactants = ['[CH]=C[CH]C(C=O)CC=O(3063)'],
    products = ['S(754)(753)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R5HJ_2;Cd_rad_out_singleH;Cs_H_out_H/CO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = '177',
    isomers = [
        'S(754)(753)',
    ],
    reactants = [
        ('HCO(14)(15)', 'C6H7O(506)(505)'),
        ('CHCHO(45)(45)', 'C5H6O(217)(216)'),
        ('CHCHO(45)(45)', 'C5H6O(219)(218)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '177',
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

