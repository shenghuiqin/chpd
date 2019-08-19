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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.319136,0.089632,-9.48754e-05,4.97278e-08,-9.46921e-12,84244.6,28.8601], Tmin=(100,'K'), Tmax=(948.213,'K')), NASAPolynomial(coeffs=[18.251,0.0228745,-7.58846e-06,1.2375e-09,-8.01996e-14,80202.3,-62.4874], Tmin=(948.213,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(699.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_P)"""),
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
    label = '[CH]=C(C=C)C1CC1=[CH](28233)',
    structure = SMILES('[CH]=C(C=C)C1CC1=[CH]'),
    E0 = (803.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0760444,0.0753784,-5.30727e-05,8.09982e-09,4.26092e-12,96735,26.1351], Tmin=(100,'K'), Tmax=(984.787,'K')), NASAPolynomial(coeffs=[18.7664,0.0226048,-7.93948e-06,1.40927e-09,-9.85205e-14,91931.7,-69.446], Tmin=(984.787,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(803.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = 'C#CC([CH2])C1=CC1[CH2](28234)',
    structure = SMILES('C#CC([CH2])C1=CC1[CH2]'),
    E0 = (802.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0118463,0.0895393,-0.000102845,6.03915e-08,-1.24611e-11,96697.8,26.959], Tmin=(100,'K'), Tmax=(785.408,'K')), NASAPolynomial(coeffs=[14.3224,0.0288025,-1.02742e-05,1.70069e-09,-1.08917e-13,94067.8,-41.1416], Tmin=(785.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(802.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropene) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH]=C1C([CH2])CC1C#C(28235)',
    structure = SMILES('[CH]=C1C([CH2])CC1C#C'),
    E0 = (734.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.7098,0.00121743,0.000103239,-1.06398e-07,2.44459e-11,87911.7,-18.5579], Tmin=(100,'K'), Tmax=(1768.54,'K')), NASAPolynomial(coeffs=[85.2364,0.0166553,-6.59147e-05,1.61936e-08,-1.19961e-12,32776.1,-501.95], Tmin=(1768.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(734.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(methylenecyclobutane) + radical(Isobutyl) + radical(Cds_P)"""),
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
    label = '[CH]=C(C=C)C(=C)C#C(28236)',
    structure = SMILES('[CH]=C(C=C)C(=C)C#C'),
    E0 = (641.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,325,375,415,465,420,450,1700,1750,180],'cm^-1')),
        HinderedRotor(inertia=(1.4883,'amu*angstrom^2'), symmetry=1, barrier=(34.219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48362,'amu*angstrom^2'), symmetry=1, barrier=(34.1113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48232,'amu*angstrom^2'), symmetry=1, barrier=(34.0814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.577618,0.0631255,-2.89282e-05,-1.17899e-08,1.00793e-11,77273.7,25.2657], Tmin=(100,'K'), Tmax=(992.133,'K')), NASAPolynomial(coeffs=[17.1672,0.0229133,-8.45758e-06,1.55196e-09,-1.10668e-13,72669.2,-61.2543], Tmin=(992.133,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(641.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)Ct) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Cds_P)"""),
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
    label = '[CH]=C(C=C)C=C(26631)',
    structure = SMILES('[CH]=C(C=C)C=C'),
    E0 = (411.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.17526,'amu*angstrom^2'), symmetry=1, barrier=(27.0215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17459,'amu*angstrom^2'), symmetry=1, barrier=(27.0061,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75262,0.0382342,7.41856e-06,-3.70507e-08,1.69951e-11,49629.8,19.502], Tmin=(100,'K'), Tmax=(974.113,'K')), NASAPolynomial(coeffs=[12.41,0.0212121,-7.54597e-06,1.37118e-09,-9.78936e-14,46284.8,-38.1462], Tmin=(974.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
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
    label = 'C#CC([CH2])C#C(27501)',
    structure = SMILES('C#CC([CH2])C#C'),
    E0 = (605.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(1.10375,'amu*angstrom^2'), symmetry=1, barrier=(60.786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.494061,'amu*angstrom^2'), symmetry=1, barrier=(27.2294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10575,'amu*angstrom^2'), symmetry=1, barrier=(60.6988,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34005,0.0508861,-4.74747e-05,2.33958e-08,-4.52454e-12,72935.7,20.5937], Tmin=(100,'K'), Tmax=(1340.28,'K')), NASAPolynomial(coeffs=[12.7012,0.0146548,-4.32414e-06,6.38341e-10,-3.82616e-14,70099,-36.7658], Tmin=(1340.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCtCsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C(C=C)[C](C)C#C(28237)',
    structure = SMILES('[CH]C(C=C)=C(C)C#C'),
    E0 = (615.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.178184,0.0817907,-6.86429e-05,3.03063e-08,-5.38513e-12,74140.4,27.8105], Tmin=(100,'K'), Tmax=(1349.59,'K')), NASAPolynomial(coeffs=[17.268,0.0300826,-1.11719e-05,1.91697e-09,-1.26247e-13,69431.3,-61.5877], Tmin=(1349.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(615.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCtCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#C[C]([CH2])C(=C)C=C(28238)',
    structure = SMILES('C#CC([CH2])=C([CH2])C=C'),
    E0 = (483.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.447834,0.0687575,-4.0266e-05,-1.6155e-09,7.44731e-12,58336.7,25.6565], Tmin=(100,'K'), Tmax=(943.612,'K')), NASAPolynomial(coeffs=[15.826,0.0262043,-8.60377e-06,1.43631e-09,-9.63458e-14,54426.8,-52.9815], Tmin=(943.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCtCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CC=CCJ) + radical(CTCC=CCJ)"""),
)

species(
    label = 'C#CC([CH2])C(=C)[C]=C(28239)',
    structure = SMILES('C#CC([CH2])C(=C)[C]=C'),
    E0 = (651.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.231271,0.0914804,-0.000102779,5.86568e-08,-1.21801e-11,78453,27.744], Tmin=(100,'K'), Tmax=(858.56,'K')), NASAPolynomial(coeffs=[16.2807,0.026133,-8.84287e-06,1.4274e-09,-9.05734e-14,75190.9,-51.8843], Tmin=(858.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJC=C) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C([C]=C)C(C)C#C(28240)',
    structure = SMILES('[CH]C(=C=C)C(C)C#C'),
    E0 = (668.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.095791,0.0868248,-8.38117e-05,4.51141e-08,-1.00852e-11,80594.6,26.488], Tmin=(100,'K'), Tmax=(1064.07,'K')), NASAPolynomial(coeffs=[12.7512,0.0392511,-1.6748e-05,3.09705e-09,-2.13372e-13,77901.3,-35.3531], Tmin=(1064.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(668.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC(=C)C([CH2])C#C(28241)',
    structure = SMILES('[CH]=CC(=C)C([CH2])C#C'),
    E0 = (699.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2175,525,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.05344,'amu*angstrom^2'), symmetry=1, barrier=(24.2206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06108,'amu*angstrom^2'), symmetry=1, barrier=(24.3962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05771,'amu*angstrom^2'), symmetry=1, barrier=(24.3188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05364,'amu*angstrom^2'), symmetry=1, barrier=(24.2254,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.319127,0.0896319,-9.4875e-05,4.97273e-08,-9.46896e-12,84244.6,28.86], Tmin=(100,'K'), Tmax=(948.208,'K')), NASAPolynomial(coeffs=[18.2509,0.0228745,-7.58849e-06,1.23751e-09,-8.02002e-14,80202.3,-62.4873], Tmin=(948.208,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(699.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = '[C]#CC(C)C(=[CH])C=C(28242)',
    structure = SMILES('[C]#CC(C)C(=[CH])C=C'),
    E0 = (831.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2175,525,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.920724,'amu*angstrom^2'), symmetry=1, barrier=(21.1692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.920081,'amu*angstrom^2'), symmetry=1, barrier=(21.1545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.920516,'amu*angstrom^2'), symmetry=1, barrier=(21.1645,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.920805,'amu*angstrom^2'), symmetry=1, barrier=(21.1711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.252118,0.0907497,-0.000102339,6.07282e-08,-1.42187e-11,100123,27.5809], Tmin=(100,'K'), Tmax=(1047.06,'K')), NASAPolynomial(coeffs=[16.7076,0.0259585,-9.51817e-06,1.62746e-09,-1.07263e-13,96571.6,-55.0193], Tmin=(1047.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(831.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(Acetyl)"""),
)

species(
    label = '[C]#CC([CH2])C(=C)C=C(28243)',
    structure = SMILES('[C]#CC([CH2])C(=C)C=C'),
    E0 = (789.163,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.662229,'amu*angstrom^2'), symmetry=1, barrier=(15.2259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.662181,'amu*angstrom^2'), symmetry=1, barrier=(15.2249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.13215,'amu*angstrom^2'), symmetry=1, barrier=(95.0062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.13885,'amu*angstrom^2'), symmetry=1, barrier=(95.1603,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.131126,0.0885105,-9.64816e-05,5.31448e-08,-1.05078e-11,95065.5,28.5287], Tmin=(100,'K'), Tmax=(873.094,'K')), NASAPolynomial(coeffs=[16.1606,0.0255955,-8.53421e-06,1.37121e-09,-8.7027e-14,91773.8,-50.4178], Tmin=(873.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(789.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=CC(=[CH])C(C)C#C(28244)',
    structure = SMILES('[CH]=CC(=[CH])C(C)C#C'),
    E0 = (741.13,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.04167,'amu*angstrom^2'), symmetry=1, barrier=(23.9501,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04277,'amu*angstrom^2'), symmetry=1, barrier=(23.9754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04381,'amu*angstrom^2'), symmetry=1, barrier=(23.9992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04399,'amu*angstrom^2'), symmetry=1, barrier=(24.0033,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.424057,0.0916201,-9.95214e-05,5.51683e-08,-1.19703e-11,89301.5,27.8582], Tmin=(100,'K'), Tmax=(1131.04,'K')), NASAPolynomial(coeffs=[19.0067,0.0229022,-8.38716e-06,1.4514e-09,-9.70089e-14,84906.1,-68.2765], Tmin=(1131.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(Cds_P)"""),
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
    label = 'C#CC([CH2])C1[CH]CC=1(28245)',
    structure = SMILES('C#CC([CH2])C1[CH]CC=1'),
    E0 = (665.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.663781,0.060553,-1.5815e-05,-2.66011e-08,1.58489e-11,80190.3,24.5905], Tmin=(100,'K'), Tmax=(959.819,'K')), NASAPolynomial(coeffs=[16.4445,0.0253352,-8.51645e-06,1.48817e-09,-1.04113e-13,75753.8,-58.2252], Tmin=(959.819,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(665.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene) + radical(Isobutyl) + radical(cyclobutene-allyl)"""),
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
    label = '[CH]=C1[CH]CCC1C#C(28246)',
    structure = SMILES('[CH]C1=CCCC1C#C'),
    E0 = (566.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.657797,0.0587742,-2.37209e-07,-4.02186e-08,1.91421e-11,68313.1,23.4987], Tmin=(100,'K'), Tmax=(993.964,'K')), NASAPolynomial(coeffs=[15.4165,0.0331881,-1.26438e-05,2.32173e-09,-1.6428e-13,63709.2,-56.0149], Tmin=(993.964,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(566.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopentene) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#CC(=C)C(=C)C=C(28247)',
    structure = SMILES('C#CC(=C)C(=C)C=C'),
    E0 = (394.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.654562,0.0595643,-1.21309e-05,-2.93606e-08,1.60785e-11,47553.9,24.5864], Tmin=(100,'K'), Tmax=(987.414,'K')), NASAPolynomial(coeffs=[16.911,0.0257735,-9.50743e-06,1.75447e-09,-1.25847e-13,42780.4,-61.5506], Tmin=(987.414,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)Ct) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = '[CH]=C(C=C)C[CH]C#C(26591)',
    structure = SMILES('[CH]=C(C=C)C[CH]C#C'),
    E0 = (673.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.4191,0.0846785,-8.48777e-05,4.40797e-08,-8.84377e-12,81151.1,30.4148], Tmin=(100,'K'), Tmax=(1330.66,'K')), NASAPolynomial(coeffs=[19.0101,0.0201252,-5.17825e-06,6.77401e-10,-3.70899e-14,76524.7,-66.8249], Tmin=(1330.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(673.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]CC#C)C=C(28248)',
    structure = SMILES('[CH]C(C=C)=CCC#C'),
    E0 = (642.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0858222,0.0759763,-5.00392e-05,9.51664e-09,2.26662e-12,77484.2,28.3844], Tmin=(100,'K'), Tmax=(1028,'K')), NASAPolynomial(coeffs=[16.1958,0.0315835,-1.19539e-05,2.12679e-09,-1.45565e-13,73205.5,-54.4826], Tmin=(1028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(642.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
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
    E0 = (699.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (803.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (802.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (752.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (758.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (870.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (743.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (982.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (912.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (797.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (891.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (921.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (809.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (1102.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (927.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (862.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (774.17,'kJ/mol'),
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
    E0 = (829.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (837.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (749.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (813.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (762.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (859.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (859.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (707.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (703.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (1037.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['CH2CHCCH(26391)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['[CH]=C(C=C)C1CC1=[CH](28233)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(103.929,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for R4_S_T;triplebond_intra_H;radadd_intra_cs2H
Exact match found for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 102.0 to 103.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['C#CC([CH2])C1=CC1[CH2](28234)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(103.848,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 14 used for R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH
Exact match found for rate rule [R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['[CH]=C1C=C(C=C)C1[CH2](28032)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['[CH]=C1C([CH2])CC1C#C(28235)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.07e+06,'s^-1'), n=1.46, Ea=(58.9944,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 343 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]=C(C=C)C(=C)C#C(28236)'],
    products = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(454.257,'m^3/(mol*s)'), n=1.51715, Ea=(17.7318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-TwoDe_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]C=C(4699)', 'CH2CHCCH(26391)'],
    products = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00294841,'m^3/(mol*s)'), n=2.48333, Ea=(17.6885,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CtH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2H(33)', '[CH]=C(C=C)C=C(26631)'],
    products = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0101073,'m^3/(mol*s)'), n=2.41967, Ea=(13.2813,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeH_Cds;CJ] for rate rule [Cds-OneDeH_Cds;CtJ_Ct]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C2H3(30)', 'C#CC([CH2])C#C(27501)'],
    products = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(189200,'cm^3/(mol*s)'), n=2.41, Ea=(20.7945,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2340 used for Ct-Cs_Ct-H;CdsJ-H
Exact match found for rate rule [Ct-Cs_Ct-H;CdsJ-H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['[CH]=C(C=C)[C](C)C#C(28237)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['C#C[C]([CH2])C(=C)C=C(28238)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['C#CC([CH2])C(=C)[C]=C(28239)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['[CH]=C([C]=C)C(C)C#C(28240)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(356.546,'s^-1'), n=2.79552, Ea=(110.814,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;Cd_H_out_doubleC] + [R4H_SSS;C_rad_out_2H;XH_out] for rate rule [R4H_SSS;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC(=C)C([CH2])C#C(28241)'],
    products = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.27529e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_singleH] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[C]#CC(C)C(=[CH])C=C(28242)'],
    products = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.39293e+07,'s^-1'), n=1.32074, Ea=(96.0416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_TSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[C]#CC([CH2])C(=C)C=C(28243)'],
    products = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(544115,'s^-1'), n=1.7475, Ea=(73.8225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cd_H_out_singleH] for rate rule [R5H_TSSD;Ct_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=CC(=[CH])C(C)C#C(28244)'],
    products = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]C=C(4699)', '[CH]=[C]C=C(4699)'],
    products = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['[CH]=C(C=C)C1[C]=CC1(28058)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(130.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['C#CC([CH2])C1[CH]CC=1(28245)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R4_D_D;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['[CH2]C1[C]=CC=C1C=C(28006)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.88e+10,'s^-1'), n=0.31, Ea=(50.6264,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleH
Exact match found for rate rule [R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['[CH]=C1[CH]CCC1C#C(28246)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.12e+09,'s^-1'), n=0.63, Ea=(114.642,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['C#CC(=C)C(=C)C=C(28247)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['[CH]=C(C=C)C[CH]C#C(26591)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['[CH]=C([CH]CC#C)C=C(28248)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    products = ['C#CC1CC=C1C=C(28249)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#C[CH]CC=[C]C=C(26590)'],
    products = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.213e+11,'s^-1'), n=0.07, Ea=(76.6885,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for 1_4_5_hexatriene
Exact match found for rate rule [1_4_5_hexatriene]
Euclidian distance = 0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(19)', '[CH]=C=CC(=[CH])C=C(28060)'],
    products = ['[CH]=C(C=C)C([CH2])C#C(26594)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/TwoDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4566',
    isomers = [
        '[CH]=C(C=C)C([CH2])C#C(26594)',
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
    label = '4566',
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

