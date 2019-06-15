species(
    label = 'C1=CC(C1)C1CC1(146)',
    structure = SMILES('C1=CC(C1)C1CC1'),
    E0 = (226.142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3725.22,'J/mol'), sigma=(6.44341,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=581.87 K, Pc=31.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00395,0.0190803,0.000102584,-1.43669e-07,5.52787e-11,27292.8,20.706], Tmin=(100,'K'), Tmax=(962.107,'K')), NASAPolynomial(coeffs=[14.8853,0.0282669,-9.55752e-06,1.81747e-09,-1.38108e-13,21910.4,-56.0326], Tmin=(962.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'cC3H5(120)',
    structure = SMILES('[CH]1CC1'),
    E0 = (280.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,894.423,894.423,894.423,894.423,894.423,894.424,894.425,3153.42],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.52414,-0.00309112,7.19158e-05,-9.14814e-08,3.49781e-11,33782.5,9.37886], Tmin=(100,'K'), Tmax=(937.435,'K')), NASAPolynomial(coeffs=[9.02695,0.0081344,-1.57947e-06,2.7853e-10,-2.51237e-14,31225.9,-24.9471], Tmin=(937.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""cC3H5""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]1C=CC1(258)',
    structure = SMILES('[CH]1C=CC1'),
    E0 = (309.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,718.425,1076.94,1076.94,1076.94,1076.94,1076.94,1076.94,1076.94,1076.94,1076.94,1821.09],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65584,-0.010946,0.00010676,-1.29424e-07,4.82651e-11,37295.1,8.17286], Tmin=(100,'K'), Tmax=(946.825,'K')), NASAPolynomial(coeffs=[9.86439,0.010387,-2.38638e-06,4.80812e-10,-4.35568e-14,33987.5,-32.6988], Tmin=(946.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutene-allyl)"""),
)

species(
    label = 'CH2(S)(28)',
    structure = SMILES('[CH2]'),
    E0 = (418.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1358.21,2621.43,3089.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19331,-0.00233105,8.15676e-06,-6.62986e-09,1.93233e-12,50366.2,-0.746734], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.13502,0.00289594,-8.16668e-07,1.13573e-10,-6.36263e-15,50504.1,4.06031], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(418.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=CC1C=CC1(265)',
    structure = SMILES('C=CC1C=CC1'),
    E0 = (219.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26734,0.01959,6.97114e-05,-1.02576e-07,4.00366e-11,26421.7,17.9956], Tmin=(100,'K'), Tmax=(961.591,'K')), NASAPolynomial(coeffs=[12.7616,0.022543,-7.59802e-06,1.42673e-09,-1.07239e-13,22248.7,-43.4265], Tmin=(961.591,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]1C[C](C1)C1CC1(266)',
    structure = SMILES('[CH]1C[C](C1)C1CC1'),
    E0 = (425.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26239,0.0235711,6.26587e-05,-8.3059e-08,2.908e-11,51206.8,23.4974], Tmin=(100,'K'), Tmax=(1035.64,'K')), NASAPolynomial(coeffs=[7.59327,0.040686,-1.67407e-05,3.20661e-09,-2.30345e-13,48080.6,-12.1699], Tmin=(1035.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Tertalkyl) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1C[CH]C1C1CC1(267)',
    structure = SMILES('[CH]1C[CH]C1C1CC1'),
    E0 = (427.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34594,0.016722,9.11e-05,-1.18079e-07,4.26676e-11,51499.1,24.4074], Tmin=(100,'K'), Tmax=(998.766,'K')), NASAPolynomial(coeffs=[9.89274,0.0367693,-1.45089e-05,2.80976e-09,-2.06585e-13,47484.2,-24.5446], Tmin=(998.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1CC(C1)[C]1CC1(268)',
    structure = SMILES('[CH]1CC(C1)[C]1CC1'),
    E0 = (425.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26239,0.0235711,6.26587e-05,-8.3059e-08,2.908e-11,51206.8,23.4974], Tmin=(100,'K'), Tmax=(1035.64,'K')), NASAPolynomial(coeffs=[7.59327,0.040686,-1.67407e-05,3.20661e-09,-2.30345e-13,48080.6,-12.1699], Tmin=(1035.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Tertalkyl) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1[CH]C(C1)C1CC1(157)',
    structure = SMILES('[CH]1[CH]C(C1)C1CC1'),
    E0 = (427.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34594,0.016722,9.11e-05,-1.18079e-07,4.26676e-11,51499.1,24.4074], Tmin=(100,'K'), Tmax=(998.766,'K')), NASAPolynomial(coeffs=[9.89274,0.0367693,-1.45089e-05,2.80976e-09,-2.06585e-13,47484.2,-24.5446], Tmin=(998.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1CC[C]1C1CC1(269)',
    structure = SMILES('[CH]1CC[C]1C1CC1'),
    E0 = (425.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26239,0.0235711,6.26588e-05,-8.30591e-08,2.908e-11,51206.8,24.1905], Tmin=(100,'K'), Tmax=(1035.64,'K')), NASAPolynomial(coeffs=[7.59326,0.040686,-1.67407e-05,3.20662e-09,-2.30345e-13,48080.6,-11.4767], Tmin=(1035.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutane) + radical(Tertalkyl)"""),
)

species(
    label = '[CH]1CCC1[C]1CC1(270)',
    structure = SMILES('[CH]1CCC1[C]1CC1'),
    E0 = (425.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26239,0.0235711,6.26587e-05,-8.3059e-08,2.908e-11,51206.8,23.4974], Tmin=(100,'K'), Tmax=(1035.64,'K')), NASAPolynomial(coeffs=[7.59327,0.040686,-1.67407e-05,3.20661e-09,-2.30345e-13,48080.6,-12.1699], Tmin=(1035.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Tertalkyl) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1CC(C1)C1[CH]C1(271)',
    structure = SMILES('[CH]1CC(C1)C1[CH]C1'),
    E0 = (465.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34594,0.016722,9.11e-05,-1.18079e-07,4.26676e-11,56078.4,24.4074], Tmin=(100,'K'), Tmax=(998.766,'K')), NASAPolynomial(coeffs=[9.89274,0.0367693,-1.45089e-05,2.80976e-09,-2.06585e-13,52063.5,-24.5446], Tmin=(998.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1CC1C1[CH]CC1(272)',
    structure = SMILES('[CH]1CC1C1[CH]CC1'),
    E0 = (465.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34594,0.016722,9.11e-05,-1.18079e-07,4.26676e-11,56078.4,24.4074], Tmin=(100,'K'), Tmax=(998.766,'K')), NASAPolynomial(coeffs=[9.89274,0.0367693,-1.45089e-05,2.80976e-09,-2.06585e-13,52063.5,-24.5446], Tmin=(998.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclopropane) + radical(cyclobutane)"""),
)

species(
    label = '[CH2]C([CH2])C1C=CC1(273)',
    structure = SMILES('[CH2]C([CH2])C1C=CC1'),
    E0 = (474.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,2950,3050,3150,900,950,1000,1050,1100,1380,1390,370,380,2900,435,600,893.017,893.017,893.017,893.017,893.017,893.017,893.017,893.017,893.017,893.017,893.017,893.017,2336.88],'cm^-1')),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48771,0.0392203,3.99592e-05,-7.78546e-08,3.26208e-11,57191.1,25.5725], Tmin=(100,'K'), Tmax=(948.084,'K')), NASAPolynomial(coeffs=[13.0689,0.0304739,-9.67025e-06,1.67214e-09,-1.17752e-13,53192.2,-39.1907], Tmin=(948.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C[CH]C1C=CC1(78)',
    structure = SMILES('[CH2]C[CH]C1C=CC1'),
    E0 = (464.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3025,407.5,1350,352.5,670.258,905.37,905.37,905.37,905.37,905.37,905.37,905.37,905.37,905.37,905.37,905.37,905.37,2427.95,2427.95],'cm^-1')),
        HinderedRotor(inertia=(0.00366254,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00366254,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00366254,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77537,0.0318221,5.39347e-05,-8.39388e-08,3.16566e-11,55915.7,26.4999], Tmin=(100,'K'), Tmax=(1007.23,'K')), NASAPolynomial(coeffs=[11.6271,0.0350252,-1.38707e-05,2.66208e-09,-1.93799e-13,51784.1,-31.7583], Tmin=(1007.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(464.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(RCCJ) + radical(Cs_S)"""),
)

species(
    label = '[CH]=CC([CH2])C1CC1(274)',
    structure = SMILES('[CH]=CC([CH2])C1CC1'),
    E0 = (529.38,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,3000,3100,440,815,1455,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,1380,1390,370,380,2900,435,600,831.352,831.352,831.352,831.352,831.352,831.352,831.352,831.352,831.352,2391.1],'cm^-1')),
        HinderedRotor(inertia=(0.00664932,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00664932,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00664932,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28042,0.0443502,2.66989e-05,-6.37255e-08,2.70623e-11,63781.5,25.9729], Tmin=(100,'K'), Tmax=(970.743,'K')), NASAPolynomial(coeffs=[13.6751,0.0309276,-1.07385e-05,1.93916e-09,-1.38206e-13,59601.1,-42.5937], Tmin=(970.743,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC[CH]C1CC1(133)',
    structure = SMILES('[CH]=CC[CH]C1CC1'),
    E0 = (516.767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,383.996,981.78,981.78,981.78,981.78,981.78,981.78,981.78,981.78,981.78,981.78,2287.33],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56103,0.037257,4.12569e-05,-7.29286e-08,2.83113e-11,62254.9,26.7717], Tmin=(100,'K'), Tmax=(1007.19,'K')), NASAPolynomial(coeffs=[12.4171,0.0340365,-1.33605e-05,2.54948e-09,-1.84976e-13,58044.6,-35.7257], Tmin=(1007.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(Cds_P)"""),
)

species(
    label = '[C]1CCC1C1CC1(275)',
    structure = SMILES('[C]1CCC1C1CC1'),
    E0 = (461.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92373,0.0213567,9.68216e-05,-1.37703e-07,5.30858e-11,55570,21.5028], Tmin=(100,'K'), Tmax=(964.906,'K')), NASAPolynomial(coeffs=[14.9395,0.0287086,-9.91474e-06,1.89137e-09,-1.43174e-13,50204.2,-55.6149], Tmin=(964.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[C]1CC(C1)C1CC1(276)',
    structure = SMILES('[C]1CC(C1)C1CC1'),
    E0 = (451.551,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89631,0.0219127,9.59076e-05,-1.37259e-07,5.30599e-11,54406.8,20.4973], Tmin=(100,'K'), Tmax=(964.159,'K')), NASAPolynomial(coeffs=[15.1292,0.0284212,-9.75405e-06,1.85864e-09,-1.40836e-13,49000.8,-57.6627], Tmin=(964.159,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: FFCM1(-)"""),
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
    E0 = (590.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (645.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (448.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (450.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (488.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (490.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (488.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (433.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (473.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (473.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (482.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (471.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (537.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (524.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (552.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (542.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction261',
    reactants = ['cC3H5(120)', '[CH]1C=CC1(258)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(6.5e+14,'cm^3/(mol*s)','*|/',2), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [C_sec_rad;C_rad/H/NonDeC] for rate rule [C_rad/H/CdCs;C_rad/H/NonDeC]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction268',
    reactants = ['CH2(S)(28)', 'C=CC1C=CC1(265)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.56622e+08,'m^3/(mol*s)'), n=-0.441667, Ea=(7.08269,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;mb_db]
Euclidian distance = 0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction269',
    reactants = ['[CH]1C[C](C1)C1CC1(266)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.0206e+11,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction270',
    reactants = ['[CH]1C[CH]C1C1CC1(267)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.0206e+11,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction271',
    reactants = ['[CH]1CC(C1)[C]1CC1(268)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction272',
    reactants = ['[CH]1[CH]C(C1)C1CC1(157)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction273',
    reactants = ['[CH]1CC[C]1C1CC1(269)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction274',
    reactants = ['[CH]1CCC1[C]1CC1(270)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction275',
    reactants = ['[CH]1CC(C1)C1[CH]C1(271)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.104e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction276',
    reactants = ['[CH]1CC1C1[CH]CC1(272)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction277',
    reactants = ['[CH2]C([CH2])C1C=CC1(273)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction278',
    reactants = ['[CH2]C[CH]C1C=CC1(78)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction279',
    reactants = ['[CH]=CC([CH2])C1CC1(274)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction280',
    reactants = ['[CH]=CC[CH]C1CC1(133)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction281',
    reactants = ['[C]1CCC1C1CC1(275)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.83662e+12,'s^-1'), n=0.345439, Ea=(91.3359,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction282',
    reactants = ['[C]1CC(C1)C1CC1(276)'],
    products = ['C1=CC(C1)C1CC1(146)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.67324e+12,'s^-1'), n=0.345439, Ea=(91.3359,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

network(
    label = '41',
    isomers = [
        'C1=CC(C1)C1CC1(146)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '41',
    Tmin = (300,'K'),
    Tmax = (3000,'K'),
    Tcount = 8,
    Tlist = ([302.617,324.619,374.997,470.374,649.057,1000.02,1706.11,2761.25],'K'),
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

