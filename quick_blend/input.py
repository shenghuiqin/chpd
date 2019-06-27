#Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],  # 'primaryThermoLibrary' 'BurkeH2O2','DFT_QCI_thermo'],#['BurkeH2O2','FFCM1(-)','primaryThermoLibrary','DFT_QCI_thermo','CBS_QB3_1dHR'],
    reactionLibraries = [],  #[('FFCM1(-)',False),('2005_Senosiain_OH_C2H2',False)],
    seedMechanisms = [],    #['BurkeH2O2inN2'],
    kineticsDepositories = ['training'], 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    maximumRadicalElectrons = 2,
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=10,
    allowSingletO2 = False,
)
#this is a test#
#test local vs#
# List of species
species(
    label='C7H10',
    reactive=True,
    structure=SMILES("C1C=CCCCC=1"),
)
species(
    label='C7H16',
    reactive=True,
    structure=SMILES("CCCCCCC"),
)
species(
    label='O2',
    reactive=True,
        structure=adjacencyList(
        """
        1 O u1 p2 {2,S}
        2 O u1 p2 {1,S}
        """),
)
species(
    label='N2',
    reactive=False,
    structure=SMILES("N#N")
)

#species choosen from "file:///C:/Users/qinsh/Desktop/chpd/quick_chpd/output.html", based on their concentration value.

"""
species(
    label='C7H10(146)',
    reactive=True,
    structure=SMILES("C=CC=CCC=C")
)
species(
    label='C7H10(153)',
    reactive=True,
    structure=SMILES("C=CC=CC1CC1")
)
species(
    label='C7H9(5)',
    reactive=True,
    structure=SMILES("[CH]1C=CCCC=C1")
)
species(
    label='C5H7(210)',
    reactive=True,
    structure=SMILES("[CH2]C=CC=C")
)
species(
    label='C2H4(115)',
    reactive=True,
    structure=SMILES("C=C")
)
species(
    label='C7H11(22)',
    reactive=True,
    structure=SMILES("[CH]1C=CCCCC1")
)
species(
    label='C7H10(21)',
    reactive=True,
    structure=SMILES("C1=CC2CCCC12")
)
species(
    label='C7H10(65)',
    reactive=True,
    structure=SMILES("C1=CCCC=CC1")
)
species(
    label='C7H10(18)',
    reactive=True,
    structure=SMILES("[CH2]C=CC=CC[CH2]")
)
species(
    label='C7H10(238)',
    reactive=True,
    structure=SMILES("C=CCC1C=CC1")
)
species(
    label='S(66)',
    reactive=True,
    structure=SMILES("[O]OC1C=CCCCC1")
)
species(
    label='C7H11(46)',
    reactive=True,
    structure=SMILES("[CH2]CCC=CC=C")
)
species(
    label='C7H10(187)',
    reactive=True,
    structure=SMILES("C1=CC(C1)C1CC1")
)
species(
    label='Cycloheptene', # (19)
    reactive=True,
    structure=adjacencyList(
        ""
        1 C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
        2 C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
        3 C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
        4 C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
        5 C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
        6 C u0 p0 c0 {4,S} {7,D} {18,S}
        7 C u0 p0 c0 {5,S} {6,D} {19,S}
        8 H u0 p0 c0 {1,S}
        9 H u0 p0 c0 {1,S}
        10 H u0 p0 c0 {2,S}
        11 H u0 p0 c0 {2,S}
        12 H u0 p0 c0 {3,S}
        13 H u0 p0 c0 {3,S}
        14 H u0 p0 c0 {4,S}
        15 H u0 p0 c0 {4,S}
        16 H u0 p0 c0 {5,S}
        17 H u0 p0 c0 {5,S}
        18 H u0 p0 c0 {6,S}
        19 H u0 p0 c0 {7,S}
        ""),
)
species(
    label='C7H8(58)',
    reactive=True,
    structure=SMILES("C1=CC=CCC=C1")
)
species(
    label='C5H8(336)',
    reactive=True,
    structure=SMILES("C=CC=CC")
)
species(
    label='S(388)',
    reactive=True,
    structure=SMILES("C=CC=CCCC=CC=C")
)
species(
    label='C7H9O2(68)',
    reactive=True,
    structure=SMILES("[O]OC1C=CC=CCC1")
)
species(
    label='s(399)',
    reactive=True,
    structure=SMILES("C=CC=CCC(C=C)C=C")
)


species(
    label='OH',
    reactive=True,
    structure=SMILES("[OH]")
)
species(
    label='H2O',
    reactive=True,
    structure=SMILES("O")
)

species(
    label='Ar',
    reactive=False,
    structure=SMILES("[Ar]")
)

species(
    label='Ne',
    reactive=False,
    structure=SMILES("[Ne]")
)
species(
    label='He',
    reactive=False,
    structure=SMILES("[He]")
)
"""



# Reaction systems
simpleReactor(
    temperature=[(700,'K'),(1500,'K')],
    pressure=[(10.0,'bar'),(40.0,'bar')],
    nSims=8,
    initialMoleFractions={
        "C7H10": 1,
#	    "C7H16": 0.00463,
        "O2":   19, # phi=1 means 9.5 O2 per C7H10
        "N2":    154, # 8.1 times as much N2 as O2
        "C7H16":1,
      
        
        
        
    },
    terminationRateRatio=0.01,
    #  terminationTime=(1e0,'s'),
    sensitivity=['C7H10','O2'],
    sensitivityThreshold=0.001,
)

simulator(
    atol=1e-16,
    rtol=1e-8,
    sens_atol=1e-6,
    sens_rtol=1e-4,
)

model(
    #toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.05,
    toleranceInterruptSimulation=0.05,
    maximumEdgeSpecies=100000,
    filterReactions=True,
    maxNumObjsPerIter=2,
)

uncertainty(
    localAnalysis=False,
    globalAnalysis=False,
    uncorrelated=True,
    correlated=True,
    localNumber=10,
    globalNumber=5,
    terminationTime=None,
    pceRunTime=1800,
    logx=True
)

options(
    units='si',
    # saveRestartPeriod=(1,'hour'),
    generateOutputHTML=True,
    generatePlots=True,
    saveSimulationProfiles=True,
    saveEdgeSpecies=True,
)

#quantumMechanics(
#    software='mopac',
#    method='pm3',
    # fileStore='QMfiles', # relative to where you run it from. Defaults to inside the output folder if not defined.
#    scratchDirectory = None, # not currently used
#    onlyCyclics = True,
#    maxRadicalNumber = 0,
#    )

#pressureDependence(
#    method='modified strong collision',
#    maximumGrainSize=(0.5,'kcal/mol'),
#    minimumNumberOfGrains=250,
#    temperatures=(300,2000,'K',8),
#    pressures=(0.01,100,'bar',5),
#    interpolation=('Chebyshev', 6, 4),
#)


