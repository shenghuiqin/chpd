#Data sources
database(
    thermoLibraries = ['primaryThermoLibrary', 'BurkeH2O2','DFT_QCI_thermo'],#['BurkeH2O2','FFCM1(-)','primaryThermoLibrary','DFT_QCI_thermo','CBS_QB3_1dHR'],
    reactionLibraries =[]  #[('FFCM1(-)',False),('2005_Senosiain_OH_C2H2',False)],
    seedMechanisms = []    #['BurkeH2O2inN2'],
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
    label='CO',
    reactive=True,
    structure=SMILES("[C-]#[O+]"),
)
species(
    label='CO2',
    reactive=True,
    structure=SMILES("C(=O)=O"),
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
    structure=SMILES("N#N"),
)   
)
species(
    label = "CH4",
    reactive=True,
    structure = SMILES("C"))

species(
    label = "C2H6",
    reactive=True,
    structure = SMILES("CC"))

species(
    label = "C3H8",
    reactive=True,
    structure = SMILES("CCC"))

species(
    label = "C4H10",
    structure = SMILES("CCCC"))

species(
    label = "C5H12",
    reactive=True,
    structure = SMILES("CCCCC"))

species(
    label = "C6H14",
    reactive=True,
    structure = SMILES("CCCCCC"))

species(
    label = "C2H4",
    reactive=True,
    structure = SMILES("C=C"))

species(
    label = "C3H6",
    reactive=True,
    structure = SMILES("CC=C"))

species(
    label = "C4H8",
    reactive=True,
    structure = SMILES("CCC=C"))

species(
    label = "C5H10",
    reactive=True,
    structure = SMILES("CCCC=C"))

species(
    label = "C6H12",
    reactive=True,
    structure = SMILES("CCCCC=C"))
species(
    label = "METHYL",
    reactive=True,
    structure = SMILES("[CH3]"))

species(
    label = "ETHYL",
    reactive=True,
    structure = SMILES("C[CH2]"))

species(
    label = "PROPYL",
    reactive=True,
    structure = SMILES("CC[CH2]"))

species(
    label = "BUTYL",
    reactive=True,
    structure = SMILES("CCC[CH2]"))

species(
    label = "PENTYL",
    reactive=True,
    structure = SMILES("CCCC[CH2]"))

species(
    label = "HEXYL",
    reactive=True,
    structure = SMILES("CCCCC[CH2]"))

species(
    label = "BENZYL",
    reactive=True,
    structure = SMILES("[CH2]c1ccccc1"))

#species(
    label='C7H16',
    reactive=True,
    structure=SMILES("CCCCCCC"),
)



# Reaction systems
simpleReactor(
    temperature=[(700,'K'),(1500,'K')],
    pressure=[(10.0,'bar'),(40.0,'bar')],
    nSims=8,
    initialMoleFractions={
        "C7H10": 0.00463,
#	    "C7H16": 0.00463,
        "O2":    0.10880, # phi=1 means 9.5 O2 per C7H10
        "N2":     0.88035, # 8.1 times as much N2 as O2
        
    },
    terminationRateRatio=0.01,
    #  terminationTime=(1e0,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
    sens_atol=1e-6,
    sens_rtol=1e-4,
)

model(
    #toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
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


