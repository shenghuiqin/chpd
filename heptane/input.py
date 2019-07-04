#Data sources
database(
    thermoLibraries =[], # 'FFCM1(-)','primaryThermoLibrary', 'BurkeH2O2','DFT_QCI_thermo','CBS_QB3_1dHR'
    reactionLibraries = [], # ('FFCM1(-)',False),('2005_Senosiain_OH_C2H2',False)
    seedMechanisms = [], #'BurkeH2O2inN2'
    kineticsDepositories = ['training'], 
    kineticsFamilies = ['default','!surface', '!surface_development']
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    maximumRadicalElectrons = 2,
  # allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=8,
)
#this is a test#
#test local vs#
# List of species
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

# Reaction systems
simpleReactor(
    temperature=[(700,'K'),(1000,'K')],
    pressure=[(20.0,'bar'),(40.0,'bar')],
    nSims=8,
    initialMoleFractions={
        #"C7H10": 1,
        "O2":   11, # phi=1 means 9.5 O2 per C7H10
        "N2":    89.1, # 8.1 times as much N2 as O2
        "C7H16":1,  
    },
    terminationTime=(1e0,'s'), 
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)


model(
    toleranceKeepInEdge=0.0001,
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=0.5,
    maximumEdgeSpecies=1000
)

options(
    units='si',
    # saveRestartPeriod=(1,'hour'),
    generateOutputHTML=False,
    generatePlots=True,
    saveSimulationProfiles=True,
    saveEdgeSpecies=False,
)


"""
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
""" 

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


