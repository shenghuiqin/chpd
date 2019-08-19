#Data sources
database(
    thermoLibraries =['BurkeH2O2','Klippenstein_Glarborg2016','thermo_DFT_CCSDTF12_BAC','CBS_QB3_1dHR','DFT_QCI_thermo','primaryThermoLibrary'], # 'FFCM1(-)','primaryThermoLibrary', 'BurkeH2O2','DFT_QCI_thermo','CBS_QB3_1dHR'
    reactionLibraries = [('2005_Senosiain_OH_C2H2',False),('Glarborg/C3', False)], # 
    seedMechanisms = ['BurkeH2O2inN2','Klippenstein_Glarborg2016',], #
    kineticsDepositories = ['training'], 
    kineticsFamilies = ['default'],
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    maximumRadicalElectrons = 2,
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=9,
    maximumOxygenAtoms=4,
    allowSingletO2=True,
)

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
    temperature=[(600,'K'),(1500,'K')],
    pressure=[(10.0,'bar'),(40.0,'bar')],
    nSims=8,
    initialMoleFractions={
        #"C7H10": 1,
        "O2":   0.10885, # phi=1 means 9.5 O2 per C7H10
        "N2":   0.881684808, # 8.1 times as much N2 as O2
        "C7H16":0.00946, #  
           
    },
    terminationTime = (10.0, 's'),
    terminationConversion={
                'O2': 0.5,
        },
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)


model(
    toleranceKeepInEdge=0, # No pruning to start
    toleranceMoveToCore=0.4,
    toleranceInterruptSimulation=1,
    maxNumObjsPerIter=3, 
    terminateAtMaxObjects=True,
    maxNumSpecies=200, # first stage is until core reaches 100 species
    filterReactions=True, # should speed up model generation
)

model(
    toleranceMoveToCore=0.25,
    toleranceInterruptSimulation=1e8,
    toleranceKeepInEdge=0.01, # Pruning enabled for stage 2
    maximumEdgeSpecies=200000,
    minCoreSizeForPrune=100,
    minSpeciesExistIterationsForPrune=2,
    filterReactions=True,
    )

options(
    units='si',
    # saveRestartPeriod=(1,'hour'),
    generateOutputHTML=False,
    generatePlots=True,
    saveSimulationProfiles=True,
    saveEdgeSpecies=False,
)


pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2000,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
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




