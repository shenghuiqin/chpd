#Data sources
database(
    thermoLibraries =['primaryThermoLibrary', 'BurkeH2O2','DFT_QCI_thermo'], #['BurkeH2O2','FFCM1(-)','primaryThermoLibrary','DFT_QCI_thermo','CBS_QB3_1dHR'],
    reactionLibraries = [('FFCM1(-)',False),('2005_Senosiain_OH_C2H2',False)],
    seedMechanisms = ['BurkeH2O2inN2'],
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


# Reaction systems
simpleReactor(
    temperature=[(1000,'K')],
    pressure=[(40.0,'bar')],
    initialMoleFractions={
        "C7H10": 1,
        "O2":   19, # phi=1 means 9.5 O2 per C7H10
        "N2":    154, # 8.1 times as much N2 as O2
        "C7H16":1,
       
    },
    #terminationRateRatio=0.01,
    terminationTime=(0.1e0,'s'),
    #terminationConversion={
     #   'C7H16': 0.9,
      #  },
    
    #sensitivity=['C7H10','O2'],
    #sensitivityThreshold=0.001,
)

simulator(
    atol=1e-16,
    rtol=1e-8,
    sens_atol=1e-6,
    sens_rtol=1e-4,
)

model(
        toleranceKeepInEdge=0.0,
        toleranceMoveToCore=0.5,
        toleranceInterruptSimulation=0.5,
        maximumEdgeSpecies=100000,
        maxNumSpecies=100
)
model(
        toleranceKeepInEdge=0.0,
        toleranceMoveToCore=0.4,
        toleranceInterruptSimulation=0.4,
        maximumEdgeSpecies=100000,
)




options(
    units='si',
    # saveRestartPeriod=(1,'hour'),
    generateOutputHTML=False,
    generatePlots=False,
    saveSimulationProfiles=True,
    saveEdgeSpecies=False,
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


