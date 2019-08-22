#Data sources
database(
    thermoLibraries =['primaryThermoLibrary'], # 'FFCM1(-)','primaryThermoLibrary', 'BurkeH2O2','DFT_QCI_thermo','CBS_QB3_1dHR'
    reactionLibraries = [], # 
    seedMechanisms = [], #
    kineticsDepositories = ['training'], 
    kineticsFamilies = ['Disproportionation','R_Recombination','intra_H_migration','H_abstraction'],#,'Cyclic_Ether_Formation','HO2_Elimination_from_PeroxyRadical','Intra_Disproportionation','intra_H_migration','ketoenol'
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    maximumRadicalElectrons = 2,
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=9,
    maximumOxygenAtoms=4,
)

# List of species
species(
    label='N2(1)',
    reactive=True,
    structure=adjacencyList(
        """
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
        """),
)


species(
    label='Ar(2)',
    reactive=True,
    structure=adjacencyList(
        """
1 Ar u0 p4 c0
        """),
)


species(
    label='He(3)',
    reactive=True,
    structure=adjacencyList(
        """
1 He u0 p1 c0
        """),
)


species(
    label='Ne(4)',
    reactive=True,
    structure=adjacencyList(
        """
1 Ne u0 p4 c0
        """),
)


species(
    label='CHPD(1)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,D} {14,S}
5  C u0 p0 c0 {3,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {7,S} {16,S}
7  C u0 p0 c0 {4,D} {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='O2(2)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
        """),
)


species(
    label='HO2(3)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H9(5)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {6,D} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {7,S} {15,S}
6  C u0 p0 c0 {3,D} {7,S} {16,S}
7  C u1 p0 c0 {5,S} {6,S} {14,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C7H10(9)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,D} {15,S}
5  C u0 p0 c0 {1,S} {6,D} {14,S}
6  C u0 p0 c0 {3,S} {5,D} {16,S}
7  C u0 p0 c0 {3,S} {4,D} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='OO(8)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
        """),
)


species(
    label='H(12)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 H u1 p0 c0
        """),
)


species(
    label='C7H9O2(20)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {9,D} {16,S}
7  C u0 p0 c0 {5,S} {8,D} {15,S}
8  C u0 p0 c0 {7,D} {9,S} {18,S}
9  C u0 p0 c0 {6,D} {8,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C7H9O2(19)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {14,S}
6  C u0 p0 c0 {4,S} {9,D} {16,S}
7  C u0 p0 c0 {3,S} {8,D} {15,S}
8  C u0 p0 c0 {5,S} {7,D} {17,S}
9  C u0 p0 c0 {5,S} {6,D} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(36)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {19,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {9,D} {16,S}
7  C u0 p0 c0 {5,S} {8,D} {15,S}
8  C u0 p0 c0 {7,D} {9,S} {18,S}
9  C u0 p0 c0 {6,D} {8,S} {17,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(35)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {19,S}
3  C u0 p0 c0 {1,S} {8,S} {9,S} {14,S}
4  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {9,D} {16,S}
7  C u0 p0 c0 {4,S} {8,D} {15,S}
8  C u0 p0 c0 {3,S} {7,D} {17,S}
9  C u0 p0 c0 {3,S} {6,D} {18,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C7H9O(53)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {13,S}
5  C u0 p0 c0 {3,S} {8,D} {15,S}
6  C u0 p0 c0 {2,S} {7,D} {14,S}
7  C u0 p0 c0 {4,S} {6,D} {16,S}
8  C u0 p0 c0 {4,S} {5,D} {17,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='OH(D)(44)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H10O(67)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {18,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {8,D} {15,S}
6  C u0 p0 c0 {4,S} {7,D} {14,S}
7  C u0 p0 c0 {6,D} {8,S} {17,S}
8  C u0 p0 c0 {5,D} {7,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H9O(69)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {17,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {4,D} {8,S} {14,S}
7  C u0 p0 c0 {5,D} {8,S} {16,S}
8  C u1 p0 c0 {6,S} {7,S} {15,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H9O3(82)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {15,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {17,S}
9  C u0 p0 c0 {5,S} {10,D} {16,S}
10 C u0 p0 c0 {6,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C7H9O3(83)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {13,S}
6  C u0 p0 c0 {4,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {5,S} {10,D} {16,S}
9  C u0 p0 c0 {7,D} {10,S} {17,S}
10 C u0 p0 c0 {8,D} {9,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C7H9O3(84)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {19,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {10,D} {16,S}
8  C u0 p0 c0 {6,S} {9,D} {15,S}
9  C u0 p0 c0 {8,D} {10,S} {18,S}
10 C u0 p0 c0 {7,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(94)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {4,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {10,D} {16,S}
8  C u0 p0 c0 {6,S} {9,D} {15,S}
9  C u0 p0 c0 {8,D} {10,S} {18,S}
10 C u0 p0 c0 {7,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(93)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {10,D} {16,S}
9  C u0 p0 c0 {7,D} {10,S} {17,S}
10 C u0 p0 c0 {8,D} {9,S} {18,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(92)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {8,S} {10,S} {15,S}
6  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {17,S}
9  C u0 p0 c0 {6,S} {10,D} {16,S}
10 C u0 p0 c0 {5,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(133)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {18,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {14,S}
6  C u0 p0 c0 {1,S} {3,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {16,S}
8  C u0 p0 c0 {4,S} {9,D} {15,S}
9  C u0 p0 c0 {5,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(98)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {6,S} {19,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {12,S}
5  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u0 p0 c0 {4,S} {9,D} {15,S}
8  C u0 p0 c0 {6,D} {9,S} {16,S}
9  C u0 p0 c0 {7,D} {8,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(97)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {6,S} {19,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {14,S}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {16,S}
8  C u0 p0 c0 {4,S} {9,D} {15,S}
9  C u0 p0 c0 {5,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(99)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {18,S}
2  O u0 p2 c0 {3,S} {19,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {9,D} {15,S}
7  C u0 p0 c0 {5,S} {8,D} {14,S}
8  C u0 p0 c0 {7,D} {9,S} {17,S}
9  C u0 p0 c0 {6,D} {8,S} {16,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C7H9O(45)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {8,D} {15,S}
6  C u0 p0 c0 {4,S} {7,D} {14,S}
7  C u0 p0 c0 {6,D} {8,S} {17,S}
8  C u0 p0 c0 {5,D} {7,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H8O(75)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {4,D} {8,S} {14,S}
7  C u0 p0 c0 {5,D} {8,S} {16,S}
8  C u1 p0 c0 {6,S} {7,S} {15,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H8O3(85)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  O u1 p2 c0 {7,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {15,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {17,S}
9  C u0 p0 c0 {5,S} {10,D} {16,S}
10 C u0 p0 c0 {6,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(101)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {13,S}
6  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {10,S} {16,S}
8  C u0 p0 c0 {2,D} {6,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {17,S}
10 C u0 p0 c0 {7,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(179)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {15,S}
6  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {2,D} {4,S} {5,S}
8  C u1 p0 c0 {5,S} {10,S} {17,S}
9  C u0 p0 c0 {6,S} {10,D} {16,S}
10 C u0 p0 c0 {8,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(180)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {5,S} {16,S}
8  C u0 p0 c0 {2,D} {6,S} {10,S}
9  C u0 p0 c0 {5,S} {10,D} {17,S}
10 C u0 p0 c0 {8,S} {9,D} {18,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(108)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {4,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {10,D} {16,S}
8  C u0 p0 c0 {6,S} {9,D} {15,S}
9  C u0 p0 c0 {8,D} {10,S} {18,S}
10 C u0 p0 c0 {7,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(183)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {19,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {7,S} {9,S} {15,S}
5  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {5,S} {16,S}
8  C u0 p0 c0 {3,D} {6,S} {10,S}
9  C u0 p0 c0 {4,S} {10,D} {17,S}
10 C u0 p0 c0 {8,S} {9,D} {18,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(126)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {19,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u1 p0 c0 {4,S} {9,S} {16,S}
8  C u0 p0 c0 {3,D} {6,S} {10,S}
9  C u0 p0 c0 {7,S} {10,D} {17,S}
10 C u0 p0 c0 {8,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(134)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {19,S}
3  O u1 p2 c0 {9,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {10,S} {14,S} {15,S}
7  C u0 p0 c0 {4,S} {9,D} {16,S}
8  C u0 p0 c0 {4,S} {10,D} {17,S}
9  C u0 p0 c0 {3,S} {5,S} {7,D}
10 C u0 p0 c0 {6,S} {8,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(221)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,S} {20,S}
3  O u0 p2 c0 {1,S} {21,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {15,S}
8  C u0 p0 c0 {6,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {7,S} {11,D} {18,S}
10 C u0 p0 c0 {4,D} {8,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(222)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,S} {20,S}
3  O u0 p2 c0 {1,S} {21,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {9,S} {11,S} {17,S}
9  C u0 p0 c0 {4,D} {7,S} {8,S}
10 C u0 p0 c0 {5,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(234)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {19,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {2,S} {4,S} {9,S} {13,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
8  C u0 p0 c0 {3,D} {5,S} {7,S}
9  C u0 p0 c0 {6,S} {10,D} {17,S}
10 C u0 p0 c0 {7,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(226)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {19,S}
2  O u1 p2 c0 {4,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {10,D} {17,S}
9  C u0 p0 c0 {3,D} {7,S} {10,S}
10 C u0 p0 c0 {8,D} {9,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(225)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,S} {20,S}
3  O u0 p2 c0 {1,S} {21,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {8,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {7,S} {11,S} {16,S} {17,S}
9  C u0 p0 c0 {4,D} {6,S} {7,S}
10 C u0 p0 c0 {5,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(255)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {19,S}
2  O u1 p2 c0 {5,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
5  C u0 p0 c0 {2,S} {4,S} {9,S} {12,S}
6  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {10,S} {15,S} {16,S}
8  C u0 p0 c0 {3,D} {4,S} {6,S}
9  C u0 p0 c0 {5,S} {10,D} {17,S}
10 C u0 p0 c0 {7,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(245)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {9,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {2,D} {4,S} {5,S}
8  C u1 p0 c0 {4,S} {10,S} {16,S}
9  C u0 p0 c0 {3,S} {6,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(273)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {6,S} {20,S}
3  O u0 p2 c0 {10,D}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
6  C u0 p0 c0 {2,S} {9,S} {11,S} {17,S}
7  C u0 p0 c0 {8,S} {10,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {11,S} {15,S} {16,S}
9  C u1 p0 c0 {5,S} {6,S} {18,S}
10 C u0 p0 c0 {3,D} {5,S} {7,S}
11 C u0 p0 c0 {4,D} {6,S} {8,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(184)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {14,S}
5  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {2,D} {3,S} {4,S}
7  C u1 p0 c0 {4,S} {9,S} {16,S}
8  C u0 p0 c0 {5,S} {9,D} {15,S}
9  C u0 p0 c0 {7,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(191)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {8,S} {20,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {9,S} {11,S} {17,S}
9  C u0 p0 c0 {3,D} {7,S} {8,S}
10 C u0 p0 c0 {6,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(188)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
7  C u0 p0 c0 {8,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {7,S} {11,S} {16,S} {17,S}
9  C u0 p0 c0 {3,D} {5,S} {7,S}
10 C u0 p0 c0 {6,S} {11,D} {19,S}
11 C u0 p0 c0 {8,S} {10,D} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(151)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {17,S}
2  O u0 p2 c0 {5,S} {18,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {4,S} {8,D}
6  C u0 p0 c0 {1,S} {3,S} {7,D}
7  C u0 p0 c0 {6,D} {9,S} {14,S}
8  C u0 p0 c0 {5,D} {9,S} {16,S}
9  C u1 p0 c0 {7,S} {8,S} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(279)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {20,S}
2  O u0 p2 c0 {6,S} {21,S}
3  O u0 p2 c0 {10,D}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
6  C u0 p0 c0 {2,S} {7,S} {11,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
8  C u0 p0 c0 {9,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {8,S} {11,S} {18,S} {19,S}
10 C u0 p0 c0 {3,D} {5,S} {8,S}
11 C u0 p0 c0 {4,D} {6,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {1,S}
21 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(312)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u0 p2 c0 {1,S} {21,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u0 p0 c0 {5,S} {11,D} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(310)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u0 p0 c0 {5,S} {11,D} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {18,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(311)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u0 p2 c0 {1,S} {21,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {10,S} {11,S} {16,S}
8  C u0 p0 c0 {3,S} {6,S} {11,D}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 C u0 p0 c0 {7,S} {8,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(309)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {10,S} {11,S} {16,S}
8  C u0 p0 c0 {3,S} {6,S} {11,D}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 C u0 p0 c0 {7,S} {8,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(330)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,S} {18,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {9,S} {10,S} {15,S}
7  C u0 p0 c0 {2,S} {5,S} {10,D}
8  C u0 p0 c0 {1,S} {4,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {16,S}
10 C u0 p0 c0 {6,S} {7,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(314)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {4,S} {19,S}
3  O u0 p2 c0 {7,S} {20,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {10,D} {15,S}
9  C u0 p0 c0 {7,D} {10,S} {16,S}
10 C u0 p0 c0 {8,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(315)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {10,D} {15,S}
9  C u0 p0 c0 {7,D} {10,S} {16,S}
10 C u0 p0 c0 {8,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(313)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,S} {18,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {7,S} {20,S}
4  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {9,S} {10,S} {15,S}
7  C u0 p0 c0 {3,S} {5,S} {10,D}
8  C u0 p0 c0 {2,S} {4,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {16,S}
10 C u0 p0 c0 {6,S} {7,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(301)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {17,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {5,D} {9,S} {14,S}
8  C u0 p0 c0 {6,D} {9,S} {16,S}
9  C u1 p0 c0 {7,S} {8,S} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(362)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u0 p2 c0 {8,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
8  C u0 p0 c0 {3,D} {6,S} {7,S}
9  C u0 p0 c0 {2,S} {5,S} {11,D}
10 C u1 p0 c0 {7,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(363)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {7,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u1 p0 c0 {2,S} {5,S} {6,S}
9  C u0 p0 c0 {3,D} {7,S} {11,S}
10 C u0 p0 c0 {6,S} {11,D} {17,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(335)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u1 p2 c0 {10,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {10,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {9,S} {11,S} {16,S}
8  C u0 p0 c0 {2,S} {5,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {17,S}
10 C u0 p0 c0 {3,S} {6,S} {11,D}
11 C u0 p0 c0 {7,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(324)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u1 p0 c0 {5,S} {11,S} {16,S}
9  C u0 p0 c0 {3,D} {7,S} {10,S}
10 C u0 p0 c0 {9,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(317)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u1 p0 c0 {5,S} {10,S} {16,S}
9  C u0 p0 c0 {4,D} {7,S} {11,S}
10 C u0 p0 c0 {8,S} {11,D} {17,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(366)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
6  C u0 p0 c0 {7,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u1 p0 c0 {2,S} {5,S} {6,S}
9  C u0 p0 c0 {4,D} {7,S} {11,S}
10 C u0 p0 c0 {5,S} {11,D} {17,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(323)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u1 p2 c0 {6,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {7,S} {10,D}
9  C u0 p0 c0 {6,S} {11,D} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(331)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
6  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {11,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {6,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {17,S}
10 C u0 p0 c0 {5,S} {11,D} {18,S}
11 C u0 p0 c0 {4,S} {7,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(365)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
6  C u0 p0 c0 {7,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {4,D} {5,S} {7,S}
9  C u0 p0 c0 {2,S} {6,S} {11,D}
10 C u1 p0 c0 {5,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(355)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {18,S}
2  O u1 p2 c0 {4,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u1 p0 c0 {4,S} {10,S} {15,S}
8  C u0 p0 c0 {3,D} {6,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {16,S}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(316)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {7,S} {10,D}
9  C u0 p0 c0 {5,S} {11,D} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {18,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(339)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {7,S} {18,S}
2  O u1 p2 c0 {6,S}
3  O u1 p2 c0 {9,S}
4  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {8,S} {10,S} {15,S}
7  C u0 p0 c0 {1,S} {4,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {16,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u0 p0 c0 {6,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(404)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {10,S} {20,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
6  C u0 p0 c0 {7,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {5,S} {11,S} {17,S}
9  C u0 p0 c0 {4,D} {5,S} {7,S}
10 C u0 p0 c0 {2,S} {6,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(401)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {7,S} {20,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {7,S} {11,D} {17,S}
10 C u0 p0 c0 {4,D} {8,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(274)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {20,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u1 p2 c0 {6,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {3,S} {5,S} {10,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {9,S} {11,S} {14,S}
9  C u0 p0 c0 {4,D} {7,S} {8,S}
10 C u0 p0 c0 {6,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(388)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {8,S} {18,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {7,S} {9,S} {15,S}
7  C u0 p0 c0 {3,D} {5,S} {6,S}
8  C u0 p0 c0 {1,S} {4,S} {10,D}
9  C u1 p0 c0 {6,S} {10,S} {16,S}
10 C u0 p0 c0 {8,D} {9,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(427)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {9,S} {11,S} {16,S}
9  C u0 p0 c0 {4,D} {7,S} {8,S}
10 C u0 p0 c0 {5,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(426)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {10,S} {20,S}
3  O u1 p2 c0 {6,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {9,S} {17,S}
7  C u0 p0 c0 {8,S} {10,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {4,D} {6,S} {8,S}
10 C u0 p0 c0 {2,S} {7,S} {11,D}
11 C u0 p0 c0 {5,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(436)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {7,S} {18,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u1 p2 c0 {10,S}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {11,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {16,S}
8  C u0 p0 c0 {2,S} {5,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {17,S}
10 C u0 p0 c0 {3,S} {7,S} {11,D}
11 C u0 p0 c0 {4,S} {6,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)

# Reaction systems
simpleReactor(
    temperature=[(600,'K'),(1000,'K')],
    pressure=[(10.0,'bar'),(40.0,'bar')],
    nSims=6,
    initialMoleFractions={
        #"C7H10": 1,
        "O2(2)":   0.108803, # phi=1 means 9.5 O2 per C7H10
        "N2(1)":   0.36, # 8.1 times as much N2 as O2
        "CHPD(1)":0.00989, # Cycloheptadiene C7H10    
"Ar(2)": 0,
"He(3)": 0,
"Ne(4)": 0,
"HO2(3)": 0,
"C7H9(5)": 0,
"C7H10(9)": 0,
"OO(8)": 0,
"H(12)": 0,
"C7H9O2(20)": 0,
"C7H9O2(19)": 0,
"S(36)": 0,
"S(35)": 0,
"C7H9O(53)": 0,
"OH(D)(44)": 0,
"C7H10O(67)": 0,
"C7H9O(69)": 0,
"C7H9O3(82)": 0,
"C7H9O3(83)": 0,
"C7H9O3(84)": 0,
"S(94)": 0,
"S(93)": 0,
"S(92)": 0,
"S(133)": 0,
"S(98)": 0,
"S(97)": 0,
"S(99)": 0,
"C7H9O(45)": 0,
"C7H8O(75)": 0,
"C7H8O3(85)": 0,
"S(101)": 0,
"S(179)": 0,
"S(180)": 0,
"S(108)": 0,
"S(183)": 0,
"S(126)": 0,
"S(134)": 0,
"S(221)": 0,
"S(222)": 0,
"S(234)": 0,
"S(226)": 0,
"S(225)": 0,
"S(255)": 0,
"S(245)": 0,
"S(273)": 0,
"S(184)": 0,
"S(191)": 0,
"S(188)": 0,
"S(151)": 0,
"S(279)": 0,
"S(312)": 0,
"S(310)": 0,
"S(311)": 0,
"S(309)": 0,
"S(330)": 0,
"S(314)": 0,
"S(315)": 0,
"S(313)": 0,
"S(301)": 0,
"S(362)": 0,
"S(363)": 0,
"S(335)": 0,
"S(324)": 0,
"S(317)": 0,
"S(366)": 0,
"S(323)": 0,
"S(331)": 0,
"S(365)": 0,
"S(355)": 0,
"S(316)": 0,
"S(339)": 0,
"S(404)": 0,
"S(401)": 0,
"S(274)": 0,
"S(388)": 0,
"S(427)": 0,
"S(426)": 0,
"S(436)": 0,
        
    },
    terminationTime = (5.0, 's'),
    terminationRateRatio = 0.01,
    terminationConversion={
                'CHPD(1)': 0.99,
    },
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)


model(
    toleranceKeepInEdge=0, # No pruning to start
    toleranceMoveToCore=0.5,
    toleranceInterruptSimulation=1,
    maxNumObjsPerIter=3,      #
    terminateAtMaxObjects=True,
    maxNumSpecies=100, # first stage is until core reaches 100 species
    filterReactions=True, # should speed up model generation
    filterThreshold=2e8,
)

model(
    toleranceMoveToCore=0.8,
    toleranceInterruptSimulation=1e8,
    toleranceKeepInEdge=0.01, # Pruning enabled for stage 2
    maximumEdgeSpecies=100000,
    minCoreSizeForPrune=100,
    minSpeciesExistIterationsForPrune=2,
    filterReactions=True,
    filterThreshold=2e8,
    maxNumObjsPerIter=4,
    terminateAtMaxObjects=True,
    )

options(
    units='si',
    # saveRestartPeriod=(1,'hour'),
    generateOutputHTML=True,
    generatePlots=True,
    saveSimulationProfiles=True,
    saveEdgeSpecies=False,
    generateSeedEachIteration=False,
)



