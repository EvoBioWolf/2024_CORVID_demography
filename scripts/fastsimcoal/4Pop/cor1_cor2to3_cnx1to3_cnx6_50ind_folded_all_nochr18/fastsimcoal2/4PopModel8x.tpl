//Parameters for the coalescence simulation program : fastsimcoal.exe
4 samples to simulate :
//Population effective sizes (number of genes)
NPOP1
NPOP2
NPOP3
NPOP4
//Samples sizes and samples age 
30
30
30
10
//Growth rates	: negative growth implies population expansion
0
GROW1
GROW2
0
//Number of migration matrices : 0 implies no migration between demes
4
//Migration matrix 0
0 MIG01R 0 0 
MIG10R 0 MIG12R 0
0 MIG21R 0 MIG23R
0 0 MIG32R 0
//Migration matrix 1
0 0 0 0
0 0 MIG12R 0
0 MIG21R 0 MIG23R
0 0 MIG32R 0
//Migration matrix 2
0 0 0 0
0 0 0 0
0 0 0 MIG23R
0 0 MIG32R 0
//Migration matrix 3
0 0 0 0
0 0 0 0  
0 0 0 0
0 0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
9 historical event
TEGR1 1 1 0 1 0 0
TMRMG 0 0 0 1 0 1
TMRMG 1 1 0 1 0 1
TMRMG 2 2 0 1 GROW2 1
TMRMG 3 3 0 1 0 1
TEGR2 2 2 0 1 0 1
TDIV1 1 2 1 NANC1 0 2 absoluteResize
TDIV2 2 3 1 NANC2 0 3 absoluteResize
TDIV3 3 0 1 NANC3 0 3 absoluteResize
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 3.18e-9 OUTEXP
