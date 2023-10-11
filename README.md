# Maximum Biplex Search 

> Accelerating maximum biplex search over large bipartite graphs



## Algorithm
- MBS, BK branching search algorithm.
- CSS, Symmetric BK branching search algorithm.
- CMBS-MBS, CMBS framework + MBS algorithm.
- CMBS-CSS, CMBS framework + CSS algorithm.
- CMBS$^+$-MBS, CMBS framework with bidirectional search strategy + MBS algorithm.
- CMBS$^+$-CSS, CMBS framework with bidirectional search strategy + CSS algorithm(the best algorithm in out paper).
- HFS, heuristic algorithm.
- FPS, FPGA-based parallel search algorithm for HFS.
- PBIE-FastBB, the state-of-art algorithm proposed in SIGMOD 23.

## Folder/File description
- MBS, BK branching search algorithm.
- CSS, Symmetric BK branching search algorithm.
- CMBS_MBS, source code of CMBS-MBS algorithm.
- CMBS_CSS, source code of CMBS-CSS algorithm.
- CMBS_PLUS_MBS, source code of CMBS$^+$-MBS algorithm.
- CMBS_PLUS_CSS, source code of CMBS$^+$-CSS algorithm.
- PBIE-FastBB, the state-of-art algorithm proposed in SIGMOD 23.
- Example, executables for all algorithms except the FPS algorithm (since the FPS algorithm requires an FPGA environment) as well as a sample diagram and a readme file.
