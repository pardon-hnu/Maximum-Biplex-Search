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

## Folder/File description
- MBS: source code of MBS algorithm.
- CSS: source code of CSS algorithm.
- CMBS_MBS: source code of CMBS-MBS algorithm.
- CMBS_CSS: source code of CMBS-CSS algorithm.
- CMBS_PLUS_MBS: source code of CMBS$^+$-MBS algorithm.
- CMBS_PLUS_CSS: source code of CMBS$^+$-CSS algorithm.
- HFS_FPS(C++): source code of HFS algorithm which is also the C++ version source code of FPS algorithm
- Datasets: some example datasets

## Compile command
- MBS: 
  ```shell
  g++ -O3 biplex.cpp biplex1.cpp core.cpp mbs.cpp -o mbs
  ```
- CSS: source code of CSS algorithm.
  ```shell
  g++ -O3 css.cpp -o css
  ```
- CMBS_MBS: source code of CMBS-MBS algorithm.
  ```shell
  g++ -O3 biplex.cpp biplex1.cpp core.cpp CMBS_MBS.cpp -o cmbs_mbs
  ```
- CMBS_CSS: source code of CMBS-CSS algorithm.
  ```shell
  g++ -O3 CMBS_CSS.cpp -o cmbs_css
  ```
- CMBS_PLUS_MBS: source code of CMBS$^+$-MBS algorithm.
  ```shell
  g++ -O3 main.cpp -o cmbs_plus_mbs
  ```
- CMBS_PLUS_CSS: source code of CMBS$^+$-CSS algorithm.
  ```shell
  g++ -O3 main.cpp -o cmbs_plus_css
  ```
- HFS_FPS(C++): source code of HFS algorithm which is also the C++ version source code of FPS algorithm
  ```shell
  g++ -O3 greedy.cpp -o hfs
  ```

## run
### css&cmbs_plus_mbs&cmbs_plus_css
```shell
    ./execfilename -f PathToBenchmark -u 3 -v 3 -k 1
```
```shell
  OPTIONS:
      -h, --help                        Display this help menu
      -f[benchmark], --file=[benchmark] Path to benchmark
      -u[Lower Bound], --u=[Lower
      Bound]                            The lower bound of the size of k-biplex
      -v[Lower Bound], --v=[Lower
      Bound]                            The lower bound of the size of k-biplex
      -k[para k], --k=[para k]          The parameter k
```
### mbs&cmbs_mbs&cmbs_css&hfs
```shell
    ./execfilename -f PathToBenchmark -s 3 -k 1
```
```shell
  options:
    -f, --file    filename (string)
    -k, --k       parameter k (int [=1])
    -s, --s       minimum size (int [=7])
    -?, --help    print this message
```