PartialDecomposition Program(R) Version 5.0 2018/4/26

This is an updated version on Vertion 4. In this version, all demands are combined together to perform the LP, so it can only run up to the parameter of (4, 8), which takes about 9 hours on my laptop (CPU: i7-6700HQ, Memory: 16GB). 

When finished, go to folder and find the plot function. The figure will be generated for the paper.

### Usage
---------------------

Function "func_partialDecomposition(N, K, t)" takes three arguments as input,

* N: number of files
* K: number of users
* t: integer from [0:K]

Function outputs two figures

* The 1st figure is the (M,R)-tradeoff line of Tian-Chen scheme and Yu scheme.
* The 2nd figure is a zoomed-in version of the 1st figure, focucing on the region where K = t.
* The traddoff line for each demand 'd' is plotted. For all file requested demands 'd', the corner point is labeled with decomposision pattern. e.g. (N,K,t) = (3,4,2), for demand (2, 1, 1), there is one corner point (4/3, 5/6), with label

Note: files such as 'A', 'B', 'C', etc are denoted as 1, 2, 3, etc... to not to be limited to only 26 files in the caching system.

<center><img src='https://github.com/kzhang14/A-New-Coding-Scheme-Uncoded-Caching-to-Coded-Caching-/blob/master/resources/plot34-page-001.jpg'  width="400" height="300"></center>

A new corner point (inner bound) of the caching system (3,4). In fact, the new point (4/3, 5/6) also locates on the outer bound, meaning it is an optimal point.

[1 3]

[1 2]

[1], [2 3]

meaning the transmission types: A + C and A + B are kept as they are, and only transmission type A + B + C will be decomposed into A, B+C.

When a corner point is get from more than one decomposition patterns, e.g. (N,K,t) = (3,5,2), for demand (2, 2, 1), there is one corner point with label

[2 3]

[1 2]

[1 3]

[1 2]

[1],[2],[3]

'------'

[2 3]

[1], [2]

[1 3]

[1], [2]

[1],[2],[3]  

meaning two decomposition patterns should be mixed to get this corner point. '-----' are used to separate each decomposition pattern.

The red star '*' indicates an achievable corner point by all demands.

Program can be reached at:
Email:kaizhang@tamu.edu
