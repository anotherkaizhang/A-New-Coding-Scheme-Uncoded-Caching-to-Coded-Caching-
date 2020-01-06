PartialDecomposition Program(R) Version 1.0 2017/12/18

Update 2018/04/26£º
This program is very inefficiency and might be wrong, for the reason that when finding new corner points, it samples along the M axis, and find the lowest R, this could potentially miss corner points, and the number of linear programing is quite large depending the number of sample points. 

Also, in the paper the LP is actually to combine all demands, but not for each demand, though there is no difference. 

General Usage Notes
---------------------

- Function "func_partialDecomposition(N, K, t)" takes three arguments as input,
  -- N: number of files
  -- K: number of users
  -- t: integer from [0:K]

- Function outputs two figures
  -- The 1st figure is the (M,R)-tradeoff line of Tian-Chen scheme and Yu scheme.
  -- The 2nd figure is a "zoom-in" version of the 1st figure, focucing on the region of the 't'-value.
  -- The traddoff line for each demand 'd' is plotted. For all file requested demands 'd', the corner point is labeled with decomposision pattern. e.g. (N,K,t) = (3,4,2), for demand (2, 1, 1), there is one corner point (4/3, 5/6), with label
	----
	[1 3]
	[1 2]
	[1], [2 3]
     meaning only for transmission type including all three files are decomposed, into A, B+C.

     When a corner point is get from more than one decomposition patterns, e.g. (N,K,t) = (3,5,2), for demand (2, 2, 1), there is one corner point with label
	-----
	[2 3]
	[1 2]
	[1 3]
	[1 2]
	[1],[2],[3]
	-----
	[2 3]
	[1], [2]
	[1 3]
	[1], [2]
	[1],[2],[3]  
     meaning two decomposition patterns should be mixed to get this corner point. '-----' are used to separate each decomposition pattern.
  -- The red star '*' indicates an achievable corner point by all demands.
============================================================================================================================================
Program can be reached at:
E-mail:kaizhang@tamu.edu