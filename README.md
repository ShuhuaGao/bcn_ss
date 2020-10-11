# Set Invariance and Optimal Set Stabilization of Boolean Control Networks: A Graphical Approach
This repository contains the code for the paper: Shuhua Gao, Cheng Xiang, and Tong Heng Lee. *Set Invariance and Optimal Set Stabilization of Boolean Control Networks: A Graphical Approach* 

**Organization**

+ The core algorithm implementation is in the folder *src/algorithm*.  Specifically, the proposed algorithms is implemented in *proposed.py*, and the existing algebraic approach is implemented in *related_work.py*. Please refer to the comments therein for more details.
+ The network transition matrix (i.e., *L*) of the $\lambda$ switch network is stored in *src/networks/lambda_switch.txt*, and the network transition matrix of the *ara* operon network is presented in  *src/networks/ara_operon.assr*, and the one for the T-LGL network is given in *src/networks/T_LGL.assr*.
+ Programs to reproduce the examples in the paper are given in four files as follows,
  + *src/example_lambda_switch.py*:  Example 1, 3, and 4
  + *src/example_lambda_switch_GYQ.py*: results of Example 1, 3 with an existing algebraic approach for verification purpose
  + *src/example_ara_operon.py*: reproduce the results for the *ara* operon network using our approach in Section VII
  + *src/benchmark_ara_operon.py*: running time comparison with the ara operon network in Section VII
  + src/benchmark_TLGL.py:  running time comparison with the TLGL network in Section VII
  
+ Full-sized colored figures in the paper are provided in the folder *fig* in PDF format.

**Requirement**

Python 3.6 or higher.

Packages:

+ [networkx]( https://pypi.org/project/networkx/ ) required for our graphical approach
+ [numpy]( https://pypi.org/project/numpy/ ) required for the algebraic approach
+ [pygraphviz](https://pypi.org/project/pygraphviz/ ) (optional) required for visualization of graphs and trees as shown in our paper. This is optional: if you don't want visualization, just set the variable `to_visualize` to `False` in related program files. By default, the variable `to_visualize` is `False`.

**How to run**

+ Download or clone this repository to your local computer.

+ In the command line  (such as *cmd*, *PowerShell* on Windows or *terminal* on Ubuntu), go into the *src* folder

+ To run an example, say *example_ara_operon.py*, just type `python ./example_ara_operon.py`. (Of course, you can also use any IDEs like PyCharm or Visual Studio Code as you want.)

**How to cite this work**
```
@ARTICLE{9209076,  
	author={S. {Gao} and C. {Xiang} and T. H. {Lee}},  
	journal={IEEE Transactions on Control of Network Systems},   
	title={Set Invariance and Optimal Set Stabilization of Boolean Control Networks: A Graphical Approach},   
	year={2020},  volume={},  number={},  pages={1-1},}
```
