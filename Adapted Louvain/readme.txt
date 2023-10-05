The code can be run like any other python file.
First install the requirements using pip install -r /path/to/requirements.txt
If python is installed on the machine, run the following command in a terminal: Python.exe [path_to_file]
We advise using an IDE like visual studio code, which does this automatically.
This script was built for python 3 or later. If previous versions are used, the clock time is measured instead of the cpu time.

The parameters for the runs can be edited at the top of the script:
Runs : Int : the number of runs per scenario.
Nodes : Int Array : The number of nodes to be used in the networks for each bach of runs.
Fractions : Int Array : The fraction of nodes that need to be of degree one.
AvgK: float : The average degree of nodes in the network.
Debug : Bool : prints additional details during run if set to true.

To reproduce the experiments, use the default parameters or the parameters as seen in the paper.