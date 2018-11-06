# ips
imaging point sources

# to build the source code, in MATLAB.
```
cd femm
make
cd ../lbfgsc/src
make
```

# to run the code
```
load ./data/options/homogeneous4
ps = pointSolver(opt)
ps.reconstruction_lbfgs()
```