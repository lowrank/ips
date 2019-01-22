# ips
imaging point sources

# Build the source code
```
cd femm
make
cd ../lbfgsc/src
make
```

# Run the code
```
load ./data/options/homogeneous4
ps = pointSolver(opt)
ps.reconstruction_lbfgs()
```