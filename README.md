# TS-Aggregator

Source Code for Multi-armed Bandit Algorithms on System-on-Chip: Go Frequentist or Bayesian?

- To replicate the various MAB algorithms, synthesize the .cpp files in their respective folders inside the 'Hardware' folder using Vivado HLS with the following parameters:
    1) Board: ZC706
    2) Clock period: 20ns

- To replicate the Verilog IPs, synthesize the .v file in the inside the 'Verilog Hardware' folder using Vivado with the following parameters:
    1) Board: ZC706

- To replicate Dynamic Partial Reconfiguration, run the provided TCL files in Vivado, while following 'Tutorial.pdf' given in the 'Tutorial' folder.
