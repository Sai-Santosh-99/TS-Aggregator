read_checkpoint -cell design_1_i/MAB1/inst/u1 Synth/UCB/ucb_synth.dcp
read_checkpoint -cell design_1_i/MAB2/inst/u1 Synth/UCB/ucb_synth.dcp
read_checkpoint -cell design_1_i/MAB3/inst/u1 Synth/UCB/ucb_synth.dcp

set_property HD.RECONFIGURABLE 1 [get_cells design_1_i/MAB1/inst/u1]
set_property HD.RECONFIGURABLE 1 [get_cells design_1_i/MAB2/inst/u1]
set_property HD.RECONFIGURABLE 1 [get_cells design_1_i/MAB3/inst/u1]

write_checkpoint -force Checkpoint/top_link_ucb.dcp