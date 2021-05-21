open_checkpoint Checkpoint/static_route_design.dcp

read_checkpoint -cell design_1_i/MAB1/inst/u1 Synth/TS_1/ts_synth.dcp
read_checkpoint -cell design_1_i/MAB2/inst/u1 Synth/TS_2/ts_synth.dcp
read_checkpoint -cell design_1_i/MAB3/inst/u1 Synth/TS_3/ts_synth.dcp

opt_design
place_design
route_design

write_checkpoint -force Implement/TS/top_ts_synth.dcp
close_project