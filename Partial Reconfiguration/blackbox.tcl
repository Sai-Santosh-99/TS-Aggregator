update_design -cell design_1_i/MAB1/inst/u1 -black_box
update_design -cell design_1_i/MAB2/inst/u1 -black_box
update_design -cell design_1_i/MAB3/inst/u1 -black_box

lock_design -level routing

write_checkpoint -force Checkpoint/static_route_design.dcp
close_project