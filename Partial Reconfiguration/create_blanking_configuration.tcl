open_checkpoint Checkpoint/static_route_design.dcp

update_design -buffer_ports -cell design_1_i/MAB1/inst/u1
update_design -buffer_ports -cell design_1_i/MAB2/inst/u1
update_design -buffer_ports -cell design_1_i/MAB3/inst/u1
u
write_checkpoint -force Checkpoint/top_link_blank.dcp

place_design
route_design

write_checkpoint -force Implement/BLANK/top_blank_synth.dcp
close_project