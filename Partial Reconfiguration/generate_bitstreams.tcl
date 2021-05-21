open_checkpoint Implement/TS/top_ts_synth.dcp
write_bitstream -file Bitstreams/TS.bit 
close_project 

open_checkpoint Implement/UCB/top_ucb_synth.dcp
write_bitstream -file Bitstreams/UCB.bit 
close_project 

open_checkpoint Implement/BLANK/top_ucb_synth.dcp
write_bitstream -file Bitstreams/top_blank_synth.bit 
close_project 