// `$readmemh` fills a memory from a file, one of several system
// tasks that take an argument by reference and write through it
// (`$deposit`, `$cast`, `$sformat`, `$swrite`, `$q_add`, the PLA
// loaders). Skipping such a call is not like skipping a `$display`:
// the target keeps whatever it held, and here that means a memory
// the design goes on to read is never filled.
module system_task_output_arg (
    input logic       clk,
    input logic [1:0] idx
);

  logic [7:0] rom[0:3];
  logic [7:0] q;

  initial $readmemh("rom.hex", rom);

  always_ff @(posedge clk) q <= rom[idx];

  assert property (@(posedge clk) q != 8'hff);

endmodule
