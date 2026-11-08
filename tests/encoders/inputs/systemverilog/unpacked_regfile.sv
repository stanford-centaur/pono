// A true SV "memory" -- an *unpacked* array of registers
// (`logic [7:0] mem[0:15];`), unlike the packed-array tests in
// test_systemverilog_types.cpp, which are bit ranges of one wide
// vector. This becomes an smt::ARRAY state var, so the element write
// is a Store and the read a Select.
//
// The read-after-write invariant below is 1-inductive, so BMC can
// only exhaust its bound; UnpackedRegfileMemoryProvable proves it,
// and unpacked_regfile_fails.sv is the paired refutation that rules
// out a vacuously over-constrained encoding.
module unpacked_regfile (
    input logic clk,
    input logic rst,
    input logic [3:0] waddr,
    input logic [7:0] wdata
);
  logic [7:0] mem[0:15];
  always_ff @(posedge clk) begin
    if (rst) mem[0] <= 8'd0;
    else mem[waddr] <= wdata;
  end
  // Read-after-write: once unpacked arrays are supported, a write from
  // the prior cycle should be visible at its address this cycle.
  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> mem[$past(waddr)] == $past(wdata));
endmodule
