// Out-of-scope: a true SV "memory" -- an *unpacked* array of
// registers (`logic [7:0] mem [0:15];`), unlike the packed-array
// tests in test_systemverilog_types.cpp.  No smt::ARRAY sort is ever
// constructed anywhere in systemverilog_encoder.cpp (confirmed by
// inspection), and an unpacked array's type is not integral, so
// type_to_sort() throws.  This is a common real-RTL construct
// (register files, small memories), not just a verification-only
// feature -- worth flagging as a gap to close, not a permanent
// non-goal.
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
  assert property (@(posedge clk) 1'b1);
endmodule
