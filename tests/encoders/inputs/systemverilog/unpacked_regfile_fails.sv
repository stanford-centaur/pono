// The falsifying counterpart to unpacked_regfile.sv, and the guard
// against a vacuous encoding: if the array were over-constrained (or
// the write never modelled), a read-after-write invariant could
// "prove" for the wrong reason.  Here the assertion reads the address
// *next* to the one written, which a real memory does not preserve, so
// a correct encoding must find a counterexample.
module unpacked_regfile_fails (
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

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> mem[$past(waddr) + 4'd1] == $past(wdata));

endmodule
