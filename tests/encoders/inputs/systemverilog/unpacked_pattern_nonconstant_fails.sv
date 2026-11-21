// The reversed reading of unpacked_pattern_nonconstant.sv's
// ascending array: filling by position rather than declared index
// would put v at the far end. Must be refuted.
module unpacked_pattern_nonconstant_fails (
    input logic clk,
    input logic rst,
    input logic [7:0] v,
    input logic [7:0] w
);

  logic [7:0] m[0:3];

  always_ff @(posedge clk) m <= '{v, w, w, w};

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> m[3] == $past(v));

endmodule
