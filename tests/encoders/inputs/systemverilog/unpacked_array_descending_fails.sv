// The reversed reading of unpacked_array_descending.sv's pattern,
// which is what filling by position rather than declared index
// produces. Must be refuted.
module unpacked_array_descending_fails (
    input logic clk,
    input logic rst
);

  logic [7:0] m[3:0];

  always_ff @(posedge clk) m <= '{8'd10, 8'd20, 8'd30, 8'd40};

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (m[3] == 8'd40 && m[0] == 8'd10));

endmodule
