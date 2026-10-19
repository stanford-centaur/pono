// A variable whose only driver is an `initial` block. It was never
// classified as state, so it stayed a free input var and the
// initial assignment reached the core as a constraint over
// something that is not a state variable -- an error that said
// nothing about the initial block it came from.
//
// Classifying it is only half of it: a state var with no update is
// treated as an input and would drift, so one that nothing else
// drives is given an update that holds. A clocked driver still
// takes over where there is one, which `r` checks against `seen`.
//
// Elements and scalars both, and an element index that is itself a
// state variable -- the index only has to be something the design
// knows at time 0.
module initial_only_driver (
    input logic clk,
    input logic rst,
    input logic [7:0] v
);

  logic [7:0] m[0:3];
  logic [7:0] s;
  logic [7:0] r;
  logic seen;
  logic [1:0] i;

  initial begin
    m[1] = 8'd7;
    s = 8'd9;
    r = 8'd7;
  end

  always_ff @(posedge clk) begin
    r <= v;
    seen <= 1'b1;
    i <= i + 2'd1;
  end

  assert property (@(posedge clk)
                   m[1] == 8'd7 && s == 8'd9 && (seen || r == 8'd7));

endmodule
