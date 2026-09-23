// An abort over a temporal operand. The true-padded truncation at the
// abort forgives the operand's eventualities from that point on while
// still holding it to its safety part before, which asks for a split
// of the property into those two parts.
//
// `abrt` is tied high, so every attempt is waived and the eventuality
// beneath it never has to be met.
module property_accept_on_temporal (
    input logic clk,
    input logic a
);

  logic abrt;
  assign abrt = 1'b1;

  assert property (@(posedge clk) accept_on (abrt) (s_eventually a));

endmodule
