// A recursive named property. slang pre-expands a named assertion's
// body, and a recursive one has no point at which the expansion
// stops, so there is nothing to recurse into -- a real
// implementation needs a fixpoint in the tableau rather than
// expansion.
//
// `p_rec` says "b holds until a does". `a` fires at cnt==3 and `b`
// holds until then, so intended: provable.
module recursive_named_property (input logic clk);

  logic [3:0] cnt;
  initial cnt = 4'd0;
  always_ff @(posedge clk) if (cnt != 4'd15) cnt <= cnt + 4'd1;

  logic a, b;
  assign a = (cnt == 4'd3);
  assign b = (cnt < 4'd3);

  property p_rec;
    a or (b and nexttime p_rec);
  endproperty

  assert property (@(posedge clk) p_rec);

endmodule
