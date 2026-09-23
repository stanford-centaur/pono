// A `force` whose extent is given by a delay: held for one time
// unit, then released. There is no time between clock edges in a
// cycle-based model, so there is no number of cycles `#1` stands
// for -- the same reason `#delay` has no meaning anywhere here.
//
// The clock-scoped spellings are a different matter and are modelled
// elsewhere; this is the form that cannot be.
module force_release_delay (
    input logic clk
);

  logic [3:0] x;
  always_ff @(posedge clk) x <= x + 4'd1;

  initial begin
    force x = 4'd5;
    #1 release x;
  end

  assert property (@(posedge clk) 1'b1);

endmodule
