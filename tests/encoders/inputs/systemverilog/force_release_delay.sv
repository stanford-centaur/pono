// A procedural continuous assignment whose extent is given by
// delays: applied one time unit in, released one time unit after
// that. There is no number of cycles `#1` stands for, which is the
// same reason `#delay` has no meaning anywhere here.
//
// The delay is what puts this out of scope, not the force or the
// release -- both are on delays here so that neither reads as the
// cause on its own. The clock-scoped spellings are modelled.
module force_release_delay (
    input logic clk
);

  logic [3:0] x;
  always_ff @(posedge clk) x <= x + 4'd1;

  initial begin
    #1 force x = 4'd5;
    #1 release x;
  end

  assert property (@(posedge clk) 1'b1);

endmodule
