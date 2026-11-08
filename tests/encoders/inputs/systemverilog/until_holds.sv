// The genuinely-*holding* side of an until-family property (the
// falsifying side is already covered by until_live.sv).  SVA has no
// direct "release" keyword -- proving `busy s_until done` holds
// requires negating it internally, which for an Until-family operator
// is exactly the Release gadget (ltl_make_R).  `cnt` deterministically
// saturates at 3 (no free inputs at all), latching `done` permanently
// true from cycle 4 onward, so the strong-until obligation is
// genuinely, unconditionally discharged -- BMC must exhaust the bound
// with UNKNOWN rather than manufacture a counterexample.
module until_holds (input logic clk, input logic rst);

  logic       busy;
  logic       done;
  logic [3:0] cnt;

  assign busy = 1'b1;
  assign done = (cnt == 4'd3);

  always_ff @(posedge clk) begin
    if (rst) cnt <= 4'd0;
    else if (cnt < 4'd3) cnt <= cnt + 4'd1;
  end

  assert property (@(posedge clk) busy s_until done);

endmodule
