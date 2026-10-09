// A `.v` variable pattern matches anything and names what it
// matched, so the arm below can read it. Two things follow that a
// plain `case` never has to deal with.
//
// Items are first-match-wins: the catch-all matches `x == 1` too, so
// without guarding each arm against the earlier ones both would run
// and the second would overwrite the first.
//
// And an irrefutable pattern makes the case exhaustive with no
// `default` in sight, so `y` is combinational here rather than the
// latch an incompletely assigned always_comb otherwise infers.
module pattern_case_first_match (
    input logic clk,
    input logic [3:0] x
);

  logic [3:0] y;

  always_comb begin
    case (x) matches
      4'd1: y = 4'd10;
      .v:   y = v;
    endcase
  end

  assert property (@(posedge clk)
                   ((x == 4'd1) |-> (y == 4'd10))
                   and ((x != 4'd1) |-> (y == x)));

endmodule
