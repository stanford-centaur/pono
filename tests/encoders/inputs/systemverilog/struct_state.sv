// A packed struct used as a *state* variable (not just a wire), with
// cross-field dependence: `flag` latches whether `cnt` was 3 on the
// previous cycle.  Exercises struct-typed always_ff registers plus
// per-field nonblocking assignment (`p.cnt <= ...`, `p.flag <= ...`)
// reading the *pre-update* value of a sibling field.
//
// cnt(k) = k-1 for k >= 1 (plain incrementer).  flag(k+1) = (cnt(k)
// == 3), so flag first becomes 1 at cycle 5, at which point cnt(5) is
// 4, not 5 -- the (deliberately wrong) property below expects 5 and
// is falsified exactly there.
typedef struct packed {
  logic [3:0] cnt;
  logic       flag;
} pair_t;

module struct_state (input logic clk, input logic rst);

  pair_t p;

  always_ff @(posedge clk) begin
    if (rst) begin
      p.cnt  <= 4'd0;
      p.flag <= 1'b0;
    end else begin
      p.flag <= (p.cnt == 4'd3);
      p.cnt  <= p.cnt + 4'd1;
    end
  end

  assert property (@(posedge clk) p.flag |-> p.cnt == 4'd5);

endmodule
