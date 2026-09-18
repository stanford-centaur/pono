// `casex`/`casez` wildcard matching (`?` don't-care bits in a case
// item), checked as a correctness invariant rather than one hand-
// picked value: whenever bits 3 and 0 of the (same-cycle-paired,
// registered) input are both set, the wildcard pattern 4'b1??1 must
// match under both casex and casez, regardless of bits 1 and 2.
// Genuinely holds for every `a` if wildcard matching is implemented;
// a plain-equality mis-implementation (treating `?` literally) would
// never match and falsify this immediately.
module casex_casez (input logic clk, input logic rst, input logic [3:0] a);

  logic [3:0] a_prev;
  logic       hit_x, hit_z;

  always_ff @(posedge clk) begin
    if (rst) begin
      a_prev <= 4'd0;
      hit_x  <= 1'b0;
      hit_z  <= 1'b0;
    end else begin
      a_prev <= a;
      casex (a)
        4'b1??1: hit_x <= 1'b1;
        default: hit_x <= 1'b0;
      endcase
      casez (a)
        4'b1??1: hit_z <= 1'b1;
        default: hit_z <= 1'b0;
      endcase
    end
  end

  assert property (@(posedge clk)
                    !(a_prev[3] && a_prev[0]) || (hit_x && hit_z));

endmodule
