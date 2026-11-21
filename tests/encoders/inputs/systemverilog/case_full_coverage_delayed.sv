// The guard against the register reading of case_full_coverage.sv
// returning. `y` is combinational, so it tracks this cycle's inputs
// and this claim about last cycle's is false -- and has to be
// refuted. Inferring a latch here would make it a register, and a
// register satisfies it exactly, so the property would be *proved*:
// a false proof of something the design does not do, which is worse
// than the false alarm the same mistake causes in the twin.
module case_full_coverage_delayed (
    input logic clk,
    input logic s,
    input logic [3:0] p,
    input logic [3:0] q
);

  logic [3:0] y;
  logic started;

  initial started = 1'b0;

  always_ff @(posedge clk) started <= 1'b1;

  always_comb begin
    case (s)
      1'b0: y = p;
      1'b1: y = q;
    endcase
  end

  assert property (@(posedge clk)
                   started |-> y == ($past(s) ? $past(q) : $past(p)));

endmodule
