// weak() over a sequence shape that leading_condition() could not
// name the start of, so weak_seq_bool() threw instead of modelling
// it. Each of these is the same sequence as `x ##1 y` written
// through a different operator, so all of them must be refuted at
// exactly the depth weak_seq_fails.sv is -- a wrong leading
// condition moves that depth or loses the refutation entirely.
module weak_seq_within (
    input logic clk,
    input logic x,
    input logic y
);

  assert property (@(posedge clk) weak((x ##1 y) within (x ##1 y)));

endmodule
