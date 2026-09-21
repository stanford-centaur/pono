// sequence_repetition.sv's expansion shifted one cycle late, which
// is the error convolving offset vectors invites. Must be refuted.
module sequence_repetition_fails (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b
);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2)
                                   && !$past(rst, 3) && !$past(rst, 4))
                   |-> (((a ##1 b)[*2])
                        |-> $past(b && $past(a) && $past(b, 2)
                                  && $past(a, 3))));

endmodule
