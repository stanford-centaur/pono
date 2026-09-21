// The twin of pattern_case_tagged.sv, claiming the `B` arm always
// wins. That only holds if the tag is never consulted and the payload
// bits alone decide the match -- which is exactly the wrong answer a
// tagged pattern exists to prevent, so this must be refuted.
typedef union tagged packed {
  logic [3:0] A;
  logic [7:0] B;
} pattern_case_tagged_ignored_u;

module pattern_case_tagged_ignored (
    input logic clk,
    input logic [8:0] raw
);

  pattern_case_tagged_ignored_u u;
  logic [7:0] y;

  always_comb begin
    u = raw;
    case (u) matches
      tagged A .n: y = {4'd0, n};
      tagged B .n: y = n;
      default: y = 8'd0;
    endcase
  end

  assert property (@(posedge clk) y == raw[7:0]);

endmodule
