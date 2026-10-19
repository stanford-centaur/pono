// A structure pattern matches a packed struct field by field,
// each field's own pattern applying to that field's bits -- found
// the same way a packed-struct member read finds them. Here one
// field is pinned to a constant and the other binds a name the arm
// then reads, so the match is only right if both the test and the
// slicing are.
typedef struct packed {
  logic [3:0] a;
  logic [3:0] b;
} p1_pair;

module pattern_case_structure (input logic clk, input logic [7:0] raw);
  p1_pair s;
  logic [3:0] y;
  always_comb begin
    s = raw;
    case (s) matches
      '{a: 4'd1, b: .bb} : y = bb;
      default            : y = 4'd0;
    endcase
  end
  assert property (@(posedge clk) raw[7:4] == 4'd1 |-> y == raw[3:0]);
endmodule
