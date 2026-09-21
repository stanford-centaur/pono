// A `tagged` pattern matches on a tagged union's discriminant,
// which this encoder does not model -- see the unpacked-union
// exclusion. Refused rather than matched on the payload bits
// alone, which would answer as though the tag always agreed.
typedef union tagged packed {
  logic [6:0] Small;
  logic [6:0] Large;
} p3_u;

module pattern_case_tagged (input logic clk, input logic [7:0] raw);
  p3_u u;
  logic [6:0] y;
  always_comb begin
    u = raw;
    case (u) matches
      tagged Small .n : y = n;
      tagged Large .n : y = n;
    endcase
  end
  assert property (@(posedge clk) 1'b1);
endmodule
