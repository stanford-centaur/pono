// A `void` member, which LRM 7.3.2 allows precisely so the tag can
// carry all the information -- the valid-bit idiom. `tagged Invalid`
// contributes no payload at all, so the whole value below the tag is
// the undefined region, and the match has only the tag to go on.
typedef union tagged packed {
  void        Invalid;
  logic [7:0] Valid;
} tagged_union_void_t;

module tagged_union_void (
    input logic clk,
    input logic [7:0] v,
    input logic sel
);

  tagged_union_void_t u;
  logic [7:0] y;

  always_comb begin
    u = sel ? tagged_union_void_t'(tagged Valid (v))
            : tagged_union_void_t'(tagged Invalid);
    case (u) matches
      tagged Valid .n: y = n;
      tagged Invalid: y = 8'hFF;
      default: y = 8'd0;
    endcase
  end

  assert property (@(posedge clk)
                   (sel |-> (y == v)) and ((!sel) |-> (y == 8'hFF)));

endmodule
