// Packed-union construction via an assignment-pattern literal
// (`'{default: ...}`), distinct from union_overlap.sv's member-access
// through `.b`/`.parts`. expr_to_term()'s StructuredAssignmentPattern
// case only builds a PackedStructType target, so a union canonical
// type throws "unsupported assignment pattern target type".
typedef union packed {
  logic [7:0] b;
  struct packed {
    logic [3:0] hi;
    logic [3:0] lo;
  } parts;
} u8_t;

module union_literal (input logic clk, input logic rst);

  u8_t u;

  always_ff @(posedge clk) begin
    if (rst) u <= '{default: 8'h00};
    else u <= '{default: 8'h44};
  end

  assert property (@(posedge clk) rst || u.b == 8'h44);

endmodule
