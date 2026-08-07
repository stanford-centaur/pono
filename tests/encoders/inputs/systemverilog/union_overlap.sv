// A packed union aliasing a plain byte with a struct of two nibbles,
// to check that member-access bit offsets are computed correctly for
// union members (all starting at bit 0, unlike struct members which
// are packed end-to-end).  `u.b` is written by an always_ff counter
// stepping by 0x11 each cycle; `u.parts.hi`/`u.parts.lo` are read back
// through the *other* union member.  0x11 * 4 = 0x44 = 4'h4 :: 4'h4,
// so if the aliasing is correct the property holds at every cycle,
// including the one where u.b first equals 8'h44 (cycle 5).
typedef union packed {
  logic [7:0] b;
  struct packed {
    logic [3:0] hi;
    logic [3:0] lo;
  } parts;
} u8_t;

module union_overlap (input logic clk, input logic rst);

  u8_t u;

  always_ff @(posedge clk) begin
    if (rst) u.b <= 8'h00;
    else u.b <= u.b + 8'h11;
  end

  assert property (@(posedge clk)
                    !(u.b == 8'h44 && (u.parts.hi != 4'h4 || u.parts.lo != 4'h4)));

endmodule
