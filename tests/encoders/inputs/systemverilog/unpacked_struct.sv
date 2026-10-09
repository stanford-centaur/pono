// An unpacked struct. It has no bit width of its own, which is why
// it had no sort and was rejected everywhere -- not only as an
// array element, despite how the gap was recorded. It does have a
// *selectable* width, though, and that is the space each field's
// own bitOffset is measured in, so laying the struct out flat there
// makes a field a bit range: exactly how a packed struct is already
// read and written.
//
// The fields have different widths deliberately, and each is
// written from its own input, so an offset that overlapped two of
// them would show. `m` covers an unpacked struct as an array
// element, and `t <= s` a whole-struct copy.
typedef struct {
  logic [3:0] a;
  logic [7:0] b;
  logic       c;
} su_t;

module unpacked_struct (
    input logic clk,
    input logic rst,
    input logic [3:0] da,
    input logic [7:0] db,
    input logic dc
);

  su_t s;
  su_t t;
  su_t m[0:1];

  always_ff @(posedge clk) begin
    s.a <= da;
    s.b <= db;
    s.c <= dc;
    t   <= s;
    m[0].a <= da;
    m[1].b <= db;
  end

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2))
                   |-> (s.a == $past(da) && s.b == $past(db)
                        && s.c == $past(dc)
                        && t.a == $past(s.a) && t.b == $past(s.b)
                        && m[0].a == $past(da) && m[1].b == $past(db)));

endmodule
