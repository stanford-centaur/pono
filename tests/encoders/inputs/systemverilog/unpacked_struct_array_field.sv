// An unpacked struct one of whose fields is itself an unpacked array.
// Unpacked structs and unpacked arrays are each modelled; this is the
// two composed, and a field with no bit width is what a flat layout
// has nowhere to put.
//
// Writing the scalar field and an element of the array field must
// each land where it was aimed.
typedef struct {
  logic [3:0] a;
  logic [3:0] v[0:1];
} sn_t;

module unpacked_struct_array_field (
    input logic clk,
    input logic [3:0] d
);

  sn_t s;
  logic started;
  initial started = 1'b0;

  always_ff @(posedge clk) begin
    started <= 1'b1;
    s.a <= d;
    s.v[0] <= ~d;
  end

  assert property (@(posedge clk) started
                   |-> (s.a == $past(d) && s.v[0] == ~$past(d)));

endmodule
