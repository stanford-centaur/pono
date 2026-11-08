// Writing a bit range *inside* an unpacked-array element, in the
// three shapes RTL uses: a bit select, a part select (byte enables),
// and a packed-struct field. None of these is a bit range of the
// array itself, so each is a Store over a Select-and-splice.
//
// Every write is at the same dynamic address, and the property reads
// that address back one cycle later, so a dropped write leaves the
// array unconstrained and is refuted rather than proved.
module unpacked_array_subelement (
    input logic clk,
    input logic rst,
    input logic [1:0] a,
    input logic [3:0] v
);

  typedef struct packed {
    logic [3:0] hi;
    logic [3:0] lo;
  } pair_t;

  logic [7:0] bitsel[0:3];
  logic [7:0] partsel[0:3];
  pair_t field[0:3];

  always_ff @(posedge clk) begin
    bitsel[a][5] <= 1'b1;
    partsel[a][3:0] <= v;
    field[a].hi <= v;
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (bitsel[$past(a)][5] == 1'b1
                        && partsel[$past(a)][3:0] == $past(v)
                        && field[$past(a)].hi == $past(v)));

endmodule
