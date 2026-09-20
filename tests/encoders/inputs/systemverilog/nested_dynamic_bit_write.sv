// A runtime-indexed write whose base is itself a select (`p[2][j]`).
// The inner select resolves to a fixed bit range of p, so the write
// is a dynamic splice at that range's offset -- element 2's bits,
// never element 1's or 3's.
module nested_dynamic_bit_write (
    input logic clk,
    input logic rst,
    input logic [2:0] j,
    input logic val
);

  logic [3:0][7:0] p;

  always_ff @(posedge clk) p[2][j] <= val;

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (p[2][$past(j)] == $past(val)
                        && p[1] == $past(p[1]) && p[3] == $past(p[3])));

endmodule
