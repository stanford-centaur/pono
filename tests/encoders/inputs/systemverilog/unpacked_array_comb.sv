// Unpacked arrays built outside a clocked block: a combinational
// lookup table, a dynamic-index write over a whole-array default,
// and a whole-array assignment on its own. Each is read back at a
// runtime index, so an element that never received its write is
// refuted rather than left quietly unconstrained.
//
// A combinational block's writes compose in order and are pinned by
// one constraint per array once the block ends, the same way a
// non-wire scalar target is: a constraint per write would bind the
// same array term several times over.
module unpacked_array_comb (
    input logic clk,
    input logic [1:0] i,
    input logic [7:0] d
);

  logic [7:0] lut[0:3];
  logic [7:0] m[0:3];
  logic [7:0] w[0:3];

  always_comb begin
    lut[0] = 8'd10;
    lut[1] = 8'd11;
    lut[2] = 8'd12;
    lut[3] = 8'd13;
  end

  always_comb begin
    m = '{default: 0};
    m[i] = d;
  end

  always_comb w = '{default: 9};

  assert property (@(posedge clk) (lut[i] == (8'd10 + {6'd0, i})
                                   && m[i] == d
                                   && w[i] == 8'd9));

endmodule
