// A `ref` argument, which aliases the caller's variable for the whole
// call rather than being copied in and back out, so the body's reads
// and writes interleave with the caller's.
//
// The same design with the call written out by hand -- `t = d;
// t = t ^ 8'hFF; q = t;` -- holds, and so must this. XOR rather than
// `+ 1` keeps the assertion from comparing an 8-bit result against a
// 32-bit expression.
module subroutine_ref_argument (input logic clk, input logic [7:0] d,
                                output logic [7:0] q);

  function automatic void invert(ref logic [7:0] x);
    x = x ^ 8'hFF;
  endfunction

  logic [7:0] t;
  always_comb begin
    t = d;
    invert(t);
    q = t;
  end

  assert property (@(posedge clk) q == (d ^ 8'hFF));

endmodule
