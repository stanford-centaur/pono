// Wildcard patterns wider than 64 bits, in both spellings that use
// them: `==?` and a casex arm. The mask and value are built as bit
// strings, so width is no longer what decides whether a pattern can
// be masked at all -- it used to fall back to plain equality past 64
// bits, which handed a pattern with unknown bits to the ordinary
// literal path and aborted the solver on an X in a decimal string.
//
// Each compare is asserted equivalent to comparing just the bits the
// pattern actually pins, which is what a wildcard means. The low and
// high variants together pin the bit order: a mask built the wrong
// way round could not satisfy both.
module wildcard_wide_pattern (
    input logic clk,
    input logic [71:0] d
);

  logic hit;

  always_comb begin
    hit = 1'b0;
    casex (d)
      72'hxxxxxxxxxxxxxxxx0F: hit = 1'b1;
      default: hit = 1'b0;
    endcase
  end

  assert property (@(posedge clk)
                   ((d ==? 72'hxxxxxxxxxxxxxxxx0F) == (d[7:0] == 8'h0F)
                    && (d ==? 72'h0Fxxxxxxxxxxxxxxxx) == (d[71:64] == 8'h0F)
                    && (d !=? 72'hxxxxxxxxxxxxxxxx0F) == (d[7:0] != 8'h0F)
                    && hit == (d[7:0] == 8'h0F)));

endmodule
