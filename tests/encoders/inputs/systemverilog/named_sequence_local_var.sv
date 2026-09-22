// A named sequence with a local variable. `saved` captures `d` at
// the cycle the attempt starts and is read back one cycle later,
// which is what a local variable is for -- carrying a value along a
// single match attempt. slang expands a named assertion's body with
// arguments substituted, but a local variable needs a binding
// environment per attempt that the expansion does not supply.
//
// Intended: provable. `held` is `d` delayed one cycle, so it equals
// whatever `saved` captured. (Confirmed provable in that form: see
// the `$past` spelling of the same invariant, which passes today.)
module named_sequence_local_var (input logic clk, input logic [7:0] d);

  logic [7:0] held;
  always_ff @(posedge clk) held <= d;

  sequence capture;
    logic [7:0] saved;
    (1'b1, saved = d) ##1 (held == saved);
  endsequence

  assert property (@(posedge clk) capture);

endmodule
