// A local read on a path that never wrote it is not a temporary: it
// keeps its previous value, which is what synthesis infers storage
// for. Both shapes appear here.
//
// `held` is written only under `c`, so when !c it keeps what it had,
// and since y is assigned from it every cycle, y follows a when c
// and keeps its own previous value otherwise. `never` is written
// nowhere, so it is one unknown value that never changes, making z
// constant. Modelling either as a fresh unconstrained value each
// cycle refutes both conjuncts.
module procedural_temp_hold (
    input logic clk,
    input logic rst,
    input logic c,
    input logic [7:0] a
);

  logic [7:0] y;
  logic [7:0] z;

  always_ff @(posedge clk) begin
    logic [7:0] held;
    logic [7:0] never;
    if (c) held = a;
    y <= held;
    z <= never;
  end

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2))
                   |-> (y == ($past(c) ? $past(a) : $past(y))
                        && z == $past(z)));

endmodule
