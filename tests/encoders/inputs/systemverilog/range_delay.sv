// Sequence range delay (`##[1:3]`), as opposed to the fixed `##N`
// already covered by sequence_delay.sv: if `arm` holds at a cycle,
// `data == 10` must hold at *some* cycle within the next 1-3 cycles.
// Both `arm` and `data` are free; the exact earliest falsifying bound
// is verified against the actual BMC run rather than hand-derived,
// since the "for some cycle in a window" semantics interacts with the
// per-cycle antecedent in a way that's easy to get subtly wrong by
// hand.
module range_delay (
    input logic clk,
    input logic rst,
    input logic arm,
    input logic [3:0] data
);

  assert property (@(posedge clk) arm |-> ##[1:3] data == 4'b1010);

endmodule
