// An initial value that depends on an input. An initial block
// describes the state at time 0, before any input has a value, so
// there is nothing for the index to name -- a core restriction
// rather than anything about the index being computed: the same
// write indexed by a state variable works.
module initial_value_from_input (
    input logic clk,
    input logic [1:0] i
);

  logic [3:0] arr;
  initial arr[i] = 1'b1;

  assert property (@(posedge clk) arr == arr);

endmodule
