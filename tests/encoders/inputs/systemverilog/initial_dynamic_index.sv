// A runtime-indexed write in an initial block: constrain_init() takes
// a fixed slice, and a runtime index names none.
module initial_dynamic_index (
    input logic clk,
    input logic [1:0] i
);

  logic [3:0] arr;

  initial arr[i] = 1'b1;

  assert property (@(posedge clk) arr == arr);

endmodule
