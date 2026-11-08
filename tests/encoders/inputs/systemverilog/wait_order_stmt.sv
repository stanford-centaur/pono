// `wait_order (a, b)` blocks until the named events arrive, and fails
// if they arrive out of order (LRM 15.5.4). Declaring the events is
// fine on its own; naming them here is what needs a notion of an
// event arriving, and of a process parked until it does.
module wait_order_stmt (
    input logic clk
);

  event a, b;
  logic [3:0] x;

  initial begin
    x = 4'd0;
    wait_order (a, b);
    x = 4'd5;
  end

  assert property (@(posedge clk) x == 4'd5);

endmodule
