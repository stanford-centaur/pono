// `default:` must apply only when no other case item matches: `s == 0`
// lets x reach 1 at the first post-reset transition, falsifying the
// property at cycle 2. A case implementation that applied `default`
// unconditionally would instead always drive x to 9, and the property
// would hold forever.
module case_default (input logic clk, input logic rst, input logic [1:0] s);

  logic [3:0] x;

  always_ff @(posedge clk) begin
    if (rst) x <= 4'd0;
    else
      case (s)
        2'd0: x <= 4'd1;
        2'd1: x <= 4'd2;
        default: x <= 4'd9;
      endcase
  end

  assert property (@(posedge clk) x != 4'd1);

endmodule
