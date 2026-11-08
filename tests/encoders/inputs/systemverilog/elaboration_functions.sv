// `$clog2`/`$bits` as elaboration-time functions sizing an RTL
// register. IDX_W = $clog2(5) = 3 sizes the `idx` register; BITW =
// $bits(logic [N-1:0]) should equal N. idx counts up from 0 and first
// reaches N-1 == 4 at cycle 5.
module elaboration_functions #(
    parameter int N = 5
) (
    input logic clk,
    input logic rst
);

  localparam int IDX_W = $clog2(N);
  localparam int BITW = $bits(logic [N-1:0]);

  logic [IDX_W-1:0] idx;

  always_ff @(posedge clk) begin
    if (rst) idx <= '0;
    else if (idx == IDX_W'(N - 1)) idx <= '0;
    else idx <= idx + 1'b1;
  end

  assert property (@(posedge clk) BITW == N && idx != IDX_W'(N - 1));

endmodule
