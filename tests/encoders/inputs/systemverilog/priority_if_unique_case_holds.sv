// Same shape as priority_if_unique_case.sv, but `sel` is hardwired
// low so `x` can never reach 4'd5 -- the property genuinely holds no
// matter what `s`/`y` do.
module priority_if_unique_case_holds (input logic clk, input logic rst, input logic [1:0] s);

  logic [3:0] x, y;
  logic       sel;

  assign sel = 1'b0;

  always_ff @(posedge clk) begin
    if (rst) begin
      x <= 4'd0;
      y <= 4'd0;
    end else begin
      priority if (sel) x <= 4'd5;
      else x <= 4'd3;

      unique case (s)
        2'd0: y <= 4'd1;
        2'd1: y <= 4'd2;
        2'd2: y <= 4'd4;
        2'd3: y <= 4'd8;
      endcase
    end
  end

  assert property (@(posedge clk) x != 4'd5 || y != 4'd8);

endmodule
