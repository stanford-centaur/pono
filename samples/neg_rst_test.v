module neg_rst_test(clk, rst, cnt);
   input clk, rst;
   output reg [15:0] cnt;

   always @(posedge clk) begin
      if (!rst)
        cnt <= 16'd0;
      else if (cnt < 16'd10)
        cnt <= cnt + 16'd1;
      else
        cnt <= 16'd0;
   end

   always @* begin
      le_10: assert (cnt <= 16'd10);
   end

endmodule // neg_rst_test
