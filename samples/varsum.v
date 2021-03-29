module varsum(clk, pick);
   input clk;
   input [2:0] pick;

   reg [7:0]   x1 = 8'd1;
   reg [7:0]   x2 = 8'd1;
   reg [7:0]   x3 = 8'd1;
   reg [7:0]   x4 = 8'd1;
   reg [7:0]   x5 = 8'd1;

   always @(posedge clk) begin
      if (pick == 3'b000) begin
         x1 <= x1 + x2; x2 <= 8'd0;
      end
      else if (pick == 3'b001) begin
         x1 <= 8'd0; x2 <= x2 + x1;
      end
      else if (pick == 3'b010) begin
         x1 <= x1 + x3; x3 <= 8'd0;
      end
      else if (pick == 3'b011) begin
         x1 <= 8'd0; x3 <= x3 + x1;
      end
      else if (pick == 3'b100) begin
         x1 <= x1 + x4; x4 <= 8'd0;
      end
      else if (pick == 3'b101) begin
         x1 <= 8'd0; x4 <= x4 + x1;
      end
      else if (pick == 3'b110) begin
         x1 <= x1 + x5; x5 <= 8'd0;
      end
      else begin
         x1 <= 8'd0; x5 <= x5 + x1;
      end
   end

   always @(posedge clk) begin
      assert property (x1 <= 8'd5);
   end

endmodule // varsum
