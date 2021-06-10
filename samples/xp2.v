module xp2(clk);
   input clk;
   reg [5:0] x = 6'd0;
   reg [5:0] y = 6'd2;

   wire      stop;
   assign stop = (x == 6'd61);

   always @(posedge clk) begin
      if (!stop) begin
         x <= x + 1;
         y <= y + 1;
      end
   end

   always @* begin
      x_lt_y: assert(x < y);
   end

endmodule // xp2
