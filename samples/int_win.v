// a property that interpolation should beat induction on
module int_win(clk, a, b, out);
   input clk;
   input [15:0] a, b;
   output [15:0] out;

   reg [3:0] op;
   reg    cfg = 1'b1;
   reg [15:0] counter = 0;

   always @(posedge clk) begin
      counter <= counter + 16'd1;
   end

   always @(posedge clk) begin
      cfg <= 1'b0;
   end

   always @(posedge clk) begin
      if (cfg) begin
         op <= 4'd0;
      end
   end

   assign out = (op == 4'd0) ? a + b : a - b;

   always @* begin
      if (counter > 16'd0) begin
         assert(out == a + b);
      end
   end

endmodule // int_win
