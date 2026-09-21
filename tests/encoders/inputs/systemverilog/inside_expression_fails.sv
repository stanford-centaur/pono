// inside_expression.sv's range widened by one at the top end, and
// its wildcard member read as a plain equality. Both readings must
// be refuted, or the equivalences there would hold for an `inside`
// that ignored its bounds or its don't-care bits. Conjoined rather
// than disjoined because no single selector falsifies both: the
// range is wrong only at 4'd4 and the wildcard only at 4'b1101.
module inside_expression_fails (
    input logic clk,
    input logic [3:0] s
);

  assert property (@(posedge clk)
                   ((s inside {[4'd1:4'd3]}) == (s >= 4'd1 && s <= 4'd4))
                   && ((s inside {4'b1?01}) == (s == 4'b1001)));

endmodule
