// Compile-time-constant identities over replication and every sized-
// literal base/signedness form (binary, hex, decimal, negative signed),
// latched into a register so the check exercises a real next-state
// path rather than a bare literal comparison.  All are true by
// construction, so this genuinely holds forever once reset releases.
module replication_sized_literals (input logic clk, input logic rst);

  logic ok;

  always_ff @(posedge clk) begin
    if (rst) begin
      ok <= 1'b1;
    end else begin
      ok <= ({2{2'b10}} == 4'b1010) && (8'hFF == {8{1'b1}})
            && (-4'sd1 == 4'b1111) && (4'sd2 + (-4'sd1) == 4'sd1);
    end
  end

  assert property (@(posedge clk) ok);

endmodule
