// Three narrow unpacked-array forms that each used to throw.
//
// `neg` has a negative declared lower bound, so element 0 of the SMT
// array sits at index -1 and the index expression is signed. The two
// cells asserted sit either side of zero, so an offset computed the
// wrong way round either collides them or lands outside the array.
//
// `bits` is written at a runtime bit position inside a runtime-
// indexed element -- a dynamic splice within the element rather than
// a fixed range of it.
//
// `step` is incremented and decremented in place. An element is not
// a bit range of its base, so those never reached the write path
// that handles `++`/`--` for everything else.
module unpacked_array_small_forms (
    input logic clk,
    input logic rst,
    input logic [1:0] i,
    input logic [2:0] j,
    input logic [7:0] d
);

  logic [7:0] neg[-1:2];
  logic [7:0] bits[0:3];
  logic [7:0] step[0:3];

  always_ff @(posedge clk) begin
    if (rst) begin
      neg  <= '{default: 0};
      bits <= '{default: 0};
      step <= '{default: 0};
    end else begin
      neg[-1]   <= d;
      neg[2]    <= d + 8'd1;
      bits[i][j] <= 1'b1;
      step[1]++;
      step[2]--;
    end
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (neg[-1] == $past(d)
                        && neg[2] == $past(d) + 8'd1
                        && bits[$past(i)][$past(j)] == 1'b1
                        && step[1] == $past(step[1]) + 8'd1
                        && step[2] == $past(step[2]) - 8'd1));

endmodule
