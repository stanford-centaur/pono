// A function whose body is more than a single return: locals, an
// assignment to the function's own name, a conditional overwrite of
// it, and a compile-time-unrolled loop. The call is inlined where it
// appears, so all of that has to run through the ordinary statement
// path with the formals and the return variable bound as locals.
//
// Each call is asserted equal to the expression it stands for, so a
// formal bound to the wrong actual, or a return value taken before
// the conditional overwrite, is refuted rather than merely encoded.
module user_function_body (
    input logic clk,
    input logic [3:0] a,
    input logic [3:0] b
);

  function automatic logic [4:0] add_clamp(logic [3:0] x, logic [3:0] y);
    logic [4:0] t;
    t = x + y;
    add_clamp = t;
    if (t > 5'd20) add_clamp = 5'd20;
  endfunction

  function automatic logic [2:0] popcount(logic [3:0] x);
    popcount = 3'd0;
    for (int i = 0; i < 4; i++) popcount = popcount + {2'd0, x[i]};
  endfunction

  assert property (@(posedge clk)
                   (add_clamp(a, b) == (((a + b) > 5'd20) ? 5'd20 : (a + b))
                    && popcount(a) == ({2'd0, a[0]} + {2'd0, a[1]}
                                       + {2'd0, a[2]} + {2'd0, a[3]})));

endmodule
