// An empty slot in a port list. It connects to nothing inside the
// instance and has no type, so reading its type reported "unsupported
// type kind" -- a throw that named neither the port nor the port
// list. There is nothing to declare; the rest of the module encodes.
module null_port_slot (
    clk,
    ,
    d,
    q
);

  input clk;
  input d;
  output q;

  assign q = d;

  assert property (@(posedge clk) q == d);

endmodule
