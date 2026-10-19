// `intersect` used directly as the whole property (as opposed to
// inside a bounded sequence match feeding `|->`/`|=>`, which
// offsets_ending_now() does model -- see SeqIntersect in
// test_systemverilog_sva.cpp). ltl_to_sat()'s Binary dispatch has no
// gadget for a top-level Intersect/Within/Throughout connective.
module property_intersect_toplevel (
    input logic clk,
    input logic a,
    input logic b
);
  assert property (@(posedge clk) a intersect b);
endmodule
