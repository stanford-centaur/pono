// Package providing a typedef and a localparam, imported by
// pkg_macro_top.sv across a --sv-filelist boundary (this file is
// named only by pkg_filelist.f, not by pkg_macro_top.sv directly).
package common_pkg;
  typedef logic [3:0] nibble_t;
  localparam nibble_t STEP = 4'd2;
endpackage
