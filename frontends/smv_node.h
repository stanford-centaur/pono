#pragma once

#include <stdio.h>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>

#include "assert.h"
#include "exceptions.h"
#include "smt-switch/smt.h"

/* node to help construct syntax tree in the parser */

struct node {
      smt::Term tm;
      smt::Sort st;
      std::string type;

      node (smt::Sort s){
        type = "SORT";
        st = s;
      }
      node (smt::Term t){
        type = "TERM";
        tm = t;
      }

      smt::Sort getSort(){
            return st;
      }
      smt::Term getTerm(){
            return tm;
      }

      std::string getType() { return type; }
};