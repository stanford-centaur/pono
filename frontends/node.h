#ifndef NODE_H
#define NODE_H
#pragma once
#include <stdio.h>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>

#include "assert.h"
#include "exceptions.h"
#include "rts.h"
#include "smt-switch/smt.h"


struct node {
      smt::Term tm;
      smt::Sort st;

      node (smt::Sort s){
            st = s;
      }
      node (smt::Term t){
            tm = t;
      }

      smt::Sort getSort(){
            return st;
      }
      smt::Term getTerm(){
            return tm;
      }
};



#endif