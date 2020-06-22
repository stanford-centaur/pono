%{
    #include<cstdio>
    #include<iostream>
    #include<string>
    #include"stdlib.h"
    #include"smv_encoder.h"
    #include"smv_node.h"
    using namespace std;
    int case_start = 0;
    bool case_true = false;
%}

%code requires{
  #include "smv_node.h"
  #include <string>
  namespace cosa{
    class SMVEncoder;
    class SMVscanner;
  }
}

%skeleton "lalr1.cc"
%require "3.2"
%locations 
%defines 
%lex-param {SMVEncoder &enc}
%parse-param {SMVscanner &smvscanner}
%parse-param {SMVEncoder &enc}

%code{
  #include "smvscanner.h"
  #undef yylex
  #define yylex smvscanner.yylex
}

%define parse.trace
%define parse.error verbose
%define api.token.constructor
%define api.value.type variant
%define api.namespace{cosa}
%define api.parser.class{smvparser}

%token MODULE IVAR INVAR VAR FROZENVAR INVARSPEC
%token INIT TRANS READ WRITE ASSIGN CONSTARRAY CONSTANTS FUN DEFINE TOK_CASE TOK_ESAC TOK_INIT
%token TOK_NEXT signed_word unsigned_word arrayword arrayinteger tok_array
%token pi ABS MAX MIN SIN COS EXP TAN ln of word1
%token tok_bool tok_toint tok_count swconst uwconst tok_sizeof tok_floor extend resize tok_typeof 
%token tok_unsigned tok_signed tok_word tok_set OP_IN time_type
%token TO ASSIGNSYM IF_ELSE 
%token ENDL

%token <std::string> integer_val real_val fraction_prefix exponential_prefix
%token bool_type integer_type real_type set_tok array_tok 
%token <std::string> word_index1 word_index2
%token <std::string> tok_name
%token <bool> TOK_TRUE TOK_FALSE
%token 
END 0
LPARE "("
RPARE ")"
COMMA  ","
UNDER "_" 
COLON ":"
SEMICOLON ";"
LBRA "["
RBRA "]"
LBRACE "{"
RBRACE "}" 
OP_PLUS "+"
OP_MINUS "-"
OP_MUL "*"
OP_DIV "/"
DOT ".";

%right OP_IMPLY
%left OP_BI
%left IF_ELSE
%left OP_OR OP_XOR OP_XNOR
%left OP_AND
%left OP_EQ OP_NEQ OP_LT OP_GT OP_LTE OP_GTE
%left OP_IN
%left UNION
%left OP_SHIFTR OP_SHIFTL
%left "+" "-" 
%left "*" "/" OP_MOD
%precedence UMINUS
%left OP_CON
%left OP_NOT
%left "[" ":" "]"

%nterm <SMVnode*> word_value basic_expr next_expr case_expr constant simple_expr 
%nterm <type_node*> type_identifier word_type array_type
%nterm <int> sizev 
%nterm <bool> boolean_constant
%nterm <std::string> integer_constant real_constant float_number fractional_number exponential_number
%nterm <std::string> complex_identifier 
%nterm <element_node*> module_element
%nterm <type_node*> module_type_identifier
%nterm <std::vector<string> > module_parameters
%nterm <std::unordered_map<SMVnode::NodeMtype,element_node*> > module_body
%nterm <std::vector<SMVnode*> > parameter_list case_body
%%

header:
    module_decl
    | header module_decl

module_decl:
    MODULE complex_identifier {
      if(!enc.module_flat){
      enc.module_list[$2] = new module_node($2);
      }
    }
  | MODULE complex_identifier module_body {
     if(!enc.module_flat){
      enc.module_list[$2] = new module_node($2,$3);
      enc.define_list_.clear();
      enc.ivar_list_.clear();
      enc.var_list_.clear();
      enc.init_list_.clear();
      enc.frozenvar_list_.clear();
      enc.trans_list_.clear();
      enc.invar_list_.clear();
      enc.invarspec_list_.clear();
     }
  }
  | MODULE complex_identifier "(" module_parameters ")" {
   if(!enc.module_flat){
    enc.module_list[$2] = new module_node($2,$4);
   }
  }
  | MODULE complex_identifier "(" module_parameters ")" module_body { 
    if(!enc.module_flat){
      enc.module_list[$2] = new module_node($2,$4,$6);
            enc.define_list_.clear();
      enc.ivar_list_.clear();
      enc.var_list_.clear();
      enc.init_list_.clear();
      enc.frozenvar_list_.clear();
      enc.trans_list_.clear();
      enc.invar_list_.clear();
      enc.invarspec_list_.clear();
    }
  }

module_parameters:
  complex_identifier {    
  if(!enc.module_flat){
    std::vector<string> id_list;
    id_list.push_back($1);
    $$ = id_list;
    }
  }
  | module_parameters "," complex_identifier{
  if(!enc.module_flat){
    std::vector<string> id_list;
    id_list = $1;
    id_list.push_back($3);
    $$ = id_list;
  }
  }

module_body:
  module_element{
    if(!enc.module_flat){
    element_node *a = $1;
    std::unordered_map<SMVnode::NodeMtype,element_node*> decl_map;
    decl_map[a->getNodeMType()] = $1;
    $$ = decl_map;
    }
  }
| module_body module_element{
  if(!enc.module_flat){
    std::unordered_map<SMVnode::NodeMtype,element_node*> decl_map = $1;
    element_node *a = $2;
    decl_map[a->getNodeMType()] = $2;
    $$ = decl_map;
  }
}

module_element:
     define_decl {
      if(!enc.module_flat){
       $$ = new define_node(enc.define_list_,SMVnode::DEFINE);
      }
     }
    | assign_decl {
      $$ = new assign_node(enc.assign_list_,SMVnode::ASSIGN);
    }
    | ivar_test {
      if(!enc.module_flat){
       $$ = new ivar_node(enc.ivar_list_,SMVnode::IVAR);
      }
    }
    | var_test {
      if(!enc.module_flat){
      $$ = new var_node(enc.var_list_,SMVnode::VAR);
      }
    }
    | frozenvar_test {
      if(!enc.module_flat){
      $$ = new frozenvar_node(enc.frozenvar_list_,SMVnode::FROZENVAR);
      }
    }
    | init_constraint{
     if(!enc.module_flat){
      $$ = new init_node(enc.init_list_,SMVnode::INIT);
     }
    }
    | trans_constraint{
    if(!enc.module_flat){
      $$ = new trans_node(enc.trans_list_,SMVnode::TRANS);
      }
    }
    | invar_constraint{
    if(!enc.module_flat){
      $$ = new invar_node(enc.invar_list_,SMVnode::INVAR);
    }
    }
    | invarspec_test{
   if(!enc.module_flat){
      $$ = new invarspec_node(enc.invarspec_list_,SMVnode::INVARSPEC);
    }
    }


define_decl:
    DEFINE define_body
    | define_decl define_body;

define_body: complex_identifier ASSIGNSYM basic_expr ";" { 
  if(enc.module_flat){
              SMVnode *a = $3;
              smt::Term define_var = a->getTerm();
              enc.terms_[$1] = define_var; 
              if (a->getBVType() == SMVnode::Unsigned){
                enc.unsignedbv_[$1] = define_var; 
              }else if (a->getBVType() == SMVnode::Signed){
                enc.signedbv_[$1] = define_var; 
              }    
  }else{
      enc.define_list_.push_back(new define_node_c($1,$3));
  }         
};

assign_decl: ASSIGN assign_list;

assign_list: assign_test ";" 
            | assign_list assign_test ";";

assign_test: complex_identifier ASSIGNSYM simple_expr {
  if(enc.module_flat){
          SMVnode *a = $3; 
          smt::Term e = enc.terms_[$1];
          smt::Term init = enc.solver_->make_term(smt::Equal, e, a->getTerm());
          enc.rts_.constrain_init(init);
          smt::Term n = enc.solver_->make_term(smt::Equal, e, a->getTerm());
          smt::Term next = enc.rts_.next(e);         
          e = enc.solver_->make_term(smt::Equal, next, n);
          enc.rts_.constrain_trans(e);
  }else{
      enc.assign_list_.push_back(new assign_node_c("",$1,$3));
  }   
}
        | TOK_INIT "(" complex_identifier ")" ASSIGNSYM simple_expr{
        if(enc.module_flat){
          SMVnode *a = $6; 
          smt::Term init = enc.terms_[$3];
          smt::Term e = enc.solver_->make_term(smt::Equal, init, a->getTerm());
          enc.rts_.constrain_init(e);
        }else{
         enc.assign_list_.push_back(new assign_node_c("init",$3,$6));
        }  
        }
        | TOK_NEXT "("complex_identifier ")" ASSIGNSYM next_expr {
        if(enc.module_flat){
          SMVnode *a = $6; 
          smt::Term state = enc.terms_[$3];
          smt::Term e = enc.solver_->make_term(smt::Equal, state, a->getTerm());
          enc.rts_.constrain_trans(e);
          }else{
          enc.assign_list_.push_back(new assign_node_c("next",$3,$6));
          }  
        };

ivar_test:
    IVAR ivar_list
    | ivar_test ivar_list;


ivar_list:
    complex_identifier ":" type_identifier ";" {
        if(enc.module_flat){
         SMVnode *a = $3;
         smt::Term input = enc.rts_.make_inputvar($1, a->getSort());
         if (a->getBVType() == SMVnode::Unsigned){
           enc.unsignedbv_[$1] = input; 
         }else if (a->getBVType() == SMVnode::Signed){
           enc.signedbv_[$1] = input; 
         }
         enc.terms_[$1] = input;
        }else{
          SMVnode *a = new ivar_node_c($1,$3);
          enc.ivar_list_.push_back(a);
        }
    };

var_test:
    VAR var_list
   | var_test var_list;

var_list:
    complex_identifier ":" type_identifier ";"{
      if(enc.module_flat){
         SMVnode *a = $3;
         smt::Term state = enc.rts_.make_statevar($1, a->getSort());
         enc.terms_[$1] = state;
         if (a->getBVType() == SMVnode::Unsigned){
            enc.unsignedbv_[$1] = state; 
         } else if(a->getBVType() == SMVnode::Signed){
           enc.signedbv_[$1] = state;
         } 
      }else{
          SMVnode *a = new var_node_c($1,$3,SMVnode::BasicT);
          enc.var_list_.push_back(new var_node_c($1,$3,SMVnode::BasicT));
      }
    }
    | complex_identifier ":" module_type_identifier ";"{
      if(enc.module_flat){
        throw CosaException("module preprocess error");
      }else{
         enc.var_list_.push_back(new var_node_c($1,$3,SMVnode::ModuleT));
      }
    }
    ;
  
  module_type_identifier:
    complex_identifier{
      $$ = new type_node($1);
    }
    | complex_identifier "(" parameter_list ")" {
      $$ = new type_node($1,$3);
    }
  
  parameter_list:
    simple_expr {
      std::vector<SMVnode*> s;
      s.push_back($1);
      $$ = s;
    }
    | parameter_list "," simple_expr {
      std::vector<SMVnode*> s= $1;
      s.push_back($3);
      $$ = s;
    }

frozenvar_test:
  FROZENVAR frozenvar_list
  | frozenvar_test frozenvar_list ;

frozenvar_list:
  complex_identifier ":" type_identifier ";" {
    if(enc.module_flat){
      SMVnode *a = $3;
      smt::Term state = enc.rts_.make_statevar($1, a->getSort());
      enc.terms_[$1] = state;
      if (a->getBVType() == SMVnode::Unsigned){
            enc.unsignedbv_[$1] = state; 
         } else if(a->getBVType() == SMVnode::Signed){
           enc.signedbv_[$1] = state;
         } 
      smt::Term n = enc.rts_.next(state);
      smt::Term e = enc.solver_->make_term(smt::Equal, n, state);
      enc.rts_.constrain_trans(e);
      enc.transterm_.push_back(make_pair(enc.loc.end.line,e));
    }else{
      SMVnode *a = new frozenvar_node_c($1,$3);
      enc.frozenvar_list_.push_back(a);
    }
  };

init_constraint: INIT init_list
                | init_constraint init_list;

init_list: simple_expr ";"{
  if(enc.module_flat){
        SMVnode *a = $1;
        enc.rts_.constrain_init(a->getTerm());
  }else{
    SMVnode *a = new init_node_c($1);
    enc.init_list_.push_back(a);
  }
      };

trans_constraint: TRANS trans_list
                | trans_constraint trans_list;

trans_list: basic_expr ";"{
  if(enc.module_flat){
            SMVnode *a = $1;
            if(!case_true){
              enc.rts_.constrain_trans(a->getTerm());
            }else{
              enc.casestore_[case_start] = a->getTerm();
            }
            case_true = false;
  }else{
    SMVnode *a = new trans_node_c($1);
    enc.trans_list_.push_back(a);
  }
};

invar_constraint: INVAR invar_list
                | invar_constraint invar_list;

invar_list: simple_expr ";"{
  if(enc.module_flat){
            SMVnode *a = $1;
            enc.rts_.constrain_trans(a->getTerm());
            enc.transterm_.push_back(make_pair(enc.loc.end.line,a->getTerm()));
  }else{
     SMVnode *a = new invar_node_c($1);
    enc.invar_list_.push_back(a);
  }
};


invarspec_test: INVARSPEC invarspec_list;

invarspec_list: basic_expr ";" {
  if(enc.module_flat){
                SMVnode *a = $1;
                smt::Term prop = a->getTerm();
                enc.propvec_.push_back(prop);
  }else{
    SMVnode *a = new invarspec_node_c($1);
    enc.invarspec_list_.push_back(a);
  }
};

constant: boolean_constant {
  if(enc.module_flat){
      smt::Term con = enc.solver_->make_term($1);
      $$ = new SMVnode(con);
  }else{
    if($1) $$ = new constant("TRUE");
    else $$ = new constant("FALSE");
  }
}
          | integer_constant {
            if(enc.module_flat){
            smt::Sort sort_ = enc.solver_->make_sort(smt::INT);
            smt::Term con = enc.solver_->make_term($1,sort_);
            $$ = new SMVnode(con,SMVnode::Integer);
            }else{
          $$ = new constant($1);
          }
}
          | real_constant{
            if(enc.module_flat){
            smt::Sort sort_ = enc.solver_->make_sort(smt::REAL);
            smt::Term con = enc.solver_->make_term($1,sort_);
            $$ = new SMVnode(con);
            }else{
            $$ = new constant($1);
            }
          }
          | word_value {
            if(enc.module_flat){
              $$ = $1;
            }else{
              SMVnode *wv = $1;
              string n = wv ->getName();
              $$ = new type_node(n);
            }
          }
          | range_constant{
            throw CosaException("No range constant now");
          };

word_value: word_index1 integer_val "_" integer_val {
      if(enc.module_flat){
          smt::Sort sort_ = enc.solver_->make_sort(smt::BV, stoi($2));
          std::string temp = $1;
          int base = 2;
          switch (temp[1]){
            case 'b':
              base = 2;
              break;
            case 'd':
              base = 10;
              break;
            case 'h':
              base = 16;
              break;
            default:
              base = 2;
          }
          smt::Term num = enc.solver_->make_term($4, sort_, base);
          $$ = new SMVnode(num, SMVnode::Unsigned); }
          else{
            string n = $1 + $2 + "_" + $4;
            $$ = new type_node(n);
        } }
        | word_index2 integer_val "_" integer_val {
        if(enc.module_flat){
          smt::Sort sort_ = enc.solver_->make_sort(smt::BV, stoi($2));
          std::string temp = $1;
          SMVnode::BVtype bvt;
          int base = 2;
          switch (temp[1]){
            case 'u':
              bvt = SMVnode::Unsigned;
              break;
            case 's':
              bvt = SMVnode::Signed;
              break;
            default:
              bvt = SMVnode::Unsigned;
          }
          switch (temp[2]){
            case 'b':
              base = 2;
              break;
            case 'd':
              base = 10;
              break;
            case 'h':
              base = 16;
              break;
            default:
              base = 2;
          }
          smt::Term num = enc.solver_->make_term($4, sort_, base);
          $$ = new SMVnode(num,bvt); 
        }else{
            string n = $1 + $2 + "_" + $4;
            $$ = new type_node(n);
        }
   };


boolean_constant: TOK_TRUE{
                $$ = true;
          }
                | TOK_FALSE{
                $$ = false;
          };

integer_constant: integer_val{ $$ = $1; };

real_constant: real_val{
  $$ = $1;
}| float_number { $$ = $1; }
| fractional_number{  $$ = $1; }
| exponential_number{ $$ = $1; };

float_number: integer_val "." integer_val{ $$ = $1 + "." + $3; }

fractional_number: fraction_prefix "'" integer_val "/" integer_val{ $$ = $1 + "'" + $3 + "/" + $5; }

exponential_number: integer_val exponential_prefix "-" integer_val{  $$ =  $1 + $2 + "-" + $4; }
| integer_val exponential_prefix integer_val{  $$ =  $1 + $2 + $3; }
| float_number exponential_prefix "-" integer_val {  $$ =  $1 + $2 + "-" + $4; };

range_constant: integer_val TO integer_val{
  throw CosaException("No range constant now");
};
basic_expr: simple_expr{ $$ = $1;}
          | next_expr{ $$ = $1;}

simple_expr: constant {
            $$ = $1; 
          }
            | complex_identifier {
            if(enc.module_flat){
              smt::Term tok = enc.terms_.at($1);
              if (enc.unsignedbv_.find($1) != enc.unsignedbv_.end() ) {
                $$ = new SMVnode(tok, SMVnode::Unsigned);
              } else if(enc.signedbv_.find($1) != enc.signedbv_.end()){
                $$ = new SMVnode(tok, SMVnode::Signed);
              } else{
                $$ = new SMVnode(tok);
              }
            }else{
              $$ = new identifier($1);
              }
            }
            | "(" basic_expr ")"{
              if(enc.module_flat){
              $$ = $2;
              }else{
              $$ = new par_expr($2);
              }
            }
            | OP_NOT basic_expr {
            if(enc.module_flat){
              SMVnode *a = $2;
              SMVnode::BVtype bvs_a = a->getBVType();
              smt::Term e = enc.solver_->make_term(smt::BVNot, a->getTerm());
              $$ = new SMVnode(e,bvs_a);
            }else{
              $$ = new not_expr($2);
              }
            }
            | basic_expr OP_AND basic_expr {
            if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              SMVnode::BVtype bvs_a = a->getBVType();
              SMVnode::BVtype bvs_b = b->getBVType();
              if(bvs_a != bvs_b){
                 throw CosaException("Unsigned/Signed unmatch");
              } else{
                smt::Term e = enc.solver_->make_term(smt::BVAnd, a->getTerm(), b->getTerm());
                $$ = new SMVnode(e,bvs_a);
              }
            }else{
              $$ = new and_expr($1,$3);
              }
            }
            | basic_expr OP_OR basic_expr{
            if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              SMVnode::BVtype bvs_a = a->getBVType();
              SMVnode::BVtype bvs_b = b->getBVType();
              if(bvs_a != bvs_b){
                 throw CosaException("Unsigned/Signed unmatch");
              } else{
                smt::Term e = enc.solver_->make_term(smt::BVOr, a->getTerm(), b->getTerm());
                $$ = new SMVnode(e,bvs_a);
              }
              }else{
              $$ = new or_expr($1,$3);
              }
            }
            | basic_expr OP_XOR basic_expr {
            if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              SMVnode::BVtype bvs_a = a->getBVType();
              SMVnode::BVtype bvs_b = b->getBVType();
              if(bvs_a != bvs_b){
                 throw CosaException("Unsigned/Signed unmatch");
              } else{
                smt::Term e = enc.solver_->make_term(smt::BVXor, a->getTerm(), b->getTerm());
                $$ = new SMVnode(e,bvs_a);
              }
              }else{
              $$ = new xor_expr($1,$3);
              }
            }
            | basic_expr OP_XNOR basic_expr{
              if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              SMVnode::BVtype bvs_a = a->getBVType();
              SMVnode::BVtype bvs_b = b->getBVType();
              if(bvs_a != bvs_b){
                 throw CosaException("Unsigned/Signed unmatch");
              } else{
              smt::Term e = enc.solver_->make_term(smt::BVXnor, a->getTerm(), b->getTerm());
              $$ = new SMVnode(e,bvs_a);
              }
              }else{
              $$ = new xor_expr($1,$3);
              }
            }
            | basic_expr OP_IMPLY basic_expr{
              if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term e = enc.solver_->make_term(smt::Implies, a->getTerm(), b->getTerm());
              $$ = new SMVnode(e);
              }else{
              $$ = new imp_expr($1,$3);
              }
            }
            | basic_expr OP_BI basic_expr{
              if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term e = enc.solver_->make_term(smt::Iff, a->getTerm(), b->getTerm());
              $$ = new SMVnode(e);
              }else{
              $$ = new iff_expr($1,$3);
              }
            }
            | basic_expr OP_EQ basic_expr {
              if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              SMVnode::BVtype bvs_a = a->getBVType();
              SMVnode::BVtype bvs_b = b->getBVType();
              if(bvs_a != bvs_b){
                 throw CosaException("Unsigned/Signed unmatch");
              } else{
              smt::Term e = enc.solver_->make_term(smt::Equal, a->getTerm(), b->getTerm());
              $$ = new SMVnode(e);
                }
              }else{
              $$ = new eq_expr($1,$3);
              }
            }
            | basic_expr OP_NEQ basic_expr {
              if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              SMVnode::BVtype bvs_a = a->getBVType();
              SMVnode::BVtype bvs_b = b->getBVType();
              if(bvs_a != bvs_b){
                 throw CosaException("Unsigned/Signed unmatch");
              } else{
                smt::Term e = enc.solver_->make_term(smt::Distinct, a->getTerm(), b->getTerm());
                $$ = new SMVnode(e);
              }
              }else{
              $$ = new neq_expr($1,$3);
              }
            }
            | basic_expr OP_LT basic_expr  {
            if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::SortKind kind_ = a->getTerm()->get_sort()->get_sort_kind();
              if ( (kind_ == smt::INT) || (kind_ == smt::REAL) ){
                  smt::Term res = enc.solver_->make_term(smt::Lt, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }else{
                  SMVnode::BVtype bvs_a = a->getBVType();
                  SMVnode::BVtype bvs_b = b->getBVType();
                  if (bvs_a == bvs_b == SMVnode::Unsigned){
                    smt::Term res = enc.solver_->make_term(smt::BVUlt, a->getTerm(), b->getTerm());
                    $$ = new SMVnode(res);
                  } else if (bvs_a == bvs_b == SMVnode::Signed){
                    smt::Term res = enc.solver_->make_term(smt::BVSlt, a->getTerm(), b->getTerm());
                    $$ = new SMVnode(res);
                  } else{
                    throw CosaException ("Unsigned/Signed unmatch");
                  } 
              }
              }else{
              $$ = new lt_expr($1,$3);
              }
            }
            | basic_expr OP_GT basic_expr {
            if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::SortKind kind_ = a->getTerm()->get_sort()->get_sort_kind();
              if ( (kind_ == smt::INT) || (kind_ == smt::REAL) ){
                  smt::Term res = enc.solver_->make_term(smt::Gt, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }else{
                  SMVnode::BVtype bvs_a = a->getBVType();
                  SMVnode::BVtype bvs_b = b->getBVType();
                  if (bvs_a == bvs_b == SMVnode::Unsigned){
                    smt::Term res = enc.solver_->make_term(smt::BVUgt, a->getTerm(), b->getTerm());
                    $$ = new SMVnode(res);
                  } else if (bvs_a == bvs_b == SMVnode::Signed){
                    smt::Term res = enc.solver_->make_term(smt::BVSgt, a->getTerm(), b->getTerm());
                    $$ = new SMVnode(res);
                  } else{
                    throw CosaException ("Unsigned/Signed unmatch");
                  } 
              }
              }else{
              $$ = new gt_expr($1,$3);
              }
            }
            | basic_expr OP_LTE basic_expr{
            if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::SortKind kind_ = a->getTerm()->get_sort()->get_sort_kind();
              if ((kind_ == smt::INT) || (kind_ == smt::REAL) ){
                  smt::Term res = enc.solver_->make_term(smt::Le, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }else{
                  SMVnode::BVtype bvs_a = a->getBVType();
                  SMVnode::BVtype bvs_b = b->getBVType();
                  if (bvs_a == bvs_b == SMVnode::Unsigned){
                    smt::Term res = enc.solver_->make_term(smt::BVUle, a->getTerm(), b->getTerm());
                    $$ = new SMVnode(res);
                  } else if (bvs_a == bvs_b == SMVnode::Signed){
                    smt::Term res = enc.solver_->make_term(smt::BVSle, a->getTerm(), b->getTerm());
                    $$ = new SMVnode(res);
                  } else{
                    throw CosaException ("Unsigned/Signed unmatch");
                  } 
              } 
              }else{
              $$ = new lte_expr($1,$3);
              }
            }
            | basic_expr OP_GTE basic_expr{
            if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::SortKind kind_ = a->getTerm()->get_sort()->get_sort_kind();
              if ((kind_ == smt::INT) || (kind_ == smt::REAL) ){
                  smt::Term res = enc.solver_->make_term(smt::Ge, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }else{
                  SMVnode::BVtype bvs_a = a->getBVType();
                  SMVnode::BVtype bvs_b = b->getBVType();
                  if (bvs_a == bvs_b == SMVnode::Unsigned){
                    smt::Term res = enc.solver_->make_term(smt::BVUge, a->getTerm(), b->getTerm());
                    $$ = new SMVnode(res);
                  } else if (bvs_a == bvs_b == SMVnode::Signed){
                    smt::Term res = enc.solver_->make_term(smt::BVSge, a->getTerm(), b->getTerm());
                    $$ = new SMVnode(res);
                  } else{
                    throw CosaException ("Unsigned/Signed unmatch");
                  } 
              }
              }else{
              $$ = new gte_expr($1,$3);
              }
            }
            | OP_MINUS basic_expr %prec UMINUS{
            if(enc.module_flat){
                SMVnode *a = $2;
                smt::Term res = enc.solver_->make_term(smt::BVNeg, a->getTerm());
                $$ = new SMVnode(res,a->getBVType());
              }else{
              $$ = new uminus_expr($2);
              }
            }
            | basic_expr "+" basic_expr{
            if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::SortKind kind_ = a->getTerm()->get_sort()->get_sort_kind();
              if ( (kind_ == smt::INT) || (kind_ == smt::REAL) ){
                  smt::Term res = enc.solver_->make_term(smt::Plus, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }else{
                  SMVnode::BVtype bvs_a = a->getBVType();
                  SMVnode::BVtype bvs_b = b->getBVType();
                  if(bvs_a != bvs_b){
                   throw CosaException("Unsigned/Signed unmatch");
                  } else{
                  smt::Term res = enc.solver_->make_term(smt::BVAdd, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res,bvs_a);
                }
              }
              }else{
              $$ = new add_expr($1,$3);
              }
            }
            | basic_expr "-" basic_expr{
            if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::SortKind kind_ = a->getTerm()->get_sort()->get_sort_kind();
              if ( (kind_ == smt::INT) || (kind_ == smt::REAL) ){
                  smt::Term res = enc.solver_->make_term(smt::Minus, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }else{
                  SMVnode::BVtype bvs_a = a->getBVType();
                  SMVnode::BVtype bvs_b = b->getBVType();
                  if(bvs_a != bvs_b){
                   throw CosaException("Unsigned/Signed unmatch");
                  } else{
                  smt::Term res = enc.solver_->make_term(smt::BVSub, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res,bvs_a);
                  }
              }
              }else{
              $$ = new sub_expr($1,$3);
              }
            }
            | basic_expr "*" basic_expr{
              if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::SortKind kind_ = a->getTerm()->get_sort()->get_sort_kind();
              if ( (kind_ == smt::INT) || (kind_ == smt::REAL) ){
                  smt::Term res = enc.solver_->make_term(smt::Mult, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }else{
                  SMVnode::BVtype bvs_a = a->getBVType();
                  SMVnode::BVtype bvs_b = b->getBVType();
                  if(bvs_a != bvs_b){
                   throw CosaException("Unsigned/Signed unmatch");
                  } else{
                  smt::Term res = enc.solver_->make_term(smt::BVMul, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res,bvs_a);
                  } 
              }
              }else{
              $$ = new mul_expr($1,$3);
              }
            }
            | basic_expr "/" basic_expr{
              if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::SortKind kind_ = a->getTerm()->get_sort()->get_sort_kind();
              if ( (kind_ == smt::INT) || (kind_ == smt::REAL) ){
                  smt::Term res = enc.solver_->make_term(smt::Div, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }else{
                  SMVnode::BVtype bvs_a = a->getBVType();
                  SMVnode::BVtype bvs_b = b->getBVType();
                  if (bvs_a == bvs_b == SMVnode::Unsigned){
                    smt::Term res = enc.solver_->make_term(smt::BVUdiv, a->getTerm(), b->getTerm());
                    $$ = new SMVnode(res,SMVnode::Unsigned);
                  } else if (bvs_a == bvs_b == SMVnode::Signed){
                    smt::Term res = enc.solver_->make_term(smt::BVSdiv, a->getTerm(), b->getTerm());
                    $$ = new SMVnode(res,SMVnode::Signed);
                  } else{
                    throw CosaException ("Unsigned/Signed unmatch");
                  } 
              }
              }else{
              $$ = new div_expr($1,$3);
              }
            }
            | basic_expr OP_MOD basic_expr{
              if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::SortKind kind_ = a->getTerm()->get_sort()->get_sort_kind();
              if ( (kind_ == smt::INT) || (kind_ == smt::REAL) ){
                  smt::Term res = enc.solver_->make_term(smt::Mod, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }else{
                  SMVnode::BVtype bvs_a = a->getBVType();
                  SMVnode::BVtype bvs_b = b->getBVType();
                  if (bvs_a == bvs_b == SMVnode::Unsigned){
                    smt::Term res = enc.solver_->make_term(smt::BVUrem, a->getTerm(), b->getTerm());
                    $$ = new SMVnode(res,SMVnode::Unsigned);
                  } else if (bvs_a == bvs_b == SMVnode::Signed){
                    smt::Term res = enc.solver_->make_term(smt::BVSmod, a->getTerm(), b->getTerm());
                    $$ = new SMVnode(res,SMVnode::Signed);
                  } else{
                    throw CosaException ("Unsigned/Signed unmatch");
                  } 
              }
              }else{
              $$ = new mod_expr($1,$3);
              }
            }
            | basic_expr OP_SHIFTR basic_expr{
              if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVLshr, a->getTerm(), b->getTerm());
              $$ = new SMVnode(res,a->getBVType());
              }else{
              $$ = new sr_expr($1,$3);
              }
            }
            | basic_expr OP_SHIFTL basic_expr{
              if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVShl, a->getTerm(), b->getTerm());
              $$ = new SMVnode(res,a->getBVType());
              }else{
              $$ = new sl_expr($1,$3);
              }
            }
            | basic_expr OP_CON basic_expr  {
              if(enc.module_flat){
              SMVnode *a = $1;
              SMVnode *b = $3;
              SMVnode::BVtype bvs_a = a->getBVType();
              SMVnode::BVtype bvs_b = b->getBVType();
              if(bvs_a != bvs_b){
                throw CosaException("Unsigned/Signed unmatch");
              } else{
                smt::Term res = enc.solver_->make_term(smt::Concat, a->getTerm(), b->getTerm());
                $$ = new SMVnode(res,SMVnode::Unsigned);
              }               
              }else{
                $$ = new con_expr($1,$3);
              }
            }
            | basic_expr sizev {
                throw CosaException("No index Subscript");
            }
            | basic_expr "[" basic_expr ":" basic_expr "]"{
              if(enc.module_flat){
                SMVnode *a = $1;
                SMVnode *b = $3;
                SMVnode *c = $5;
                SMVnode::BVtype bvs_a = a->getBVType();
                SMVnode::BVtype bvs_b = b->getBVType();
                SMVnode::BVtype bvs_c = c->getBVType();
                if(bvs_a == SMVnode::BVnot || bvs_b != SMVnode::Integer || bvs_c != SMVnode::Integer){
                  smt::Term res = enc.solver_->make_term(smt::Extract, a->getTerm(), b->getTerm(),c->getTerm());
                  $$ = new SMVnode(res,SMVnode::Unsigned);
                }else{
                  throw CosaException("Bit selection type is uncompatible");
                }
                }else{
              $$ = new sel_expr($1,$3,$5);
              }
            }
            | word1 "(" basic_expr ")" {
              throw CosaException("No word1");
            }
            | tok_bool "(" basic_expr ")"{
              throw CosaException("No type convert");
            }
            | tok_toint "(" basic_expr ")"{
              throw CosaException("No type convert");
            }
            | tok_count "(" basic_expr_list ")"{
              throw CosaException("No type convert");
            }
            | swconst "(" basic_expr ")"{
              throw CosaException("No type convert");
            }
            | uwconst "(" basic_expr ")"{
              throw CosaException("No type convert");
            }
            | tok_signed "(" basic_expr ")"{
                if(enc.module_flat){
                $$ = $3;
                }else{
                $$ = new signed_expr($3);
                }
            }    //unsigned word convert to 
            | tok_unsigned "(" basic_expr ")"{
              if(enc.module_flat){
                $$ = $3;
                }else{
              $$ = new unsigned_expr($3);
              }
            }
            | tok_sizeof "(" basic_expr ")"{
              throw CosaException("No array size");
            }
            | tok_floor "(" basic_expr ")"{
              throw CosaException("No floor operation");
            }
            | extend "(" basic_expr ")"{
              throw CosaException("No extend now");
            }
            | resize "(" basic_expr ")" {
              throw CosaException("No resize");
            }
            | signed_word sizev "(" basic_expr ")"{
               throw CosaException("No resize");
            }
            | unsigned_word sizev "(" basic_expr ")"{
               throw CosaException("No resize");
            }
            | basic_expr UNION  basic_expr {
               throw CosaException("No union");
            }
            |"{" set_body_expr "}" {
               throw CosaException("No set");
            }
            | basic_expr OP_IN basic_expr {
               throw CosaException("No array");
            }
            | basic_expr IF_ELSE basic_expr ":" basic_expr  {
              if(enc.module_flat){
                 SMVnode *a = $1;
                 SMVnode *b = $3;
                 SMVnode *c = $5;
                 smt::Term e = enc.solver_->make_term(smt::Ite, a->getTerm(),b->getTerm(),c->getTerm());
                 $$ = new SMVnode(e);
              }else{
                $$ = new ite_expr($1,$3,$5);
              }
            }        
          | WRITE "(" basic_expr "," basic_expr "," basic_expr ")"{
            if(enc.module_flat){
            SMVnode *a = $3;
            SMVnode *b = $5;
            SMVnode *c = $7;
            smt::Term write_r =  enc.solver_->make_term(smt::Store, a->getTerm(),b->getTerm(),c->getTerm());
            $$ = new SMVnode(write_r);
            }else{
              $$ = new write_expr($3,$5,$7);
            }
          }
          | READ "(" basic_expr "," basic_expr ")"{
            if(enc.module_flat){
            SMVnode *a = $3;
            SMVnode *b = $5;
            smt::Term read_r =  enc.solver_->make_term(smt::Select, a->getTerm(),b->getTerm());
            $$ = new SMVnode(read_r);
            }else{
              $$ = new read_expr($3,$5);
            }
          }
          | CONSTARRAY "(" tok_typeof "(" complex_identifier ")" "," basic_expr ")" {
             throw CosaException("No constarray");
          }
          | CONSTARRAY "(" arrayword sizev of type_identifier "," basic_expr ")" {
            throw CosaException("No constarray");
          }
          | CONSTARRAY "(" arrayinteger of type_identifier "," basic_expr ")" {
            throw CosaException("No constarray");
          }
          | case_expr {
            $$ = $1;
          };

next_expr: TOK_NEXT "(" basic_expr ")"{
  if(enc.module_flat){
          SMVnode *a = $3;
          smt::Term n = enc.rts_.next(a->getTerm());
          $$ = new SMVnode(n,a->getBVType());
  }else{
    $$ = new next_expr($3);
  }
};

case_expr: TOK_CASE case_body TOK_ESAC {
  if(enc.module_flat){
          std::pair<smt::Term,smt::Term> term_last = enc.caseterm_.back();
          smt::Term cond = term_last.first;
          smt::Term final_term = term_last.second;
          enc.caseterm_.pop_back();
          int total_num = enc.caseterm_.size();
          for (int i = 0; i < total_num; i++){
            std::pair<smt::Term,smt::Term> term_pair = enc.caseterm_.back();
            enc.caseterm_.pop_back();
            case_true = true;
            smt::Term e = enc.solver_->make_term(smt::Ite, term_pair.first, term_pair.second, final_term);
            cond = enc.solver_->make_term(smt::BVOr,cond, term_pair.first);
            final_term = e;
          }
          enc.casestore_[case_start] = final_term;
          enc.casecheck_[case_start] = cond;
          $$ = new SMVnode(final_term);
  }else{
    $$ = new case_expr($2);
  }
}

case_body: basic_expr ":" basic_expr ";"{
  if(enc.module_flat){
      case_start = enc.loc.end.line;
      SMVnode *a = $1;
      SMVnode *b = $3;
      enc.caseterm_.clear();
      enc.caseterm_.push_back(make_pair(a->getTerm(),b->getTerm()));
  }else{
    vector<SMVnode*> body;
    body.push_back(new case_body_ex($1,$3));
    $$ = body;
  }
} | case_body basic_expr ":" basic_expr ";" {
  if(enc.module_flat){
      SMVnode *a = $2;
      SMVnode *b = $4;
      enc.caseterm_.push_back(make_pair(a->getTerm(),b->getTerm()));
  }else{
    vector<SMVnode*> body = $1;
    body.push_back(new case_body_ex($2,$4));
    $$ = body;
  }
};

basic_expr_list: basic_expr
                | basic_expr_list "," basic_expr;

set_body_expr: basic_expr
                | set_body_expr "," basic_expr;

complex_identifier: tok_name{
                $$ = $1;
 }
              | complex_identifier "." tok_name{
                $$ = $1 + "." + $3;
 }
              | complex_identifier "." integer_val{
                $$ = $1 + "." + $3;
 };

type_identifier: real_type{
        if(enc.module_flat){
                smt::Sort sort_ = enc.solver_->make_sort(smt::REAL);
                $$ =  new type_node (sort_,SMVnode::BVnot);
        }else{
          $$ = new type_node("real");
        }
                }
                | integer_type{
                  if(enc.module_flat){
                  smt::Sort sort_ = enc.solver_->make_sort(smt::INT);
                  $$ =  new type_node (sort_,SMVnode::BVnot);
                  }else{
                   $$ = new type_node("integer");
                  }
                }
                | bool_type {
                  if(enc.module_flat){
                  smt::Sort sort_ = enc.solver_->make_sort(smt::BOOL);
                  $$ =  new type_node(sort_,SMVnode::BVnot);
                  }else{
                    $$ = new type_node("boolean");
                  }
                }
                | array_type{
                  $$ = $1 ;
                }
                | word_type {
                $$ = $1;
                }           
                | integer_val TO integer_val{
                  throw CosaException("No range now");
                };

word_type: signed_word sizev {
  if(enc.module_flat){
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
        $$ =  new type_node (sort_,SMVnode::Signed);
  }else{
    string n = "signed word [" + std::to_string($2) + "]";
    $$ = new type_node(n);
  }
}
          | unsigned_word sizev{
    if(enc.module_flat){
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
        $$ =  new type_node (sort_,SMVnode::Unsigned);
    }else{
        string n = "unsigned word [" + std::to_string($2) + "]";
         $$ = new type_node(n);
    }
}
          | tok_word sizev{
    if(enc.module_flat){
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
        $$ =  new type_node (sort_,SMVnode::Unsigned);
    }else{
        string n = "word [" + std::to_string($2) + "]";
        $$ = new type_node(n);
    }
};

array_type: arrayword sizev of type_identifier{
  if(enc.module_flat){
            smt::Sort arraysort = enc.solver_->make_sort(smt::BV,$2);
            SMVnode *a = $4;
            smt::Sort sort_ = enc.solver_->make_sort(smt::ARRAY, arraysort,a->getSort());
            $$ = new type_node(sort_);
  }else{
    SMVnode *temp = $4;
    string n = "array word" + std::to_string($2) + "of " + temp->getName();
    $$ = new type_node(n);
  }
          }
          | arrayinteger of type_identifier{
            throw CosaException("no array integer type now");
          }
          | array_tok of type_identifier{
            throw CosaException("No other array now");
          };

sizev:
    "[" integer_val "]"{
        $$  = stoi($2);
    };
%%

void cosa::smvparser::error (const location &loc, const std::string& m)
{ 
  std::cerr << loc << " " <<  m << '\n';
   exit(-1);
}