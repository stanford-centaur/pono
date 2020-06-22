%require "3.4.2"
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

%token MODULE MAIN IVAR INVAR VAR FROZENVAR INVARSPEC
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

%nterm <SMVnode*> type_identifier word_type array_type word_value basic_expr next_expr case_expr constant simple_expr
%nterm <int> sizev 
%nterm <bool> boolean_constant
%nterm <std::string> integer_constant real_constant float_number fractional_number exponential_number
%nterm <std::string> complex_identifier 

%%

header:
     define_decl
    | constants_decl
    | assign_decl
    | ivar_test
    | var_test
    | frozenvar_test
    | init_constraint
    | trans_constraint
    | invar_constraint
    | invarspec_test
    | header define_decl
    | header constants_decl
    | header assign_decl
    | header ivar_test
    | header var_test
    | header frozenvar_test
    | header init_constraint
    | header trans_constraint
    | header invar_constraint
    | header invarspec_test
    | module_header define_decl
    | module_header constants_decl
    | module_header assign_decl
    | module_header ivar_test
    | module_header var_test
    | module_header frozenvar_test
    | module_header init_constraint
    | module_header trans_constraint
    | module_header invar_constraint
    | module_header invarspec_test;

module_header:
    MODULE MAIN;

define_decl:
    DEFINE define_body
    | define_decl define_body;

define_body: complex_identifier ASSIGNSYM basic_expr ";" { 
              SMVnode *a = $3;
              smt::Term define_var = a->getTerm();
              enc.terms_[$1] = define_var; 
              if (a->getBVType() == SMVnode::Unsigned){
                enc.unsignedbv_[$1] = define_var; 
              }else if (a->getBVType() == SMVnode::Signed){
                enc.signedbv_[$1] = define_var; 
              }             
};

constants_decl:
    CONSTANTS constants_body ;

constants_body: complex_identifier
            | constants_body "," complex_identifier ";" ;

assign_decl: ASSIGN assign_list;

assign_list: assign_test ";" 
            | assign_list assign_test ";";

assign_test: complex_identifier ASSIGNSYM simple_expr {
          SMVnode *a = $3; 
          //smt::Term e = enc.solver_->make_term(smt::Equal, init, a->getTerm());
          smt::Term state = a->getTerm();
          enc.rts_.constrain_trans(state);
}
        | TOK_INIT "(" complex_identifier ")" ASSIGNSYM simple_expr{
          SMVnode *a = $6; 
          smt::Term init = enc.terms_[$3];
          smt::Term e = enc.solver_->make_term(smt::Equal, init, a->getTerm());
          enc.rts_.constrain_init(e);
        }
        | TOK_NEXT "("complex_identifier ")" ASSIGNSYM next_expr {
          SMVnode *a = $6; 
          smt::Term state = enc.terms_[$3];
          smt::Term e = enc.solver_->make_term(smt::Equal, state, a->getTerm());
          enc.rts_.constrain_trans(e);
        };

ivar_test:
    IVAR ivar_list
    | ivar_test ivar_list;


ivar_list:
    complex_identifier ":" type_identifier ";" {
         SMVnode *a = $3;
<<<<<<< HEAD
         smt::Term input = enc.rts_.make_input($1, a->getSort());
=======
         smt::Term input = enc.rts_.make_inputvar($1, a->getSort());
>>>>>>> fix a conflict bug
         if (a->getBVType() == SMVnode::Unsigned){
           enc.unsignedbv_[$1] = input; 
         }else if (a->getBVType() == SMVnode::Signed){
           enc.signedbv_[$1] = input; 
         }
         enc.terms_[$1] = input;
    };

var_test:
    VAR var_list
   | var_test var_list;

var_list:
    complex_identifier ":" type_identifier ";"{
         SMVnode *a = $3;
         smt::Term state = enc.rts_.make_statevar($1, a->getSort());
         enc.terms_[$1] = state;
         if (a->getBVType() == SMVnode::Unsigned){
            enc.unsignedbv_[$1] = state; 
         } else if(a->getBVType() == SMVnode::Signed){
           enc.signedbv_[$1] = state;
         } 
    };

frozenvar_test:
  FROZENVAR frozenvar_list
  | frozenvar_test frozenvar_list ;

frozenvar_list:
  complex_identifier ":" type_identifier ";" {
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
      //cout<<"find a frozen var" <<endl;
  };

init_constraint: INIT init_list
                | init_constraint init_list;

init_list: simple_expr ";"{
        SMVnode *a = $1;
        enc.rts_.constrain_init(a->getTerm());
      };

trans_constraint: TRANS trans_list
                | trans_constraint trans_list;

trans_list: basic_expr ";"{
            SMVnode *a = $1;
            if(!case_true){
              enc.rts_.constrain_trans(a->getTerm());
            }else{
              enc.casestore_[case_start] = a->getTerm();
            }
            case_true = false;
            //cout <<"find a trans"<<endl;
};

invar_constraint: INVAR invar_list
                | invar_constraint invar_list;

invar_list: simple_expr ";"{
            SMVnode *a = $1;
            enc.rts_.constrain_trans(a->getTerm());
            enc.transterm_.push_back(make_pair(enc.loc.end.line,a->getTerm()));
};


invarspec_test: INVARSPEC invarspec_list;

invarspec_list: basic_expr ";" {
                //cout<<"find an invarspec" <<endl;
                SMVnode *a = $1;
                smt::Term prop = a->getTerm();
                enc.propvec_.push_back(prop);
};

constant: boolean_constant {
      smt::Term con = enc.solver_->make_term($1);
      $$ = new SMVnode(con);
}
          | integer_constant {
            smt::Sort sort_ = enc.solver_->make_sort(smt::INT);
            smt::Term con = enc.solver_->make_term($1,sort_);
            $$ = new SMVnode(con,SMVnode::Integer);
}
          | real_constant{
            smt::Sort sort_ = enc.solver_->make_sort(smt::REAL);
            smt::Term con = enc.solver_->make_term($1,sort_);
            $$ = new SMVnode(con);
          }
          | word_value {
           $$ = $1;
          }
          | range_constant{
            throw CosaException("No range constant now");
          };

word_value: word_index1 integer_val "_" integer_val {
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
        | word_index2 integer_val "_" integer_val {
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
   };


boolean_constant: TOK_TRUE{
                $$ = true;
          }
                | TOK_FALSE{
                $$ = false;
          };

//define_identifier: complex_identifier{
//                $$ = $1;
//};

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
              smt::Term tok = enc.terms_.at($1);
              if (enc.unsignedbv_.find($1) != enc.unsignedbv_.end() ) {
                $$ = new SMVnode(tok, SMVnode::Unsigned);
              } else if(enc.signedbv_.find($1) != enc.signedbv_.end()){
                $$ = new SMVnode(tok, SMVnode::Signed);
              } else{
                $$ = new SMVnode(tok);
              }
            }
            //| pi
            //| ABS '(' basic_expr ')'
            //| MAX '(' basic_expr ',' basic_expr ')'
            //| MIN '(' basic_expr ',' basic_expr ')'
            //| SIN '(' basic_expr ')'
            //| COS '(' basic_expr ')'
            //| EXP '(' basic_expr ')'
            //| TAN '(' basic_expr ')'
            //| ln '(' basic_expr ')'
            | "(" basic_expr ")"{
              $$ = $2;
            }
            | OP_NOT basic_expr {
              SMVnode *a = $2;
              SMVnode::BVtype bvs_a = a->getBVType();
              smt::Term e = enc.solver_->make_term(smt::BVNot, a->getTerm());
              $$ = new SMVnode(e,bvs_a);
            }
            | basic_expr OP_AND basic_expr {
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
            }
            | basic_expr OP_OR basic_expr{
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
            }
            | basic_expr OP_XOR basic_expr {
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
            }
            | basic_expr OP_XNOR basic_expr{
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
            }
            | basic_expr OP_IMPLY basic_expr{
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term e = enc.solver_->make_term(smt::Implies, a->getTerm(), b->getTerm());
              $$ = new SMVnode(e);
            }
            | basic_expr OP_BI basic_expr{
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term e = enc.solver_->make_term(smt::Iff, a->getTerm(), b->getTerm());
              $$ = new SMVnode(e);
            }
            | basic_expr OP_EQ basic_expr {
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
            }
            | basic_expr OP_NEQ basic_expr {
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
            }
            | basic_expr OP_LT basic_expr  {
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
            }
            | basic_expr OP_GT basic_expr {
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
            }
            | basic_expr OP_LTE basic_expr{
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
            }
            | basic_expr OP_GTE basic_expr{
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
            }
            | OP_MINUS basic_expr %prec UMINUS{
                SMVnode *a = $2;
                smt::Term res = enc.solver_->make_term(smt::BVNeg, a->getTerm());
                $$ = new SMVnode(res,a->getBVType());
            }
            | basic_expr "+" basic_expr{
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
            }
            | basic_expr "-" basic_expr{
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
            }
            | basic_expr "*" basic_expr{
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
            }
            | basic_expr "/" basic_expr{
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
            }
            | basic_expr OP_MOD basic_expr{
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
            }
            | basic_expr OP_SHIFTR basic_expr{
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVLshr, a->getTerm(), b->getTerm());
              $$ = new SMVnode(res,a->getBVType());
            }
            | basic_expr OP_SHIFTL basic_expr{
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVShl, a->getTerm(), b->getTerm());
              $$ = new SMVnode(res,a->getBVType());
            }
            | basic_expr sizev {
              throw CosaException("No word1");
            }
            | basic_expr OP_CON basic_expr  {
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
            }
            | basic_expr "[" integer_val "]" {
                throw CosaException("No index Subscript");
            }
            | basic_expr "[" basic_expr ":" basic_expr "]"{
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
                $$ = $3;
            }    //unsigned word convert to 
            | tok_unsigned "(" basic_expr ")"{
                $$ = $3;
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
                 SMVnode *a = $1;
                 SMVnode *b = $3;
                 SMVnode *c = $5;
                 smt::Term e = enc.solver_->make_term(smt::Ite, a->getTerm(),b->getTerm(),c->getTerm());
                 $$ = new SMVnode(e);
            }        
          | WRITE "(" basic_expr "," basic_expr "," basic_expr ")"{
            SMVnode *a = $3;
            SMVnode *b = $5;
            SMVnode *c = $7;
            smt::Term write_r =  enc.solver_->make_term(smt::Store, a->getTerm(),b->getTerm(),c->getTerm());
            $$ = new SMVnode(write_r);
          }
          | READ "(" basic_expr "," basic_expr ")"{
            SMVnode *a = $3;
            SMVnode *b = $5;
            smt::Term read_r =  enc.solver_->make_term(smt::Select, a->getTerm(),b->getTerm());
            $$ = new SMVnode(read_r);
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
          SMVnode *a = $3;
          smt::Term n = enc.rts_.next(a->getTerm());
          $$ = new SMVnode(n,a->getBVType());
};

case_expr: TOK_CASE case_body TOK_ESAC {
          std::pair<smt::Term,smt::Term> term_last = enc.caseterm_.back();
          smt::Term cond = term_last.first;
          smt::Term final_term = term_last.second;
          enc.caseterm_.pop_back();
          int total_num = enc.caseterm_.size();
          for (int i = 0; i < total_num; i++){
            //<<"find a case" <<std::endl;
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
}

case_body: basic_expr ":" basic_expr ";"{
      case_start = enc.loc.end.line;
      SMVnode *a = $1;
      SMVnode *b = $3;
      enc.caseterm_.clear();
      enc.caseterm_.push_back(make_pair(a->getTerm(),b->getTerm()));
} | case_body basic_expr ":" basic_expr ";" {
      SMVnode *a = $2;
      SMVnode *b = $4;
      enc.caseterm_.push_back(make_pair(a->getTerm(),b->getTerm()));
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
                smt::Sort sort_ = enc.solver_->make_sort(smt::REAL);
                $$ =  new SMVnode (sort_,SMVnode::BVnot);
                //throw CosaException("No real type now in boolector");
                }
                | integer_type{
                  smt::Sort sort_ = enc.solver_->make_sort(smt::INT);
                  $$ =  new SMVnode (sort_,SMVnode::BVnot);
                  //throw CosaException("No integer type now in boolector");  
                }
                | bool_type {
                smt::Sort sort_ = enc.solver_->make_sort(smt::BOOL);
                $$ =  new SMVnode (sort_,SMVnode::BVnot);
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
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
        $$ =  new SMVnode (sort_,SMVnode::Signed);
}
          | unsigned_word sizev{
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
        $$ =  new SMVnode (sort_,SMVnode::Unsigned);
}
          | tok_word sizev{
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
        $$ =  new SMVnode (sort_,SMVnode::Unsigned);
};

array_type: arrayword sizev of type_identifier{
            smt::Sort arraysort = enc.solver_->make_sort(smt::BV,$2);
            SMVnode *a = $4;
            smt::Sort sort_ = enc.solver_->make_sort(smt::ARRAY, arraysort,a->getSort());
            $$ = new SMVnode(sort_);
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