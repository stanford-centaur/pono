%require "3.4.2"
%{
    #include<cstdio>
    #include<iostream>
    #include<string>
    #include"stdlib.h"
    #include"smv_encoder.h"
    #include"smv_node.h"
    #include "smt-switch/smt.h"
    using namespace std;
    extern FILE *yyin;
    extern int yylineno;
%}

%code requires{
  #include "smv_node.h"
  #include <string>
  namespace cosa{
    class SMVEncoder;
  }
}

%skeleton "lalr1.cc"
%require "3.2"
%language "c++"
%defines 
%lex-param {SMVEncoder &enc}
%parse-param {SMVEncoder &enc}

%define parse.trace
%define parse.error verbose
%define api.token.constructor
%define api.value.type variant
%define api.namespace{cosa}
%define api.parser.class{smvparser}
%code {
  #include "smv_encoder.h"
}

%token MODULE tok_main IVAR INVAR VAR FROZENVAR INVARSPEC
%token INIT TRANS READ WRITE ASSIGN CONSTARRAY CONSTANTS FUN DEFINE
%token tok_next tok_init signed_word unsigned_word arrayword arrayinteger tok_array tok_case tok_esac
%token pi ABS MAX MIN SIN COS EXP TAN ln of word1
%token tok_bool tok_toint tok_count swconst uwconst tok_sizeof tok_floor extend resize tok_typeof 
%token tok_unsigned tok_signed tok_word tok_set OP_in time_type
%token TO ASSIGNSYM IF_ELSE 
%token ENDL

%token <std::string> integer_val real_val
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
%left OP_in
%left UNION
%left OP_SHIFTR OP_SHIFTL
%left "+" "-" 
%left "*" "/" OP_MOD
%precedence UMINUS
%left OP_CON
%left OP_NOT
%left "[" ":" "]"

%type <SMVnode*> type_identifier word_type array_type word_value basic_expr next_expr constant 
%type <int> sizev 
%type <bool> boolean_constant
%type <std::string> integer_constant real_constant
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
    MODULE tok_main;

define_decl:
    DEFINE define_body
    | define_decl define_body;

define_body: complex_identifier ASSIGNSYM basic_expr ";" { 
              SMVnode *a = $3;
              smt::Term define_var = a->getTerm();
              //cout << "find a define" <<$1 <<endl;
              enc.terms_[$1] = define_var;              
};

constants_decl:
    CONSTANTS constants_body ;

constants_body: complex_identifier
            | constants_body "," complex_identifier ";" ;

assign_decl: ASSIGN assign_list;

assign_list: assign_test ";" 
            | assign_list assign_test ";";

assign_test: complex_identifier ASSIGNSYM basic_expr {
          SMVnode *a = $3; 
          //smt::Term e = enc.solver_->make_term(smt::Equal, init, a->getTerm());
          smt::Term state = a->getTerm();
}
        | tok_init "(" complex_identifier ")" ASSIGNSYM basic_expr{
          SMVnode *a = $6; 
          smt::Term init = enc.terms_[$3];
          smt::Term e = enc.solver_->make_term(smt::Equal, init, a->getTerm());
          enc.rts_.constrain_init(e);
        }
        | tok_next "("complex_identifier ")" ASSIGNSYM next_expr {
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
         //cout <<"find an ivar"<<endl;
         SMVnode *a = $3;
         smt::Term input = enc.rts_.make_input($1, a->getSort());
         enc.terms_[$1] = input;
    };

var_test:
    VAR var_list
   | var_test var_list;

var_list:
    complex_identifier ":" type_identifier ";"{
         SMVnode *a = $3;
         smt::Term state = enc.rts_.make_state($1, a->getSort());
         enc.terms_[$1] = state;
         //cout<<"find a var" <<endl;
    };

frozenvar_test:
  FROZENVAR frozenvar_list
  | frozenvar_test frozenvar_list ;

frozenvar_list:
  complex_identifier ":" type_identifier ";" {
      SMVnode *a = $3;
      smt::Term state = enc.rts_.make_state($1, a->getSort());
      enc.terms_[$1] = state;
      smt::Term n = enc.rts_.next(state);
      smt::Term e = enc.solver_->make_term(smt::Equal, n, state);
      enc.rts_.constrain_trans(e);
      //cout<<"find a frozen var" <<endl;
  };

init_constraint: INIT init_list
                | init_constraint init_list;

init_list: basic_expr ";"{
        SMVnode *a = $1;
        enc.rts_.constrain_init(a->getTerm());
      };

trans_constraint: TRANS trans_list
                | trans_constraint trans_list;

trans_list: basic_expr ";"{
            SMVnode *a = $1;
            enc.rts_.constrain_trans(a->getTerm());
            //cout <<"find a trans"<<endl;
};

invar_constraint: INVAR invar_list
                | invar_constraint invar_list;

invar_list: basic_expr ";"{
            SMVnode *a = $1;
            enc.rts_.constrain_trans(a->getTerm());
            //cout <<"find a invar"<<endl;
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
            $$ = new SMVnode(con);
}
          | real_constant{
            smt::Sort sort_ = enc.solver_->make_sort(smt::REAL);
            smt::Term con = enc.solver_->make_term($1,sort_);
            $$ = new SMVnode(con);
          }
          | word_value {
           $$ = $1;
          }
//          | symbolic_constant{
//            throw CosaException("No symbolic constant now");
//          }
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
          $$ = new SMVnode(num); }
        | word_index2 integer_val "_" integer_val {
          smt::Sort sort_ = enc.solver_->make_sort(smt::BV, stoi($2));
          std::string temp = $1;
          int base = 2;
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
          $$ = new SMVnode(num);
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

integer_constant: integer_val{
  $$ = $1;
  //throw CosaException("No integer constant now");
};
real_constant: real_val{
  $$ = $1;
  //throw CosaException("No real constant now");
};
range_constant: integer_val TO integer_val{
  throw CosaException("No range constant now");
};
basic_expr: constant {
            $$ = $1;
}
            | complex_identifier {
              smt::Term tok = enc.terms_.at($1);
              $$ = new SMVnode(tok);
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
              smt::Term e = enc.solver_->make_term(smt::Not, a->getTerm());
              $$ = new SMVnode(e);
            }
            | basic_expr OP_AND basic_expr {
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term e = enc.solver_->make_term(smt::BVAnd, a->getTerm(), b->getTerm());
              $$ = new SMVnode(e);
            }
            | basic_expr OP_OR basic_expr{
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term e = enc.solver_->make_term(smt::BVOr, a->getTerm(), b->getTerm());
              $$ = new SMVnode(e);
            }
            | basic_expr OP_XOR basic_expr {
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term e = enc.solver_->make_term(smt::BVXor, a->getTerm(), b->getTerm());
              $$ = new SMVnode(e);
            }
            | basic_expr OP_XNOR basic_expr{
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term e = enc.solver_->make_term(smt::BVXnor, a->getTerm(), b->getTerm());
              $$ = new SMVnode(e);
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
              smt::Term e = enc.solver_->make_term(smt::Equal, a->getTerm(), b->getTerm());
              $$ = new SMVnode(e);
            }
            | basic_expr OP_NEQ basic_expr {
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term e = enc.solver_->make_term(smt::Distinct, a->getTerm(), b->getTerm());
              $$ = new SMVnode(e);
            }
            | basic_expr OP_LT basic_expr  {
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::SortKind kind_ = a->getTerm()->get_sort()->get_sort_kind();
              if ( (kind_ == smt::INT) || (kind_ == smt::REAL) ){
                  smt::Term res = enc.solver_->make_term(smt::Lt, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }else{
                  smt::Term res = enc.solver_->make_term(smt::BVUlt, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
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
                  smt::Term res = enc.solver_->make_term(smt::BVUgt, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
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
                  smt::Term res = enc.solver_->make_term(smt::BVUle, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
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
                  smt::Term res = enc.solver_->make_term(smt::BVUge, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }
            }
            | OP_MINUS basic_expr %prec UMINUS{
                SMVnode *a = $2;
                smt::Term res = enc.solver_->make_term(smt::BVNeg, a->getTerm());
                $$ = new SMVnode(res);
            }
            | basic_expr "+" basic_expr{
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::SortKind kind_ = a->getTerm()->get_sort()->get_sort_kind();
              if ( (kind_ == smt::INT) || (kind_ == smt::REAL) ){
                  smt::Term res = enc.solver_->make_term(smt::Plus, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }else{
                  smt::Term res = enc.solver_->make_term(smt::BVAdd, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
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
                  smt::Term res = enc.solver_->make_term(smt::BVSub, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
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
                  smt::Term res = enc.solver_->make_term(smt::BVMul, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
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
                  smt::Term res = enc.solver_->make_term(smt::BVUdiv, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
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
                  smt::Term res = enc.solver_->make_term(smt::BVSmod, a->getTerm(), b->getTerm());
                  $$ = new SMVnode(res);
              }
            }
            | basic_expr OP_SHIFTR basic_expr{
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVLshr, a->getTerm(), b->getTerm());
              $$ = new SMVnode(res);
            }
            | basic_expr OP_SHIFTL basic_expr{
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVShl, a->getTerm(), b->getTerm());
              $$ = new SMVnode(res);
            }
            | basic_expr sizev {
              throw CosaException("No word1");
            }
            | basic_expr OP_CON basic_expr  {
              SMVnode *a = $1;
              SMVnode *b = $3;
              smt::Term res = enc.solver_->make_term(smt::Concat, a->getTerm(), b->getTerm());
              $$ = new SMVnode(res);
            }
            | basic_expr "[" basic_expr ":" basic_expr "]"{
              throw CosaException("No word1");
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
            }
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
            | basic_expr OP_in basic_expr {
               throw CosaException("No array");
            }
            | basic_expr IF_ELSE basic_expr ":" basic_expr  {
               throw CosaException("No case now");
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
          | next_expr{
            $$ = $1;
          };

next_expr: tok_next "(" basic_expr ")"{
          SMVnode *a = $3;
          smt::Term n = enc.rts_.next(a->getTerm());
          $$ = new SMVnode(n);
};

basic_expr_list: basic_expr
                | basic_expr_list "," basic_expr;

set_body_expr: basic_expr
                | set_body_expr "," basic_expr;
        
//symbolic_constant: complex_identifier;

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
                $$ =  new SMVnode (sort_);
                //throw CosaException("No real type now in boolector");
                }
                | integer_type{
                  smt::Sort sort_ = enc.solver_->make_sort(smt::INT);
                  $$ =  new SMVnode (sort_);
                  //throw CosaException("No integer type now in boolector");  
                }
                | bool_type {
                smt::Sort sort_ = enc.solver_->make_sort(smt::BOOL);
                $$ =  new SMVnode (sort_);
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
        $$ =  new SMVnode (sort_);
}
          | unsigned_word sizev{
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
        $$ =  new SMVnode (sort_);
}
          | tok_word sizev{
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
        $$ =  new SMVnode (sort_);
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

void cosa::SMVEncoder::parse(std::string filename){
    FILE *myfile = fopen(filename.c_str(), "r");
    std::string file = filename; 
  // make sure file is valid:
  if (!myfile) {
    std::cout << "NO input file!" << std::endl;
    exit(-1);
  }
  yyin = myfile;
  cosa::smvparser parse (*this);
  parse();
}

void
cosa::smvparser::error (const std::string& m)
{
  std::cerr << yylineno << ":" << m << '\n';
   exit(-1);
}
