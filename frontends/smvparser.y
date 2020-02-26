%{
    #include<cstdio>
    #include<iostream>
    #include<string>
    #include"encoder.h"
    #include"node.h"
    #include "smt-switch/smt.h"
    using namespace std;
    extern FILE *yyin;
    extern int yylex();
    extern int yyparse(smvEncoder &enc);
    extern int line_num;  
    extern void yyerror(smvEncoder &enc, const char*s){cout<< "error at" << line_num << s<<endl; exit(-1);};

%}

%code requires{
  #include "node.h"
  class smvEncoder;
}

%parse-param {smvEncoder &enc}

%define parse.error verbose
%locations


%union{
  node *n;
  char *str;
  int num;
  bool bl;
}

%token MODULE tok_main IVAR INVAR VAR FROZENVAR INVARSPEC
%token INIT TRANS READ WRITE ASSIGN CONSTARRAY CONSTANTS FUN DEFINE
%token INPUT OUTPUT
%token tok_next tok_init signed_word unsigned_word arrayword arrayinteger tok_array
%token pi ABS MAX MIN SIN COS EXP TAN ln of word1
%token tok_bool tok_toint tok_count swconst uwconst tok_sizeof tok_floor extend resize tok_typeof 
%token tok_unsigned tok_signed tok_word tok_set in time_type
%token TO ASSIGNSYM IF_ELSE 
%token ENDL

%token <num> integer_val real_val
%token <str> bool_type integer_type real_type clock_type set_tok array_tok word_index1 word_index2
%token <str> tok_name
%token <bl> TOK_TRUE TOK_FALSE

%left '!'
%left OP_CON
%left UMINUS
%left '*' '/' OP_MOD
%left '+' '-' 
%left OP_SHIFTR OP_SHIFTL
%left UNION
%left OP_EQ OP_NEQ '<' '>' OP_LTE OP_GTE
%left '&'
%left OP_OR OP_XOR OP_XNOR
%left OP_BI
%right OP_IMPLY

%type <n> type_identifier word_type array_type word_value basic_expr next_expr next_formula constant
%type <num> sizev 
%type <bl> boolean_constant
%type <str> complex_identifier define_identifier

%%
header: 
    ENDL
    | define_decl
    | constants_decl
    | assign_decl
    | ivar_test 
    | var_test
    | init_test
    | trans_test
    | frozenvar_test
    | invarspec_test
//   | function_test
    | header ENDL
    | header define_decl
    | header constants_decl
    | header assign_decl
    | header ivar_test
    | header var_test
    | header init_test
    | header trans_test
    | header frozenvar_test
    | header invarspec_test
//    | header function_test
    | module_header ENDL
    | module_header define_decl
    | module_header constants_decl
    | module_header assign_decl
    | module_header ivar_test
    | module_header var_test
    | module_header init_test
    | module_header frozenvar_test
    | module_header trans_test
    | module_header invarspec_test
//    | module_header function_test;

module_header:
    MODULE tok_main;

define_decl:
    DEFINE define_body 
    | DEFINE ENDL define_body;

define_body: complex_identifier ASSIGNSYM basic_expr ';' ENDL{
              node *a = $3;
              enc.terms_[$1] = a-> getTerm();
}
            | define_body complex_identifier ASSIGNSYM basic_expr ';' ENDL {
              node *a = $4;
              enc.terms_[$2] = a-> getTerm();
            };

constants_decl:
    CONSTANTS constants_body ;

constants_body: complex_identifier
            | constants_body ',' complex_identifier ';' ENDL;

assign_decl: ASSIGN assign_list 
            | ASSIGN ENDL assign_list;

assign_list: assign_test ';' ENDL
            | assign_list assign_test ';' ENDL;

assign_test: complex_identifier ASSIGNSYM basic_expr {
          node *a = $3; 
          smt::Term state = a->getTerm();
}
        | tok_init '('complex_identifier ')' ASSIGNSYM basic_expr{
          node *a = $6; 
          smt::Term init = enc.terms_[$3];
          smt::Term e = enc.solver_->make_term(smt::Equal, init, a->getTerm());
          enc.rts_.constrain_init(e);
        }
        | tok_next '('complex_identifier ')' ASSIGNSYM next_expr {
          node *a = $6; 
          smt::Term state = enc.terms_[$3];
          smt::Term e = enc.solver_->make_term(smt::Equal, state, a->getTerm());
          enc.rts_.constrain_trans(e);
        };

ivar_test:
    IVAR ENDL ivar_list ENDL 
    | IVAR ENDL ivar_list
    | ivar_test ivar_list ENDL
    | ivar_test ivar_list;


ivar_list:
    complex_identifier ':' type_identifier ';' ENDL{
         //cout <<"find an ivar"<<endl;
         node *a = $3;
         smt::Term input = enc.rts_.make_input($1, a->getSort());
         enc.terms_[$1] = input;
         //free($1);
    };

var_test:
    VAR ENDL var_list ENDL
   | var_test var_list ENDL
   | VAR ENDL var_list
   | var_test var_list;

var_list:
    complex_identifier ':' type_identifier ';' ENDL{
         node *a = $3;
         smt::Term state = enc.rts_.make_state($1, a->getSort());
         enc.terms_[$1] = state;
         //cout<<"find a var" <<endl;
         //free($1);
    };

frozenvar_test:
  FROZENVAR ENDL frozenvar_list ENDL
  | frozenvar_test frozenvar_list ENDL
  | FROZENVAR ENDL frozenvar_list
  | frozenvar_test frozenvar_test;

frozenvar_list:
  complex_identifier ':' type_identifier ';' ENDL{
      node *a = $3;
      smt::Term state = enc.rts_.make_state($1, a->getSort());
      enc.terms_[$1] = state;
      smt::Term n = enc.rts_.next(state);
      smt::Term e = enc.solver_->make_term(smt::Equal, n, state);
      enc.rts_.constrain_trans(e);
      //cout<<"find a frozen var" <<endl;
      //free($1);
  };

init_test: INIT basic_expr ';' ENDL ENDL{
        node *a = $2;
        enc.rts_.constrain_init(a->getTerm());
}
      | INIT basic_expr ';' ENDL{
        node *a = $2;
        enc.rts_.constrain_init(a->getTerm());
      };

trans_test: TRANS next_formula ';' ENDL ENDL{
            node *a = $2;
            enc.rts_.constrain_trans(a->getTerm());
            //cout <<"find a trans"<<endl;
}
            | TRANS next_formula ';' ENDL{
            node *a = $2;
            enc.rts_.constrain_trans(a->getTerm());
            //cout <<"find a trans"<<endl;
};
            

invarspec_test: INVARSPEC basic_expr ';' ENDL ENDL{
                //cout<<"find an invarspec" <<endl;
                node *a = $2;
                smt::Term prop = a->getTerm();
                enc.propvec_.push_back(prop);
}
                | INVARSPEC basic_expr ';' ENDL {
                //cout<<"find an invarspec" <<endl;
                node *a = $2;
                smt::Term prop = a->getTerm();
                enc.propvec_.push_back(prop);
                }
                | INVARSPEC basic_expr ';'{
                //cout<<"find an invarspec" <<endl;
                node *a = $2;
                smt::Term prop = a->getTerm();
                enc.propvec_.push_back(prop);
              };

//function_test: FUN function_list ';' ENDL ENDL
//              | FUN ENDL function_list ';' ENDL ;

//function_list: function_declaration
//              | function_list function_declaration
//function_declaration: complex_identifier ':' function_type_specifier
//function_type_specifier: function_args_type_specifier OP_IMPLY type_identifier
//function_args_type_specifier: type_identifier
//                | function_args_type_specifier '*' type_identifier
//fun_test: FUN ENDL basic_expr ':' integer_type '*' integer_type OP_IMPLY word_type ';' ENDL next_line;
//          | FUN ENDL basic_expr ':' integer_type '*' integer_type OP_IMPLY word_type ';' ENDL ;

next_formula: basic_expr {
              $$ = $1;
};

constant: boolean_constant {
      smt::Term con = enc.solver_->make_term($1);
      $$ = new node(con);
}
          | integer_constant 
          | real_constant
          | word_value {
           $$ = $1;
          }
          | clock_constant
          | symbolic_constant
          | range_constant;

word_value: word_index1 integer_val '_' integer_val {
          smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
          char* temp = $1;
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
          std::string int_val = std::to_string($4);
          smt::Term num = enc.solver_->make_term(int_val, sort_, base);
          $$ = new node(num); }
        | word_index2 integer_val '_' integer_val {
          smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
          char* temp = $1;
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
          std::string int_val = std::to_string($4);
          smt::Term num = enc.solver_->make_term(int_val, sort_, base);
          $$ = new node(num);
   };


boolean_constant: TOK_TRUE{
                $$ = true;
          }
                | TOK_FALSE{
                $$ = false;    
          };

define_identifier: complex_identifier{
                $$ = $1;
};

integer_constant: integer_val{
  throw CosaException("No integer constant now");
};
real_constant: real_val{
  throw CosaException("No real constant now");
};
range_constant: integer_val TO integer_val{
  throw CosaException("No range constant now");
};
clock_constant: time_type{
  throw CosaException("No clock constant now");
};
basic_expr: constant {
            $$ = $1;
}
            | complex_identifier {
              smt::Term tok = enc.terms_.at($1);
              $$ = new node(tok);
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
            | '(' basic_expr ')'{
              $$ = $2;
            }
            | '!' basic_expr
            | basic_expr '&' basic_expr {
              node *a = $1;
              node *b = $3;
              smt::Term e = enc.solver_->make_term(smt::BVAnd, a->getTerm(), b->getTerm());
              $$ = new node(e);
            }
            | basic_expr '|' basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term e = enc.solver_->make_term(smt::BVOr, a->getTerm(), b->getTerm());
              $$ = new node(e);
            }
            | basic_expr OP_XOR basic_expr {
              node *a = $1;
              node *b = $3;
              smt::Term e = enc.solver_->make_term(smt::BVXor, a->getTerm(), b->getTerm());
              $$ = new node(e);
            }
            | basic_expr OP_XNOR basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term e = enc.solver_->make_term(smt::BVXnor, a->getTerm(), b->getTerm());
              $$ = new node(e);
            }
            | basic_expr OP_IMPLY basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term e = enc.solver_->make_term(smt::Implies, a->getTerm(), b->getTerm());
              $$ = new node(e);
            }
            | basic_expr OP_BI basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term e = enc.solver_->make_term(smt::Iff, a->getTerm(), b->getTerm());
              $$ = new node(e);
            }
            | basic_expr '=' basic_expr {
              node *a = $1;
              node *b = $3;
              smt::Term e = enc.solver_->make_term(smt::Equal, a->getTerm(), b->getTerm());
              $$ = new node(e);
            }
            | basic_expr OP_NEQ basic_expr {
              node *a = $1;
              node *b = $3;
              smt::Term e = enc.solver_->make_term(smt::Distinct, a->getTerm(), b->getTerm());
              $$ = new node(e);
            }
            | basic_expr '<' basic_expr  {
              node *a = $1;
              node *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVUlt, a->getTerm(), b->getTerm());
              $$ = new node(res);
            }
            | basic_expr '>' basic_expr {
              node *a = $1;
              node *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVUgt, a->getTerm(), b->getTerm());
              $$ = new node(res);
            }
            | basic_expr OP_LTE basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVUle, a->getTerm(), b->getTerm());
              $$ = new node(res);
            }
            | basic_expr OP_GTE basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVUge, a->getTerm(), b->getTerm());
              $$ = new node(res);
            }
            | '-' basic_expr {

            }
            | basic_expr '+' basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVAdd, a->getTerm(), b->getTerm());
              $$ = new node(res);
            }
            | basic_expr '-' basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVSub, a->getTerm(), b->getTerm());
              $$ = new node(res);
            }
            | basic_expr '*' basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVMul, a->getTerm(), b->getTerm());
              $$ = new node(res);
            }
            | basic_expr '/' basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVUdiv, a->getTerm(), b->getTerm());
              $$ = new node(res);
            }
            | basic_expr OP_MOD basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVSmod, a->getTerm(), b->getTerm());
              $$ = new node(res);
            }
            | basic_expr OP_SHIFTR basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVLshr, a->getTerm(), b->getTerm());
              $$ = new node(res);
            }
            | basic_expr OP_SHIFTL basic_expr{
              node *a = $1;
              node *b = $3;
              smt::Term res = enc.solver_->make_term(smt::BVShl, a->getTerm(), b->getTerm());
              $$ = new node(res);
            }
            | basic_expr '[' integer_val ']'
            | basic_expr '[' basic_expr ':' basic_expr ']'
            | basic_expr OP_CON basic_expr  {
              node *a = $1;
              node *b = $3;
              smt::Term res = enc.solver_->make_term(smt::Concat, a->getTerm(), b->getTerm());
              $$ = new node(res);
            }
            | word1 '(' basic_expr ')'
            | tok_bool '(' basic_expr ')'
            | tok_toint '(' basic_expr ')'
            | tok_count '(' basic_expr_list ')'
            | swconst '(' basic_expr ')'
            | uwconst '(' basic_expr ')'
            | tok_signed '(' basic_expr ')'
            | tok_unsigned '(' basic_expr ')'
            | tok_sizeof '(' basic_expr ')'
            | tok_floor '(' basic_expr ')'
            | extend '(' basic_expr ')'
            | resize '(' basic_expr ')'
            | signed_word sizev '(' basic_expr ')'
            | unsigned_word sizev '(' basic_expr ')'
            | basic_expr UNION '(' basic_expr ')'
            |'{' set_body_expr '}'
            | basic_expr in basic_expr
            | basic_expr IF_ELSE basic_expr ':' basic_expr          
          | WRITE'('basic_expr','basic_expr','basic_expr')'{
            node *a = $3;
            node *b = $5;
            node *c = $7;
            smt::Term write_r =  enc.solver_->make_term(smt::Store, a->getTerm(),b->getTerm(),c->getTerm());
            $$ = new node(write_r);
          }
          | READ '(' basic_expr ',' basic_expr ')'{
            node *a = $3;
            node *b = $5;
            smt::Term read_r =  enc.solver_->make_term(smt::Select, a->getTerm(),b->getTerm());
            $$ = new node(read_r);
          }
          | CONSTARRAY '(' tok_typeof '(' complex_identifier ')' ',' basic_expr ')' 
          | CONSTARRAY '(' arrayword sizev of type_identifier ',' basic_expr ')' 
          | CONSTARRAY '(' arrayinteger of type_identifier ',' basic_expr ')' 
          | next_expr{
            $$ = $1;
          };

next_expr: tok_next '(' basic_expr ')'{
          node *a = $3;
          smt::Term n = enc.rts_.next(a->getTerm());
          $$ = new node(n);
};

basic_expr_list: basic_expr
                | basic_expr_list ',' basic_expr;

set_body_expr: basic_expr
                | set_body_expr ',' basic_expr;
        
symbolic_constant: complex_identifier;

complex_identifier: tok_name{
                $$ = $1;
 }
              | complex_identifier '.' tok_name{
                throw CosaException("No module component access now"); 
               //$$ = strcat($1,$3);
 }
              | complex_identifier '.' integer_val{
                throw CosaException("No module component access now"); 
                //$$ = $1 + $3;
 };

type_identifier: real_type{
                smt::Sort sort_ = enc.solver_->make_sort(smt::REAL);
                $$ =  new node (sort_);
                throw CosaException("No real type now in boolector");
                }
                | integer_type{
                  smt::Sort sort_ = enc.solver_->make_sort(smt::INT);
                  $$ =  new node (sort_);
                  throw CosaException("No integer type now in boolector");  
                }
                | bool_type {
                smt::Sort sort_ = enc.solver_->make_sort(smt::BOOL);
                $$ =  new node (sort_);
                }
                | clock_type{
                  throw CosaException("No clock type now");
                }
                | array_type{
                  $$ = $1;
                }
                | word_type {
                $$ = $1;
}           
                | integer_val TO integer_val;

word_type: signed_word sizev {
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
        $$ =  new node (sort_);
}
          | unsigned_word sizev{
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
        $$ =  new node (sort_);
}
          | tok_word sizev{
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, $2);
        $$ =  new node (sort_);
};

array_type: arrayword sizev of type_identifier{
            smt::Sort arraysort = enc.solver_->make_sort(smt::BV,$2);
            node *a = $4;
            smt::Sort sort_ = enc.solver_->make_sort(smt::ARRAY, arraysort,a->getSort());
            $$ = new node(sort_);
          }
          | arrayinteger of type_identifier{
            throw CosaException("no array integer type now");
          }
          | array_tok of type_identifier;

sizev:
    '[' integer_val ']'{
        $$  = $2;
    };
%%

void smvEncoder::parse(std::string filename){
    FILE *myfile = fopen(filename.c_str(), "r");
  // make sure ist's valid:
  if (!myfile) {
    std::cout << "NO input file!" << std::endl;
  }
  yyin = myfile;
  yyparse(*this); 
}
