%{
  #include"smv_node.h"
  #include "smv_encoder.h"
  #include "smvparser.h" 
  using namespace std;
  int line_num = 1;
  #include <string>
  #undef yyterminate
  #define yyterminate()  return(cosa::smvparser::make_END())
%}

%option yylineno
%option noyywrap
%x COMMENT
%x WORD_INIT WORD_BASE WORD_VALUE
%%

[ \t]     ;
MODULE {return cosa::smvparser::make_MODULE();}
main  {return cosa::smvparser::make_tok_main();}
IVAR   {return cosa::smvparser::make_IVAR();}
VAR  {return cosa::smvparser::make_VAR();}
INVAR   {return cosa::smvparser::make_INVAR();}
FROZENVAR   {return cosa::smvparser::make_FROZENVAR();}
INVARSPEC   {return cosa::smvparser::make_INVARSPEC();}
INIT   {return cosa::smvparser::make_INIT();}
TRANS   {return cosa::smvparser::make_TRANS();}
READ   {return cosa::smvparser::make_READ();}
WRITE   {return cosa::smvparser::make_WRITE();}
ASSIGN   {return cosa::smvparser::make_ASSIGN();}
CONSTARRAY   {return cosa::smvparser::make_CONSTARRAY();}
CONSTANTS   {return cosa::smvparser::make_CONSTANTS();}
FUN   {return cosa::smvparser::make_FUN();}
DEFINE {return cosa::smvparser::make_DEFINE();}

next  {return cosa::smvparser::make_tok_next();}
init {return cosa::smvparser::make_tok_init();}
case {return cosa::smvparser::make_tok_case();}
esac {return cosa::smvparser::make_tok_esac();}

TRUE   {return cosa::smvparser::make_TOK_TRUE(true);}
FALSE   {return cosa::smvparser::make_TOK_FALSE(false);}

"signed word" {return cosa::smvparser::make_signed_word();}
"unsigned word" {return cosa::smvparser::make_unsigned_word();}
"array word"   {return cosa::smvparser::make_arrayword();}
"array integer"   {return cosa::smvparser::make_arrayinteger();}
"array" {return cosa::smvparser::make_tok_array();}

pi {return cosa::smvparser::make_pi();}
abs {return cosa::smvparser::make_ABS();}
max {return cosa::smvparser::make_MAX();}
min {return cosa::smvparser::make_MIN();}
sin {return cosa::smvparser::make_SIN();}
cos {return cosa::smvparser::make_COS();}
exp {return cosa::smvparser::make_EXP();}
tan {return cosa::smvparser::make_TAN();}
ln {return cosa::smvparser::make_ln();}
of  {return cosa::smvparser::make_of();}
word1 {return cosa::smvparser::make_word1();}
bool {return cosa::smvparser::make_tok_bool();}
toint   {return cosa::smvparser::make_tok_toint();}
count {return cosa::smvparser::make_tok_count();}
swconst {return cosa::smvparser::make_swconst();}
uwconst {return cosa::smvparser::make_uwconst();}
sizeof {return cosa::smvparser::make_tok_sizeof();}
floor {return cosa::smvparser::make_tok_floor();}
extend {return cosa::smvparser::make_extend();}
resize {return cosa::smvparser::make_resize();}
typeof {return cosa::smvparser::make_tok_typeof();}

unsigned   {return cosa::smvparser::make_tok_unsigned();}
signed   {return cosa::smvparser::make_tok_signed();}
boolean   {return cosa::smvparser::make_bool_type();}
real  {return cosa::smvparser::make_real_type();}
integer   {return cosa::smvparser::make_integer_type();}
word   {return cosa::smvparser::make_tok_word();}
set {return cosa::smvparser::make_tok_set();}
in {return cosa::smvparser::make_OP_in();}
time {return cosa::smvparser::make_time_type();}

"->"  {return cosa::smvparser::make_OP_IMPLY();}
"<->"  {return cosa::smvparser::make_OP_BI();}
"!"  {return cosa::smvparser::make_OP_NOT();}
"::"  {return cosa::smvparser::make_OP_CON();}
"-"  {return cosa::smvparser::make_OP_MINUS();}
"+"  {return cosa::smvparser::make_OP_PLUS();}
"*"  {return cosa::smvparser::make_OP_MUL();}
"/"  {return cosa::smvparser::make_OP_DIV(); }
"mod"  {return cosa::smvparser::make_OP_MOD();}
">>"  {return cosa::smvparser::make_OP_SHIFTR();}
"<<"  {return cosa::smvparser::make_OP_SHIFTL();}
"union"  {return cosa::smvparser::make_UNION();}
"="  {return cosa::smvparser::make_OP_EQ();}
"!="  {return cosa::smvparser::make_OP_NEQ();}
">"  {return cosa::smvparser::make_OP_GT();}
"<"  {return cosa::smvparser::make_OP_LT();}
"<="  {return cosa::smvparser::make_OP_LTE();}
">="  {return cosa::smvparser::make_OP_GTE();}
"&"  {return cosa::smvparser::make_OP_AND();}
"|"  {return cosa::smvparser::make_OP_OR();}
"xor"  {return cosa::smvparser::make_OP_XOR();}
"xnor"  {return cosa::smvparser::make_OP_XNOR();}

".." {return cosa::smvparser::make_TO();}
"(" {return cosa::smvparser::make_LPARE();}
")" {return cosa::smvparser::make_RPARE();}
"[" {return cosa::smvparser::make_LBRA();}
"]" {return cosa::smvparser::make_RBRA();}
"{" {return cosa::smvparser::make_LBRACE();}
"}" {return cosa::smvparser::make_RBRACE();}
";" {return cosa::smvparser::make_SEMICOLON();}
":" {return cosa::smvparser::make_COLON();}
"," {return cosa::smvparser::make_COMMA();}
":=" {return cosa::smvparser::make_ASSIGNSYM();}
"--".*  {}
"?"  {return cosa::smvparser::make_IF_ELSE();}
"_" {return cosa::smvparser::make_UNDER();}

0(b|d|o|h)/[0-9]+_[0-9]+ {BEGIN(WORD_INIT); return cosa::smvparser::make_word_index1(yytext);}
0(u|s)(b|d|o|h)/[0-9]+_[0-9]+ {BEGIN(WORD_INIT);  return cosa::smvparser::make_word_index2(yytext);}
<WORD_INIT>[0-9]+/_[0-9]+ {BEGIN(WORD_BASE); return cosa::smvparser::make_integer_val(yytext);}
<WORD_BASE>_/[0-9]+ {BEGIN(WORD_VALUE); return cosa::smvparser::make_UNDER();}
<WORD_VALUE>[0-9]+ {BEGIN(INITIAL); return cosa::smvparser::make_integer_val(yytext);}

[0-9]+           {return cosa::smvparser::make_integer_val(yytext);}
[a-zA-Z_][A-Za-z0-9_\$#-]*  {return cosa::smvparser::make_tok_name(yytext);}
\n   {}
.                {return cosa::smvparser::make_DOT();};
<<EOF>>    return cosa::smvparser::make_END();

%%
