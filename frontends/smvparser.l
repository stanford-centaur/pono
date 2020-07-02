%{
  #include"smv_node.h"
  #include "smv_encoder.h"
  #include "smvparser.h" 
  using namespace std;
  int line_num = 1;
  #include <string>
  #undef yyterminate
  #define yyterminate()  return(pono::smvparser::make_END())
%}

%option yylineno
%option noyywrap
%x COMMENT
%x WORD_INIT WORD_BASE WORD_VALUE
%%

[ \t]     ;
MODULE {return pono::smvparser::make_MODULE();}
main  {return pono::smvparser::make_tok_main();}
IVAR   {return pono::smvparser::make_IVAR();}
VAR  {return pono::smvparser::make_VAR();}
INVAR   {return pono::smvparser::make_INVAR();}
FROZENVAR   {return pono::smvparser::make_FROZENVAR();}
INVARSPEC   {return pono::smvparser::make_INVARSPEC();}
INIT   {return pono::smvparser::make_INIT();}
TRANS   {return pono::smvparser::make_TRANS();}
READ   {return pono::smvparser::make_READ();}
WRITE   {return pono::smvparser::make_WRITE();}
ASSIGN   {return pono::smvparser::make_ASSIGN();}
CONSTARRAY   {return pono::smvparser::make_CONSTARRAY();}
CONSTANTS   {return pono::smvparser::make_CONSTANTS();}
FUN   {return pono::smvparser::make_FUN();}
DEFINE {return pono::smvparser::make_DEFINE();}

next  {return pono::smvparser::make_tok_next();}
init {return pono::smvparser::make_tok_init();}
case {return pono::smvparser::make_tok_case();}
esac {return pono::smvparser::make_tok_esac();}

TRUE   {return pono::smvparser::make_TOK_TRUE(true);}
FALSE   {return pono::smvparser::make_TOK_FALSE(false);}

"signed word" {return pono::smvparser::make_signed_word();}
"unsigned word" {return pono::smvparser::make_unsigned_word();}
"array word"   {return pono::smvparser::make_arrayword();}
"array integer"   {return pono::smvparser::make_arrayinteger();}
"array" {return pono::smvparser::make_tok_array();}

pi {return pono::smvparser::make_pi();}
abs {return pono::smvparser::make_ABS();}
max {return pono::smvparser::make_MAX();}
min {return pono::smvparser::make_MIN();}
sin {return pono::smvparser::make_SIN();}
cos {return pono::smvparser::make_COS();}
exp {return pono::smvparser::make_EXP();}
tan {return pono::smvparser::make_TAN();}
ln {return pono::smvparser::make_ln();}
of  {return pono::smvparser::make_of();}
word1 {return pono::smvparser::make_word1();}
bool {return pono::smvparser::make_tok_bool();}
toint   {return pono::smvparser::make_tok_toint();}
count {return pono::smvparser::make_tok_count();}
swconst {return pono::smvparser::make_swconst();}
uwconst {return pono::smvparser::make_uwconst();}
sizeof {return pono::smvparser::make_tok_sizeof();}
floor {return pono::smvparser::make_tok_floor();}
extend {return pono::smvparser::make_extend();}
resize {return pono::smvparser::make_resize();}
typeof {return pono::smvparser::make_tok_typeof();}

unsigned   {return pono::smvparser::make_tok_unsigned();}
signed   {return pono::smvparser::make_tok_signed();}
boolean   {return pono::smvparser::make_bool_type();}
real  {return pono::smvparser::make_real_type();}
integer   {return pono::smvparser::make_integer_type();}
word   {return pono::smvparser::make_tok_word();}
set {return pono::smvparser::make_tok_set();}
in {return pono::smvparser::make_OP_in();}
time {return pono::smvparser::make_time_type();}

"->"  {return pono::smvparser::make_OP_IMPLY();}
"<->"  {return pono::smvparser::make_OP_BI();}
"!"  {return pono::smvparser::make_OP_NOT();}
"::"  {return pono::smvparser::make_OP_CON();}
"-"  {return pono::smvparser::make_OP_MINUS();}
"+"  {return pono::smvparser::make_OP_PLUS();}
"*"  {return pono::smvparser::make_OP_MUL();}
"/"  {return pono::smvparser::make_OP_DIV(); }
"mod"  {return pono::smvparser::make_OP_MOD();}
">>"  {return pono::smvparser::make_OP_SHIFTR();}
"<<"  {return pono::smvparser::make_OP_SHIFTL();}
"union"  {return pono::smvparser::make_UNION();}
"="  {return pono::smvparser::make_OP_EQ();}
"!="  {return pono::smvparser::make_OP_NEQ();}
">"  {return pono::smvparser::make_OP_GT();}
"<"  {return pono::smvparser::make_OP_LT();}
"<="  {return pono::smvparser::make_OP_LTE();}
">="  {return pono::smvparser::make_OP_GTE();}
"&"  {return pono::smvparser::make_OP_AND();}
"|"  {return pono::smvparser::make_OP_OR();}
"xor"  {return pono::smvparser::make_OP_XOR();}
"xnor"  {return pono::smvparser::make_OP_XNOR();}

".." {return pono::smvparser::make_TO();}
"(" {return pono::smvparser::make_LPARE();}
")" {return pono::smvparser::make_RPARE();}
"[" {return pono::smvparser::make_LBRA();}
"]" {return pono::smvparser::make_RBRA();}
"{" {return pono::smvparser::make_LBRACE();}
"}" {return pono::smvparser::make_RBRACE();}
";" {return pono::smvparser::make_SEMICOLON();}
":" {return pono::smvparser::make_COLON();}
"," {return pono::smvparser::make_COMMA();}
":=" {return pono::smvparser::make_ASSIGNSYM();}
"--".*  {}
"?"  {return pono::smvparser::make_IF_ELSE();}
"_" {return pono::smvparser::make_UNDER();}

0(b|d|o|h)/[0-9]+_[0-9]+ {BEGIN(WORD_INIT); return pono::smvparser::make_word_index1(yytext);}
0(u|s)(b|d|o|h)/[0-9]+_[0-9]+ {BEGIN(WORD_INIT);  return pono::smvparser::make_word_index2(yytext);}
<WORD_INIT>[0-9]+/_[0-9]+ {BEGIN(WORD_BASE); return pono::smvparser::make_integer_val(yytext);}
<WORD_BASE>_/[0-9]+ {BEGIN(WORD_VALUE); return pono::smvparser::make_UNDER();}
<WORD_VALUE>[0-9]+ {BEGIN(INITIAL); return pono::smvparser::make_integer_val(yytext);}

[0-9]+           {return pono::smvparser::make_integer_val(yytext);}
[a-zA-Z_][A-Za-z0-9_\$#-]*  {return pono::smvparser::make_tok_name(yytext);}
\n   {}
.                {return pono::smvparser::make_DOT();};
<<EOF>>    return pono::smvparser::make_END();

%%
