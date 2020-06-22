%{
  #include "smv_node.h"
  #include "smv_encoder.h"
  #include "smvparser.h" 
  #include "smvscanner.h"
  #include "location.hh"
  #include <string>
  #undef YY_DECL
  #define YY_DECL cosa::smvparser::symbol_type cosa::SMVscanner::yylex(SMVEncoder &encoder)
  using namespace std;
  # define YY_USER_ACTION  _encoder.loc.step();
  #undef yyterminate
  #define yyterminate()  return(cosa::smvparser::make_END())
%}
%option c++
%option yyclass="SMVscanner"
%option yylineno
%option noyywrap

%x COMMENT
%x WORD_INIT WORD_BASE WORD_VALUE
%%

[ \t]     ;
MODULE {return cosa::smvparser::make_MODULE(_encoder.loc);}
IVAR   {return cosa::smvparser::make_IVAR(_encoder.loc);}
VAR  {return cosa::smvparser::make_VAR(_encoder.loc);}
INVAR   {return cosa::smvparser::make_INVAR(_encoder.loc);}
FROZENVAR   {return cosa::smvparser::make_FROZENVAR(_encoder.loc);}
INVARSPEC   {return cosa::smvparser::make_INVARSPEC(_encoder.loc);}
INIT   {return cosa::smvparser::make_INIT(_encoder.loc);}
TRANS   {return cosa::smvparser::make_TRANS(_encoder.loc);}
READ   {return cosa::smvparser::make_READ(_encoder.loc);}
WRITE   {return cosa::smvparser::make_WRITE(_encoder.loc);}
ASSIGN   {return cosa::smvparser::make_ASSIGN(_encoder.loc);}
CONSTARRAY   {return cosa::smvparser::make_CONSTARRAY(_encoder.loc);}
CONSTANTS   {return cosa::smvparser::make_CONSTANTS(_encoder.loc);}
FUN   {return cosa::smvparser::make_FUN(_encoder.loc);}
DEFINE {return cosa::smvparser::make_DEFINE(_encoder.loc);}

next  {return cosa::smvparser::make_TOK_NEXT(_encoder.loc);}
case {return cosa::smvparser::make_TOK_CASE(_encoder.loc);}
esac {return cosa::smvparser::make_TOK_ESAC(_encoder.loc);}
init {return cosa::smvparser::make_TOK_INIT(_encoder.loc);}

TRUE   {return cosa::smvparser::make_TOK_TRUE(true,_encoder.loc);}
FALSE   {return cosa::smvparser::make_TOK_FALSE(false,_encoder.loc);}

"signed word" {return cosa::smvparser::make_signed_word(_encoder.loc);}
"unsigned word" {return cosa::smvparser::make_unsigned_word(_encoder.loc);}
"array word"   {return cosa::smvparser::make_arrayword(_encoder.loc);}
"array integer"   {return cosa::smvparser::make_arrayinteger(_encoder.loc);}
"array" {return cosa::smvparser::make_tok_array(_encoder.loc);}

pi {return cosa::smvparser::make_pi(_encoder.loc);}
abs {return cosa::smvparser::make_ABS(_encoder.loc);}
max {return cosa::smvparser::make_MAX(_encoder.loc);}
min {return cosa::smvparser::make_MIN(_encoder.loc);}
sin {return cosa::smvparser::make_SIN(_encoder.loc);}
cos {return cosa::smvparser::make_COS(_encoder.loc);}
exp {return cosa::smvparser::make_EXP(_encoder.loc);}
tan {return cosa::smvparser::make_TAN(_encoder.loc);}
ln {return cosa::smvparser::make_ln(_encoder.loc);}
of  {return cosa::smvparser::make_of(_encoder.loc);}
word1 {return cosa::smvparser::make_word1(_encoder.loc);}
bool {return cosa::smvparser::make_tok_bool(_encoder.loc);}
toint   {return cosa::smvparser::make_tok_toint(_encoder.loc);}
count {return cosa::smvparser::make_tok_count(_encoder.loc);}
swconst {return cosa::smvparser::make_swconst(_encoder.loc);}
uwconst {return cosa::smvparser::make_uwconst(_encoder.loc);}
sizeof {return cosa::smvparser::make_tok_sizeof(_encoder.loc);}
floor {return cosa::smvparser::make_tok_floor(_encoder.loc);}
extend {return cosa::smvparser::make_extend(_encoder.loc);}
resize {return cosa::smvparser::make_resize(_encoder.loc);}
typeof {return cosa::smvparser::make_tok_typeof(_encoder.loc);}

unsigned   {return cosa::smvparser::make_tok_unsigned(_encoder.loc);}
signed   {return cosa::smvparser::make_tok_signed(_encoder.loc);}
boolean   {return cosa::smvparser::make_bool_type(_encoder.loc);}
real  {return cosa::smvparser::make_real_type(_encoder.loc);}
integer   {return cosa::smvparser::make_integer_type(_encoder.loc);}
word   {return cosa::smvparser::make_tok_word(_encoder.loc);}
set {return cosa::smvparser::make_tok_set(_encoder.loc);}
time {return cosa::smvparser::make_time_type(_encoder.loc);}
in {return cosa::smvparser::make_OP_IN(_encoder.loc);}

[fF] {return cosa::smvparser::make_fraction_prefix(yytext,_encoder.loc);}
"->"  {return cosa::smvparser::make_OP_IMPLY(_encoder.loc);}
"<->"  {return cosa::smvparser::make_OP_BI(_encoder.loc);}
"!"  {return cosa::smvparser::make_OP_NOT(_encoder.loc);}
"::"  {return cosa::smvparser::make_OP_CON(_encoder.loc);}
"-"  {return cosa::smvparser::make_OP_MINUS(_encoder.loc);}
"+"  {return cosa::smvparser::make_OP_PLUS(_encoder.loc);}
"*"  {return cosa::smvparser::make_OP_MUL(_encoder.loc);}
"/"  {return cosa::smvparser::make_OP_DIV(_encoder.loc); }
"mod"  {return cosa::smvparser::make_OP_MOD(_encoder.loc);}
">>"  {return cosa::smvparser::make_OP_SHIFTR(_encoder.loc);}
"<<"  {return cosa::smvparser::make_OP_SHIFTL(_encoder.loc);}
"union"  {return cosa::smvparser::make_UNION(_encoder.loc);}
"="  {return cosa::smvparser::make_OP_EQ(_encoder.loc);}
"!="  {return cosa::smvparser::make_OP_NEQ(_encoder.loc);}
">"  {return cosa::smvparser::make_OP_GT(_encoder.loc);}
"<"  {return cosa::smvparser::make_OP_LT(_encoder.loc);}
"<="  {return cosa::smvparser::make_OP_LTE(_encoder.loc);}
">="  {return cosa::smvparser::make_OP_GTE(_encoder.loc);}
"&"  {return cosa::smvparser::make_OP_AND(_encoder.loc);}
"|"  {return cosa::smvparser::make_OP_OR(_encoder.loc);}
"xor"  {return cosa::smvparser::make_OP_XOR(_encoder.loc);}
"xnor"  {return cosa::smvparser::make_OP_XNOR(_encoder.loc);}

\.    {return cosa::smvparser::make_DOT(_encoder.loc);};
".." {return cosa::smvparser::make_TO(_encoder.loc);}
"(" {return cosa::smvparser::make_LPARE(_encoder.loc);}
")" {return cosa::smvparser::make_RPARE(_encoder.loc);}
"[" {return cosa::smvparser::make_LBRA(_encoder.loc);}
"]" {return cosa::smvparser::make_RBRA(_encoder.loc);}
"{" {return cosa::smvparser::make_LBRACE(_encoder.loc);}
"}" {return cosa::smvparser::make_RBRACE(_encoder.loc);}
";" {return cosa::smvparser::make_SEMICOLON(_encoder.loc);}
":" {return cosa::smvparser::make_COLON(_encoder.loc);}
"," {return cosa::smvparser::make_COMMA(_encoder.loc);}
":=" {return cosa::smvparser::make_ASSIGNSYM(_encoder.loc);}
"--".*  {}
"?"  {return cosa::smvparser::make_IF_ELSE(_encoder.loc);}
"_" {return cosa::smvparser::make_UNDER(_encoder.loc);}


0(b|d|o|h)/[0-9]+_[0-9a-fA-F]+ {BEGIN(WORD_INIT); return cosa::smvparser::make_word_index1(yytext,_encoder.loc);}
0(u|s)(b|d|o|h)/[0-9]+_[0-9a-fA-F]+ {BEGIN(WORD_INIT);  return cosa::smvparser::make_word_index2(yytext,_encoder.loc);}
<WORD_INIT>[0-9]+/_[0-9a-fA-F]+ {BEGIN(WORD_BASE); return cosa::smvparser::make_integer_val(yytext,_encoder.loc);}
<WORD_BASE>_/[0-9a-fA-F]+ {BEGIN(WORD_VALUE); return cosa::smvparser::make_UNDER(_encoder.loc);}
<WORD_VALUE>[0-9a-fA-F]+ {BEGIN(INITIAL); return cosa::smvparser::make_integer_val(yytext,_encoder.loc);}

[0-9]+           {return cosa::smvparser::make_integer_val(yytext,_encoder.loc);}
[a-zA-Z_][A-Za-z0-9_\$#-]*  {return cosa::smvparser::make_tok_name(yytext,_encoder.loc);}
\n   {_encoder.loc.lines(yyleng);}
.    {}
<<EOF>>    return cosa::smvparser::make_END(_encoder.loc);

%%
