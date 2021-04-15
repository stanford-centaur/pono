%{
  #include "frontends/smv_node.h"
  #include "frontends/smv_encoder.h"
  #include "smvparser.h"
  #include "frontends/smvscanner.h"
  #include "location.hh"
  #include <string>
  #undef YY_DECL
  #define YY_DECL pono::smvparser::symbol_type pono::SMVscanner::yylex(SMVEncoder &encoder)
  using namespace std;
  # define YY_USER_ACTION  _encoder.loc.step();
  #undef yyterminate
  #define yyterminate()  return(pono::smvparser::make_END())
%}
%option c++
%option yyclass="SMVscanner"
%option yylineno
%option noyywrap
/** Uncomment to enable debug tracing */
/** %option debug */

%x COMMENT
%x WORD_INIT WORD_BASE WORD_VALUE
%%

[ \t]     ;
MODULE {yylineno = 1; return pono::smvparser::make_MODULE(_encoder.loc);}
IVAR   {return pono::smvparser::make_IVAR(_encoder.loc);}
VAR  {return pono::smvparser::make_VAR(_encoder.loc);}
INVAR   {return pono::smvparser::make_INVAR(_encoder.loc);}
FROZENVAR   {return pono::smvparser::make_FROZENVAR(_encoder.loc);}
INVARSPEC   {return pono::smvparser::make_INVARSPEC(_encoder.loc);}
INIT   {return pono::smvparser::make_INIT(_encoder.loc);}
TRANS   {return pono::smvparser::make_TRANS(_encoder.loc);}
READ   {return pono::smvparser::make_READ(_encoder.loc);}
WRITE   {return pono::smvparser::make_WRITE(_encoder.loc);}
ASSIGN   {return pono::smvparser::make_ASSIGN(_encoder.loc);}
CONSTARRAY   {return pono::smvparser::make_CONSTARRAY(_encoder.loc);}
CONSTANTS   {return pono::smvparser::make_CONSTANTS(_encoder.loc);}
FUN   {return pono::smvparser::make_FUN(_encoder.loc);}
DEFINE {return pono::smvparser::make_DEFINE(_encoder.loc);}

next  {return pono::smvparser::make_TOK_NEXT(_encoder.loc);}
case {return pono::smvparser::make_TOK_CASE(_encoder.loc);}
esac {return pono::smvparser::make_TOK_ESAC(_encoder.loc);}
init {return pono::smvparser::make_TOK_INIT(_encoder.loc);}

TRUE   {return pono::smvparser::make_TOK_TRUE(true,_encoder.loc);}
FALSE   {return pono::smvparser::make_TOK_FALSE(false,_encoder.loc);}

"signed word" {return pono::smvparser::make_signed_word(_encoder.loc);}
"unsigned word" {return pono::smvparser::make_unsigned_word(_encoder.loc);}
"array word"   {return pono::smvparser::make_arrayword(_encoder.loc);}
"array integer"   {return pono::smvparser::make_arrayinteger(_encoder.loc);}
"array" {return pono::smvparser::make_tok_array(_encoder.loc);}

pi {return pono::smvparser::make_pi(_encoder.loc);}
abs {return pono::smvparser::make_ABS(_encoder.loc);}
max {return pono::smvparser::make_MAX(_encoder.loc);}
min {return pono::smvparser::make_MIN(_encoder.loc);}
sin {return pono::smvparser::make_SIN(_encoder.loc);}
cos {return pono::smvparser::make_COS(_encoder.loc);}
exp {return pono::smvparser::make_EXP(_encoder.loc);}
tan {return pono::smvparser::make_TAN(_encoder.loc);}
ln {return pono::smvparser::make_ln(_encoder.loc);}
of  {return pono::smvparser::make_of(_encoder.loc);}
word1 {return pono::smvparser::make_word1(_encoder.loc);}
bool {return pono::smvparser::make_tok_bool(_encoder.loc);}
toint   {return pono::smvparser::make_tok_toint(_encoder.loc);}
count {return pono::smvparser::make_tok_count(_encoder.loc);}
swconst {return pono::smvparser::make_swconst(_encoder.loc);}
uwconst {return pono::smvparser::make_uwconst(_encoder.loc);}
sizeof {return pono::smvparser::make_tok_sizeof(_encoder.loc);}
floor {return pono::smvparser::make_tok_floor(_encoder.loc);}
extend {return pono::smvparser::make_extend(_encoder.loc);}
resize {return pono::smvparser::make_resize(_encoder.loc);}
typeof {return pono::smvparser::make_tok_typeof(_encoder.loc);}

unsigned   {return pono::smvparser::make_tok_unsigned(_encoder.loc);}
signed   {return pono::smvparser::make_tok_signed(_encoder.loc);}
boolean   {return pono::smvparser::make_bool_type(_encoder.loc);}
real  {return pono::smvparser::make_real_type(_encoder.loc);}
integer   {return pono::smvparser::make_integer_type(_encoder.loc);}
word   {return pono::smvparser::make_tok_word(_encoder.loc);}
clock {return pono::smvparser::make_time_type(_encoder.loc);}
in {return pono::smvparser::make_OP_IN(_encoder.loc);}

[fF]' {return pono::smvparser::make_fraction_prefix(yytext,_encoder.loc);}
[0-9]+           {return pono::smvparser::make_integer_val(yytext,_encoder.loc);}
-[0-9]+           {return pono::smvparser::make_neg_integer_val(yytext,_encoder.loc);}
"->"  {return pono::smvparser::make_OP_IMPLY(_encoder.loc);}
"<->"  {return pono::smvparser::make_OP_BI(_encoder.loc);}
"!"  {return pono::smvparser::make_OP_NOT(_encoder.loc);}
"::"  {return pono::smvparser::make_OP_CON(_encoder.loc);}
"-"  {return pono::smvparser::make_OP_MINUS(_encoder.loc);}
"+"  {return pono::smvparser::make_OP_PLUS(_encoder.loc);}
"*"  {return pono::smvparser::make_OP_MUL(_encoder.loc);}
"/"  {return pono::smvparser::make_OP_DIV(_encoder.loc); }
"mod"  {return pono::smvparser::make_OP_MOD(_encoder.loc);}
">>"  {return pono::smvparser::make_OP_SHIFTR(_encoder.loc);}
"<<"  {return pono::smvparser::make_OP_SHIFTL(_encoder.loc);}
"union"  {return pono::smvparser::make_UNION(_encoder.loc);}
"="  {return pono::smvparser::make_OP_EQ(_encoder.loc);}
"!="  {return pono::smvparser::make_OP_NEQ(_encoder.loc);}
">"  {return pono::smvparser::make_OP_GT(_encoder.loc);}
"<"  {return pono::smvparser::make_OP_LT(_encoder.loc);}
"<="  {return pono::smvparser::make_OP_LTE(_encoder.loc);}
">="  {return pono::smvparser::make_OP_GTE(_encoder.loc);}
"&"  {return pono::smvparser::make_OP_AND(_encoder.loc);}
"|"  {return pono::smvparser::make_OP_OR(_encoder.loc);}
"xor"  {return pono::smvparser::make_OP_XOR(_encoder.loc);}
"xnor"  {return pono::smvparser::make_OP_XNOR(_encoder.loc);}

\.    {return pono::smvparser::make_DOT(_encoder.loc);};
".." {return pono::smvparser::make_TO(_encoder.loc);}
"(" {return pono::smvparser::make_LPARE(_encoder.loc);}
")" {return pono::smvparser::make_RPARE(_encoder.loc);}
"[" {return pono::smvparser::make_LBRA(_encoder.loc);}
"]" {return pono::smvparser::make_RBRA(_encoder.loc);}
"{" {return pono::smvparser::make_LBRACE(_encoder.loc);}
"}" {return pono::smvparser::make_RBRACE(_encoder.loc);}
";" {return pono::smvparser::make_SEMICOLON(_encoder.loc);}
":" {return pono::smvparser::make_COLON(_encoder.loc);}
"," {return pono::smvparser::make_COMMA(_encoder.loc);}
":=" {return pono::smvparser::make_ASSIGNSYM(_encoder.loc);}
"--".*  {}
"?"  {return pono::smvparser::make_IF_ELSE(_encoder.loc);}
"_" {return pono::smvparser::make_UNDER(_encoder.loc);}


0(b|d|o|h)/[0-9]+_[0-9a-fA-F]+ {BEGIN(WORD_INIT); return pono::smvparser::make_word_index1(yytext,_encoder.loc);}
0(u|s)(b|d|o|h)/[0-9]+_[0-9a-fA-F]+ {BEGIN(WORD_INIT);  return pono::smvparser::make_word_index2(yytext,_encoder.loc);}
<WORD_INIT>[0-9]+/_[0-9a-fA-F]+ {BEGIN(WORD_BASE); return pono::smvparser::make_integer_val(yytext,_encoder.loc);}
<WORD_BASE>_/[0-9a-fA-F]+ {BEGIN(WORD_VALUE); return pono::smvparser::make_UNDER(_encoder.loc);}
<WORD_VALUE>[0-9a-fA-F]+ {BEGIN(INITIAL); return pono::smvparser::make_integer_val(yytext,_encoder.loc);}

[a-zA-Z_][A-Za-z0-9_\$#-]*  {return pono::smvparser::make_tok_name(yytext,_encoder.loc);}
\"(\\.|[^\"\\])*\"    {
                        for (const char * c=yytext; *c; ++c)
                        {
                          if(*c == '\n')
                          {
                            _encoder.loc.lines(1);
                          }
                        }
                        _encoder.loc.step();
                        return pono::smvparser::make_tok_name(yytext, _encoder.loc);
                      }
\n   {_encoder.loc.lines(yyleng);}
.    {}
<<EOF>>    return pono::smvparser::make_END(_encoder.loc);

%%
