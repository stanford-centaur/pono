%{
  #include "smvparser.h"
  #include"node.h"
  #include "smv_encoder.h"
  using namespace std;
  int line_num = 1;

%}
%option yylineno
%option noyywrap
%x COMMENT
%x WORD_INIT WORD_BASE WORD_VALUE
%%

[ \t]     ;
MODULE {return MODULE;}
main  {return tok_main;}
IVAR   {return IVAR;}
VAR  {return VAR;}
INVAR   {return INVAR;}
FROZENVAR   {return FROZENVAR;}
INVARSPEC   {return INVARSPEC;}
INIT   {return INIT;}
TRANS   {return TRANS;}
READ   {return READ;}
WRITE   {return WRITE;}
ASSIGN   {return ASSIGN;}
CONSTARRAY   {return CONSTARRAY;}
CONSTANTS   {return CONSTANTS;}
FUN   {return FUN;}
DEFINE {return DEFINE;}

next  {return tok_next;}
init {return tok_init;}

TRUE   {return TOK_TRUE;}
FALSE   {return TOK_FALSE;}

"signed word" {return signed_word;}
"unsigned word" {return unsigned_word;}
"array word"   {return arrayword;}
"array integer"   {return arrayinteger;}
"array" {return tok_array;}

pi {return pi;}
abs {return ABS;}
max {return MAX;}
min {return MIN;}
sin {return SIN;}
cos {return COS;}
exp {return EXP;}
tan {return TAN;}
ln {return ln;}
of  {return of;}
word1 {return word1;}
bool {return tok_bool;}
toint   {return tok_toint;}
count {return tok_count;}
swconst {return swconst;}
uwconst {return uwconst;}
sizeof {return tok_sizeof;}
floor {return tok_floor;}
extend {return extend;}
resize {return resize;}
typeof {return tok_typeof;}

unsigned   {return tok_unsigned;}
signed   {return tok_signed;}
boolean   {return bool_type;}
real  {return real_type;}
integer   {return integer_type;}
clock   {return clock_type;}
word   {return tok_word;}
set {return tok_set;}
in {return in;}
time {return time_type;}

"->"  {return OP_IMPLY;}
"<->"  {return OP_BI;}
"!"  {return '!';}
"::"  {return OP_CON;}
"-"  {return '-';}
"+"  {return '+';}
"*"  {return '*';}
"/"  {return '/'; }
"mod"  {return OP_MOD;}
">>"  {return OP_SHIFTR;}
"<<"  {return OP_SHIFTL;}
"union"  {return UNION;}
"="  {return '=';}
"!="  {return OP_NEQ;}
">"  {return '>';}
"<"  {return '<';}
"<="  {return OP_LTE;}
">="  {return OP_GTE;}
"&"  {return '&';}
"|"  {return '|';}
"xor"  {return OP_XOR;}
"xnor"  {return OP_XNOR;}

".." {return TO;}
"(" {return '(';}
")" {return ')';}
"[" {return '[';}
"]" {return ']';}
"{" {return '{';}
"}" {return '}';}
";" {return ';';}
":" {return ':';}
"," {return ',';}
":=" {return ASSIGNSYM;}
"--".*  {}
"?"  {return IF_ELSE;}
"_" {return '_';}

0(b|d|o|h)/[0-9]+_[0-9]+ {BEGIN(WORD_INIT); yylval.str = strdup(yytext); return word_index1;}
0(u|s)(b|d|o|h)/[0-9]+_[0-9]+ {BEGIN(WORD_INIT); yylval.str = strdup(yytext); return word_index2;}
<WORD_INIT>[0-9]+/_[0-9]+ {BEGIN(WORD_BASE);yylval.num = atoi(yytext);return integer_val;}
<WORD_BASE>_/[0-9]+ {BEGIN(WORD_VALUE);return '_';}
<WORD_VALUE>[0-9]+ {BEGIN(INITIAL); yylval.num = atoi(yytext); return integer_val;}

[0-9]+           {yylval.num = atoi(yytext);return integer_val;}
[a-zA-Z][A-Za-z0-9_\$#-]*  {yylval.str = strdup(yytext); return tok_name;}
\n   {line_num++; return ENDL;}
.                {return '.';};

%%
