/* A Bison parser, made by GNU Bison 3.4.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2019 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Undocumented macros, especially those whose name start with YY_,
   are private implementation details.  Do not rely on them.  */

#ifndef YY_YY_SMVPARSER_TAB_H_INCLUDED
# define YY_YY_SMVPARSER_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif
/* "%code requires" blocks.  */
#line 17 "smvparser.y"

  #include "node.h"
  class smvEncoder;

#line 53 "smvparser.tab.h"

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    MODULE = 258,
    tok_main = 259,
    IVAR = 260,
    INVAR = 261,
    VAR = 262,
    FROZENVAR = 263,
    INVARSPEC = 264,
    INIT = 265,
    TRANS = 266,
    READ = 267,
    WRITE = 268,
    ASSIGN = 269,
    CONSTARRAY = 270,
    CONSTANTS = 271,
    FUN = 272,
    DEFINE = 273,
    INPUT = 274,
    OUTPUT = 275,
    tok_next = 276,
    tok_init = 277,
    signed_word = 278,
    unsigned_word = 279,
    arrayword = 280,
    arrayinteger = 281,
    tok_array = 282,
    pi = 283,
    ABS = 284,
    MAX = 285,
    MIN = 286,
    SIN = 287,
    COS = 288,
    EXP = 289,
    TAN = 290,
    ln = 291,
    of = 292,
    word1 = 293,
    tok_bool = 294,
    tok_toint = 295,
    tok_count = 296,
    swconst = 297,
    uwconst = 298,
    tok_sizeof = 299,
    tok_floor = 300,
    extend = 301,
    resize = 302,
    tok_typeof = 303,
    tok_unsigned = 304,
    tok_signed = 305,
    tok_word = 306,
    tok_set = 307,
    in = 308,
    time_type = 309,
    TO = 310,
    ASSIGNSYM = 311,
    IF_ELSE = 312,
    ENDL = 313,
    integer_val = 314,
    real_val = 315,
    bool_type = 316,
    integer_type = 317,
    real_type = 318,
    clock_type = 319,
    set_tok = 320,
    array_tok = 321,
    word_index1 = 322,
    word_index2 = 323,
    tok_name = 324,
    TOK_TRUE = 325,
    TOK_FALSE = 326,
    OP_CON = 327,
    UMINUS = 328,
    OP_MOD = 329,
    OP_SHIFTR = 330,
    OP_SHIFTL = 331,
    UNION = 332,
    OP_EQ = 333,
    OP_NEQ = 334,
    OP_LTE = 335,
    OP_GTE = 336,
    OP_OR = 337,
    OP_XOR = 338,
    OP_XNOR = 339,
    OP_BI = 340,
    OP_IMPLY = 341
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 28 "smvparser.y"

  node *n;
  char *str;
  int num;
  bool bl;

#line 158 "smvparser.tab.h"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif

/* Location type.  */
#if ! defined YYLTYPE && ! defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE YYLTYPE;
struct YYLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
};
# define YYLTYPE_IS_DECLARED 1
# define YYLTYPE_IS_TRIVIAL 1
#endif


extern YYSTYPE yylval;
extern YYLTYPE yylloc;
int yyparse (smvEncoder &enc);

#endif /* !YY_YY_SMVPARSER_TAB_H_INCLUDED  */
