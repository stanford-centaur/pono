/* A Bison parser, made by GNU Bison 3.4.2.  */

/* Bison implementation for Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Undocumented macros, especially those whose name start with YY_,
   are private implementation details.  Do not rely on them.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.4.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* First part of user prologue.  */
#line 1 "smvparser.y"

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


#line 86 "smvparser.tab.c"

# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 1
#endif

/* Use api.header.include to #include this header
   instead of duplicating it here.  */
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

#line 125 "smvparser.tab.c"

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

#line 230 "smvparser.tab.c"

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



#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && ! defined __ICC && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL \
             && defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
  YYLTYPE yyls_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE) + sizeof (YYLTYPE)) \
      + 2 * YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  82
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1917

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  108
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  37
/* YYNRULES -- Number of rules.  */
#define YYNRULES  158
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  360

#define YYUNDEFTOK  2
#define YYMAXUTOK   341

/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                                \
  ((unsigned) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    72,     2,     2,     2,     2,    89,     2,
      97,    98,    75,    78,    96,    79,   107,    76,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    99,    95,
      85,   102,    86,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   103,     2,   104,     2,   100,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   105,   101,   106,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    73,    74,    77,
      80,    81,    82,    83,    84,    87,    88,    90,    91,    92,
      93,    94
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    70,    70,    71,    72,    73,    74,    75,    76,    77,
      78,    79,    81,    82,    83,    84,    85,    86,    87,    88,
      89,    90,    92,    93,    94,    95,    96,    97,    98,    99,
     100,   101,   105,   108,   109,   111,   115,   121,   123,   124,
     126,   127,   129,   130,   132,   136,   142,   150,   151,   152,
     153,   157,   166,   167,   168,   169,   172,   181,   182,   183,
     184,   187,   198,   202,   207,   212,   219,   225,   231,   250,
     254,   258,   259,   260,   263,   264,   265,   267,   287,   310,
     313,   321,   324,   327,   330,   333,   336,   349,   352,   353,
     359,   365,   371,   377,   383,   389,   395,   401,   407,   413,
     419,   425,   428,   434,   440,   446,   452,   458,   464,   470,
     471,   472,   478,   479,   480,   481,   482,   483,   484,   485,
     486,   487,   488,   489,   490,   491,   492,   493,   494,   495,
     496,   503,   509,   510,   511,   512,   516,   522,   523,   525,
     526,   528,   530,   533,   536,   542,   547,   552,   556,   559,
     562,   565,   567,   571,   575,   580,   586,   589,   592
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "MODULE", "tok_main", "IVAR", "INVAR",
  "VAR", "FROZENVAR", "INVARSPEC", "INIT", "TRANS", "READ", "WRITE",
  "ASSIGN", "CONSTARRAY", "CONSTANTS", "FUN", "DEFINE", "INPUT", "OUTPUT",
  "tok_next", "tok_init", "signed_word", "unsigned_word", "arrayword",
  "arrayinteger", "tok_array", "pi", "ABS", "MAX", "MIN", "SIN", "COS",
  "EXP", "TAN", "ln", "of", "word1", "tok_bool", "tok_toint", "tok_count",
  "swconst", "uwconst", "tok_sizeof", "tok_floor", "extend", "resize",
  "tok_typeof", "tok_unsigned", "tok_signed", "tok_word", "tok_set", "in",
  "time_type", "TO", "ASSIGNSYM", "IF_ELSE", "ENDL", "integer_val",
  "real_val", "bool_type", "integer_type", "real_type", "clock_type",
  "set_tok", "array_tok", "word_index1", "word_index2", "tok_name",
  "TOK_TRUE", "TOK_FALSE", "'!'", "OP_CON", "UMINUS", "'*'", "'/'",
  "OP_MOD", "'+'", "'-'", "OP_SHIFTR", "OP_SHIFTL", "UNION", "OP_EQ",
  "OP_NEQ", "'<'", "'>'", "OP_LTE", "OP_GTE", "'&'", "OP_OR", "OP_XOR",
  "OP_XNOR", "OP_BI", "OP_IMPLY", "';'", "','", "'('", "')'", "':'", "'_'",
  "'|'", "'='", "'['", "']'", "'{'", "'}'", "'.'", "$accept", "header",
  "module_header", "define_decl", "define_body", "constants_decl",
  "constants_body", "assign_decl", "assign_list", "assign_test",
  "ivar_test", "ivar_list", "var_test", "var_list", "frozenvar_test",
  "frozenvar_list", "init_test", "trans_test", "invarspec_test",
  "next_formula", "constant", "word_value", "boolean_constant",
  "integer_constant", "real_constant", "range_constant", "clock_constant",
  "basic_expr", "next_expr", "basic_expr_list", "set_body_expr",
  "symbolic_constant", "complex_identifier", "type_identifier",
  "word_type", "array_type", "sizev", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,    33,   327,   328,    42,    47,   329,    43,    45,
     330,   331,   332,   333,   334,    60,    62,   335,   336,    38,
     337,   338,   339,   340,   341,    59,    44,    40,    41,    58,
      95,   124,    61,    91,    93,   123,   125,    46
};
# endif

#define YYPACT_NINF -129

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-129)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     190,     7,     0,    14,    28,   257,   257,   257,   213,   -18,
     -34,  -129,    18,   204,  -129,  -129,  -129,   -18,   -18,    -2,
    -129,  -129,  -129,  -129,   -18,   -18,   -18,     9,    12,    29,
      65,    10,    10,    69,    75,    76,    77,    81,    82,    83,
      84,    86,    87,    89,    90,  -129,    41,  -129,    71,    96,
    -129,  -129,  -129,   257,   257,   257,   257,  -129,  -129,  -129,
    -129,  -129,  -129,  -129,   400,  -129,  -129,    85,   443,   101,
    1666,   105,   106,    16,    16,   115,   -48,   136,    85,   -18,
     -18,   -47,  -129,  -129,  -129,  -129,  -129,   -18,   -18,    -2,
    -129,  -129,  -129,  -129,  -129,  -129,  -129,   -18,   -18,    -2,
    -129,  -129,  -129,   158,   -26,   175,    -9,    -2,   178,    -8,
     179,   183,   184,   257,   257,    15,   257,   185,   146,   149,
     257,   257,   257,   257,   257,   257,   257,   257,   257,   257,
     257,   257,   188,   150,   151,  1666,  1792,   486,  1666,   -52,
     257,   257,   257,   257,   257,   257,   257,   257,   257,   257,
     152,   257,   257,   257,   257,   257,   257,   257,   257,   257,
     257,   194,   257,   257,   327,    -7,   195,   196,   -18,   -18,
      16,   160,   198,   257,   -18,   -18,   -46,   257,  -129,   359,
    -129,   359,  -129,   359,  -129,  -129,  -129,   528,   570,    10,
     222,   163,   612,   157,   257,   257,   654,   696,   738,  1666,
      -1,   780,   822,   864,   906,   948,   990,  1032,  1074,  -129,
     214,   215,  -129,   257,  -129,  1666,  1116,  1708,  1750,  1750,
    1750,  1792,  1792,  1814,  1814,   257,   137,   137,   137,   137,
     137,   -14,   164,   164,    74,    74,   218,  1666,  1666,   -40,
    1159,  -129,  -129,   219,   221,   -76,   -42,   225,  -129,  1666,
     -65,   257,  1202,    10,    10,    10,   247,    10,   230,  -129,
    -129,  -129,  -129,   249,   192,  -129,  -129,   199,   210,   257,
     257,   252,   359,   -18,  -129,  -129,  1245,  1287,  -129,  -129,
    -129,   257,  -129,  -129,  -129,  -129,  -129,  -129,  -129,  -129,
    -129,  -129,  -129,  1666,   257,  1329,  -129,  -129,   257,  -129,
    -129,   234,   235,  -129,   250,  1371,   251,  -129,  -129,   256,
     359,  -129,   253,   359,   255,   260,   261,  1414,  1456,   359,
     224,   -41,  -129,  -129,  1666,  1666,  -129,   358,   289,   257,
    -129,   263,  -129,   359,  -129,  -129,  -129,  -129,  -129,  -129,
    -129,   257,   226,   257,   236,  -129,  -129,  1666,  -129,  -129,
    1498,   257,  1540,   257,  -129,  1582,  -129,  1624,  -129,  -129
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     2,     0,     0,     3,     4,     5,     6,     7,    10,
       8,     9,    11,    32,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    84,    81,    82,     0,     0,
     142,    79,    80,     0,     0,     0,     0,    85,    73,    70,
      71,    72,    76,    74,     0,   135,    75,    86,     0,     0,
      69,     0,     0,     0,    40,     0,     0,    37,    38,     0,
      33,     0,     1,    12,    13,    14,    15,    16,    17,    20,
      18,    19,    21,    22,    23,    24,    25,    26,    27,    29,
      28,    30,    31,    50,     0,    55,     0,    60,     0,     0,
      48,    54,    59,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    88,   101,     0,   139,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    68,     0,     0,     0,     0,     0,     0,     0,     0,
      41,     0,     0,     0,     0,    34,     0,     0,    49,     0,
      53,     0,    58,     0,    47,    52,    57,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   137,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    83,
       0,     0,    87,     0,   127,   128,     0,   111,   104,   105,
     106,   102,   103,   107,   108,     0,    96,    97,    98,    99,
     100,    89,    91,    92,    94,    93,    67,    90,    95,    81,
       0,   144,   143,    63,    65,     0,     0,     0,    42,    44,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   147,
     146,   145,   148,     0,     0,   150,   149,     0,     0,     0,
       0,     0,     0,     0,   136,   158,     0,     0,   112,   113,
     114,     0,   115,   116,   117,   120,   121,   122,   123,   119,
     118,    77,    78,   140,     0,     0,    66,   109,     0,    62,
      64,     0,     0,    43,     0,     0,     0,   152,   153,     0,
       0,   154,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   124,   125,   138,   129,   126,     0,     0,     0,
      39,     0,    35,     0,   156,   151,   157,    51,    56,    61,
     131,     0,     0,     0,     0,   110,    46,    45,    36,   155,
       0,     0,     0,     0,   130,     0,   134,     0,   133,   132
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -129,  -129,  -129,     4,   254,    58,  -129,    88,   241,   -67,
      92,   291,    99,   305,    33,   308,   116,   120,   148,  -129,
    -129,  -129,  -129,  -129,  -129,  -129,  -129,    -6,     3,  -129,
    -129,  -129,    -5,  -128,  -129,  -129,   -30
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    12,    13,    14,    80,    15,    77,    16,    74,    75,
      17,   103,    18,   105,   107,   108,    20,    21,    22,    69,
      57,    58,    59,    60,    61,    62,    63,    64,    65,   200,
     139,    66,    67,   264,   265,   266,   118
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint16 yytable[] =
{
      68,    70,   119,    76,    78,    81,     4,   171,   173,   177,
     251,    23,   104,   106,   109,   132,    84,    94,    82,   104,
     106,   109,   301,     2,    79,     3,     4,     5,     6,     7,
     304,   165,     8,    19,     9,    50,    10,    71,    72,   140,
     189,   190,   165,   141,   213,    89,    99,   135,   136,   137,
     138,    50,   241,   267,   214,   268,   302,   344,    24,   165,
     165,   165,   242,   191,   297,   165,   165,    50,    76,    76,
      85,    95,    25,   179,    81,   176,    83,   157,   158,   159,
     160,   165,   104,   106,   109,    50,    26,   162,   163,   164,
     181,   183,   104,   106,   109,   281,   132,   282,   165,   165,
      86,    96,   109,   171,    87,    97,   113,   187,   188,   114,
     192,    88,    98,   117,   196,   197,   198,   199,   201,   202,
     203,   204,   205,   206,   207,   208,   115,   140,    90,   100,
     133,   141,    91,   101,   215,   216,   217,   218,   219,   220,
     221,   222,   223,   224,   320,   226,   227,   228,   229,   230,
     231,   232,   233,   234,   235,   134,   237,   238,   240,   271,
      92,   102,   116,   245,   246,    76,   120,   249,   160,   250,
     176,   252,   121,   122,   123,   162,   163,   164,   124,   125,
     126,   127,   334,   128,   129,   336,   130,   131,   276,   277,
     140,   342,   165,     1,   141,     2,   167,     3,     4,     5,
       6,     7,   168,   169,     8,   349,     9,   293,    10,     2,
     172,     3,     4,     5,     6,     7,   178,   140,     8,   295,
       9,   141,    10,   307,   308,   309,   156,   311,   157,   158,
     159,   160,   174,   180,    71,    72,   182,   184,   162,   163,
     164,   185,   186,   194,   193,   305,   195,   209,    11,   225,
     210,   211,   236,   243,   244,   247,   248,   159,   160,   272,
     273,   275,    93,   317,   318,   162,   163,   164,   321,    27,
      28,    73,    29,   291,   292,   324,   296,   299,    30,   300,
      31,    32,    50,   303,   310,   312,   313,   314,   325,   319,
     328,   329,   327,   333,   315,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,   316,    43,    44,   330,   332,
      30,    45,   335,   337,   170,   110,    46,    47,   338,   339,
     343,   348,   351,   347,    48,    49,    50,    51,    52,    53,
     111,   346,   353,   175,   112,   350,    54,   352,     0,    27,
      28,     0,    29,     0,     0,   355,     0,   357,    30,     0,
      31,    32,     0,     0,    55,     0,     0,     0,     0,     0,
       0,     0,    56,     0,     0,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,     0,    43,    44,     0,     0,
       0,    45,   253,   254,   255,   256,   239,    47,     0,     0,
       0,     0,     0,     0,    48,    49,    50,    51,    52,    53,
       0,     0,     0,     0,     0,     0,    54,     0,     0,     0,
     257,   140,     0,     0,     0,   141,     0,     0,   258,     0,
     259,   260,   261,   262,    55,   263,     0,     0,     0,     0,
       0,   142,    56,   143,   144,   145,   146,   147,   148,   149,
     150,     0,   151,   152,   153,   154,   155,   156,     0,   157,
     158,   159,   160,   140,     0,     0,     0,   141,     0,   162,
     163,   164,   345,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   142,     0,   143,   144,   145,   146,   147,
     148,   149,   150,     0,   151,   152,   153,   154,   155,   156,
       0,   157,   158,   159,   160,   161,   140,     0,     0,     0,
     141,   162,   163,   164,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   142,     0,   143,   144,
     145,   146,   147,   148,   149,   150,     0,   151,   152,   153,
     154,   155,   156,     0,   157,   158,   159,   160,   166,   140,
       0,     0,     0,   141,   162,   163,   164,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   142,
       0,   143,   144,   145,   146,   147,   148,   149,   150,     0,
     151,   152,   153,   154,   155,   156,     0,   157,   158,   159,
     160,   140,     0,     0,   212,   141,     0,   162,   163,   164,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   142,     0,   143,   144,   145,   146,   147,   148,   149,
     150,     0,   151,   152,   153,   154,   155,   156,     0,   157,
     158,   159,   160,   140,   269,     0,     0,   141,     0,   162,
     163,   164,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   142,     0,   143,   144,   145,   146,   147,
     148,   149,   150,     0,   151,   152,   153,   154,   155,   156,
       0,   157,   158,   159,   160,   140,   270,     0,     0,   141,
       0,   162,   163,   164,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   142,     0,   143,   144,   145,
     146,   147,   148,   149,   150,     0,   151,   152,   153,   154,
     155,   156,     0,   157,   158,   159,   160,   140,     0,     0,
     274,   141,     0,   162,   163,   164,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   142,     0,   143,
     144,   145,   146,   147,   148,   149,   150,     0,   151,   152,
     153,   154,   155,   156,     0,   157,   158,   159,   160,   140,
       0,     0,   278,   141,     0,   162,   163,   164,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   142,
       0,   143,   144,   145,   146,   147,   148,   149,   150,     0,
     151,   152,   153,   154,   155,   156,     0,   157,   158,   159,
     160,   140,     0,     0,   279,   141,     0,   162,   163,   164,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   142,     0,   143,   144,   145,   146,   147,   148,   149,
     150,     0,   151,   152,   153,   154,   155,   156,     0,   157,
     158,   159,   160,   140,     0,     0,   280,   141,     0,   162,
     163,   164,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   142,     0,   143,   144,   145,   146,   147,
     148,   149,   150,     0,   151,   152,   153,   154,   155,   156,
       0,   157,   158,   159,   160,   140,     0,     0,   283,   141,
       0,   162,   163,   164,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   142,     0,   143,   144,   145,
     146,   147,   148,   149,   150,     0,   151,   152,   153,   154,
     155,   156,     0,   157,   158,   159,   160,   140,     0,     0,
     284,   141,     0,   162,   163,   164,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   142,     0,   143,
     144,   145,   146,   147,   148,   149,   150,     0,   151,   152,
     153,   154,   155,   156,     0,   157,   158,   159,   160,   140,
       0,     0,   285,   141,     0,   162,   163,   164,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   142,
       0,   143,   144,   145,   146,   147,   148,   149,   150,     0,
     151,   152,   153,   154,   155,   156,     0,   157,   158,   159,
     160,   140,     0,     0,   286,   141,     0,   162,   163,   164,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   142,     0,   143,   144,   145,   146,   147,   148,   149,
     150,     0,   151,   152,   153,   154,   155,   156,     0,   157,
     158,   159,   160,   140,     0,     0,   287,   141,     0,   162,
     163,   164,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   142,     0,   143,   144,   145,   146,   147,
     148,   149,   150,     0,   151,   152,   153,   154,   155,   156,
       0,   157,   158,   159,   160,   140,     0,     0,   288,   141,
       0,   162,   163,   164,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   142,     0,   143,   144,   145,
     146,   147,   148,   149,   150,     0,   151,   152,   153,   154,
     155,   156,     0,   157,   158,   159,   160,   140,     0,     0,
     289,   141,     0,   162,   163,   164,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   142,     0,   143,
     144,   145,   146,   147,   148,   149,   150,     0,   151,   152,
     153,   154,   155,   156,     0,   157,   158,   159,   160,   140,
       0,     0,   290,   141,     0,   162,   163,   164,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   142,
       0,   143,   144,   145,   146,   147,   148,   149,   150,     0,
     151,   152,   153,   154,   155,   156,     0,   157,   158,   159,
     160,     0,   140,     0,     0,   294,   141,   162,   163,   164,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   142,     0,   143,   144,   145,   146,   147,   148,
     149,   150,     0,   151,   152,   153,   154,   155,   156,     0,
     157,   158,   159,   160,     0,   140,     0,     0,   298,   141,
     162,   163,   164,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   142,     0,   143,   144,   145,
     146,   147,   148,   149,   150,     0,   151,   152,   153,   154,
     155,   156,     0,   157,   158,   159,   160,   306,   140,     0,
       0,     0,   141,   162,   163,   164,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   142,     0,
     143,   144,   145,   146,   147,   148,   149,   150,     0,   151,
     152,   153,   154,   155,   156,     0,   157,   158,   159,   160,
     140,     0,     0,   322,   141,     0,   162,   163,   164,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     142,     0,   143,   144,   145,   146,   147,   148,   149,   150,
       0,   151,   152,   153,   154,   155,   156,     0,   157,   158,
     159,   160,   140,     0,     0,   323,   141,     0,   162,   163,
     164,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   142,     0,   143,   144,   145,   146,   147,   148,
     149,   150,     0,   151,   152,   153,   154,   155,   156,     0,
     157,   158,   159,   160,   140,     0,     0,   326,   141,     0,
     162,   163,   164,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   142,     0,   143,   144,   145,   146,
     147,   148,   149,   150,     0,   151,   152,   153,   154,   155,
     156,     0,   157,   158,   159,   160,   331,   140,     0,     0,
       0,   141,   162,   163,   164,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   142,     0,   143,
     144,   145,   146,   147,   148,   149,   150,     0,   151,   152,
     153,   154,   155,   156,     0,   157,   158,   159,   160,   140,
       0,     0,   340,   141,     0,   162,   163,   164,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   142,
       0,   143,   144,   145,   146,   147,   148,   149,   150,     0,
     151,   152,   153,   154,   155,   156,     0,   157,   158,   159,
     160,   140,   341,     0,     0,   141,     0,   162,   163,   164,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   142,     0,   143,   144,   145,   146,   147,   148,   149,
     150,     0,   151,   152,   153,   154,   155,   156,     0,   157,
     158,   159,   160,   140,     0,     0,   354,   141,     0,   162,
     163,   164,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   142,     0,   143,   144,   145,   146,   147,
     148,   149,   150,     0,   151,   152,   153,   154,   155,   156,
       0,   157,   158,   159,   160,   140,     0,     0,   356,   141,
       0,   162,   163,   164,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   142,     0,   143,   144,   145,
     146,   147,   148,   149,   150,     0,   151,   152,   153,   154,
     155,   156,     0,   157,   158,   159,   160,   140,     0,     0,
     358,   141,     0,   162,   163,   164,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   142,     0,   143,
     144,   145,   146,   147,   148,   149,   150,     0,   151,   152,
     153,   154,   155,   156,     0,   157,   158,   159,   160,   140,
       0,     0,   359,   141,     0,   162,   163,   164,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   142,
       0,   143,   144,   145,   146,   147,   148,   149,   150,     0,
     151,   152,   153,   154,   155,   156,     0,   157,   158,   159,
     160,   140,     0,     0,     0,   141,     0,   162,   163,   164,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   143,   144,   145,   146,   147,   148,   149,
     150,     0,   151,   152,   153,   154,   155,   156,     0,   157,
     158,   159,   160,   140,     0,     0,     0,   141,     0,   162,
     163,   164,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   146,   147,
     148,   149,   150,     0,   151,   152,   153,   154,   155,   156,
       0,   157,   158,   159,   160,   140,     0,     0,     0,   141,
       0,   162,   163,   164,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   140,     0,     0,
       0,   141,   148,   149,   150,     0,   151,   152,   153,   154,
     155,   156,     0,   157,   158,   159,   160,     0,     0,     0,
       0,     0,     0,   162,   163,   164,   150,     0,   151,   152,
     153,   154,   155,   156,     0,   157,   158,   159,   160,     0,
       0,     0,     0,     0,     0,   162,   163,   164
};

static const yytype_int16 yycheck[] =
{
       6,     7,    32,     8,     9,    10,     8,    74,    56,    56,
      56,     4,    17,    18,    19,    55,    12,    13,     0,    24,
      25,    26,    98,     5,    58,     7,     8,     9,    10,    11,
      95,   107,    14,     0,    16,    69,    18,    21,    22,    53,
      25,    26,   107,    57,    96,    12,    13,    53,    54,    55,
      56,    69,    59,   181,   106,   183,    98,    98,    58,   107,
     107,   107,    69,    48,   104,   107,   107,    69,    73,    74,
      12,    13,    58,    99,    79,    80,    58,    91,    92,    93,
      94,   107,    87,    88,    89,    69,    58,   101,   102,   103,
      99,    99,    97,    98,    99,    96,    55,    98,   107,   107,
      12,    13,   107,   170,    12,    13,    97,   113,   114,    97,
     116,    12,    13,   103,   120,   121,   122,   123,   124,   125,
     126,   127,   128,   129,   130,   131,    97,    53,    12,    13,
      59,    57,    12,    13,   140,   141,   142,   143,   144,   145,
     146,   147,   148,   149,   272,   151,   152,   153,   154,   155,
     156,   157,   158,   159,   160,    59,   162,   163,   164,   189,
      12,    13,    97,   168,   169,   170,    97,   173,    94,   174,
     175,   177,    97,    97,    97,   101,   102,   103,    97,    97,
      97,    97,   310,    97,    97,   313,    97,    97,   194,   195,
      53,   319,   107,     3,    57,     5,    95,     7,     8,     9,
      10,    11,    97,    97,    14,   333,    16,   213,    18,     5,
      95,     7,     8,     9,    10,    11,    58,    53,    14,   225,
      16,    57,    18,   253,   254,   255,    89,   257,    91,    92,
      93,    94,    96,    58,    21,    22,    58,    58,   101,   102,
     103,    58,    58,    97,    59,   251,    97,    59,    58,    97,
     100,   100,    58,    58,    58,    95,    58,    93,    94,    37,
      97,   104,    58,   269,   270,   101,   102,   103,   273,    12,
      13,    58,    15,    59,    59,   281,    58,    58,    21,    58,
      23,    24,    69,    58,    37,    55,    37,    95,   294,    37,
      56,    56,   298,    37,    95,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    95,    49,    50,    58,    58,
      21,    54,    59,    58,    73,    24,    59,    60,    58,    58,
      96,    58,    96,   329,    67,    68,    69,    70,    71,    72,
      25,   328,    96,    79,    26,   341,    79,   343,    -1,    12,
      13,    -1,    15,    -1,    -1,   351,    -1,   353,    21,    -1,
      23,    24,    -1,    -1,    97,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   105,    -1,    -1,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    -1,    49,    50,    -1,    -1,
      -1,    54,    23,    24,    25,    26,    59,    60,    -1,    -1,
      -1,    -1,    -1,    -1,    67,    68,    69,    70,    71,    72,
      -1,    -1,    -1,    -1,    -1,    -1,    79,    -1,    -1,    -1,
      51,    53,    -1,    -1,    -1,    57,    -1,    -1,    59,    -1,
      61,    62,    63,    64,    97,    66,    -1,    -1,    -1,    -1,
      -1,    73,   105,    75,    76,    77,    78,    79,    80,    81,
      82,    -1,    84,    85,    86,    87,    88,    89,    -1,    91,
      92,    93,    94,    53,    -1,    -1,    -1,    57,    -1,   101,
     102,   103,   104,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    73,    -1,    75,    76,    77,    78,    79,
      80,    81,    82,    -1,    84,    85,    86,    87,    88,    89,
      -1,    91,    92,    93,    94,    95,    53,    -1,    -1,    -1,
      57,   101,   102,   103,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    73,    -1,    75,    76,
      77,    78,    79,    80,    81,    82,    -1,    84,    85,    86,
      87,    88,    89,    -1,    91,    92,    93,    94,    95,    53,
      -1,    -1,    -1,    57,   101,   102,   103,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,
      -1,    75,    76,    77,    78,    79,    80,    81,    82,    -1,
      84,    85,    86,    87,    88,    89,    -1,    91,    92,    93,
      94,    53,    -1,    -1,    98,    57,    -1,   101,   102,   103,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    73,    -1,    75,    76,    77,    78,    79,    80,    81,
      82,    -1,    84,    85,    86,    87,    88,    89,    -1,    91,
      92,    93,    94,    53,    96,    -1,    -1,    57,    -1,   101,
     102,   103,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    73,    -1,    75,    76,    77,    78,    79,
      80,    81,    82,    -1,    84,    85,    86,    87,    88,    89,
      -1,    91,    92,    93,    94,    53,    96,    -1,    -1,    57,
      -1,   101,   102,   103,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    73,    -1,    75,    76,    77,
      78,    79,    80,    81,    82,    -1,    84,    85,    86,    87,
      88,    89,    -1,    91,    92,    93,    94,    53,    -1,    -1,
      98,    57,    -1,   101,   102,   103,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,    -1,    75,
      76,    77,    78,    79,    80,    81,    82,    -1,    84,    85,
      86,    87,    88,    89,    -1,    91,    92,    93,    94,    53,
      -1,    -1,    98,    57,    -1,   101,   102,   103,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,
      -1,    75,    76,    77,    78,    79,    80,    81,    82,    -1,
      84,    85,    86,    87,    88,    89,    -1,    91,    92,    93,
      94,    53,    -1,    -1,    98,    57,    -1,   101,   102,   103,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    73,    -1,    75,    76,    77,    78,    79,    80,    81,
      82,    -1,    84,    85,    86,    87,    88,    89,    -1,    91,
      92,    93,    94,    53,    -1,    -1,    98,    57,    -1,   101,
     102,   103,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    73,    -1,    75,    76,    77,    78,    79,
      80,    81,    82,    -1,    84,    85,    86,    87,    88,    89,
      -1,    91,    92,    93,    94,    53,    -1,    -1,    98,    57,
      -1,   101,   102,   103,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    73,    -1,    75,    76,    77,
      78,    79,    80,    81,    82,    -1,    84,    85,    86,    87,
      88,    89,    -1,    91,    92,    93,    94,    53,    -1,    -1,
      98,    57,    -1,   101,   102,   103,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,    -1,    75,
      76,    77,    78,    79,    80,    81,    82,    -1,    84,    85,
      86,    87,    88,    89,    -1,    91,    92,    93,    94,    53,
      -1,    -1,    98,    57,    -1,   101,   102,   103,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,
      -1,    75,    76,    77,    78,    79,    80,    81,    82,    -1,
      84,    85,    86,    87,    88,    89,    -1,    91,    92,    93,
      94,    53,    -1,    -1,    98,    57,    -1,   101,   102,   103,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    73,    -1,    75,    76,    77,    78,    79,    80,    81,
      82,    -1,    84,    85,    86,    87,    88,    89,    -1,    91,
      92,    93,    94,    53,    -1,    -1,    98,    57,    -1,   101,
     102,   103,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    73,    -1,    75,    76,    77,    78,    79,
      80,    81,    82,    -1,    84,    85,    86,    87,    88,    89,
      -1,    91,    92,    93,    94,    53,    -1,    -1,    98,    57,
      -1,   101,   102,   103,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    73,    -1,    75,    76,    77,
      78,    79,    80,    81,    82,    -1,    84,    85,    86,    87,
      88,    89,    -1,    91,    92,    93,    94,    53,    -1,    -1,
      98,    57,    -1,   101,   102,   103,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,    -1,    75,
      76,    77,    78,    79,    80,    81,    82,    -1,    84,    85,
      86,    87,    88,    89,    -1,    91,    92,    93,    94,    53,
      -1,    -1,    98,    57,    -1,   101,   102,   103,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,
      -1,    75,    76,    77,    78,    79,    80,    81,    82,    -1,
      84,    85,    86,    87,    88,    89,    -1,    91,    92,    93,
      94,    -1,    53,    -1,    -1,    99,    57,   101,   102,   103,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    73,    -1,    75,    76,    77,    78,    79,    80,
      81,    82,    -1,    84,    85,    86,    87,    88,    89,    -1,
      91,    92,    93,    94,    -1,    53,    -1,    -1,    99,    57,
     101,   102,   103,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    73,    -1,    75,    76,    77,
      78,    79,    80,    81,    82,    -1,    84,    85,    86,    87,
      88,    89,    -1,    91,    92,    93,    94,    95,    53,    -1,
      -1,    -1,    57,   101,   102,   103,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,    -1,
      75,    76,    77,    78,    79,    80,    81,    82,    -1,    84,
      85,    86,    87,    88,    89,    -1,    91,    92,    93,    94,
      53,    -1,    -1,    98,    57,    -1,   101,   102,   103,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      73,    -1,    75,    76,    77,    78,    79,    80,    81,    82,
      -1,    84,    85,    86,    87,    88,    89,    -1,    91,    92,
      93,    94,    53,    -1,    -1,    98,    57,    -1,   101,   102,
     103,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    73,    -1,    75,    76,    77,    78,    79,    80,
      81,    82,    -1,    84,    85,    86,    87,    88,    89,    -1,
      91,    92,    93,    94,    53,    -1,    -1,    98,    57,    -1,
     101,   102,   103,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    73,    -1,    75,    76,    77,    78,
      79,    80,    81,    82,    -1,    84,    85,    86,    87,    88,
      89,    -1,    91,    92,    93,    94,    95,    53,    -1,    -1,
      -1,    57,   101,   102,   103,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,    -1,    75,
      76,    77,    78,    79,    80,    81,    82,    -1,    84,    85,
      86,    87,    88,    89,    -1,    91,    92,    93,    94,    53,
      -1,    -1,    98,    57,    -1,   101,   102,   103,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,
      -1,    75,    76,    77,    78,    79,    80,    81,    82,    -1,
      84,    85,    86,    87,    88,    89,    -1,    91,    92,    93,
      94,    53,    96,    -1,    -1,    57,    -1,   101,   102,   103,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    73,    -1,    75,    76,    77,    78,    79,    80,    81,
      82,    -1,    84,    85,    86,    87,    88,    89,    -1,    91,
      92,    93,    94,    53,    -1,    -1,    98,    57,    -1,   101,
     102,   103,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    73,    -1,    75,    76,    77,    78,    79,
      80,    81,    82,    -1,    84,    85,    86,    87,    88,    89,
      -1,    91,    92,    93,    94,    53,    -1,    -1,    98,    57,
      -1,   101,   102,   103,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    73,    -1,    75,    76,    77,
      78,    79,    80,    81,    82,    -1,    84,    85,    86,    87,
      88,    89,    -1,    91,    92,    93,    94,    53,    -1,    -1,
      98,    57,    -1,   101,   102,   103,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,    -1,    75,
      76,    77,    78,    79,    80,    81,    82,    -1,    84,    85,
      86,    87,    88,    89,    -1,    91,    92,    93,    94,    53,
      -1,    -1,    98,    57,    -1,   101,   102,   103,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,
      -1,    75,    76,    77,    78,    79,    80,    81,    82,    -1,
      84,    85,    86,    87,    88,    89,    -1,    91,    92,    93,
      94,    53,    -1,    -1,    -1,    57,    -1,   101,   102,   103,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    75,    76,    77,    78,    79,    80,    81,
      82,    -1,    84,    85,    86,    87,    88,    89,    -1,    91,
      92,    93,    94,    53,    -1,    -1,    -1,    57,    -1,   101,
     102,   103,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    78,    79,
      80,    81,    82,    -1,    84,    85,    86,    87,    88,    89,
      -1,    91,    92,    93,    94,    53,    -1,    -1,    -1,    57,
      -1,   101,   102,   103,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    53,    -1,    -1,
      -1,    57,    80,    81,    82,    -1,    84,    85,    86,    87,
      88,    89,    -1,    91,    92,    93,    94,    -1,    -1,    -1,
      -1,    -1,    -1,   101,   102,   103,    82,    -1,    84,    85,
      86,    87,    88,    89,    -1,    91,    92,    93,    94,    -1,
      -1,    -1,    -1,    -1,    -1,   101,   102,   103
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,     5,     7,     8,     9,    10,    11,    14,    16,
      18,    58,   109,   110,   111,   113,   115,   118,   120,   122,
     124,   125,   126,     4,    58,    58,    58,    12,    13,    15,
      21,    23,    24,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    49,    50,    54,    59,    60,    67,    68,
      69,    70,    71,    72,    79,    97,   105,   128,   129,   130,
     131,   132,   133,   134,   135,   136,   139,   140,   135,   127,
     135,    21,    22,    58,   116,   117,   140,   114,   140,    58,
     112,   140,     0,    58,   111,   113,   115,   118,   120,   122,
     124,   125,   126,    58,   111,   113,   115,   118,   120,   122,
     124,   125,   126,   119,   140,   121,   140,   122,   123,   140,
     119,   121,   123,    97,    97,    97,    97,   103,   144,   144,
      97,    97,    97,    97,    97,    97,    97,    97,    97,    97,
      97,    97,    55,    59,    59,   135,   135,   135,   135,   138,
      53,    57,    73,    75,    76,    77,    78,    79,    80,    81,
      82,    84,    85,    86,    87,    88,    89,    91,    92,    93,
      94,    95,   101,   102,   103,   107,    95,    95,    97,    97,
     116,   117,    95,    56,    96,   112,   140,    56,    58,    99,
      58,    99,    58,    99,    58,    58,    58,   135,   135,    25,
      26,    48,   135,    59,    97,    97,   135,   135,   135,   135,
     137,   135,   135,   135,   135,   135,   135,   135,   135,    59,
     100,   100,    98,    96,   106,   135,   135,   135,   135,   135,
     135,   135,   135,   135,   135,    97,   135,   135,   135,   135,
     135,   135,   135,   135,   135,   135,    58,   135,   135,    59,
     135,    59,    69,    58,    58,   140,   140,    95,    58,   135,
     140,    56,   135,    23,    24,    25,    26,    51,    59,    61,
      62,    63,    64,    66,   141,   142,   143,   141,   141,    96,
      96,   144,    37,    97,    98,   104,   135,   135,    98,    98,
      98,    96,    98,    98,    98,    98,    98,    98,    98,    98,
      98,    59,    59,   135,    99,   135,    58,   104,    99,    58,
      58,    98,    98,    58,    95,   135,    95,   144,   144,   144,
      37,   144,    55,    37,    95,    95,    95,   135,   135,    37,
     141,   140,    98,    98,   135,   135,    98,   135,    56,    56,
      58,    95,    58,    37,   141,    59,   141,    58,    58,    58,
      98,    96,   141,    96,    98,   104,   136,   135,    58,   141,
     135,    96,   135,    96,    98,   135,    98,   135,    98,    98
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,   108,   109,   109,   109,   109,   109,   109,   109,   109,
     109,   109,   109,   109,   109,   109,   109,   109,   109,   109,
     109,   109,   109,   109,   109,   109,   109,   109,   109,   109,
     109,   109,   110,   111,   111,   112,   112,   113,   114,   114,
     115,   115,   116,   116,   117,   117,   117,   118,   118,   118,
     118,   119,   120,   120,   120,   120,   121,   122,   122,   122,
     122,   123,   124,   124,   125,   125,   126,   126,   126,   127,
     128,   128,   128,   128,   128,   128,   128,   129,   129,   130,
     130,   131,   132,   133,   134,   135,   135,   135,   135,   135,
     135,   135,   135,   135,   135,   135,   135,   135,   135,   135,
     135,   135,   135,   135,   135,   135,   135,   135,   135,   135,
     135,   135,   135,   135,   135,   135,   135,   135,   135,   135,
     135,   135,   135,   135,   135,   135,   135,   135,   135,   135,
     135,   135,   135,   135,   135,   135,   136,   137,   137,   138,
     138,   139,   140,   140,   140,   141,   141,   141,   141,   141,
     141,   141,   142,   142,   142,   143,   143,   143,   144
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     3,     5,     6,     2,     1,     5,
       2,     3,     3,     4,     3,     6,     6,     4,     3,     3,
       2,     5,     4,     3,     3,     2,     5,     4,     3,     3,
       2,     5,     5,     4,     5,     4,     5,     4,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     4,     4,     1,
       1,     1,     1,     3,     1,     1,     1,     3,     2,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     2,     3,     3,     3,     3,     3,     3,     3,     4,
       6,     3,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     5,     5,     5,     3,     3,     5,
       8,     6,     9,     9,     8,     1,     4,     1,     3,     1,
       3,     1,     1,     3,     3,     1,     1,     1,     1,     1,
       1,     3,     2,     2,     2,     4,     3,     3,     3
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (enc, YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)                                \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;        \
          (Current).first_column = YYRHSLOC (Rhs, 1).first_column;      \
          (Current).last_line    = YYRHSLOC (Rhs, N).last_line;         \
          (Current).last_column  = YYRHSLOC (Rhs, N).last_column;       \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).first_line   = (Current).last_line   =              \
            YYRHSLOC (Rhs, 0).last_line;                                \
          (Current).first_column = (Current).last_column =              \
            YYRHSLOC (Rhs, 0).last_column;                              \
        }                                                               \
    while (0)
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K])


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL

/* Print *YYLOCP on YYO.  Private, do not rely on its existence. */

YY_ATTRIBUTE_UNUSED
static int
yy_location_print_ (FILE *yyo, YYLTYPE const * const yylocp)
{
  int res = 0;
  int end_col = 0 != yylocp->last_column ? yylocp->last_column - 1 : 0;
  if (0 <= yylocp->first_line)
    {
      res += YYFPRINTF (yyo, "%d", yylocp->first_line);
      if (0 <= yylocp->first_column)
        res += YYFPRINTF (yyo, ".%d", yylocp->first_column);
    }
  if (0 <= yylocp->last_line)
    {
      if (yylocp->first_line < yylocp->last_line)
        {
          res += YYFPRINTF (yyo, "-%d", yylocp->last_line);
          if (0 <= end_col)
            res += YYFPRINTF (yyo, ".%d", end_col);
        }
      else if (0 <= end_col && yylocp->first_column < end_col)
        res += YYFPRINTF (yyo, "-%d", end_col);
    }
  return res;
 }

#  define YY_LOCATION_PRINT(File, Loc)          \
  yy_location_print_ (File, &(Loc))

# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, Location, enc); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, smvEncoder &enc)
{
  FILE *yyoutput = yyo;
  YYUSE (yyoutput);
  YYUSE (yylocationp);
  YYUSE (enc);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyo, yytoknum[yytype], *yyvaluep);
# endif
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, smvEncoder &enc)
{
  YYFPRINTF (yyo, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  YY_LOCATION_PRINT (yyo, *yylocationp);
  YYFPRINTF (yyo, ": ");
  yy_symbol_value_print (yyo, yytype, yyvaluep, yylocationp, enc);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, YYLTYPE *yylsp, int yyrule, smvEncoder &enc)
{
  unsigned long yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &yyvsp[(yyi + 1) - (yynrhs)]
                       , &(yylsp[(yyi + 1) - (yynrhs)])                       , enc);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, yylsp, Rule, enc); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            else
              goto append;

          append:
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return (YYSIZE_T) (yystpcpy (yyres, yystr) - yyres);
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM)
                    yysize = yysize1;
                  else
                    return 2;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
    default: /* Avoid compiler warnings. */
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM)
      yysize = yysize1;
    else
      return 2;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, YYLTYPE *yylocationp, smvEncoder &enc)
{
  YYUSE (yyvaluep);
  YYUSE (yylocationp);
  YYUSE (enc);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Location data for the lookahead symbol.  */
YYLTYPE yylloc
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
  = { 1, 1, 1, 1 }
# endif
;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (smvEncoder &enc)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.
       'yyls': related to locations.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    /* The location stack.  */
    YYLTYPE yylsa[YYINITDEPTH];
    YYLTYPE *yyls;
    YYLTYPE *yylsp;

    /* The locations where the error started and ended.  */
    YYLTYPE yyerror_range[3];

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
  YYLTYPE yyloc;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N), yylsp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yylsp = yyls = yylsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  yylsp[0] = yylloc;
  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yynewstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  *yyssp = (yytype_int16) yystate;

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    goto yyexhaustedlab;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = (YYSIZE_T) (yyssp - yyss + 1);

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;
        YYLTYPE *yyls1 = yyls;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yyls1, yysize * sizeof (*yylsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
        yyls = yyls1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
        YYSTACK_RELOCATE (yyls_alloc, yyls);
# undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
      yylsp = yyls + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END
  *++yylsp = yylloc;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

  /* Default location. */
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
  yyerror_range[1] = yyloc;
  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 35:
#line 111 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              enc.terms_[(yyvsp[-4].str)] = a-> getTerm();
}
#line 2009 "smvparser.tab.c"
    break;

  case 36:
#line 115 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              enc.terms_[(yyvsp[-4].str)] = a-> getTerm();
            }
#line 2018 "smvparser.tab.c"
    break;

  case 44:
#line 132 "smvparser.y"
    {
          node *a = (yyvsp[0].n); 
          smt::Term state = a->getTerm();
}
#line 2027 "smvparser.tab.c"
    break;

  case 45:
#line 136 "smvparser.y"
    {
          node *a = (yyvsp[0].n); 
          smt::Term init = enc.terms_[(yyvsp[-3].str)];
          smt::Term e = enc.solver_->make_term(smt::Equal, init, a->getTerm());
          enc.rts_.constrain_init(e);
        }
#line 2038 "smvparser.tab.c"
    break;

  case 46:
#line 142 "smvparser.y"
    {
          node *a = (yyvsp[0].n); 
          smt::Term state = enc.terms_[(yyvsp[-3].str)];
          smt::Term e = enc.solver_->make_term(smt::Equal, state, a->getTerm());
          enc.rts_.constrain_trans(e);
        }
#line 2049 "smvparser.tab.c"
    break;

  case 51:
#line 157 "smvparser.y"
    {
         //cout <<"find an ivar"<<endl;
         node *a = (yyvsp[-2].n);
         smt::Term input = enc.rts_.make_input((yyvsp[-4].str), a->getSort());
         enc.terms_[(yyvsp[-4].str)] = input;
         //free($1);
    }
#line 2061 "smvparser.tab.c"
    break;

  case 56:
#line 172 "smvparser.y"
    {
         node *a = (yyvsp[-2].n);
         smt::Term state = enc.rts_.make_state((yyvsp[-4].str), a->getSort());
         enc.terms_[(yyvsp[-4].str)] = state;
         //cout<<"find a var" <<endl;
         //free($1);
    }
#line 2073 "smvparser.tab.c"
    break;

  case 61:
#line 187 "smvparser.y"
    {
      node *a = (yyvsp[-2].n);
      smt::Term state = enc.rts_.make_state((yyvsp[-4].str), a->getSort());
      enc.terms_[(yyvsp[-4].str)] = state;
      smt::Term n = enc.rts_.next(state);
      smt::Term e = enc.solver_->make_term(smt::Equal, n, state);
      enc.rts_.constrain_trans(e);
      //cout<<"find a frozen var" <<endl;
      //free($1);
  }
#line 2088 "smvparser.tab.c"
    break;

  case 62:
#line 198 "smvparser.y"
    {
        node *a = (yyvsp[-3].n);
        enc.rts_.constrain_init(a->getTerm());
}
#line 2097 "smvparser.tab.c"
    break;

  case 63:
#line 202 "smvparser.y"
    {
        node *a = (yyvsp[-2].n);
        enc.rts_.constrain_init(a->getTerm());
      }
#line 2106 "smvparser.tab.c"
    break;

  case 64:
#line 207 "smvparser.y"
    {
            node *a = (yyvsp[-3].n);
            enc.rts_.constrain_trans(a->getTerm());
            //cout <<"find a trans"<<endl;
}
#line 2116 "smvparser.tab.c"
    break;

  case 65:
#line 212 "smvparser.y"
    {
            node *a = (yyvsp[-2].n);
            enc.rts_.constrain_trans(a->getTerm());
            //cout <<"find a trans"<<endl;
}
#line 2126 "smvparser.tab.c"
    break;

  case 66:
#line 219 "smvparser.y"
    {
                //cout<<"find an invarspec" <<endl;
                node *a = (yyvsp[-3].n);
                smt::Term prop = a->getTerm();
                enc.propvec_.push_back(prop);
}
#line 2137 "smvparser.tab.c"
    break;

  case 67:
#line 225 "smvparser.y"
    {
                //cout<<"find an invarspec" <<endl;
                node *a = (yyvsp[-2].n);
                smt::Term prop = a->getTerm();
                enc.propvec_.push_back(prop);
                }
#line 2148 "smvparser.tab.c"
    break;

  case 68:
#line 231 "smvparser.y"
    {
                //cout<<"find an invarspec" <<endl;
                node *a = (yyvsp[-1].n);
                smt::Term prop = a->getTerm();
                enc.propvec_.push_back(prop);
              }
#line 2159 "smvparser.tab.c"
    break;

  case 69:
#line 250 "smvparser.y"
    {
              (yyval.n) = (yyvsp[0].n);
}
#line 2167 "smvparser.tab.c"
    break;

  case 70:
#line 254 "smvparser.y"
    {
      smt::Term con = enc.solver_->make_term((yyvsp[0].bl));
      (yyval.n) = new node(con);
}
#line 2176 "smvparser.tab.c"
    break;

  case 73:
#line 260 "smvparser.y"
    {
           (yyval.n) = (yyvsp[0].n);
          }
#line 2184 "smvparser.tab.c"
    break;

  case 77:
#line 267 "smvparser.y"
    {
          smt::Sort sort_ = enc.solver_->make_sort(smt::BV, (yyvsp[-2].num));
          char* temp = (yyvsp[-3].str);
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
          std::string int_val = std::to_string((yyvsp[0].num));
          smt::Term num = enc.solver_->make_term(int_val, sort_, base);
          (yyval.n) = new node(num); }
#line 2209 "smvparser.tab.c"
    break;

  case 78:
#line 287 "smvparser.y"
    {
          smt::Sort sort_ = enc.solver_->make_sort(smt::BV, (yyvsp[-2].num));
          char* temp = (yyvsp[-3].str);
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
          std::string int_val = std::to_string((yyvsp[0].num));
          smt::Term num = enc.solver_->make_term(int_val, sort_, base);
          (yyval.n) = new node(num);
   }
#line 2235 "smvparser.tab.c"
    break;

  case 79:
#line 310 "smvparser.y"
    {
                (yyval.bl) = true;
          }
#line 2243 "smvparser.tab.c"
    break;

  case 80:
#line 313 "smvparser.y"
    {
                (yyval.bl) = false;    
          }
#line 2251 "smvparser.tab.c"
    break;

  case 81:
#line 321 "smvparser.y"
    {
  throw CosaException("No integer constant now");
}
#line 2259 "smvparser.tab.c"
    break;

  case 82:
#line 324 "smvparser.y"
    {
  throw CosaException("No real constant now");
}
#line 2267 "smvparser.tab.c"
    break;

  case 83:
#line 327 "smvparser.y"
    {
  throw CosaException("No range constant now");
}
#line 2275 "smvparser.tab.c"
    break;

  case 84:
#line 330 "smvparser.y"
    {
  throw CosaException("No clock constant now");
}
#line 2283 "smvparser.tab.c"
    break;

  case 85:
#line 333 "smvparser.y"
    {
            (yyval.n) = (yyvsp[0].n);
}
#line 2291 "smvparser.tab.c"
    break;

  case 86:
#line 336 "smvparser.y"
    {
              smt::Term tok = enc.terms_.at((yyvsp[0].str));
              (yyval.n) = new node(tok);
            }
#line 2300 "smvparser.tab.c"
    break;

  case 87:
#line 349 "smvparser.y"
    {
              (yyval.n) = (yyvsp[-1].n);
            }
#line 2308 "smvparser.tab.c"
    break;

  case 89:
#line 353 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term e = enc.solver_->make_term(smt::BVAnd, a->getTerm(), b->getTerm());
              (yyval.n) = new node(e);
            }
#line 2319 "smvparser.tab.c"
    break;

  case 90:
#line 359 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term e = enc.solver_->make_term(smt::BVOr, a->getTerm(), b->getTerm());
              (yyval.n) = new node(e);
            }
#line 2330 "smvparser.tab.c"
    break;

  case 91:
#line 365 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term e = enc.solver_->make_term(smt::BVXor, a->getTerm(), b->getTerm());
              (yyval.n) = new node(e);
            }
#line 2341 "smvparser.tab.c"
    break;

  case 92:
#line 371 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term e = enc.solver_->make_term(smt::BVXnor, a->getTerm(), b->getTerm());
              (yyval.n) = new node(e);
            }
#line 2352 "smvparser.tab.c"
    break;

  case 93:
#line 377 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term e = enc.solver_->make_term(smt::Implies, a->getTerm(), b->getTerm());
              (yyval.n) = new node(e);
            }
#line 2363 "smvparser.tab.c"
    break;

  case 94:
#line 383 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term e = enc.solver_->make_term(smt::Iff, a->getTerm(), b->getTerm());
              (yyval.n) = new node(e);
            }
#line 2374 "smvparser.tab.c"
    break;

  case 95:
#line 389 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term e = enc.solver_->make_term(smt::Equal, a->getTerm(), b->getTerm());
              (yyval.n) = new node(e);
            }
#line 2385 "smvparser.tab.c"
    break;

  case 96:
#line 395 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term e = enc.solver_->make_term(smt::Distinct, a->getTerm(), b->getTerm());
              (yyval.n) = new node(e);
            }
#line 2396 "smvparser.tab.c"
    break;

  case 97:
#line 401 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term res = enc.solver_->make_term(smt::BVUlt, a->getTerm(), b->getTerm());
              (yyval.n) = new node(res);
            }
#line 2407 "smvparser.tab.c"
    break;

  case 98:
#line 407 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term res = enc.solver_->make_term(smt::BVUgt, a->getTerm(), b->getTerm());
              (yyval.n) = new node(res);
            }
#line 2418 "smvparser.tab.c"
    break;

  case 99:
#line 413 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term res = enc.solver_->make_term(smt::BVUle, a->getTerm(), b->getTerm());
              (yyval.n) = new node(res);
            }
#line 2429 "smvparser.tab.c"
    break;

  case 100:
#line 419 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term res = enc.solver_->make_term(smt::BVUge, a->getTerm(), b->getTerm());
              (yyval.n) = new node(res);
            }
#line 2440 "smvparser.tab.c"
    break;

  case 101:
#line 425 "smvparser.y"
    {

            }
#line 2448 "smvparser.tab.c"
    break;

  case 102:
#line 428 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term res = enc.solver_->make_term(smt::BVAdd, a->getTerm(), b->getTerm());
              (yyval.n) = new node(res);
            }
#line 2459 "smvparser.tab.c"
    break;

  case 103:
#line 434 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term res = enc.solver_->make_term(smt::BVSub, a->getTerm(), b->getTerm());
              (yyval.n) = new node(res);
            }
#line 2470 "smvparser.tab.c"
    break;

  case 104:
#line 440 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term res = enc.solver_->make_term(smt::BVMul, a->getTerm(), b->getTerm());
              (yyval.n) = new node(res);
            }
#line 2481 "smvparser.tab.c"
    break;

  case 105:
#line 446 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term res = enc.solver_->make_term(smt::BVUdiv, a->getTerm(), b->getTerm());
              (yyval.n) = new node(res);
            }
#line 2492 "smvparser.tab.c"
    break;

  case 106:
#line 452 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term res = enc.solver_->make_term(smt::BVSmod, a->getTerm(), b->getTerm());
              (yyval.n) = new node(res);
            }
#line 2503 "smvparser.tab.c"
    break;

  case 107:
#line 458 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term res = enc.solver_->make_term(smt::BVLshr, a->getTerm(), b->getTerm());
              (yyval.n) = new node(res);
            }
#line 2514 "smvparser.tab.c"
    break;

  case 108:
#line 464 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term res = enc.solver_->make_term(smt::BVShl, a->getTerm(), b->getTerm());
              (yyval.n) = new node(res);
            }
#line 2525 "smvparser.tab.c"
    break;

  case 111:
#line 472 "smvparser.y"
    {
              node *a = (yyvsp[-2].n);
              node *b = (yyvsp[0].n);
              smt::Term res = enc.solver_->make_term(smt::Concat, a->getTerm(), b->getTerm());
              (yyval.n) = new node(res);
            }
#line 2536 "smvparser.tab.c"
    break;

  case 130:
#line 496 "smvparser.y"
    {
            node *a = (yyvsp[-5].n);
            node *b = (yyvsp[-3].n);
            node *c = (yyvsp[-1].n);
            smt::Term write_r =  enc.solver_->make_term(smt::Store, a->getTerm(),b->getTerm(),c->getTerm());
            (yyval.n) = new node(write_r);
          }
#line 2548 "smvparser.tab.c"
    break;

  case 131:
#line 503 "smvparser.y"
    {
            node *a = (yyvsp[-3].n);
            node *b = (yyvsp[-1].n);
            smt::Term read_r =  enc.solver_->make_term(smt::Select, a->getTerm(),b->getTerm());
            (yyval.n) = new node(read_r);
          }
#line 2559 "smvparser.tab.c"
    break;

  case 135:
#line 512 "smvparser.y"
    {
            (yyval.n) = (yyvsp[0].n);
          }
#line 2567 "smvparser.tab.c"
    break;

  case 136:
#line 516 "smvparser.y"
    {
          node *a = (yyvsp[-1].n);
          smt::Term n = enc.rts_.next(a->getTerm());
          (yyval.n) = new node(n);
}
#line 2577 "smvparser.tab.c"
    break;

  case 142:
#line 530 "smvparser.y"
    {
                (yyval.str) = (yyvsp[0].str);
 }
#line 2585 "smvparser.tab.c"
    break;

  case 143:
#line 533 "smvparser.y"
    {
               (yyval.str) = strcat((yyvsp[-2].str),(yyvsp[0].str));
 }
#line 2593 "smvparser.tab.c"
    break;

  case 144:
#line 536 "smvparser.y"
    {
                std::string temp = (yyvsp[-2].str);
                cout << "find a string" <<(yyvsp[-2].str)<< endl;
                (yyval.str) = (yyvsp[-2].str) + (yyvsp[0].num);
 }
#line 2603 "smvparser.tab.c"
    break;

  case 145:
#line 542 "smvparser.y"
    {
                smt::Sort sort_ = enc.solver_->make_sort(smt::REAL);
                (yyval.n) =  new node (sort_);
                throw CosaException("No real type now in boolector");
                }
#line 2613 "smvparser.tab.c"
    break;

  case 146:
#line 547 "smvparser.y"
    {
                  smt::Sort sort_ = enc.solver_->make_sort(smt::INT);
                  (yyval.n) =  new node (sort_);
                  throw CosaException("No integer type now in boolector");  
                }
#line 2623 "smvparser.tab.c"
    break;

  case 147:
#line 552 "smvparser.y"
    {
                smt::Sort sort_ = enc.solver_->make_sort(smt::BOOL);
                (yyval.n) =  new node (sort_);
                }
#line 2632 "smvparser.tab.c"
    break;

  case 148:
#line 556 "smvparser.y"
    {
                  throw CosaException("No clock type now");
                }
#line 2640 "smvparser.tab.c"
    break;

  case 149:
#line 559 "smvparser.y"
    {
                  (yyval.n) = (yyvsp[0].n);
                }
#line 2648 "smvparser.tab.c"
    break;

  case 150:
#line 562 "smvparser.y"
    {
                (yyval.n) = (yyvsp[0].n);
}
#line 2656 "smvparser.tab.c"
    break;

  case 152:
#line 567 "smvparser.y"
    {
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, (yyvsp[0].num));
        (yyval.n) =  new node (sort_);
}
#line 2665 "smvparser.tab.c"
    break;

  case 153:
#line 571 "smvparser.y"
    {
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, (yyvsp[0].num));
        (yyval.n) =  new node (sort_);
}
#line 2674 "smvparser.tab.c"
    break;

  case 154:
#line 575 "smvparser.y"
    {
        smt::Sort sort_ = enc.solver_->make_sort(smt::BV, (yyvsp[0].num));
        (yyval.n) =  new node (sort_);
}
#line 2683 "smvparser.tab.c"
    break;

  case 155:
#line 580 "smvparser.y"
    {
            smt::Sort arraysort = enc.solver_->make_sort(smt::BV,(yyvsp[-2].num));
            node *a = (yyvsp[0].n);
            smt::Sort sort_ = enc.solver_->make_sort(smt::ARRAY, arraysort,a->getSort());
            (yyval.n) = new node(sort_);
          }
#line 2694 "smvparser.tab.c"
    break;

  case 156:
#line 586 "smvparser.y"
    {
            throw CosaException("no array integer type now");
          }
#line 2702 "smvparser.tab.c"
    break;

  case 158:
#line 592 "smvparser.y"
    {
        (yyval.num)  = (yyvsp[-1].num);
    }
#line 2710 "smvparser.tab.c"
    break;


#line 2714 "smvparser.tab.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;
  *++yylsp = yyloc;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (enc, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (enc, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }

  yyerror_range[1] = yylloc;

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, &yylloc, enc);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;

      yyerror_range[1] = *yylsp;
      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp, yylsp, enc);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  yyerror_range[2] = yylloc;
  /* Using YYLLOC is tempting, but would change the location of
     the lookahead.  YYLOC is available though.  */
  YYLLOC_DEFAULT (yyloc, yyerror_range, 2);
  *++yylsp = yyloc;

  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;


#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (enc, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif


/*-----------------------------------------------------.
| yyreturn -- parsing is finished, return the result.  |
`-----------------------------------------------------*/
yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, &yylloc, enc);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, yylsp, enc);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 595 "smvparser.y"


void smvEncoder::parse(std::string filename){
    FILE *myfile = fopen(filename.c_str(), "r");
  // make sure ist's valid:
  if (!myfile) {
    std::cout << "NO input file!" << std::endl;
  }
  yyin = myfile;
  yyparse(*this); 
}
