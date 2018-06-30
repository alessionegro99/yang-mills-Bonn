#ifndef FUNCTION_POINTERS_C
#define FUNCTION_POINTERS_C

#include"../include/function_pointers.h"
#include"../include/macro.h"
#include"../include/su2.h"
#include"../include/su2_upd.h"
#include"../include/sun.h"
#include"../include/sun_upd.h"
#include"../include/tens_prod.h"
#include"../include/u1.h"
#include"../include/u1_upd.h"

void init_function_pointers(void)
  {
  #if NCOLOR == 1

  one  = &one_U1;
  zero = &zero_U1;
  equal     = &equal_U1;
  equal_dag = &equal_dag_U1;
  plus_equal     = &plus_equal_U1;
  plus_equal_dag = &plus_equal_dag_U1;
  minus_equal     = &minus_equal_U1;
  minus_equal_times_real  = &minus_equal_times_real_U1;
  minus_equal_dag = &minus_equal_dag_U1;
  lin_comb       = &lin_comb_U1;
  lin_comb_dag1  = &lin_comb_dag1_U1;
  lin_comb_dag2  = &lin_comb_dag2_U1;
  lin_comb_dag12 = &lin_comb_dag12_U1;
  times_equal_real = &times_equal_real_U1;
  times_equal     = &times_equal_U1;
  times_equal_dag = &times_equal_dag_U1;
  times       = &times_U1;
  times_dag1  = &times_dag1_U1;
  times_dag2  = &times_dag2_U1;
  times_dag12 = &times_dag12_U1;
  rand_matrix = &rand_matrix_U1;
  norm = &norm_U1;
  retr = &retr_U1;
  imtr = &imtr_U1;
  unitarize = &unitarize_U1;
  ta = &ta_U1;
  taexp = &taexp_U1;
  print_on_screen = &print_on_screen_U1;
  print_on_file   = &print_on_file_U1;
  print_on_binary_file_bigen   = &print_on_binary_file_bigen_U1;
  read_from_file   = &read_from_file_U1;
  read_from_binary_file_bigen   = &read_from_binary_file_bigen_U1;
  TensProd_init=&TensProd_init_U1;
  single_heatbath = &single_heatbath_U1;
  single_overrelaxation = &single_overrelaxation_U1;
  cool = &cool_U1;

  #elif NCOLOR == 2

  one  = &one_Su2;
  zero = &zero_Su2;
  equal     = &equal_Su2;
  equal_dag = &equal_dag_Su2;
  plus_equal     = &plus_equal_Su2;
  plus_equal_dag = &plus_equal_dag_Su2;
  minus_equal     = &minus_equal_Su2;
  minus_equal_times_real  = &minus_equal_times_real_Su2;
  minus_equal_dag = &minus_equal_dag_Su2;
  lin_comb       = &lin_comb_Su2;
  lin_comb_dag1  = &lin_comb_dag1_Su2;
  lin_comb_dag2  = &lin_comb_dag2_Su2;
  lin_comb_dag12 = &lin_comb_dag12_Su2;
  times_equal_real = &times_equal_real_Su2;
  times_equal     = &times_equal_Su2;
  times_equal_dag = &times_equal_dag_Su2;
  times       = &times_Su2;
  times_dag1  = &times_dag1_Su2;
  times_dag2  = &times_dag2_Su2;
  times_dag12 = &times_dag12_Su2;
  rand_matrix = &rand_matrix_Su2;
  norm = &norm_Su2;
  retr = &retr_Su2;
  imtr = &imtr_Su2;
  unitarize = &unitarize_Su2;
  ta = &ta_Su2;
  taexp = &taexp_Su2;
  print_on_screen = &print_on_screen_Su2;
  print_on_file   = &print_on_file_Su2;
  print_on_binary_file_bigen   = &print_on_binary_file_bigen_Su2;
  read_from_file   = &read_from_file_Su2;
  read_from_binary_file_bigen   = &read_from_binary_file_bigen_Su2;
  TensProd_init=&TensProd_init_Su2;
  single_heatbath = &single_heatbath_Su2;
  single_overrelaxation = &single_overrelaxation_Su2;
  cool = &cool_Su2;

  #else

  one  = &one_SuN;
  zero = &zero_SuN;
  equal     = &equal_SuN;
  equal_dag = &equal_dag_SuN;
  plus_equal     = &plus_equal_SuN;
  plus_equal_dag = &plus_equal_dag_SuN;
  minus_equal     = &minus_equal_SuN;
  minus_equal_times_real= &minus_equal_times_real_SuN;
  minus_equal_dag = &minus_equal_dag_SuN;
  lin_comb       = &lin_comb_SuN;
  lin_comb_dag1  = &lin_comb_dag1_SuN;
  lin_comb_dag2  = &lin_comb_dag2_SuN;
  lin_comb_dag12 = &lin_comb_dag12_SuN;
  times_equal_real = &times_equal_real_SuN;
  times_equal     = &times_equal_SuN;
  times_equal_dag = &times_equal_dag_SuN;
  times       = &times_SuN;
  times_dag1  = &times_dag1_SuN;
  times_dag2  = &times_dag2_SuN;
  times_dag12 = &times_dag12_SuN;
  rand_matrix = &rand_matrix_SuN;
  norm = &norm_SuN;
  retr = &retr_SuN;
  imtr = &imtr_SuN;
  unitarize = &unitarize_SuN;
  ta = &ta_SuN;
  taexp = &taexp_SuN;
  print_on_screen = &print_on_screen_SuN;
  print_on_file   = &print_on_file_SuN;
  print_on_binary_file_bigen   = &print_on_binary_file_bigen_SuN;
  read_from_file   = &read_from_file_SuN;
  read_from_binary_file_bigen   = &read_from_binary_file_bigen_SuN;
  TensProd_init=&TensProd_init_SuN;
  single_heatbath = &single_heatbath_SuN;
  single_overrelaxation = &single_overrelaxation_SuN;
  cool = &cool_SuN;

  #endif
  }


#endif
