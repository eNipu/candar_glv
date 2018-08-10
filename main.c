//
//  main.c
//  Fp16_EffSM
//
//

#include"KSS_16.h"

int count_eca_new, count_ecd_new,count_eca_BIN, count_ecd_BIN,count_eca_BIN_map, count_ecd_BIN_map, count_eca_WIN, count_ecd_WIN, count_eca_WIN_map,count_ecd_WIN_map, count_eca_NAF, count_ecd_NAF,count_eca_NAF_map,count_ecd_NAF_map, count_eca_ML, count_ecd_ML,count_eca_ML_map, count_ecd_ML_map;

int count_eca_new_pre, count_ecd_new_pre;

int Big_M, Small_m, Big_add, Sqr;

long int myNAF[5000];
long int naf_index = 0;

struct timeval tonaf_t0;
struct timeval tonaf_t1;
float elapsed_tonaf_t0;
float elapsed_tonaf_t1;
mpz_t B, A;

struct Fp C1,C1_INV;
struct Fp p_C4, p_M_C4;
struct Fp p_C8, p_M_C_C8;
struct Fp p_C_C4_C8, p_M_C_C4_C8;
struct Fp p_C16, p_M_C_C16;
struct Fp p_C4_C16, p_C_C4_C16, p_M_C_C4_C16;
struct Fp p_C4_C8_C16, p_C_C4_C8_C16, p_M_C_C_C4_C8_C16;
struct Fp p_C8_C16, p_C_C8_C16, p_M_C_C8_C16;
mpz_t p8p1dr;

struct Fp p3_C4;
struct Fp p3_M_C4;
struct Fp p3_M_C_C8;
struct Fp p3_C8;
struct Fp p3_M_C_C4_C8;
struct Fp p3_C4_C8;
struct Fp p3_M_C_C4_C16;
struct Fp p3_C4_C16;
struct Fp p3_C16;
struct Fp p3_M_C16;
struct Fp p3_C_C4_C8_C16;
struct Fp p3_M_C_C4_C8_C16;
struct Fp p3_M_C_C8_C16;
struct Fp p3_C8_C16;

struct Fp p5_C4, p5_M_C4;
struct Fp p5_C8, p5_M_C8;
struct Fp p5_C_C4_C8, p5_M_C_C4_C8;
struct Fp p5_C16, p5_M_C_C16;
struct Fp p5_C4_C16, p5_C_C4_C16, p5_M_C_C4_C16;
struct Fp p5_C4_C8_C16, p5_C_C4_C8_C16, p5_M_C_C_C4_C8_C16;
struct Fp p5_C8_C16, p5_C_C8_C16, p5_M_C_C8_C16;

struct Fp p7_C4;
struct Fp p7_M_C4;
struct Fp p7_M_C_C8;
struct Fp p7_C8;
struct Fp p7_M_C_C4_C8;
struct Fp p7_C8_C4;
struct Fp p7_M_C_C4_C16;
struct Fp p7_C4_C16;
struct Fp p7_C16;
struct Fp p7_M_C16;
struct Fp p7_C_C4_C8_C16;
struct Fp p7_M_C_C4_C8_C16;
struct Fp p7_M_C_C8_C16;
struct Fp p7_C8_C16;

struct Fp p4_C4;
struct Fp p4_C8;
struct Fp p4_C16;
struct Fp p4_C4_C8, p4_C8_C16;
struct Fp p4_C4_C16, p4_C4_C8_C16;

struct Fp p8_C4;
struct Fp p8_C8;
struct Fp p8_C16;
struct Fp p8_C4_C8, p8_C8_C16;
struct Fp p8_C4_C16, p8_C4_C8_C16;

struct Fp p2_C8, p2_M_C8;
struct Fp p2_C4, p2_C4_C8;
struct Fp p2_C16, p2_M_C16, p2_C_C16, p2_M_C_C16;
struct Fp p2_C4_C8_C16;
struct Fp p2_C8_C16, p2_C_C8_C16, p2_M_C8_C16, p2_M_C_C8_C16;
struct Fp p2_C_C4_C8_C16;

struct Fp p6_C8, p6_M_C8;
struct Fp p6_C16, p6_M_C16, p6_C_C16, p6_M_C_C16;
struct Fp p6_C_C8_C16, p6_M_C_C8_C16, p6_C8_C16, p6_M_C8_C16;

struct Fp4 z_inv2;

struct scm_cost{
    unsigned long int pre_eca;
    unsigned long int pre_ecd;
    unsigned long int loop_eca;
    unsigned long int loop_ecd;
};
struct scm_cost scm_cost;
void Init_scm_cost(struct scm_cost *cost){
    cost->pre_eca=0;
    cost->pre_ecd=0;
    cost->loop_eca=0;
    cost->loop_ecd=0;
}
void Print_scm_Cost(struct scm_cost *cost,char *str){
    printf("%s\n",str);
    printf("pre-eca:%ld,pre-ecd:%ld,loop-eca:%ld,loop-ecd:%ld",cost->pre_eca,cost->pre_ecd,cost->loop_eca,cost->loop_ecd);
}

struct mpz_Cost{
    unsigned long int mpz_mpz_mul;
    unsigned long int mpz_ui_mul;
    unsigned long int sqr;
    unsigned long int basis;
    unsigned long int mpz_mpz_add;
    unsigned long int mpz_ui_add;
    unsigned long int inv;
};
struct mpz_Cost mpz_cost;
void Init_mpz_Cost(struct mpz_Cost *cost){
    cost->mpz_mpz_mul=0;
    cost->mpz_ui_mul=0;
    cost->sqr=0;
    cost->basis=0;
    cost->mpz_mpz_add=0;
    cost->mpz_ui_add=0;
    cost->inv=0;
}
void Print_mpz_Cost(struct mpz_Cost *cost,char *str){
    printf("%s\n",str);
    printf("mpz_mpz_mul:%ld,mpz_ui_mul:%ld,sqr:%ld,basis:%ld,mpz_mpz_add:%ld,mpz_ui_add:%ld,inv:%ld\n",cost->mpz_mpz_mul,cost->mpz_ui_mul,cost->sqr,cost->basis,cost->mpz_mpz_add,cost->mpz_ui_add,cost->inv);
}


mpz_t trace4,prime4;

int fp16_pow;

void check_skew_frobenius(void);

int main(void){
    mpz_init(X);
    mpz_init(X1);
    generate_X();
    mpz_init(PRIME_P);
    mpz_init(order_r);
    mpz_init(order_EFp);
    mpz_init(trace_t);
    mpz_init(a_x);
    mpz_init(EFp_total);
    mpz_init(EFp16_total);
    mpz_init(minus1_root1);
    mpz_init(minus1_root2);
    
    KSS_16_parameters();
    pre_calculate_frob_p();
    pre_calculate_frob_p2();
    pre_calculate_frob_p3();
    pre_calculate_frob_p4();
    pre_calculate_frob_p5();
    pre_calculate_frob_p6();
    pre_calculate_frob_p7();
    pre_calculate_frob_p8();
    pre_calculate_root_minus1();
//    pre_calc_vector_final_exp();
    weil();
    
    gmp_printf("X = %Zd\n",X);
    gmp_printf("p = %Zd\n",PRIME_P);
    gmp_printf("r = %Zd\n",order_r);
    gmp_printf("t = %Zd\n",trace_t);
    gmp_printf("#E(Fp) = %Zd\n",order_EFp);
    
    printf("X = %dbit\n",(int)mpz_sizeinbase(X,2));
    printf("p = %dbit\n",(int)mpz_sizeinbase(PRIME_P,2));
    printf("r = %dbit\n",(int)mpz_sizeinbase(order_r,2));
    printf("t = %dbit\n",(int)mpz_sizeinbase(trace_t,2));
    gmp_printf("y^2 = x^3 + %Zdx\n",a_x);
    
    //Check_SCM();
//    check_skew_frobenius();
//    check_Pairing();
    //check_G1_SCM();
    check_G2_SCM();
    //check_G3_EXP();
//    check_G2_order();
    
    
    mpz_clear(PRIME_P);
    mpz_clear(order_r);
    mpz_clear(order_EFp);
    mpz_clear(trace_t);
    //    mpz_clear(b);
    mpz_clear(a_x);
    dealloc_constants();
    return 0;
}

void check_skew_frobenius(void){
    
    struct EFp16 Q, tmp;
    struct EFp4 Qt, qt_tmp;
    
    EFp16_init(&Q);
    EFp16_init(&tmp);
    EFp4_init(&Qt);
    EFp4_init(&qt_tmp);
    
    EFp16_random_set_G2(&Q);
    
    EFp16_SCM_BIN(&tmp, &Q, PRIME_P);
    printf("EFp16 SCM [p]Q=");
    EFp16_printf(&tmp);
    
    printf("\n EFp4 SCM [p]Q=");
    EFp16_to_EFp4_map(&Qt, &Q);
    EFp4_SCM_BIN_Sparse(&qt_tmp, &Qt, PRIME_P);
    EFp4_to_EFp16_map(&tmp, &qt_tmp);
    EFp16_printf(&tmp);
    printf("\n");
    
    printf("skew frobenius p times");
    EFp16_to_EFp4_map(&Qt, &Q);
    EFp4_Skew_Frobenius_p(&qt_tmp, &Qt);
//    EFp4_Skew_Frobenius_p(&qt_tmp, &qt_tmp);
//    EFp4_Skew_Frobenius_p(&qt_tmp, &qt_tmp);
//    EFp4_Skew_Frobenius_p(&qt_tmp, &qt_tmp);
//    EFp4_Skew_Frobenius_p(&qt_tmp, &qt_tmp);
//    EFp4_Skew_Frobenius_p(&qt_tmp, &qt_tmp);
//    EFp4_Skew_Frobenius_p(&qt_tmp, &qt_tmp);
    EFp4_printf(&qt_tmp);
    
//    printf("skew frobenius p^3 and p^3");
//    EFp4_Skew_Frobenius_p3(&qt_tmp, &Qt);
//    EFp4_Skew_Frobenius_p4(&qt_tmp, &qt_tmp);
//    EFp4_printf(&qt_tmp);
//
//
//    printf("skew frobenius p^7 and p");
//    EFp4_Skew_Frobenius_p7(&qt_tmp, &Qt);
//    EFp4_printf(&qt_tmp);
    
    
}

//-----------------------------------------------------------------------------------------
#pragma mark util methods
void dealloc_constants(){
    Fp_clear(&C1);
    Fp_init(&C1_INV);
    Fp_clear(&p_C16);
    Fp_clear(&p_C8);
    Fp_clear(&p_C4);
    Fp_clear(&p_M_C4);
    Fp_clear(&p_M_C_C8);
    Fp_clear(&p_M_C_C16);
    Fp_clear(&p_C8_C16);
    Fp_clear(&p_C_C8_C16);
    Fp_clear(&p_M_C_C8_C16);
    Fp_clear(&p_C4_C16);
    Fp_clear(&p_C_C4_C8);
    Fp_clear(&p_M_C_C4_C8);
    Fp_clear(&p_C4_C8_C16);
    Fp_clear(&p_M_C_C_C4_C8_C16);
    Fp_clear(&p_C_C4_C8_C16);
}
void pre_calculate_frob_p(){
    mpz_init(p8p1dr);
    Fp_init(&C1);
    Fp_init(&p_C16);
    Fp_init(&p_C8);
    Fp_init(&p_C4);
    Fp_init(&p_M_C4);
    Fp_init(&p_M_C_C8);
    Fp_init(&p_M_C_C16);
    Fp_init(&p_C8_C16);
    Fp_init(&p_C_C8_C16);
    Fp_init(&p_M_C_C8_C16);
    Fp_init(&p_C4_C16);
    Fp_init(&p_C_C4_C8);
    Fp_init(&p_M_C_C4_C8);
    Fp_init(&p_C4_C8_C16);
    Fp_init(&p_M_C_C_C4_C8_C16);
    Fp_init(&p_C_C4_C8_C16);
    Fp_init(&p_C_C4_C16);
    Fp_init(&p_M_C_C4_C16);
    
    Fp_set_ui(&C1,c1);
    
    mpz_pow_ui(p8p1dr,PRIME_P,8);
    mpz_add_ui(p8p1dr,p8p1dr,1);
    mpz_tdiv_q(p8p1dr,p8p1dr,order_r);
    
    mpz_invert(C1_INV.x0,C1.x0,PRIME_P);
    mpz_cost.inv++;
    
    mpz_sub_ui(p_C16.x0,PRIME_P,13);
    mpz_tdiv_q_ui(p_C16.x0,p_C16.x0,16);
    Fp_pow(&p_C16,&C1,p_C16.x0);
    
    mpz_sub_ui(p_C8.x0,PRIME_P,5);
    mpz_tdiv_q_ui(p_C8.x0,p_C8.x0,8);
    Fp_pow(&p_C8,&C1,p_C8.x0);
    
    mpz_sub_ui(p_C4.x0,PRIME_P,1);
    mpz_tdiv_q_ui(p_C4.x0,p_C4.x0,4);
    Fp_pow(&p_C4,&C1,p_C4.x0);
    Fp_neg(&p_M_C4, &p_C4);
    
    mpz_sub(p_M_C_C8.x0,PRIME_P,C1.x0);
    Fp_mul(&p_M_C_C8, &p_C8, &p_M_C_C8);
    
    Fp_mul(&p_M_C_C16, &p_C16, &C1);
    Fp_neg(&p_M_C_C16, &p_M_C_C16);
    
    Fp_mul(&p_C4_C16, &p_C4, &p_C8);
    Fp_mul(&p_C_C4_C8, &p_C4_C16, &C1); // c.c^(p-1)/4.c^(p-5)/8
    mpz_sub(p_M_C_C4_C8.x0,PRIME_P,p_C_C4_C8.x0);
    
    Fp_mul(&p_C8_C16,&p_C16,&p_C8);// c^(p-13)/16.c^(p-5)/8
    Fp_mul(&p_C_C8_C16,&p_C8_C16,&C1);// c.c^(p-13)/16.c^(p-5)/8
    mpz_sub(p_M_C_C8_C16.x0,PRIME_P,p_C_C8_C16.x0);
    
    Fp_mul(&p_C4_C8_C16, &p_C_C8_C16, &p_C4); // c.c^(p-1)/4.c^(p-13)/16.c^(p-5)/8
    
    Fp_mul(&p_C_C4_C8_C16, &p_C4_C8_C16, &C1); // c.c.c^(p-13)/16.c^(p-5)/8
    mpz_sub(p_M_C_C_C4_C8_C16.x0,PRIME_P,p_C_C4_C8_C16.x0);
    
    Fp_mul(&p_C_C4_C16, &p_C4, &p_C16);
    Fp_mul(&p_C_C4_C16, &p_C_C4_C16, &C1);
    Fp_neg(&p_M_C_C4_C16, &p_C_C4_C16);
}
void pre_calculate_frob_p2(){
    Fp_init(&p2_C4);
    Fp_init(&p2_C4_C8);
    Fp_init(&p2_C16);
    Fp_init(&p2_M_C16);
    Fp_init(&p2_C8);
    Fp_init(&p2_M_C8);
    Fp_init(&p2_C_C16);
    Fp_init(&p2_M_C_C16);
    Fp_init(&p2_C8_C16);
    Fp_init(&p2_M_C8_C16);
    Fp_init(&p2_C_C8_C16);
    Fp_init(&p2_M_C_C8_C16);

    mpz_t p2; mpz_init(p2);
    mpz_pow_ui(p2,PRIME_P,2);
    
    mpz_sub_ui(p2_C4.x0,p2,1);
    mpz_tdiv_q_ui(p2_C4.x0,p2_C4.x0,4);
    Fp_pow(&p2_C4,&C1,p2_C4.x0); // 2^(P^2-1)/4
    
    mpz_sub_ui(p2_C8.x0,p2,1);
    mpz_tdiv_q_ui(p2_C8.x0,p2_C8.x0,8);
    Fp_pow(&p2_C8,&C1,p2_C8.x0);
    Fp_neg(&p2_M_C8,&p2_C8);
    
    Fp_mul(&p2_C4_C8, &p2_C4, &p2_C8);
    
    mpz_sub_ui(p2_C16.x0,p2,9);
    mpz_tdiv_q_ui(p2_C16.x0,p2_C16.x0,16);
    Fp_pow(&p2_C16,&C1,p2_C16.x0);
    Fp_neg(&p2_M_C16,&p2_C16);
    
    Fp_mul(&p2_C4_C8_C16, &p2_C4_C8, &p2_C16);
    Fp_mul(&p2_C_C4_C8_C16, &C1, &p2_C4_C8_C16);

    
    Fp_mul(&p2_C_C16,&p2_C16,&C1);//c1*c1^(p^2-9/16)
    Fp_neg(&p2_M_C_C16, &p2_C_C16);
    
    Fp_mul(&p2_C8_C16, &p2_C8, &p2_C16);
    Fp_mul(&p2_C_C8_C16, &p2_C8_C16, &C1);
    Fp_neg(&p2_M_C8_C16, &p2_C8_C16);
    Fp_neg(&p2_M_C_C8_C16, &p2_C_C8_C16);
    
    mpz_clear(p2);
}
void pre_calculate_frob_p3(){
    struct Fp TMP;
    Fp_init(&TMP);
    Fp_init(&p3_C4);
    Fp_init(&p3_M_C4);
    Fp_init(&p3_M_C_C8);
    Fp_init(&p3_M_C_C4_C8);
    Fp_init(&p3_C4_C8);
    Fp_init(&p3_M_C_C4_C16);
    Fp_init(&p3_C4_C16);
    Fp_init(&p3_C16);
    Fp_init(&p3_M_C16);
    Fp_init(&p3_C4);
    Fp_init(&p3_C_C4_C8_C16);
    Fp_init(&p3_M_C_C4_C8_C16);
    Fp_init(&p3_M_C_C8_C16);
    Fp_init(&p3_C8_C16);
    
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_t prime_3;
    mpz_init(prime_3);
    mpz_pow_ui(prime_3,PRIME_P,3);
    
    mpz_sub_ui(tmp,prime_3,5);
    mpz_tdiv_q_ui(tmp,tmp,16);
    Fp_pow(&p3_C16,&C1,tmp);
    Fp_neg(&p3_M_C16, &p3_C16);
    
    mpz_sub_ui(tmp,prime_3,5);
    mpz_tdiv_q_ui(tmp,tmp,8);
    Fp_pow(&p3_C8,&C1,tmp);
    
    mpz_sub_ui(tmp,prime_3,1);
    mpz_tdiv_q_ui(tmp,tmp,4);
    Fp_pow(&p3_C4,&C1,tmp);
    Fp_printf(&p3_C4);
    
    Fp_neg(&p3_M_C4, &p3_C4);
    
    Fp_mul(&p3_M_C_C8, &p3_C8, &C1);
    Fp_neg(&p3_M_C_C8, &p3_M_C_C8);
    ;
    Fp_mul(&p3_M_C_C4_C8, &p3_M_C_C8, &p3_C4);
    
    Fp_mul(&p3_C4_C8, &p3_C8, &p3_C4);
    
    Fp_mul(&p3_C8_C16, &p3_C8, &p3_C16);
    Fp_mul(&TMP, &p3_C8_C16, &C1);
    Fp_neg(&p3_M_C_C8_C16, &TMP);
    
    Fp_mul(&p3_C_C4_C8_C16, &TMP, &p3_C4);
    Fp_neg(&p3_M_C_C4_C8_C16, &p3_C_C4_C8_C16);
    
    Fp_mul(&p3_C4_C16, &p3_C4, &p3_C16);
    Fp_mul(&p3_M_C_C4_C16, &p3_C4_C16, &C1);
    Fp_neg(&p3_M_C_C4_C16, &p3_M_C_C4_C16);
    
    Fp_clear(&TMP);
    mpz_clear(tmp);
    mpz_clear(prime_3);
}
void pre_calculate_frob_p4(){
    Fp_init(&p4_C4);
    Fp_init(&p4_C8);
    Fp_init(&p4_C16);
    Fp_init(&p4_C4_C8);
    Fp_init(&p4_C8_C16);
    Fp_init(&p4_C4_C16);
    Fp_init(&p4_C4_C8_C16);
    
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_t prime_4;
    mpz_init(prime_4);
    mpz_pow_ui(prime_4,PRIME_P,4);
    
    mpz_sub_ui(tmp,prime_4,1);
    mpz_tdiv_q_ui(tmp,tmp,16);
    Fp_pow(&p4_C16,&C1,tmp);
    
    mpz_sub_ui(tmp,prime_4,1);
    mpz_tdiv_q_ui(tmp,tmp,8);
    Fp_pow(&p4_C8,&C1,tmp);
    
    mpz_sub_ui(tmp,prime_4,1);
    mpz_tdiv_q_ui(tmp,tmp,4);
    Fp_pow(&p4_C4,&C1,tmp);
    
    Fp_mul(&p4_C4_C8, &p4_C4, &p4_C8);
    Fp_mul(&p4_C4_C8_C16, &p4_C4_C8, &p4_C16);
    
    Fp_mul(&p4_C8_C16, &p4_C16, &p4_C8);
    Fp_mul(&p4_C4_C16, &p4_C16, &p4_C4);
    
    mpz_clear(tmp);
    mpz_clear(prime_4);
}
void pre_calculate_frob_p8(){
    Fp_init(&p8_C4);
    Fp_init(&p8_C8);
    Fp_init(&p8_C16);
    Fp_init(&p8_C4_C8);
    Fp_init(&p8_C8_C16);
    Fp_init(&p8_C4_C16);
    Fp_init(&p8_C4_C8_C16);
    
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_t prime_8;
    mpz_init(prime_8);
    mpz_pow_ui(prime_8,PRIME_P,8);
    
    mpz_sub_ui(tmp,prime_8,1);
    mpz_tdiv_q_ui(tmp,tmp,16);
    Fp_pow(&p8_C16,&C1,tmp);
    
    mpz_sub_ui(tmp,prime_8,1);
    mpz_tdiv_q_ui(tmp,tmp,8);
    Fp_pow(&p8_C8,&C1,tmp);
    
    mpz_sub_ui(tmp,prime_8,1);
    mpz_tdiv_q_ui(tmp,tmp,4);
    Fp_pow(&p8_C4,&C1,tmp);
    
    Fp_mul(&p8_C4_C8, &p8_C4, &p8_C8);
    Fp_mul(&p8_C4_C8_C16, &p8_C4_C8, &p8_C16);
    
    Fp_mul(&p8_C8_C16, &p8_C16, &p8_C8);
    Fp_mul(&p8_C4_C16, &p8_C16, &p8_C4);
    
    mpz_clear(tmp);
    mpz_clear(prime_8);
}
void pre_calculate_frob_p5(){
    Fp_init(&p5_C16);
    Fp_init(&p5_C8);
    Fp_init(&p5_C4);
    Fp_init(&p5_M_C4);
    Fp_init(&p5_M_C8);
    Fp_init(&p5_M_C_C16);
    Fp_init(&p5_C8_C16);
    Fp_init(&p5_C_C8_C16);
    Fp_init(&p5_M_C_C8_C16);
    Fp_init(&p5_C4_C16);
    Fp_init(&p5_C_C4_C8);
    Fp_init(&p5_M_C_C4_C8);
    Fp_init(&p5_C4_C8_C16);
    Fp_init(&p5_M_C_C_C4_C8_C16);
    Fp_init(&p5_C_C4_C8_C16);
    Fp_init(&p5_C_C4_C16);
    Fp_init(&p5_M_C_C4_C16);
    
    mpz_t prime_5;
    mpz_init(prime_5);
    mpz_pow_ui(prime_5,PRIME_P,5);
    
    mpz_sub_ui(p5_C16.x0,prime_5,13);
    mpz_tdiv_q_ui(p5_C16.x0,p5_C16.x0,16);
    Fp_pow(&p5_C16,&C1,p5_C16.x0);
    
    mpz_sub_ui(p5_C8.x0,prime_5,5);
    mpz_tdiv_q_ui(p5_C8.x0,p5_C8.x0,8);
    Fp_pow(&p5_C8,&C1,p5_C8.x0);
    
    mpz_sub_ui(p5_C4.x0,prime_5,1);
    mpz_tdiv_q_ui(p5_C4.x0,p5_C4.x0,4);
    Fp_pow(&p5_C4,&C1,p5_C4.x0);
    Fp_neg(&p5_M_C4, &p5_C4);
    
    mpz_sub(p5_M_C8.x0,PRIME_P,C1.x0);
    Fp_mul(&p5_M_C8, &p5_C8, &p5_M_C8);
    
    Fp_mul(&p5_M_C_C16, &p5_C16, &C1);
    Fp_neg(&p5_M_C_C16, &p5_M_C_C16);
    
    Fp_mul(&p5_C4_C16, &p5_C4, &p5_C8);
    Fp_mul(&p5_C_C4_C8, &p5_C4_C16, &C1); // c.c^(p-1)/4.c^(p-5)/8
    mpz_sub(p5_M_C_C4_C8.x0,PRIME_P,p5_C_C4_C8.x0);
    
    Fp_mul(&p5_C8_C16,&p5_C16,&p5_C8);// c^(p-13)/16.c^(p-5)/8
    Fp_mul(&p5_C_C8_C16,&p5_C8_C16,&C1);// c.c^(p-13)/16.c^(p-5)/8
    mpz_sub(p5_M_C_C8_C16.x0,PRIME_P,p5_C_C8_C16.x0);
    
    Fp_mul(&p5_C4_C8_C16, &p5_C_C8_C16, &p5_C4); // c.c^(p-1)/4.c^(p-13)/16.c^(p-5)/8
    
    Fp_mul(&p5_C_C4_C8_C16, &p5_C4_C8_C16, &C1); // c.c.c^(p-13)/16.c^(p-5)/8
    Fp_neg(&p5_M_C_C_C4_C8_C16, &p5_C_C4_C8_C16);
    
    Fp_mul(&p5_C_C4_C16, &p5_C4, &p5_C16);
    Fp_mul(&p5_C_C4_C16, &p5_C_C4_C16, &C1);
    Fp_neg(&p5_M_C_C4_C16, &p5_C_C4_C16);
}
void pre_calculate_frob_p6()
{
    Fp_init(&C1);
    Fp_init(&p6_C16);
    Fp_init(&p6_M_C16);
    Fp_init(&p6_C8);
    Fp_init(&p6_M_C8);
    Fp_init(&p6_C_C16);
    Fp_init(&p6_M_C_C16);
    Fp_init(&p6_C_C8_C16);
    Fp_init(&p6_M_C_C8_C16);
    Fp_init(&p6_C8_C16);
    Fp_init(&p6_M_C8_C16);
    
    Fp_set_ui(&C1,c1);
    
    mpz_pow_ui(p6_C16.x0,PRIME_P,6);
    mpz_sub_ui(p6_C16.x0,p6_C16.x0,9);
    mpz_tdiv_q_ui(p6_C16.x0,p6_C16.x0,16);
    Fp_pow(&p6_C16,&C1,p6_C16.x0);
    Fp_neg(&p6_M_C16,&p6_C16);
    
    mpz_pow_ui(p6_C8.x0,PRIME_P,6);
    mpz_sub_ui(p6_C8.x0,p6_C8.x0,1);
    mpz_tdiv_q_ui(p6_C8.x0,p6_C8.x0,8);
    Fp_pow(&p6_C8,&C1,p6_C8.x0);
    Fp_neg(&p6_M_C8,&p6_C8);
    
    Fp_mul(&p6_C_C16,&p6_C16,&C1);//c1*c1^(p^2-9/16)
    Fp_neg(&p6_M_C_C16, &p6_C_C16);
    
    Fp_mul(&p6_C8_C16, &p6_C8, &p6_C16);
    Fp_mul(&p6_C_C8_C16, &p6_C8_C16, &C1);
    Fp_neg(&p6_M_C8_C16, &p6_C8_C16);
    Fp_neg(&p6_M_C_C8_C16, &p6_C_C8_C16);
}
void pre_calculate_frob_p7(){
    struct Fp TMP;
    Fp_init(&TMP);
    Fp_init(&p7_C4);
    Fp_init(&p7_M_C4);
    Fp_init(&p7_M_C_C8);
    Fp_init(&p7_M_C_C4_C8);
    Fp_init(&p7_C8_C4);
    Fp_init(&p7_M_C_C4_C16);
    Fp_init(&p7_C4_C16);
    Fp_init(&p7_C16);
    Fp_init(&p7_M_C16);
    Fp_init(&p7_C4);
    Fp_init(&p7_C_C4_C8_C16);
    Fp_init(&p7_M_C_C4_C8_C16);
    Fp_init(&p7_M_C_C8_C16);
    Fp_init(&p7_C8_C16);
    
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_t prime_7;
    mpz_init(prime_7);
    mpz_pow_ui(prime_7,PRIME_P,7);
    
    mpz_sub_ui(tmp,prime_7,5);
    mpz_tdiv_q_ui(tmp,tmp,16);
    Fp_pow(&p7_C16,&C1,tmp);
    Fp_neg(&p7_M_C16, &p7_C16);
    
    mpz_sub_ui(tmp,prime_7,5);
    mpz_tdiv_q_ui(tmp,tmp,8);
    Fp_pow(&p7_C8,&C1,tmp);
    
    mpz_sub_ui(tmp,prime_7,1);
    mpz_tdiv_q_ui(tmp,tmp,4);
    Fp_pow(&p7_C4,&C1,tmp);
    Fp_printf(&p7_C4);
    
    Fp_neg(&p7_M_C4, &p7_C4);
    
    Fp_mul(&p7_M_C_C8, &p7_C8, &C1);
    Fp_neg(&p7_M_C_C8, &p7_M_C_C8);
    ;
    Fp_mul(&p7_M_C_C4_C8, &p7_M_C_C8, &p7_C4);
    
    Fp_mul(&p7_C8_C4, &p7_C8, &p7_C4);
    
    Fp_mul(&p7_C8_C16, &p7_C8, &p7_C16);
    Fp_mul(&TMP, &p7_C8_C16, &C1);
    Fp_neg(&p7_M_C_C8_C16, &TMP);
    
    Fp_mul(&p7_C_C4_C8_C16, &TMP, &p7_C4);
    Fp_neg(&p7_M_C_C4_C8_C16, &p7_C_C4_C8_C16);
    
    Fp_mul(&p7_C4_C16, &p7_C4, &p7_C16);
    Fp_mul(&p7_M_C_C4_C16, &p7_C4_C16, &C1);
    Fp_neg(&p7_M_C_C4_C16, &p7_M_C_C4_C16);
    
    Fp_clear(&TMP);
    mpz_clear(tmp);
    mpz_clear(prime_7);
}
void generate_X(){
    //c1 = 2
    // 2^ -2^32-2^18+2^8+1
    X_bit_binary[35]=1;
    X_bit_binary[32]=-1;
    X_bit_binary[18]=-1;
    X_bit_binary[8]=1;
    X_bit_binary[0]=1;
    
    X1_bit_binary[35]=1;
    X1_bit_binary[32]=-1;
    X1_bit_binary[18]=-1;
    X1_bit_binary[8]=1;
    X1_bit_binary[1]=1;
    
    // -2^34-2^32-2^13-2^13+2^6+1
//    X_bit_binary[34]=-1;
//    X_bit_binary[32]=-1;
//    X_bit_binary[13]=-1;
//    X_bit_binary[11]=-1;
//    X_bit_binary[6]=1;
//    X_bit_binary[0]=1;
//
//    X1_bit_binary[34]=-1;
//    X1_bit_binary[32]=-1;
//    X1_bit_binary[13]=-1;
//    X1_bit_binary[11]=-1;
//    X1_bit_binary[6]=1;
//    X1_bit_binary[1]=1;
    
    mpz_t tmp,set_2;
    mpz_init(tmp);
    mpz_init(set_2);
    mpz_set_ui(set_2,2);
    
    int i;
    for(i=x_bit;i>=0;i--){
//        printf("%d",X_bit_binary[i]);
        if(X_bit_binary[i]==1){
            mpz_pow_ui(tmp,set_2,i);
            mpz_add(X,X,tmp);
        }else if(X_bit_binary[i]==-1){
            mpz_pow_ui(tmp,set_2,i);
            mpz_sub(X,X,tmp);
        }
    }
    printf("\n");
    for(i=x_bit;i>=0;i--){
//        printf("%d",X1_bit_binary[i]);
        if(X1_bit_binary[i]==1){
            mpz_pow_ui(tmp,set_2,i);
            mpz_add(X1,X1,tmp);
        }else if(X1_bit_binary[i]==-1){
            mpz_pow_ui(tmp,set_2,i);
            mpz_sub(X1,X1,tmp);
        }
    }
    
//    printf("\n");
//    mpz_out_str(stdout,10,X);
//    printf("\n");
//    mpz_out_str(stdout,10,X1);
//    printf("\n");
    return;
}
void pre_calculate_root_minus1(){
    struct Fp buf;
    Fp_init(&buf);
    mpz_sub_ui(buf.x0,PRIME_P,1);
    if(mpz_legendre(buf.x0,PRIME_P)==1){
        Fp_sqrt(&buf,&buf);
    }
    mpz_set(minus1_root1,buf.x0);
    mpz_sub(minus1_root2,PRIME_P,minus1_root1);
    
    Fp_clear(&buf);
}
void weil(){
    mpz_init(trace4);
    mpz_init(prime4);
    mpz_t t2,t4,t8,t16,p_exp_2,p_exp_4,p_exp_8,buf;
    mpz_init(t2);
    mpz_init(t4);
    mpz_init(t8);
    mpz_init(t16);
    mpz_init(p_exp_2);
    mpz_init(p_exp_4);
    mpz_init(p_exp_8);
    mpz_init(buf);
    
    //EFp_total
    mpz_add_ui(buf,PRIME_P,1);
    mpz_sub(EFp_total,buf,trace_t);
    //t2←α^2+β^2
    mpz_pow_ui(t2,trace_t,2);
    mpz_mul_ui(buf,PRIME_P,2);
    mpz_sub(t2,t2,buf);
    mpz_pow_ui(p_exp_2,PRIME_P,2);
    //α^4+β^4
    mpz_pow_ui(t4,t2,2);
    mpz_sub(t4,t4,p_exp_2);
    mpz_sub(t4,t4,p_exp_2);
    mpz_pow_ui(p_exp_4,p_exp_2,2);
    //α^8+β^8
    mpz_pow_ui(t8,t4,2);
    mpz_sub(t8,t8,p_exp_4);
    mpz_sub(t8,t8,p_exp_4);
    mpz_pow_ui(p_exp_8,p_exp_4,2);
    //α^16+β^16
    mpz_pow_ui(t16,t8,2);
    mpz_sub(t16,t16,p_exp_8);
    mpz_sub(t16,t16,p_exp_8);
    //EFp16_total
    mpz_pow_ui(buf,p_exp_8,2);
    mpz_sub(buf,buf,t16);
    mpz_add_ui(EFp16_total,buf,1);
    
    mpz_set(trace4,t4);
    mpz_set(prime4,p_exp_4);
    
    mpz_clear(t2);
    mpz_clear(t4);
    mpz_clear(t8);
    mpz_clear(t16);
    mpz_clear(p_exp_2);
    mpz_clear(p_exp_4);
    mpz_clear(p_exp_8);
    mpz_clear(buf);
}

//-----------------------------------------------------------------------------------------
// #pragma mark Fp methods
void Fp_init(struct Fp *A){
    mpz_init(A->x0);
}
void Fp_set(struct Fp *ANS,struct Fp *E){
    mpz_set(ANS->x0,E->x0);
}
void Fp_set_ui(struct Fp *A,signed long int B){
    mpz_set_ui(A->x0,B);
}
void Fp_set_mpz(struct Fp *A, mpz_t a)
{
    mpz_set(A->x0,a);
}
void Fp_random(struct Fp *A){
    mpz_random(A->x0,10);
    mpz_mod(A->x0,A->x0,PRIME_P);
}
void Fp_clear(struct Fp *A){
    mpz_clear(A->x0);
}
void Fp_printf(struct Fp *A){
    gmp_printf("%Zd\n",A->x0);
}
void Fp_add(struct Fp *ANS,struct Fp *A,struct Fp *B){
    mpz_add(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
    Big_add++;
    mpz_cost.mpz_mpz_add++;
}
void Fp_add_ui(struct Fp *ANS,struct Fp *A,unsigned long int B){
    mpz_add_ui(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
    Big_add++;
    mpz_cost.mpz_ui_add++;
}
void Fp_add_mpz(struct Fp *ANS,struct Fp *A,mpz_t B){
    mpz_add(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
    Big_add++;
    mpz_cost.mpz_mpz_add++;
}
void Fp_sub(struct Fp *ANS,struct Fp *A,struct Fp *B){
    mpz_sub(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
    Big_add++;
    mpz_cost.mpz_mpz_add++;
}
void Fp_sub_ui(struct Fp *ANS,struct Fp *A,unsigned long int B){
    mpz_sub_ui(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
    Big_add++;
    mpz_cost.mpz_ui_add++;
}
void Fp_sqr(struct Fp *ANS,struct Fp *A){
    mpz_mul(ANS->x0,A->x0,A->x0);//A^2: GMP mpz_mul take care of squaring
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
    Sqr++;
    mpz_cost.sqr++;
}
void Fp_mul(struct Fp *ANS,struct Fp *A,struct Fp *B){
    
    if (Fp_cmp(A, B) == 0) {
        Fp_sqr(ANS, A);
    }
    else{
        if (mpz_cmp_ui(B->x0, 2)==0) {
            mpz_mul_2exp(ANS->x0,A->x0,1);
            Small_m++;
            mpz_cost.basis++;
        }
        else{
            mpz_mul(ANS->x0,A->x0,B->x0);
            Big_M++;
            mpz_cost.mpz_mpz_mul++;
        }
        mpz_mod(ANS->x0,ANS->x0,PRIME_P);
    }
}
void Fp_mul_mpz(struct Fp *ANS,struct Fp *A,mpz_t B){
    mpz_mul(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
    Big_M++;
    if(mpz_cmp(A->x0,B)==0){
        mpz_cost.sqr++;
    }else{
        mpz_cost.mpz_mpz_mul++;
    }
}
void Fp_mul_ui(struct Fp *ANS,struct Fp *A,unsigned long int B){
    if (B == 2) {
        mpz_mul_2exp(ANS->x0,A->x0,1);
        Small_m++;
        mpz_cost.basis++;
    }
    else{
        mpz_mul_ui(ANS->x0,A->x0,B);
        Big_M++;
        mpz_cost.mpz_ui_mul++;
    }
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
}
void Fp_div(struct Fp *ANS,struct Fp *A,struct Fp *B){
    mpz_cost.inv++;
    mpz_cost.mpz_mpz_mul++;
    mpz_invert(ANS->x0,B->x0,PRIME_P);
    mpz_mul(ANS->x0,A->x0,ANS->x0);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
}
void Fp_pow(struct Fp *ANS,struct Fp *A,mpz_t j){
    struct Fp tmp;
    Fp_init(&tmp);
    Fp_set(&tmp,A);
    int i = 0;
    int bit_length = (int)mpz_sizeinbase(j,2);
    for(i = bit_length - 2; i >= 0; i--){
        if(mpz_tstbit(j,i) == 1){
            Fp_mul(&tmp,&tmp,&tmp);//a*2
            Fp_mul(&tmp,&tmp,A);//*a
        }
        else{
            Fp_mul(&tmp,&tmp,&tmp);//a*2
        }
    }
    Fp_set(ANS,&tmp);
    Fp_clear(&tmp);
}
void Fp_sqrt(struct Fp *ANS,struct Fp *A){
    struct Fp n_tmp,y_tmp,x_tmp,b_tmp,t_tmp,tmp_Fp;
    Fp_init(&n_tmp);
    Fp_init(&y_tmp);
    Fp_init(&x_tmp);
    Fp_init(&b_tmp);
    Fp_init(&t_tmp);
    Fp_init(&tmp_Fp);
    
    Fp_set(&n_tmp,A);
    
    mpz_t tmp_mpz,q_tmp,e_tmp,r_tmp,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q_tmp);
    mpz_init(e_tmp);
    mpz_init(r_tmp);
    mpz_init(set_1);
    mpz_init(set_2);
    
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    while(mpz_legendre(n_tmp.x0,PRIME_P) != -1){
        Fp_add_ui(&n_tmp,&n_tmp,1);
    }
    
    mpz_set(q_tmp,PRIME_P);
    mpz_sub_ui(q_tmp,q_tmp,1);
    mpz_set_ui(e_tmp,0);
    
    while(mpz_odd_p(q_tmp)==0){
        mpz_add_ui(e_tmp,e_tmp,1);
        mpz_div_ui(q_tmp,q_tmp,2);
    }
    
    Fp_pow(&y_tmp,&n_tmp,q_tmp);
    
    mpz_set(r_tmp,e_tmp);
    mpz_sub_ui(tmp_mpz,q_tmp,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    Fp_pow(&x_tmp,A,tmp_mpz);
    Fp_pow(&tmp_Fp,&x_tmp,set_2);
    Fp_mul(&b_tmp,&tmp_Fp,A);
    Fp_mul(&x_tmp,&x_tmp,A);
    
    int m =0;
    while(Fp_cmp_mpz(&b_tmp,set_1)==1){
        m=1;
        while(Fp_cmp_mpz(&tmp_Fp,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp_pow(&tmp_Fp,&b_tmp,tmp_mpz);
        }
        mpz_sub_ui(tmp_mpz,r_tmp,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,PRIME_P);
        Fp_pow(&t_tmp,&y_tmp,tmp_mpz);
        Fp_pow(&y_tmp,&t_tmp,set_2);
        mpz_set_ui(r_tmp,m);
        Fp_mul(&x_tmp,&x_tmp,&t_tmp);
        Fp_mul(&b_tmp,&b_tmp,&y_tmp);
        Fp_set(&tmp_Fp, &b_tmp);
    }
    Fp_set(ANS,&x_tmp);
    
    Fp_clear(&n_tmp);
    Fp_clear(&y_tmp);
    Fp_clear(&x_tmp);
    Fp_clear(&b_tmp);
    Fp_clear(&t_tmp);
    Fp_clear(&tmp_Fp);
    mpz_clear(tmp_mpz);
    mpz_clear(q_tmp);
    mpz_clear(e_tmp);
    mpz_clear(r_tmp);
    mpz_clear(set_1);
    mpz_clear(set_2);
}
void Fp_neg(struct Fp *ANS,struct Fp *A){
//    if (mpz_cmp_ui(A->x0,0) != 0) {
//        mpz_sub(ANS->x0,PRIME_P,A->x0);
//    }
//    else{
//        Fp_set(ANS,A);
//    }
    
    mpz_sub(ANS->x0,PRIME_P,A->x0);
}
int Fp_cmp(struct Fp *A,struct Fp *B){
    if(mpz_cmp(A->x0,B->x0) == 0){
        return 0;
    }
    return 1;
}
int Fp_cmp_mpz(struct Fp *A,mpz_t B){
    if(mpz_cmp(A->x0,B)==0){
        return 0;
    }
    return 1;
}

//-----------------------------------------------------------------------------------------
void Fp2_init(struct Fp2 *A){
    Fp_init(&A->x0);
    Fp_init(&A->x1);
}
void Fp2_set(struct Fp2 *ANS,struct Fp2 *A){
    Fp_set(&ANS->x0,&A->x0);
    Fp_set(&ANS->x1,&A->x1);
}
void Fp2_set_ui(struct Fp2 *A,signed long int B){
    Fp_set_ui(&A->x0,B);
    Fp_set_ui(&A->x1,B);
}
void Fp2_random(struct Fp2 *A){
    Fp_random(&A->x0);
    Fp_random(&A->x1);
}
void Fp2_clear(struct Fp2 *A){
    Fp_clear(&A->x0);
    Fp_clear(&A->x1);
}
void Fp2_printf(struct Fp2 *A){
    gmp_printf("%Zd,%Zd\n",A->x0.x0,A->x1.x0);
}
void Fp2_add(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    Fp_add(&ANS->x0,&A->x0,&B->x0);
    Fp_add(&ANS->x1,&A->x1,&B->x1);
}
void Fp2_add_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B){
    Fp_add_ui(&ANS->x0,&A->x0,B);
    Fp_add_ui(&ANS->x1,&A->x1,B);
}
void Fp2_sub(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    Fp_sub(&ANS->x0,&A->x0,&B->x0);
    Fp_sub(&ANS->x1,&A->x1,&B->x1);
}
//TODO: Optimize temp variables
void Fp2_sqr(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp tmp1,tmp2,tmp3;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp3);
    
    //Cost = 2M + 4Ap + 3m in Fp
    Fp_mul(&tmp1,&A->x0,&A->x1);//t=a0a1
    Fp_add(&tmp2,&A->x0,&A->x1);//(a0+a1)
    Fp_mul_ui(&tmp3,&A->x1,c1);//a1*i
    Fp_add(&tmp3,&tmp3,&A->x0);//(a0+a1*i)
    Fp_mul(&ANS->x0,&tmp2,&tmp3);// (a0+a1)(a0+a1*i)
    Fp_sub(&ANS->x0, &ANS->x0,&tmp1);
    Fp_mul_ui(&tmp2,&tmp1,c1);//
    
    Fp_sub(&ANS->x0,&ANS->x0,&tmp2);//
    //    Fp_mul_ui(&ANS->x1,&tmp1,c1);//
    Fp_set(&ANS->x1,&tmp2);
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp3);
}
void Fp2_mul(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    if (Fp2_cmp(A, B) == 0) {
        Fp2_sqr(ANS, A);
    }
    else{
        struct Fp tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
        Fp_init(&tmp1);
        Fp_init(&tmp2);
        Fp_init(&tmp3);
        Fp_init(&tmp4);
        Fp_init(&tmp5);
        Fp_init(&tmp6);
        
        struct Fp2 t_ans;
        Fp2_init(&t_ans);
        //3Mp+5Ap+1M\alpha
        Fp_mul(&tmp1,&A->x0,&B->x0);//a*c
        Fp_mul(&tmp2,&A->x1,&B->x1);//b*d
        Fp_mul_ui(&tmp3,&tmp2,c1);//b*d*v
        Fp_add(&t_ans.x0,&tmp1,&tmp3);//a*c+b*d*v
        Fp_add(&tmp4,&A->x0,&A->x1);//a+b
        Fp_add(&tmp5,&B->x0,&B->x1);//c+d
        Fp_mul(&tmp6,&tmp4,&tmp5);//(a+b)(c+d)
        Fp_sub(&t_ans.x1,&tmp6,&tmp1);
        Fp_sub(&t_ans.x1,&t_ans.x1,&tmp2);
        
        Fp2_set(ANS,&t_ans);
        
        Fp_clear(&tmp1);
        Fp_clear(&tmp2);
        Fp_clear(&tmp3);
        Fp_clear(&tmp4);
        Fp_clear(&tmp5);
        Fp_clear(&tmp6);
        Fp2_clear(&t_ans);
    }
}
void Fp2_mul_basis(struct Fp2 *ANS,struct Fp2 *A){
    //(a,b)(1,1)=(a-b,a+b)
    struct Fp2 tmp;
    Fp2_init(&tmp);
    
    Fp_mul_ui(&tmp.x0, &A->x1, c1);
    Fp_set(&tmp.x1,&A->x0);
    
    Fp2_set(ANS,&tmp);
    Fp2_clear(&tmp);
}
void Fp2_mul_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B){
    Fp_mul_ui(&ANS->x0,&A->x0,B);
    Fp_mul_ui(&ANS->x1,&A->x1,B);
}
void Fp2_mul_Fp(struct Fp2 *ANS,struct Fp2 *A,struct Fp *B){
    Fp_mul(&ANS->x0,&A->x0,B);
    Fp_mul(&ANS->x1,&A->x1,B);
}
void Fp2_mul_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t B){
    Fp_mul_mpz(&ANS->x0,&A->x0,B);
    Fp_mul_mpz(&ANS->x1,&A->x1,B);
}
void Fp2_invert(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    
    // tmp=A^(p)=(x0,-x1)
    Fp_set(&tmp.x0,&A->x0);
    Fp_neg(&tmp.x1,&A->x1);
    
    struct Fp c,a,b;
    Fp_init(&c);
    Fp_init(&a);
    Fp_init(&b);
    
    //1Ip+4Mp+1Ap+1m
    Fp_mul(&a,&A->x0,&A->x0); // a=x0^2
    Fp_mul(&b,&A->x1,&A->x1); // b=x1^2
    Fp_mul_ui(&b,&b,c1); // b=x1^2*v
    Fp_sub(&c,&a,&b); // c=x0^2-x1^2*v mod q
    
    mpz_invert(c.x0,c.x0,PRIME_P);
     mpz_cost.inv++;
    
    // ANS=A^{-1}=(c)^{-1}*A^(p) A which c is Fp2-element and tmp is a vector A Fp2
    Fp_mul(&tmp.x0,&tmp.x0,&c);
    Fp_mul(&tmp.x1,&tmp.x1,&c);
    
    Fp2_set(ANS,&tmp);
    
    Fp_clear(&c);
    Fp_clear(&a);
    Fp_clear(&b);
    Fp2_clear(&tmp);
}
void Fp2_div(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    struct Fp2 tmp,t_ans;
    Fp2_init(&tmp);
    Fp2_init(&t_ans);
    //1Ip^2=1Ip+4Mp+1Ap+1m + 1Mp^2=3Mp+5Ap+1M\alpha
    //1Ip+7Mp+6Ap+2m
    Fp2_invert(&tmp,B);
    Fp2_mul(&t_ans,A,&tmp);
    
    Fp2_set(ANS,&t_ans);
    
    Fp2_clear(&tmp);
    Fp2_clear(&t_ans);
}
void Fp2_pow(struct Fp2 *ANS,struct Fp2 *A,mpz_t B){
    int i;
    int r;//bit数
    r= (int)mpz_sizeinbase(B,2);
    // printf("r= %d\n",r);
    
    struct Fp2 answer_tmp;
    Fp2_init(&answer_tmp);
    Fp2_set(&answer_tmp,A);
    
    struct Fp2 in_tmp;
    Fp2_init(&in_tmp);
    Fp2_set(&in_tmp,A);
    
    for(i=r-2;i>=0;i--){
        if(mpz_tstbit(B,i)==1){
            Fp2_mul(&answer_tmp,&answer_tmp,&answer_tmp);//a*2
            Fp2_mul(&answer_tmp,&answer_tmp,&in_tmp);//*a
        }else{
            Fp2_mul(&answer_tmp,&answer_tmp,&answer_tmp);//a*2
        }
    }
    
    Fp2_set(ANS,&answer_tmp);
    
    Fp2_clear(&answer_tmp);
    Fp2_clear(&in_tmp);
}
void Fp2_sqrt(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp2 n,y,x,b,t,tmp_Fp2;
    Fp2_init(&n);
    Fp2_init(&y);
    Fp2_init(&x);
    Fp2_init(&b);
    Fp2_init(&t);
    Fp2_init(&tmp_Fp2);
    Fp2_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    
    while(Fp2_legendre(&n)!=-1){
        Fp2_random(&n);
    }
    
    mpz_pow_ui(q,PRIME_P,2);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    
    Fp2_pow(&y,&n,q);
    mpz_set(r,e);
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    Fp2_pow(&x,A,tmp_mpz);
    Fp2_pow(&tmp_Fp2,&x,set_2);
    Fp2_mul(&b,&tmp_Fp2,A);
    Fp2_mul(&x,&x,A);
    
    int m;
    
    while(Fp2_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp2_set(&tmp_Fp2,&b);
        while(Fp2_cmp_mpz(&tmp_Fp2,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp2_pow(&tmp_Fp2,&b,tmp_mpz);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,PRIME_P);
        Fp2_pow(&t,&y,tmp_mpz);
        Fp2_pow(&y,&t,set_2);
        mpz_set_ui(r,m);
        Fp2_mul(&x,&x,&t);
        Fp2_mul(&b,&b,&y);
    }
    
    Fp2_set(ANS,&x);
    
    Fp2_clear(&n);
    Fp2_clear(&y);
    Fp2_clear(&x);
    Fp2_clear(&b);
    Fp2_clear(&t);
    Fp2_clear(&tmp_Fp2);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp2_cmp(struct Fp2 *A,struct Fp2 *B){
    if(Fp_cmp(&A->x0,&B->x0)==0 && Fp_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp2_cmp_mpz(struct Fp2 *A,mpz_t B){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    if(Fp_cmp_mpz(&A->x0,B)==0 && Fp_cmp(&A->x1,&tmp.x1)==0){
        Fp2_clear(&tmp);
        return 0;
    }
    Fp2_clear(&tmp);
    return 1;
}
int Fp2_legendre(struct Fp2 *a){
    mpz_t i;
    struct Fp2 tmp;
    Fp2_init(&tmp);
    mpz_init(i);
    
    mpz_pow_ui(i,PRIME_P,2);
    mpz_sub_ui(i,i,1);
    mpz_div_ui(i,i,2);
    
    Fp2_pow(&tmp,a,i);
    
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    
    if((Fp2_cmp_mpz(&tmp,cmp))==0){
        Fp2_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp2_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
void Fp2_neg(struct Fp2 *ans,struct Fp2 *a){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    
    Fp_neg(&tmp.x0,&a->x0);
    Fp_neg(&tmp.x1,&a->x1);
    
    Fp2_set(ans,&tmp);
    
    Fp2_clear(&tmp);
}

void Fp2_frobenius_map(struct Fp2 *ANS, struct Fp2 *A){
    struct Fp2 t_ans;
    Fp2_init(&t_ans);
    
    Fp_set(&t_ans.x0,&A->x0);
    if (mpz_cmp_ui(A->x1.x0,0)==0) {
        Fp_set(&t_ans.x1,&A->x1);
    }
    else{
        mpz_sub(t_ans.x1.x0,PRIME_P,A->x1.x0);
    }
    
    
    Fp2_set(ANS,&t_ans);
    
    Fp2_clear(&t_ans);
}

//-----------------------------------------------------------------------------------------
#pragma mark Fp4 methods
void Fp4_init(struct Fp4 *A){
    Fp2_init(&A->x0);
    Fp2_init(&A->x1);
}
void Fp4_set(struct Fp4 *ANS,struct Fp4 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_set(&ANS->x1,&A->x1);
}
void Fp4_set_ui(struct Fp4 *A,signed long int B){
    Fp2_set_ui(&A->x0,B);
    Fp2_set_ui(&A->x1,B);
}
void Fp4_random(struct Fp4 *A){
    Fp2_random(&A->x0);
    Fp2_random(&A->x1);
}
void Fp4_clear(struct Fp4 *A){
    Fp2_clear(&A->x0);
    Fp2_clear(&A->x1);
}
void Fp4_printf(struct Fp4 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd\n",A->x0.x0.x0,A->x0.x1.x0,A->x1.x0.x0,A->x1.x1.x0);
}
void Fp4_add(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_add(&tmp.x0,&A->x0,&B->x0);
    Fp2_add(&tmp.x1,&A->x1,&B->x1);
    
    Fp4_set(ANS,&tmp);
    
    Fp4_clear(&tmp);
}
void Fp4_add_ui(struct Fp4 *ANS,struct Fp4 *A,unsigned long int B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_add_ui(&tmp.x0,&A->x0,B);
    Fp2_add_ui(&tmp.x1,&A->x1,B);
    
    Fp4_set(ANS,&tmp);
    
    Fp4_clear(&tmp);
}
void Fp4_sub(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_sub(&tmp.x0,&A->x0,&B->x0);
    Fp2_sub(&tmp.x1,&A->x1,&B->x1);
    
    Fp4_set(ANS,&tmp);
    
    Fp4_clear(&tmp);
}
void Fp4_sqr(struct Fp4 *ANS,struct Fp4 *A){
    
    struct Fp2 tmp1,tmp2,tmp3;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    Fp2_init(&tmp3);
    
    struct Fp4 t_ans;
    Fp4_init(&t_ans);
    
    //    printf("BIG M %d\n",Big_M);
    //Cost = 2M + 4Ap + 3m in Fp
    Fp2_mul(&tmp1,&A->x0,&A->x1);//t=a0a1
    Fp2_add(&tmp2,&A->x0,&A->x1);//(a0+a1)
    Fp2_mul_basis(&tmp3,&A->x1);//a1*i
    Fp2_add(&tmp3,&tmp3,&A->x0);//(a0+a1*i)
    Fp2_mul(&t_ans.x0,&tmp2,&tmp3);// (a0+a1)(a0+a1*i)
    Fp2_sub(&t_ans.x0, &t_ans.x0,&tmp1);
    Fp2_mul_basis(&tmp2,&tmp1);//
    Fp2_sub(&t_ans.x0,&t_ans.x0,&tmp2);//
    
    Fp2_mul_ui(&t_ans.x1,&tmp1,c1);//
    
    Fp4_set(ANS,&t_ans);
    
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
    Fp2_clear(&tmp3);
    Fp4_clear(&t_ans);
    
}
void Fp4_mul(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B){
    if (Fp4_cmp(A, B) == 0) {
        Fp4_sqr(ANS, A);
    }
    else{
        //x^2-v=0
        struct Fp2 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
        Fp2_init(&tmp1);
        Fp2_init(&tmp2);
        Fp2_init(&tmp3);
        Fp2_init(&tmp4);
        Fp2_init(&tmp5);
        Fp2_init(&tmp6);
        
        struct Fp4 t_ans;
        Fp4_init(&t_ans);
        //3Mp2+5Ap2+1M\beta
        Fp2_mul(&tmp1,&A->x0,&B->x0);//a*c
        Fp2_mul(&tmp2,&A->x1,&B->x1);//b*d
        Fp2_mul_basis(&tmp3,&tmp2);//b*d*v
        Fp2_add(&t_ans.x0,&tmp1,&tmp3);//a*c+b*d*v
        Fp2_add(&tmp4,&A->x0,&A->x1);//a+b
        Fp2_add(&tmp5,&B->x0,&B->x1);//c+d
        Fp2_mul(&tmp6,&tmp4,&tmp5);//(a+b)(c+d)
        Fp2_sub(&t_ans.x1,&tmp6,&tmp1);
        Fp2_sub(&t_ans.x1,&t_ans.x1,&tmp2);
        
        Fp4_set(ANS,&t_ans);
        
        Fp2_clear(&tmp1);
        Fp2_clear(&tmp2);
        Fp2_clear(&tmp3);
        Fp2_clear(&tmp4);
        Fp2_clear(&tmp5);
        Fp2_clear(&tmp6);
        Fp4_clear(&t_ans);
    }
}
void Fp4_mul_v(struct Fp4 *ANS,struct Fp4 *A){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_mul_basis(&tmp.x0,&A->x1);
    Fp2_set(&tmp.x1,&A->x0);
    
    Fp4_set(ANS,&tmp);
    Fp4_clear(&tmp);
}
void Fp4_mul_ui(struct Fp4 *ANS,struct Fp4 *A,unsigned long int B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_mul_ui(&tmp.x0,&A->x0,B);
    Fp2_mul_ui(&tmp.x1,&A->x1,B);
    
    Fp4_set(ANS,&tmp);
    
    Fp4_clear(&tmp);
}
void Fp4_mul_mpz(struct Fp4 *ANS,struct Fp4 *A,mpz_t B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_mul_mpz(&tmp.x0,&A->x0,B);
    Fp2_mul_mpz(&tmp.x1,&A->x1,B);
    
    Fp4_set(ANS,&tmp);
    
    Fp4_clear(&tmp);
}
void Fp4_mul_Fp(struct Fp4 *ANS,struct Fp4 *A,struct Fp *B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_mul_Fp(&tmp.x0,&A->x0,B);
    Fp2_mul_Fp(&tmp.x1,&A->x1,B);
    
    Fp4_set(ANS,&tmp);
    
    Fp4_clear(&tmp);
}
void Fp4_invert(struct Fp4 *ANS,struct Fp4 *A){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    // tmp=A^(p^2)=(x0,-x1)
    Fp2_set(&tmp.x0,&A->x0);
    Fp2_neg(&tmp.x1,&A->x1);
    
    struct Fp2 c,a,b;
    Fp2_init(&c);
    Fp2_init(&a);
    Fp2_init(&b);
    //Ip^2+ 4Mp2+1Ap2+1m\beta
    Fp2_mul(&a,&A->x0,&A->x0); // a=x0^2
    Fp2_mul(&b,&A->x1,&A->x1); // b=x1^2
    Fp2_mul_basis(&b,&b); // b=x1^2*v
    Fp2_sub(&c,&a,&b); // c=x0^2-x1^2*v mod q
    
    Fp2_invert(&c,&c);
    
    // ANS=A^{-1}=(c)^{-1}*A^(p^2) A which c is Fp4-element and tmp is a vector A Fp4
    Fp2_mul(&tmp.x0,&tmp.x0,&c);
    Fp2_mul(&tmp.x1,&tmp.x1,&c);
    
    Fp4_set(ANS,&tmp);
    
    Fp2_clear(&c);
    Fp2_clear(&a);
    Fp2_clear(&b);
    Fp4_clear(&tmp);
}
void Fp4_div(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B){
    struct Fp4 tmp,t_ans;
    Fp4_init(&tmp);
    Fp4_init(&t_ans);
    
    Fp4_invert(&tmp,B);
    Fp4_mul(&t_ans,A,&tmp);
    
    Fp4_set(ANS,&t_ans);
    
    Fp4_clear(&tmp);
    Fp4_clear(&t_ans);
}
void Fp4_pow(struct Fp4 *ANS,struct Fp4 *A,mpz_t B){
    int i,length;
    length= (int)mpz_sizeinbase(B,2);
    char B_binary[length];
    mpz_get_str(B_binary,2,B);
    struct Fp4 tmp;
    Fp4_init(&tmp);
    Fp4_set(&tmp,A);
    for(i=1;B_binary[i]!='\0';i++){
        Fp4_mul(&tmp,&tmp,&tmp);
        if(B_binary[i]=='1'){
            Fp4_mul(&tmp,&tmp,A);
        }
    }
    Fp4_set(ANS,&tmp);
    Fp4_clear(&tmp);
}
void Fp4_sqrt(struct Fp4 *ANS,struct Fp4 *A){
    struct Fp4 n,y,x,b,t,tmp_Fp4;
    Fp4_init(&n);
    Fp4_init(&y);
    Fp4_init(&x);
    Fp4_init(&b);
    Fp4_init(&t);
    Fp4_init(&tmp_Fp4);
    Fp4_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    while(Fp4_legendre(&n)!=-1){
        Fp4_random(&n);
    }
    mpz_pow_ui(q,PRIME_P,4);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    Fp4_pow(&y,&n,q);
    
    mpz_set(r,e);
    
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    
    Fp4_pow(&x,A,tmp_mpz);
    Fp4_pow(&tmp_Fp4,&x,set_2);
    Fp4_mul(&b,&tmp_Fp4,A);
    Fp4_mul(&x,&x,A);
    
    int m;
    
    while(Fp4_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp4_set(&tmp_Fp4,&b);
        while(Fp4_cmp_mpz(&tmp_Fp4,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp4_pow(&tmp_Fp4,&b,tmp_mpz);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,PRIME_P);
        // gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
        Fp4_pow(&t,&y,tmp_mpz);
        Fp4_pow(&y,&t,set_2);
        // gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
        mpz_set_ui(r,m);
        Fp4_mul(&x,&x,&t);
        Fp4_mul(&b,&b,&y);
    }
    
    Fp4_set(ANS,&x);
    
    Fp4_clear(&n);
    Fp4_clear(&y);
    Fp4_clear(&x);
    Fp4_clear(&b);
    Fp4_clear(&t);
    Fp4_clear(&tmp_Fp4);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp4_legendre(struct Fp4 *a){
    mpz_t i,cmp;
    struct Fp4 tmp;
    Fp4_init(&tmp);
    mpz_init(i);
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,PRIME_P,4);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp4_pow(&tmp,a,i);
    
    if((Fp4_cmp_mpz(&tmp,cmp))==0){
        Fp4_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp4_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
int Fp4_cmp(struct Fp4 *A,struct Fp4 *B){
    if(Fp2_cmp(&A->x0,&B->x0)==0 && Fp2_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp4_cmp_mpz(struct Fp4 *A,mpz_t B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    if(Fp2_cmp_mpz(&A->x0,B)==0 && Fp2_cmp(&A->x1,&tmp.x1)==0){
        Fp4_clear(&tmp);
        return 0;
    }
    Fp4_clear(&tmp);
    return 1;
}
void Fp4_neg(struct Fp4 *ans,struct Fp4 *a){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    Fp2_neg(&tmp.x0,&a->x0);
    Fp2_neg(&tmp.x1,&a->x1);
    Fp4_set(ans,&tmp);
    Fp4_clear(&tmp);
}

void Fp4_mul_betainv(struct Fp4 *ANS)
{
    struct Fp4 tmp;
    Fp4_init(&tmp);
    Fp4_set_ui(&tmp, 0);
    mpz_set(tmp.x1.x1.x0,a_x);
    
    //    mpz_t c_inv,c;
    //    mpz_init(c_inv);
    //    mpz_init(c);
    //    mpz_set_ui(c,c1);
    //    mpz_invert(c_inv,c,prime);
    mpz_mul(tmp.x1.x1.x0,tmp.x1.x1.x0,C1_INV.x0);
    
    Fp4_set(ANS, &tmp);
    
    //    mpz_clear(c);
    //    mpz_clear(c_inv);
    Fp4_clear(&tmp);
}

void Fp4_frobenius_map(struct Fp4 *ANS, struct Fp4 *A){
    struct Fp4 t_ans;
    Fp4_init(&t_ans);
    
    Fp2_frobenius_map(&t_ans.x0,&A->x0);
    Fp2_frobenius_map(&t_ans.x1,&A->x1);
    Fp2_mul_Fp(&t_ans.x1,&t_ans.x1,&p_C4);
    
    Fp4_set(ANS,&t_ans);
    
    Fp4_clear(&t_ans);
}

//-----------------------------------------------------------------------------------------
// #pragma mark Fp8 methods
void Fp8_init(struct Fp8 *A){
    Fp4_init(&A->x0);
    Fp4_init(&A->x1);
}
void Fp8_set(struct Fp8 *ANS,struct Fp8 *A){
    Fp4_set(&ANS->x0,&A->x0);
    Fp4_set(&ANS->x1,&A->x1);
}
void Fp8_set_ui(struct Fp8 *A,signed long int B){
    Fp4_set_ui(&A->x0,B);
    Fp4_set_ui(&A->x1,B);
}
void Fp8_random(struct Fp8 *A){
    Fp4_random(&A->x0);
    Fp4_random(&A->x1);
}
void Fp8_clear(struct Fp8 *A){
    Fp4_clear(&A->x0);
    Fp4_clear(&A->x1);
}
void Fp8_printf(struct Fp8 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd\n",A->x0.x0.x0.x0,A->x0.x0.x1.x0,A->x0.x1.x0.x0,A->x0.x1.x1.x0);
    gmp_printf("(%Zd,%Zd,%Zd,%Zd\n",A->x1.x0.x0.x0,A->x1.x0.x1.x0,A->x1.x1.x0.x0,A->x1.x1.x1.x0);
}
void Fp8_add(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_add(&tmp.x0,&A->x0,&B->x0);
    Fp4_add(&tmp.x1,&A->x1,&B->x1);
    
    Fp8_set(ANS,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_add_ui(struct Fp8 *ANS,struct Fp8 *A,unsigned long int B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_add_ui(&tmp.x0,&A->x0,B);
    Fp4_add_ui(&tmp.x1,&A->x1,B);
    
    Fp8_set(ANS,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_sub(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_sub(&tmp.x0,&A->x0,&B->x0);
    Fp4_sub(&tmp.x1,&A->x1,&B->x1);
    
    Fp8_set(ANS,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_sqr(struct Fp8 *ANS,struct Fp8 *A){
    
    struct Fp4 tmp1,tmp2,tmp3;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    
    struct Fp8 t_ans;
    Fp8_init(&t_ans);
    
    //Cost = 2M + 4Ap + 3m in Fp
    Fp4_mul(&tmp1,&A->x0,&A->x1);//t=a0a1
    Fp4_add(&tmp2,&A->x0,&A->x1);//(a0+a1)
    Fp4_mul_v(&tmp3,&A->x1);//a1*i
    Fp4_add(&tmp3,&tmp3,&A->x0);//(a0+a1*i)
    Fp4_mul(&t_ans.x0,&tmp2,&tmp3);// (a0+a1)(a0+a1*i)
    Fp4_sub(&t_ans.x0, &t_ans.x0,&tmp1);
    Fp4_mul_v(&tmp2,&tmp1);//
    Fp4_sub(&t_ans.x0,&t_ans.x0,&tmp2);//
    Fp4_mul_ui(&t_ans.x1,&tmp1,c1);//
    
    Fp8_set(ANS,&t_ans);
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp8_clear(&t_ans);
    
}

void Fp8_mul(struct Fp8 *ANS,struct Fp8 *A, struct Fp8 *B){
    if (Fp8_cmp(A, B) == 0) {
        Fp8_sqr(ANS, A);
    }
    else{
        //x^2-v=0
        struct Fp4 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
        Fp4_init(&tmp1);
        Fp4_init(&tmp2);
        Fp4_init(&tmp3);
        Fp4_init(&tmp4);
        Fp4_init(&tmp5);
        Fp4_init(&tmp6);
        
        struct Fp8 t_ans;
        Fp8_init(&t_ans);
        
        Fp4_mul(&tmp1,&A->x0,&B->x0);//a*c
        Fp4_mul(&tmp2,&A->x1,&B->x1);//b*d
        Fp4_mul_v(&tmp3,&tmp2);//b*d*v
        Fp4_add(&t_ans.x0,&tmp1,&tmp3);//a*c+b*d*v
        Fp4_add(&tmp4,&A->x0,&A->x1);//a+b
        Fp4_add(&tmp5,&B->x0,&B->x1);//c+d
        Fp4_mul(&tmp6,&tmp4,&tmp5);//(a+b)(c+d)
        Fp4_sub(&t_ans.x1,&tmp6,&tmp1);
        Fp4_sub(&t_ans.x1,&t_ans.x1,&tmp2);
        
        Fp8_set(ANS,&t_ans);
        
        Fp4_clear(&tmp1);
        Fp4_clear(&tmp2);
        Fp4_clear(&tmp3);
        Fp4_clear(&tmp4);
        Fp4_clear(&tmp5);
        Fp4_clear(&tmp6);
        Fp8_clear(&t_ans);
    }
}
void Fp8_mul_v(struct Fp8 *ANS,struct Fp8 *A){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_mul_v(&tmp.x0,&A->x1);
    Fp4_set(&tmp.x1,&A->x0);
    
    Fp8_set(ANS,&tmp);
    Fp8_clear(&tmp);
}
void Fp8_mul_ui(struct Fp8 *ANS,struct Fp8 *A,unsigned long int B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_mul_ui(&tmp.x0,&A->x0,B);
    Fp4_mul_ui(&tmp.x1,&A->x1,B);
    
    Fp8_set(ANS,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_mul_Fp(struct Fp8 *ANS,struct Fp8 *A,struct Fp *B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_mul_Fp(&tmp.x0,&A->x0,B);
    Fp4_mul_Fp(&tmp.x1,&A->x1,B);
    
    Fp8_set(ANS,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_mul_mpz(struct Fp8 *ANS,struct Fp8 *A,mpz_t B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_mul_mpz(&tmp.x0,&A->x0,B);
    Fp4_mul_mpz(&tmp.x1,&A->x1,B);
    
    Fp8_set(ANS,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_invert(struct Fp8 *ANS,struct Fp8 *A){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    // tmp=A^(q^4)=(x0,-x1)
    Fp4_set(&tmp.x0,&A->x0);
    Fp4_neg(&tmp.x1,&A->x1);
    
    struct Fp4 c,a,b;
    Fp4_init(&c);
    Fp4_init(&a);
    Fp4_init(&b);
    
    Fp4_mul(&a,&A->x0,&A->x0); // a=x0^2
    Fp4_mul(&b,&A->x1,&A->x1); // b=x1^2
    Fp4_mul_v(&b,&b); // b=x1^2*v
    Fp4_sub(&c,&a,&b); // c=x0^2-x1^2*v mod q
    
    Fp4_invert(&c,&c);
    
    // ANS=A^{-1}=(c)^{-1}*A^(p^6) A which c is Fp8-element and tmp is a vector A Fp8
    Fp4_mul(&tmp.x0,&tmp.x0,&c);
    Fp4_mul(&tmp.x1,&tmp.x1,&c);
    
    Fp8_set(ANS,&tmp);
    
    Fp4_clear(&c);
    Fp4_clear(&a);
    Fp4_clear(&b);
    Fp8_clear(&tmp);
}
void Fp8_div(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    struct Fp8 tmp,t_ans;
    Fp8_init(&tmp);
    Fp8_init(&t_ans);
    
    Fp8_invert(&tmp,B);
    Fp8_mul(&t_ans,A,&tmp);
    
    Fp8_set(ANS,&t_ans);
    
    Fp8_clear(&tmp);
    Fp8_clear(&t_ans);
}
void Fp8_pow(struct Fp8 *ANS,struct Fp8 *A,mpz_t B){
    int i,length;
    length= (int)mpz_sizeinbase(B,2);
    char B_binary[length];
    mpz_get_str(B_binary,2,B);
    struct Fp8 tmp;
    Fp8_init(&tmp);
    Fp8_set(&tmp,A);
    for(i=1;B_binary[i]!='\0';i++){
        Fp8_mul(&tmp,&tmp,&tmp);
        if(B_binary[i]=='1'){
            Fp8_mul(&tmp,&tmp,A);
        }
    }
    Fp8_set(ANS,&tmp);
    Fp8_clear(&tmp);
}
void Fp8_sqrt(struct Fp8 *ANS,struct Fp8 *A){
    struct Fp8 n,y,x,b,t,tmp_Fp4;
    Fp8_init(&n);
    Fp8_init(&y);
    Fp8_init(&x);
    Fp8_init(&b);
    Fp8_init(&t);
    Fp8_init(&tmp_Fp4);
    Fp8_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    while(Fp8_legendre(&n)!=-1){
        Fp8_random(&n);
    }
    mpz_pow_ui(q,PRIME_P,8);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    Fp8_pow(&y,&n,q);
    
    mpz_set(r,e);
    
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    
    Fp8_pow(&x,A,tmp_mpz);
    Fp8_pow(&tmp_Fp4,&x,set_2);
    Fp8_mul(&b,&tmp_Fp4,A);
    Fp8_mul(&x,&x,A);
    
    int m;
    
    while(Fp8_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp8_set(&tmp_Fp4,&b);
        while(Fp8_cmp_mpz(&tmp_Fp4,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp8_pow(&tmp_Fp4,&b,tmp_mpz);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,PRIME_P);
        // gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
        Fp8_pow(&t,&y,tmp_mpz);
        Fp8_pow(&y,&t,set_2);
        // gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
        mpz_set_ui(r,m);
        Fp8_mul(&x,&x,&t);
        Fp8_mul(&b,&b,&y);
    }
    
    Fp8_set(ANS,&x);
    
    Fp8_clear(&n);
    Fp8_clear(&y);
    Fp8_clear(&x);
    Fp8_clear(&b);
    Fp8_clear(&t);
    Fp8_clear(&tmp_Fp4);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp8_legendre(struct Fp8 *a){
    mpz_t i,cmp;
    struct Fp8 tmp;
    Fp8_init(&tmp);
    mpz_init(i);
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,PRIME_P,8);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp8_pow(&tmp,a,i);
    
    if((Fp8_cmp_mpz(&tmp,cmp))==0){
        Fp8_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp8_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
int Fp8_cmp(struct Fp8 *A,struct Fp8 *B){
    if(Fp4_cmp(&A->x0,&B->x0)==0 && Fp4_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp8_cmp_mpz(struct Fp8 *A,mpz_t B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    if(Fp4_cmp_mpz(&A->x0,B)==0 && Fp4_cmp(&A->x1,&tmp.x1)==0){
        Fp8_clear(&tmp);
        return 0;
    }
    Fp8_clear(&tmp);
    return 1;
}
void Fp8_neg(struct Fp8 *ans,struct Fp8 *a){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    Fp4_neg(&tmp.x0,&a->x0);
    Fp4_neg(&tmp.x1,&a->x1);
    Fp8_set(ans,&tmp);
    Fp8_clear(&tmp);
}
void Fp8_frobenius_map(struct Fp8 *ANS, struct Fp8 *A){
    struct Fp8 tmp_ans;
    struct Fp4 ans_tmp4;
    Fp4_init(&ans_tmp4);
    Fp8_init(&tmp_ans);
    
    Fp4_frobenius_map(&tmp_ans.x0,&A->x0);
    Fp4_frobenius_map(&tmp_ans.x1,&A->x1);
    Fp4_mul_Fp(&tmp_ans.x1,&tmp_ans.x1,&p_C8);
    
    Fp2_mul_basis(&tmp_ans.x1.x0, &tmp_ans.x1.x0);
    Fp2_mul_basis(&tmp_ans.x1.x1, &tmp_ans.x1.x1);
    
    Fp8_set(ANS,&tmp_ans);
    
    //    Fp_clear(&pm5d8);
    Fp4_clear(&ans_tmp4);
    //    Fp_clear(&set_c1);
    Fp8_clear(&tmp_ans);
}

//-----------------------------------------------------------------------------------------
// #pragma mark Fp16 methods
void Fp16_init(struct Fp16 *A){
    Fp8_init(&A->x0);
    Fp8_init(&A->x1);
}
void Fp16_set(struct Fp16 *ANS,struct Fp16 *A){
    Fp8_set(&ANS->x0,&A->x0);
    Fp8_set(&ANS->x1,&A->x1);
}
void Fp16_set_ui(struct Fp16 *A,signed long int B){
    Fp8_set_ui(&A->x0,B);
    Fp8_set_ui(&A->x1,B);
}
void Fp16_random(struct Fp16 *A){
    Fp8_random(&A->x0);
    Fp8_random(&A->x1);
}
void Fp16_clear(struct Fp16 *A){
    Fp8_clear(&A->x0);
    Fp8_clear(&A->x1);
}
void Fp16_printf(struct Fp16 *A){
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x0.x0.x0.x0.x0,A->x0.x0.x0.x1.x0,A->x0.x0.x1.x0.x0,A->x0.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x0.x1.x0.x0.x0,A->x0.x1.x0.x1.x0,A->x0.x1.x1.x0.x0,A->x0.x1.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x1.x0.x0.x0.x0,A->x1.x0.x0.x1.x0,A->x1.x0.x1.x0.x0,A->x1.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x1.x1.x0.x0.x0,A->x1.x1.x0.x1.x0,A->x1.x1.x1.x0.x0,A->x1.x1.x1.x1.x0);
}
void Fp16_add(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_add(&tmp.x0,&A->x0,&B->x0);
    Fp8_add(&tmp.x1,&A->x1,&B->x1);
    
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
}
void Fp16_add_ui(struct Fp16 *ANS,struct Fp16 *A,unsigned long int B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_add_ui(&tmp.x0,&A->x0,B);
    Fp8_add_ui(&tmp.x1,&A->x1,B);
    
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
}
void Fp16_sub(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_sub(&tmp.x0,&A->x0,&B->x0);
    Fp8_sub(&tmp.x1,&A->x1,&B->x1);
    
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
}
void Fp16_sqr(struct Fp16 *ANS,struct Fp16 *A){
    
    struct Fp8 tmp1,tmp2,tmp3;
    Fp8_init(&tmp1);
    Fp8_init(&tmp2);
    Fp8_init(&tmp3);
    
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    
    //Cost = 2M + 4Ap + 3m in Fp
    Fp8_mul(&tmp1,&A->x0,&A->x1);//t=a0a1
    Fp8_add(&tmp2,&A->x0,&A->x1);//(a0+a1)
    Fp8_mul_v(&tmp3,&A->x1);//a1*i
    Fp8_add(&tmp3,&tmp3,&A->x0);//(a0+a1*i)
    Fp8_mul(&t_ans.x0,&tmp2,&tmp3);// (a0+a1)(a0+a1*i)
    Fp8_sub(&t_ans.x0, &t_ans.x0,&tmp1);
    Fp8_mul_v(&tmp2,&tmp1);//
    Fp8_sub(&t_ans.x0,&t_ans.x0,&tmp2);//
    Fp8_mul_ui(&t_ans.x1,&tmp1,c1);//
    
    Fp16_set(ANS,&t_ans);
    
    Fp8_clear(&tmp1);
    Fp8_clear(&tmp2);
    Fp8_clear(&tmp3);
    Fp16_clear(&t_ans);
}

void Fp16_mul(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    
    if (Fp16_cmp(A, B)==0) {
        Fp16_sqr(ANS, A);
    }
    else{
        //x^2-v=0
        struct Fp8 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
        Fp8_init(&tmp1);
        Fp8_init(&tmp2);
        Fp8_init(&tmp3);
        Fp8_init(&tmp4);
        Fp8_init(&tmp5);
        Fp8_init(&tmp6);
        
        struct Fp16 t_ans;
        Fp16_init(&t_ans);
        
        //A=a0+a1 B=b0+b1
        //t1 = a0*b0
        //t2 = a1*b1
        //t3 = (t2)*gamma
        //ans.x0 = t1+t3
        //t4 = a0+a1
        //t5 = b0+b1
        //t6 = (t4*t5)
        //ans.x1 =t6-t1-t2
        
        Fp8_mul(&tmp1,&A->x0,&B->x0);//a*c
        Fp8_mul(&tmp2,&A->x1,&B->x1);//b*d
        Fp8_mul_v(&tmp3,&tmp2);//b*d*v
        Fp8_add(&t_ans.x0,&tmp1,&tmp3);//a*c+b*d*v
        Fp8_add(&tmp4,&A->x0,&A->x1);//a+b
        Fp8_add(&tmp5,&B->x0,&B->x1);//c+d
        Fp8_mul(&tmp6,&tmp4,&tmp5);//(a+b)(c+d)
        Fp8_sub(&t_ans.x1,&tmp6,&tmp1);
        Fp8_sub(&t_ans.x1,&t_ans.x1,&tmp2);
        
        Fp16_set(ANS,&t_ans);
        
        Fp8_clear(&tmp1);
        Fp8_clear(&tmp2);
        Fp8_clear(&tmp3);
        Fp8_clear(&tmp4);
        Fp8_clear(&tmp5);
        Fp8_clear(&tmp6);
        Fp16_clear(&t_ans);
    }
}
void Fp16_mul_v(struct Fp16 *ANS,struct Fp16 *A){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_mul_v(&tmp.x0,&A->x1);
    Fp8_set(&tmp.x1,&A->x0);
    
    Fp16_set(ANS,&tmp);
    Fp16_clear(&tmp);
}
void Fp16_mul_ui(struct Fp16 *ANS,struct Fp16 *A,unsigned long int B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_mul_ui(&tmp.x0,&A->x0,B);
    Fp8_mul_ui(&tmp.x1,&A->x1,B);
    
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
}
void Fp16_mul_Fp(struct Fp16 *ANS,struct Fp16 *A,struct Fp *B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_mul_Fp(&tmp.x0,&A->x0,B);
    Fp8_mul_Fp(&tmp.x1,&A->x1,B);
    
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
}
void Fp16_mul_mpz(struct Fp16 *ANS,struct Fp16 *A,mpz_t B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_mul_mpz(&tmp.x0,&A->x0,B);
    Fp8_mul_mpz(&tmp.x1,&A->x1,B);
    
    Fp16_set(ANS,&tmp);
    Fp16_clear(&tmp);
}
void Fp16_invert(struct Fp16 *ANS,struct Fp16 *A){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    // tmp=A^(p^8)=(x0,-x1)
    Fp8_set(&tmp.x0,&A->x0);
    Fp8_neg(&tmp.x1,&A->x1);
    
    struct Fp8 c,a,b;
    Fp8_init(&c);
    Fp8_init(&a);
    Fp8_init(&b);
    
    Fp8_mul(&a,&A->x0,&A->x0); // a=x0^2
    Fp8_mul(&b,&A->x1,&A->x1); // b=x1^2
    Fp8_mul_v(&b,&b); // b=x1^2*v
    Fp8_sub(&c,&a,&b); // c=x0^2-x1^2*v mod q
    
    Fp8_invert(&c,&c);
    
    // ANS=A^{-1}=(c)^{-1}*A^(p^6) A which c is Fp16-element and tmp is a vector A Fp16
    Fp8_mul(&tmp.x0,&tmp.x0,&c);
    Fp8_mul(&tmp.x1,&tmp.x1,&c);
    
    Fp16_set(ANS,&tmp);
    
    Fp8_clear(&c);
    Fp8_clear(&a);
    Fp8_clear(&b);
    Fp16_clear(&tmp);
}
void Fp16_div(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    struct Fp16 tmp,t_ans;
    Fp16_init(&tmp);
    Fp16_init(&t_ans);
    
    Fp16_invert(&tmp,B);
    Fp16_mul(&t_ans,A,&tmp);
    
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&tmp);
    Fp16_clear(&t_ans);
}
void Fp16_pow_mother_parameter_add_1(struct Fp16 *ANS,struct Fp16 *A){
    int i;
    struct Fp16 tmp,A_inv;
    Fp16_init(&tmp);
    Fp16_init(&A_inv);
//    Fp16_frobenius_map_p8(&A_inv,A);
    Fp8_set(&A_inv.x0, &A->x0);
    Fp8_neg(&A_inv.x1, &A->x1);
    
    Fp16_set(&tmp,&A_inv);
    for(i=x_bit-1; i>=0; i--){
        switch(X1_bit_binary[i]){
            case 0:
                Fp16_sqr(&tmp,&tmp);
                break;
            case 1:
                Fp16_sqr(&tmp,&tmp);
                Fp16_mul(&tmp,&tmp,A);
                break;
            case -1:
                Fp16_sqr(&tmp,&tmp);
                Fp16_mul(&tmp,&tmp,&A_inv);
                break;
            default:
                break;
        }
    }
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
    Fp16_clear(&A_inv);
}

void Fp16_pow_mother_parameter(struct Fp16 *ANS,struct Fp16 *A){
    int i;
    struct Fp16 tmp,A_inv;
    Fp16_init(&tmp);
    Fp16_init(&A_inv);
    Fp8_set(&A_inv.x0, &A->x0);
    Fp8_neg(&A_inv.x1, &A->x1);
//    Fp16_frobenius_map_p8(&A_inv,A);
    
    Fp16_set(&tmp,&A_inv);
    for(i=x_bit-1; i>=0; i--){
        switch(X_bit_binary[i]){
            case 0:
                Fp16_sqr(&tmp,&tmp);
                break;
            case 1:
                Fp16_sqr(&tmp,&tmp);
                Fp16_mul(&tmp,&tmp,A);
                break;
            case -1:
                Fp16_sqr(&tmp,&tmp);
                Fp16_mul(&tmp,&tmp,&A_inv);
                break;
            default:
                break;
        }
    }
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
    Fp16_clear(&A_inv);
}

void Fp16_pow(struct Fp16 *ANS,struct Fp16 *A,mpz_t B){
    //    printf("\n\n 1 Fp16_pow \n");
    //    clock_t start = clock();
    int i,length;
    length= (int)mpz_sizeinbase(B,2);
    char B_binary[length];
    mpz_get_str(B_binary,2,B);
    struct Fp16 tmp;
    Fp16_init(&tmp);
    Fp16_set(&tmp,A);
    //    for(i=1;B_binary[i]!='\0';i++){
    for(i=length-2; i>=0; i--){
        Fp16_mul(&tmp,&tmp,&tmp);
        if(mpz_tstbit(B,i)==1){
            //        if(B_binary[i]=='1'){
            Fp16_mul(&tmp,&tmp,A);
        }
    }
    
    fp16_pow++;
    //    clock_t stop = clock();
    //    double elapsed = (double)(stop - start) * 1000.0 / CLOCKS_PER_SEC;
    //    printf("%d -th Fp16 pow ms: %f [ms]\n",fp16_pow, elapsed);
    Fp16_set(ANS,&tmp);
    Fp16_clear(&tmp);
}
void Fp16_sqrt(struct Fp16 *ANS,struct Fp16 *A){
    struct Fp16 n,y,x,b,t,tmp_Fp4;
    Fp16_init(&n);
    Fp16_init(&y);
    Fp16_init(&x);
    Fp16_init(&b);
    Fp16_init(&t);
    Fp16_init(&tmp_Fp4);
    Fp16_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    while(Fp16_legendre(&n)!=-1){
        Fp16_random(&n);
    }
    mpz_pow_ui(q,PRIME_P,16);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    Fp16_pow(&y,&n,q);
    
    mpz_set(r,e);
    
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    
    Fp16_pow(&x,A,tmp_mpz);
    Fp16_pow(&tmp_Fp4,&x,set_2);
    Fp16_mul(&b,&tmp_Fp4,A);
    Fp16_mul(&x,&x,A);
    
    int m;
    
    while(Fp16_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp16_set(&tmp_Fp4,&b);
        while(Fp16_cmp_mpz(&tmp_Fp4,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp16_pow(&tmp_Fp4,&b,tmp_mpz);
            // Fp16_printf(&tmp_Fp4);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,PRIME_P);
        // gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
        Fp16_pow(&t,&y,tmp_mpz);
        Fp16_pow(&y,&t,set_2);
        // gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
        mpz_set_ui(r,m);
        Fp16_mul(&x,&x,&t);
        Fp16_mul(&b,&b,&y);
        
    }
    
    Fp16_set(ANS,&x);
    
    Fp16_clear(&n);
    Fp16_clear(&y);
    Fp16_clear(&x);
    Fp16_clear(&b);
    Fp16_clear(&t);
    Fp16_clear(&tmp_Fp4);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp16_legendre(struct Fp16 *a){
    mpz_t i,cmp;
    struct Fp16 tmp;
    Fp16_init(&tmp);
    mpz_init(i);
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,PRIME_P,16);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp16_pow(&tmp,a,i);
    
    if((Fp16_cmp_mpz(&tmp,cmp))==0){
        Fp16_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp16_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
int Fp16_cmp(struct Fp16 *A,struct Fp16 *B){
    if(Fp8_cmp(&A->x0,&B->x0)==0 && Fp8_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp16_cmp_mpz(struct Fp16 *A,mpz_t B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    if(Fp8_cmp_mpz(&A->x0,B)==0 && Fp8_cmp(&A->x1,&tmp.x1)==0){
        Fp16_clear(&tmp);
        return 0;
    }
    Fp16_clear(&tmp);
    return 1;
}
void Fp16_neg(struct Fp16 *ans,struct Fp16 *a){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    Fp8_neg(&tmp.x0,&a->x0);
    Fp8_neg(&tmp.x1,&a->x1);
    Fp16_set(ans,&tmp);
    Fp16_clear(&tmp);
}

void Fp16_frobenius_map(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 tmp_ans;
    Fp16_init(&tmp_ans);
    Fp16_set(&tmp_ans, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //a0-a3
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &tmp_ans.x0.x0.x0.x0);
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &tmp_ans.x0.x0.x0.x1);
    if (mpz_cmp_ui(tmp_ans.x0.x0.x0.x1.x0, 0) != 0) {
        Fp_neg(&ans_tmp8.x0.x0.x0.x1, &ans_tmp8.x0.x0.x0.x1);
    }
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &tmp_ans.x0.x0.x1.x0, &p_C4);
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &tmp_ans.x0.x0.x1.x1, &p_M_C4);
    
    //a4-a7
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &tmp_ans.x0.x1.x0.x1, &p_M_C_C8);
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &tmp_ans.x0.x1.x0.x0, &p_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &tmp_ans.x0.x1.x1.x1, &p_M_C_C4_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &tmp_ans.x0.x1.x1.x0, &p_C4_C16);
    
    //from a8-a11
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &tmp_ans.x1.x0.x1.x0, &p_C_C4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &tmp_ans.x1.x0.x1.x1, &p_M_C_C4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &tmp_ans.x1.x0.x0.x1, &p_M_C_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &tmp_ans.x1.x0.x0.x0, &p_C16);
    
    //a12-a15
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &tmp_ans.x1.x1.x1.x1, &p_M_C_C_C4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &tmp_ans.x1.x1.x1.x0, &p_C4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &tmp_ans.x1.x1.x0.x0, &p_C_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &tmp_ans.x1.x1.x0.x1, &p_M_C_C8_C16);
    
    Fp8_set(&tmp_ans.x0, &ans_tmp8.x0);
    Fp8_set(&tmp_ans.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&tmp_ans);
    
    Fp16_clear(&tmp_ans);
}

void Fp16_frobenius_map_p2(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 tmp_ans;
    Fp16_init(&tmp_ans);
    Fp16_set(&tmp_ans, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //a0-a3
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &tmp_ans.x0.x0.x0.x0);
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &tmp_ans.x0.x0.x0.x1);
    Fp_neg(&ans_tmp8.x0.x0.x1.x0, &tmp_ans.x0.x0.x1.x0);
    Fp_neg(&ans_tmp8.x0.x0.x1.x1, &tmp_ans.x0.x0.x1.x1);
    
    //a4-a7
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &tmp_ans.x0.x1.x0.x0, &p2_C8);
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &tmp_ans.x0.x1.x0.x1, &p2_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &tmp_ans.x0.x1.x1.x0, &p2_M_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &tmp_ans.x0.x1.x1.x1, &p2_M_C8);
    
    //from a8-a11
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &tmp_ans.x1.x0.x0.x1, &p2_C_C16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &tmp_ans.x1.x0.x0.x0, &p2_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &tmp_ans.x1.x0.x1.x1, &p2_M_C_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &tmp_ans.x1.x0.x1.x0, &p2_M_C16);
    
    //a12-a15
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &tmp_ans.x1.x1.x0.x1, &p2_C_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &tmp_ans.x1.x1.x0.x0, &p2_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &tmp_ans.x1.x1.x1.x1, &p2_M_C_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &tmp_ans.x1.x1.x1.x0, &p2_M_C8_C16);
    
    Fp8_set(&tmp_ans.x0, &ans_tmp8.x0);
    Fp8_set(&tmp_ans.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&tmp_ans);
    
    Fp16_clear(&tmp_ans);
}


void Fp16_frobenius_map_p3(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 TMP;
    Fp16_init(&TMP);
    Fp16_set(&TMP, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //1-4
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &TMP.x0.x0.x0.x0); //a0
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &TMP.x0.x0.x0.x1);//-a1
    if (mpz_cmp_ui(TMP.x0.x0.x0.x1.x0, 0) != 0) {
        Fp_neg(&ans_tmp8.x0.x0.x0.x1, &ans_tmp8.x0.x0.x0.x1);
    }
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &TMP.x0.x0.x1.x0, &p3_C4);//c4 a2
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &TMP.x0.x0.x1.x1, &p3_M_C4);//-c4 a3
    //5-8
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &TMP.x0.x1.x0.x1, &p3_M_C_C8); //-cc8 a5
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &TMP.x0.x1.x0.x0, &p3_C8);//c8 a4
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &TMP.x0.x1.x1.x1, &p3_M_C_C4_C8); // -c8c4c a7
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &TMP.x0.x1.x1.x0, &p3_C4_C8); //c8c4
    
    //from 9-16
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &TMP.x1.x0.x1.x1, &p3_M_C_C4_C16);//-cc4c16
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &TMP.x1.x0.x1.x0, &p3_C4_C16);//c4c16
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &TMP.x1.x0.x0.x0, &p3_C16);//c16
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &TMP.x1.x0.x0.x1, &p3_M_C16);//-c16
    
    
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &TMP.x1.x1.x1.x0, &p3_C_C4_C8_C16);//cc4c8c16
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &TMP.x1.x1.x1.x1, &p3_M_C_C4_C8_C16);//-cc4c8c16
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &TMP.x1.x1.x0.x1, &p3_M_C_C8_C16);//-cc8c16
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &TMP.x1.x1.x0.x0, &p3_C8_C16);//c8c16
    
    
    Fp8_set(&TMP.x0, &ans_tmp8.x0);
    Fp8_set(&TMP.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&TMP);
    
    Fp16_clear(&TMP);
    Fp16_clear(&ans_tmp8);
}

void Fp16_frobenius_map_p4(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 TMP;
    Fp16_init(&TMP);
    Fp16_set(&TMP, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //1-4
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &TMP.x0.x0.x0.x0);
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &TMP.x0.x0.x0.x1);
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &TMP.x0.x0.x1.x0, &p4_C4);
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &TMP.x0.x0.x1.x1, &p4_C4);
    
    //5-8
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &TMP.x0.x1.x0.x0, &p4_C8);
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &TMP.x0.x1.x0.x1, &p4_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &TMP.x0.x1.x1.x0, &p4_C4_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &TMP.x0.x1.x1.x1, &p4_C4_C8);
    
    //from 9-11
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &TMP.x1.x0.x0.x0, &p4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &TMP.x1.x0.x0.x1, &p4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &TMP.x1.x0.x1.x0, &p4_C4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &TMP.x1.x0.x1.x1, &p4_C4_C16);
    //12-16
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &TMP.x1.x1.x0.x0, &p4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &TMP.x1.x1.x0.x1, &p4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &TMP.x1.x1.x1.x0, &p4_C4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &TMP.x1.x1.x1.x1, &p4_C4_C8_C16);
    
    Fp8_set(&TMP.x0, &ans_tmp8.x0);
    Fp8_set(&TMP.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&TMP);
    
    Fp16_clear(&TMP);
    Fp16_clear(&ans_tmp8);
}

void Fp16_frobenius_map_p8(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 TMP;
    Fp16_init(&TMP);
    Fp16_set(&TMP, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //1-4
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &TMP.x0.x0.x0.x0);
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &TMP.x0.x0.x0.x1);
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &TMP.x0.x0.x1.x0, &p8_C4);
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &TMP.x0.x0.x1.x1, &p8_C4);
    
    //5-8
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &TMP.x0.x1.x0.x0, &p8_C8);
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &TMP.x0.x1.x0.x1, &p8_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &TMP.x0.x1.x1.x0, &p8_C4_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &TMP.x0.x1.x1.x1, &p8_C4_C8);
    
    //from 9-11
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &TMP.x1.x0.x0.x0, &p8_C16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &TMP.x1.x0.x0.x1, &p8_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &TMP.x1.x0.x1.x0, &p8_C4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &TMP.x1.x0.x1.x1, &p8_C4_C16);
    //12-16
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &TMP.x1.x1.x0.x0, &p8_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &TMP.x1.x1.x0.x1, &p8_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &TMP.x1.x1.x1.x0, &p8_C4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &TMP.x1.x1.x1.x1, &p8_C4_C8_C16);
    
    Fp8_set(&TMP.x0, &ans_tmp8.x0);
    Fp8_set(&TMP.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&TMP);
    
    Fp16_clear(&TMP);
    Fp16_clear(&ans_tmp8);
}

void Fp16_frobenius_map_p5(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 tmp_ans;
    Fp16_init(&tmp_ans);
    Fp16_set(&tmp_ans, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //a0-a3
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &tmp_ans.x0.x0.x0.x0);
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &tmp_ans.x0.x0.x0.x1);
    if (mpz_cmp_ui(tmp_ans.x0.x0.x0.x1.x0, 0) != 0) {
        Fp_neg(&ans_tmp8.x0.x0.x0.x1, &ans_tmp8.x0.x0.x0.x1);
    }
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &tmp_ans.x0.x0.x1.x0, &p5_C4);
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &tmp_ans.x0.x0.x1.x1, &p5_M_C4);
    
    //a4-a7
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &tmp_ans.x0.x1.x0.x1, &p5_M_C8);
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &tmp_ans.x0.x1.x0.x0, &p5_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &tmp_ans.x0.x1.x1.x1, &p5_M_C_C4_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &tmp_ans.x0.x1.x1.x0, &p5_C4_C16);
    
    //from a8-a11
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &tmp_ans.x1.x0.x1.x0, &p5_C_C4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &tmp_ans.x1.x0.x1.x1, &p5_M_C_C4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &tmp_ans.x1.x0.x0.x1, &p5_M_C_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &tmp_ans.x1.x0.x0.x0, &p5_C16);
    
    //a12-a15
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &tmp_ans.x1.x1.x1.x1, &p5_M_C_C_C4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &tmp_ans.x1.x1.x1.x0, &p5_C4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &tmp_ans.x1.x1.x0.x0, &p5_C_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &tmp_ans.x1.x1.x0.x1, &p5_M_C_C8_C16);
    
    Fp8_set(&tmp_ans.x0, &ans_tmp8.x0);
    Fp8_set(&tmp_ans.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&tmp_ans);
    
    Fp16_clear(&tmp_ans);
}


void Fp16_frobenius_map_p6(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 tmp_ans;
    Fp16_init(&tmp_ans);
    Fp16_set(&tmp_ans, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //a0-a3
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &tmp_ans.x0.x0.x0.x0);
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &tmp_ans.x0.x0.x0.x1);
    Fp_neg(&ans_tmp8.x0.x0.x1.x0, &tmp_ans.x0.x0.x1.x0);
    Fp_neg(&ans_tmp8.x0.x0.x1.x1, &tmp_ans.x0.x0.x1.x1);
    
    //a4-a7
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &tmp_ans.x0.x1.x0.x0, &p6_C8);
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &tmp_ans.x0.x1.x0.x1, &p6_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &tmp_ans.x0.x1.x1.x0, &p6_M_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &tmp_ans.x0.x1.x1.x1, &p6_M_C8);
    
    //from a8-a11
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &tmp_ans.x1.x0.x0.x1, &p6_C_C16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &tmp_ans.x1.x0.x0.x0, &p6_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &tmp_ans.x1.x0.x1.x1, &p6_M_C_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &tmp_ans.x1.x0.x1.x0, &p6_M_C16);
    
    //a12-a15
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &tmp_ans.x1.x1.x0.x1, &p6_C_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &tmp_ans.x1.x1.x0.x0, &p6_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &tmp_ans.x1.x1.x1.x1, &p6_M_C_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &tmp_ans.x1.x1.x1.x0, &p6_M_C8_C16);
    
    Fp8_set(&tmp_ans.x0, &ans_tmp8.x0);
    Fp8_set(&tmp_ans.x1, &ans_tmp8.x1);
    Fp16_set(ANS,&tmp_ans);
    
    Fp16_clear(&tmp_ans);
}

void Fp16_frobenius_map_p7(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 TMP;
    Fp16_init(&TMP);
    Fp16_set(&TMP, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //1-4
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &TMP.x0.x0.x0.x0); //a0
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &TMP.x0.x0.x0.x1);//-a1
    if (mpz_cmp_ui(TMP.x0.x0.x0.x1.x0, 0) != 0) {
        Fp_neg(&ans_tmp8.x0.x0.x0.x1, &ans_tmp8.x0.x0.x0.x1);
    }
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &TMP.x0.x0.x1.x0, &p7_C4);//c4 a2
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &TMP.x0.x0.x1.x1, &p7_M_C4);//-c4 a3
    //5-8
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &TMP.x0.x1.x0.x1, &p7_M_C_C8); //-cc8 a5
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &TMP.x0.x1.x0.x0, &p7_C8);//c8 a4
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &TMP.x0.x1.x1.x1, &p7_M_C_C4_C8); // -c8c4c a7
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &TMP.x0.x1.x1.x0, &p7_C8_C4); //c8c4
    //from 9-16
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &TMP.x1.x0.x1.x1, &p7_M_C_C4_C16);//-cc4c16
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &TMP.x1.x0.x1.x0, &p7_C4_C16);//c4c16
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &TMP.x1.x0.x0.x0, &p7_C16);//c16
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &TMP.x1.x0.x0.x1, &p7_M_C16);//-c16
    
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &TMP.x1.x1.x1.x0, &p7_C_C4_C8_C16);//cc4c8c16
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &TMP.x1.x1.x1.x1, &p7_M_C_C4_C8_C16);//-cc4c8c16
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &TMP.x1.x1.x0.x1, &p7_M_C_C8_C16);//-cc8c16
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &TMP.x1.x1.x0.x0, &p7_C8_C16);//c8c16
    
    
    Fp8_set(&TMP.x0, &ans_tmp8.x0);
    Fp8_set(&TMP.x1, &ans_tmp8.x1);
    Fp16_set(ANS,&TMP);
    
    Fp16_clear(&TMP);
    Fp16_clear(&ans_tmp8);
}

//-----------------------------------------------------------------------------------------
// #pragma mark EFp methods
void EFp_init(struct EFp *A){
    Fp_init(&A->x);
    Fp_init(&A->y);
    A->infity=FALSE;
}
void EFp_set(struct EFp *A,struct EFp *B){
    Fp_set(&A->x,&B->x);
    Fp_set(&A->y,&B->y);
    A->infity=B->infity;
}
void EFp_set_infity(struct EFp *A){
    Fp_set_ui(&A->x,0);
    Fp_set_ui(&A->y,0);
    A->infity=TRUE;
}
void EFp_set_EFp16(struct EFp *ANS,struct EFp16 *P){
    Fp_set(&ANS->x,&P->x.x0.x0.x0.x0);
    Fp_set(&ANS->y,&P->y.x0.x0.x0.x0);
    ANS->infity=FALSE;
}
void EFp_clear(struct EFp *A){
    Fp_clear(&A->x);
    Fp_clear(&A->y);
}
void EFp_printf(struct EFp *A){
    gmp_printf("(%Zd,%Zd)\n",A->x.x0,A->y.x0);
}
void EFp_SCM_BIN(struct EFp *ANS, struct EFp *P,mpz_t j){
    int i;
    int r;//bit数
    r= (int)mpz_sizeinbase(j,2);
    
    struct EFp Q;
    EFp_init(&Q);
    EFp_set(&Q,P);
    
    for(i=r-2;i>=0;i--){
        if(mpz_tstbit(j,i)==1){
            EFp_ECD(&Q,&Q);
            EFp_ECA(&Q,&Q,P);
        }else{
            EFp_ECD(&Q,&Q);
        }
    }
    
    EFp_set(ANS,&Q);
    EFp_clear(&Q);
    return;
}
//TODO: NAF and SLIDING WINDOW
void TO_NAF(mpz_t scalar){
    int size = (int)mpz_sizeinbase(scalar,2);
    myNAF[size];
    naf_index = 0;
    mpz_t zi,tmp1,tmp2,scalar_tmp;
    mpz_init(zi);
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(scalar_tmp);
    mpz_set(scalar_tmp,scalar);
    //    char x1_binary[(int)mpz_sizeinbase(scalar_tmp,2)];
    //    mpz_get_str(x1_binary,2,scalar_tmp);
    //    printf("scalar binary %c\n",x1_binary[5]);
    
    while (mpz_cmp_ui(scalar_tmp, 0) > 0) {
        
        if (mpz_odd_p(scalar_tmp)){
            mpz_mod_ui(tmp1,scalar_tmp,4);
            mpz_set_ui(tmp2,2);
            mpz_sub(zi,tmp2,tmp1);
            mpz_sub(scalar_tmp,scalar_tmp,zi);
        }
        else{
            mpz_set_ui(zi,0);
        }
        mpz_tdiv_q_2exp(scalar_tmp,scalar_tmp,1);
        myNAF[naf_index] = mpz_get_si(zi);
        naf_index++;
    }
    //    printf("To NAF index %ld\n",naf_index);
    signed long temp, end;
    end = naf_index - 1;
    for (int c = 0; c < naf_index/2; c++)
    {
        temp       = myNAF[c];
        myNAF[c]   = myNAF[end];
        myNAF[end] = temp;
        end--;
    }
}

void EFp_ECD(struct EFp *ANS, struct EFp *P){
    if(P->infity==TRUE){
        EFp_set(ANS,P);
        return;
    }
    if(mpz_sgn(P->y.x0)==0){//P.y==0
        EFp_set_infity(ANS);
        return;
    }
    
    struct Fp x,y,lambda,tmp;
    struct EFp t_ans;
    Fp_init(&x);
    Fp_init(&lambda);
    Fp_init(&tmp);
    Fp_init(&y);
    EFp_init(&t_ans);
    
    Fp_mul(&x,&P->x,&P->x);
    Fp_add(&tmp,&x,&x);
    Fp_add(&x,&tmp,&x);//3x^2+a
    //    gmp_printf("tmem A = %Zd\n",tmp_a);
    Fp_add_mpz(&x,&x,tmp_a);
    
    Fp_add(&y,&P->y,&P->y);//2y
    
    Fp_div(&lambda,&x,&y);
    Fp_mul(&tmp,&lambda,&lambda);
    Fp_add(&x,&P->x,&P->x);
    Fp_sub(&x,&tmp,&x);
    Fp_sub(&tmp,&P->x,&x);
    
    
    Fp_set(&t_ans.x,&x);
    
    Fp_mul(&tmp,&tmp,&lambda);
    Fp_sub(&t_ans.y,&tmp,&P->y);
    
    EFp_set(ANS,&t_ans);
    
    Fp_clear(&x);
    Fp_clear(&lambda);
    Fp_clear(&y);
    Fp_clear(&tmp);
    EFp_clear(&t_ans);
}
void EFp_ECA(struct EFp *ANS, struct EFp *P1, struct EFp *P2){
    if(P2->infity==TRUE){//if P2==inf
        EFp_set(ANS,P1);
        return;
    }
    else if(P1->infity==TRUE){//if P1==inf
        EFp_set(ANS,P2);
        return;
    }
    else if(Fp_cmp(&P1->x,&P2->x)==0 && Fp_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp_set_infity(ANS);
        return;
    }
    else if(EFp_cmp(P1,P2)==0){ // P=Q
        EFp_ECD(ANS,P1);
        return;
    }
    
    struct Fp x,y,lambda,tmp;
    struct EFp t_ans;
    
    Fp_init(&x);
    Fp_init(&y);
    Fp_init(&lambda);
    Fp_init(&tmp);
    EFp_init(&t_ans);
    
    Fp_sub(&x,&P2->x,&P1->x);
    Fp_sub(&y,&P2->y,&P1->y);
    Fp_div(&lambda,&y,&x);
    Fp_mul(&tmp,&lambda,&lambda);
    Fp_add(&x,&P1->x,&P2->x);
    Fp_sub(&x,&tmp,&x);
    Fp_sub(&tmp,&P1->x,&x);
    Fp_set(&t_ans.x,&x);
    Fp_mul(&tmp,&tmp,&lambda);
    Fp_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp_set(ANS,&t_ans);
    
    Fp_clear(&x);
    Fp_clear(&y);
    Fp_clear(&lambda);
    Fp_clear(&tmp);
    EFp_clear(&t_ans);
}
int EFp_cmp(struct EFp *A,struct EFp *B){
    if(Fp_cmp(&A->x,&B->x)==0 && Fp_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}
void EFp_neg(struct EFp *ANS, struct EFp *A){
    struct EFp tmp;
    EFp_init(&tmp);
    Fp_neg(&tmp.y,&A->y);
    Fp_set(&tmp.x,&A->x);
    EFp_set(ANS,&tmp);
    EFp_clear(&tmp);
}
void EFp_random_set(struct EFp *ANS){
    struct Fp a,x,tmp;
    Fp_init(&a);
    Fp_init(&x);
    Fp_init(&tmp);
    
    struct EFp P,Q;
    EFp_init(&P);
    EFp_init(&Q);
    
    
    do{
        Fp_random(&x);
        Fp_mul(&a,&x,&x);
        Fp_mul(&a,&a,&x);
        Fp_mul_mpz(&tmp, &x, a_x);
        Fp_add(&a, &a, &tmp);
        //    mpz_add(a.x0,a.x0,b);
    }while(mpz_legendre(a.x0,PRIME_P)!=1);
    Q.infity=0;
    Fp_sqrt(&P.y,&a);
    Fp_set(&P.x,&x);
    EFp_set(ANS,&P);
    
    
    Fp_clear(&a);
    Fp_clear(&x);
    EFp_clear(&P);
}
void EFp_Skew_Frobenius_p4(struct EFp *ANS, struct EFp *P){
    Fp_neg(&ANS->x,&P->x);
    Fp_mul_mpz(&ANS->y,&P->y,minus1_root2);
    ANS->infity=P->infity;
}
//-----------------------------------------------------------------------------------------
// #pragma mark EFp2 methods
void EFp2_init(struct EFp2 *A){
    Fp2_init(&A->x);
    Fp2_init(&A->y);
    A->infity=FALSE;
}
void EFp2_set(struct EFp2 *A,struct EFp2 *B){
    Fp2_set(&A->x,&B->x);
    Fp2_set(&A->y,&B->y);
    A->infity=B->infity;
}
void EFp2_set_infity(struct EFp2 *A){
    Fp2_set_ui(&A->x,0);
    Fp2_set_ui(&A->y,0);
    A->infity=TRUE;
}
void EFp2_clear(struct EFp2 *A){
    Fp2_clear(&A->x);
    Fp2_clear(&A->y);
}
void EFp2_printf(struct EFp2 *A){
    gmp_printf("(%Zd,%Zd)(%Zd,%Zd)\n",A->x.x0.x0,A->x.x1.x0,A->y.x0.x0,A->y.x1.x0);
}
void EFp2_SCM_BIN(struct EFp2 *ANS,struct EFp2 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp2 Q,R;
    EFp2_init(&Q);
    EFp2_set(&Q,P);
    EFp2_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp2_ECD(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp2_ECA(&Q,&Q,P);
        }
    }
    EFp2_set(ANS,&Q);
    
    EFp2_clear(&Q);
    EFp2_clear(&R);
    return;
}
void EFp2_ECD(struct EFp2 *ANS, struct EFp2 *P){
    if(P->infity==TRUE){
        EFp2_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp2_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp2_set_infity(ANS);
        return;
    }
    
    struct Fp2 x,y,lambda,tmp;
    struct EFp2 t_ans;
    Fp2_init(&x);
    Fp2_init(&lambda);
    Fp2_init(&tmp);
    Fp2_init(&y);
    EFp2_init(&t_ans);
    
    Fp2_mul(&x,&P->x,&P->x);
    Fp2_add(&tmp,&x,&x);
    Fp2_add(&x,&tmp,&x);
    Fp_add_mpz(&x.x0,&x.x0,a_x);
    Fp2_add(&y,&P->y,&P->y);
    Fp2_div(&lambda,&x,&y);
    Fp2_mul(&tmp,&lambda,&lambda);
    Fp2_add(&x,&P->x,&P->x);
    Fp2_sub(&x,&tmp,&x);
    Fp2_sub(&tmp,&P->x,&x);
    Fp2_set(&t_ans.x,&x);
    Fp2_mul(&tmp,&tmp,&lambda);
    Fp2_sub(&t_ans.y,&tmp,&P->y);
    
    EFp2_set(ANS,&t_ans);
    
    Fp2_clear(&x);
    Fp2_clear(&lambda);
    Fp2_clear(&y);
    Fp2_clear(&tmp);
    EFp2_clear(&t_ans);
}
void EFp2_ECA(struct EFp2 *ANS, struct EFp2 *P1, struct EFp2 *P2){
    if(P2->infity==TRUE){//if P2==inf
        EFp2_set(ANS,P1);
        return;
    }
    else if(P1->infity==TRUE){//if P1==inf
        EFp2_set(ANS,P2);
        return;
    }
    else if(Fp2_cmp(&P1->x,&P2->x)==0&&Fp2_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp2_set_infity(ANS);
        return;
    }
    else if(EFp2_cmp(P1,P2)==0){ // P=Q
        EFp2_ECD(ANS,P1);
        return;
    }
    
    struct Fp2 x,y,lambda,tmp;
    struct EFp2 t_ans;
    
    Fp2_init(&x);
    Fp2_init(&y);
    Fp2_init(&lambda);
    Fp2_init(&tmp);
    EFp2_init(&t_ans);
    
    Fp2_sub(&x,&P2->x,&P1->x);
    Fp2_sub(&y,&P2->y,&P1->y);
    Fp2_div(&lambda,&y,&x);
    Fp2_mul(&tmp,&lambda,&lambda);
    Fp2_add(&x,&P1->x,&P2->x);
    Fp2_sub(&x,&tmp,&x);
    Fp2_sub(&tmp,&P1->x,&x);
    Fp2_set(&t_ans.x,&x);
    Fp2_mul(&tmp,&tmp,&lambda);
    Fp2_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp2_set(ANS,&t_ans);
    
    Fp2_clear(&x);
    Fp2_clear(&y);
    Fp2_clear(&lambda);
    Fp2_clear(&tmp);
    EFp2_clear(&t_ans);
}
int EFp2_cmp(struct EFp2 *A,struct EFp2 *B){
    if(Fp2_cmp(&A->x,&B->x)==0 && Fp2_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}
void EFp2_random_set(struct EFp2 *ANS){
    struct EFp2 P;
    EFp2_init(&P);
    
    struct Fp2 x,a,tmp_fp;
    Fp2_init(&a);
    Fp2_init(&x);
    Fp2_init(&tmp_fp);
    
    mpz_t t2,p2,p22,tmp,r2;
    mpz_t set_3;
    mpz_init(set_3);
    mpz_set_ui(set_3,3);
    
    mpz_init(t2);
    mpz_init(p2);
    mpz_init(p22);
    mpz_init(tmp);
    mpz_init(r2);
    
    mpz_pow_ui(t2,trace_t,2);
    mpz_mul_ui(p2,PRIME_P,2);
    mpz_sub(tmp,t2,p2);
    
    mpz_pow_ui(p22,PRIME_P,2);
    mpz_add_ui(p22,p22,1);
    mpz_sub(p22,p22,tmp);
    
    do{
        Fp2_random(&x);
        //    Fp2_printf(&x);
        Fp2_pow(&a,&x,set_3);
        //    Fp2_printf(&a);
        Fp2_mul_mpz(&tmp_fp, &x, a_x);
        Fp2_add(&a, &a, &tmp_fp);
        //    mpz_add(a.x0.x0,a.x0.x0,tmp_fp.x0.x0);
        // Fp2_printf(&a);
        // printf("before crash %d\n",Fp2_legendre(&a));
        // printf("%d\n",Fp2_legendre(&a));
    }while(Fp2_legendre(&a)!=1);
    // Fp2_printf(&a);
    Fp2_sqrt(&P.y,&a);
    Fp2_set(&P.x,&x);
    
    mpz_t r12_div_r2;
    mpz_init(r12_div_r2);
    mpz_div(r12_div_r2,p22,order_r);
    mpz_div(r12_div_r2,r12_div_r2,order_r);
    
    printf("before scm bit leg====================\n");
    EFp2_SCM_BIN(ANS,&P,r12_div_r2);
    // EFp2_SCM_BIN(ANS,&P,p22);
    printf("after scm bit leg====================\n");
    EFp2_clear(&P);
    Fp2_clear(&a);
    Fp2_clear(&x);
    Fp2_clear(&tmp_fp);
    mpz_clear(tmp);
    mpz_clear(t2);
    mpz_clear(p22);
    mpz_clear(p2);
}

//-----------------------------------------------------------------------------------------
// #pragma mark EFp4 methods
void EFp4_init(struct EFp4 *A){
    Fp4_init(&A->x);
    Fp4_init(&A->y);
    A->infity=FALSE;
}
void EFp4_set(struct EFp4 *A,struct EFp4 *B){
    Fp4_set(&A->x,&B->x);
    Fp4_set(&A->y,&B->y);
    A->infity=B->infity;
}
void EFp4_set_infity(struct EFp4 *A){
    Fp4_set_ui(&A->x,0);
    Fp4_set_ui(&A->y,0);
    A->infity=TRUE;
}
void EFp4_set_EFp(struct EFp4 *ANS,struct EFp *A){
    Fp4_set_ui(&ANS->x,0);
    Fp4_set_ui(&ANS->y,0);
    
    Fp_set(&ANS->x.x0.x0,&A->x);
    Fp_set(&ANS->y.x0.x0,&A->y);
    ANS->infity=A->infity;
}
void EFp4_neg(struct EFp4 *ANS,struct EFp4 *A){
    Fp4_set(&ANS->x,&A->x);
    Fp4_neg(&ANS->y,&A->y);
    ANS->infity=A->infity;
}
void EFp4_clear(struct EFp4 *A){
    Fp4_clear(&A->x);
    Fp4_clear(&A->y);
}
void EFp4_printf(struct EFp4 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd)\n",A->x.x0.x0.x0,A->x.x0.x1.x0,A->x.x1.x0.x0, A->x.x1.x1.x0);
    gmp_printf("(%Zd,%Zd,%Zd,%Zd)\n",A->y.x0.x0.x0,A->y.x0.x1.x0,A->y.x1.x0.x0, A->y.x1.x1.x0);
}
void EFp4_SCM_BIN(struct EFp4 *ANS,struct EFp4 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp4 Q,R;
    EFp4_init(&Q);
    EFp4_set(&Q,P);
    EFp4_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp4_ECD(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp4_ECA(&Q,&Q,P);
        }
    }
    EFp4_set(ANS,&Q);
    
    EFp4_clear(&Q);
    EFp4_clear(&R);
    return;
}
void EFp4_ECD(struct EFp4 *ANS, struct EFp4 *P){
    if(P->infity==TRUE){
        EFp4_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp4_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp4_set_infity(ANS);
        return;
    }
    
    struct Fp4 x,y,lambda,tmp,a_betainv;
    struct EFp4 t_ans;
    Fp4_init(&x);
    Fp4_init(&lambda);
    Fp4_init(&tmp);
    Fp4_init(&y);
    Fp4_init(&a_betainv);
    
    EFp4_init(&t_ans);
    
    
    Fp4_mul(&x,&P->x,&P->x);
    Fp4_add(&tmp,&x,&x);
    Fp4_add(&x,&tmp,&x);
    
    //    Fp4_mul_betainv(&a_betainv);
    //    Fp4_add(&x, &x, &a_betainv);
    
    Fp_add_mpz(&x.x0.x0,&x.x0.x0,a_x);
    Fp4_add(&y,&P->y,&P->y);
    Fp4_div(&lambda,&x,&y);
    Fp4_mul(&tmp,&lambda,&lambda);
    Fp4_add(&x,&P->x,&P->x);
    Fp4_sub(&x,&tmp,&x);
    Fp4_sub(&tmp,&P->x,&x);
    Fp4_set(&t_ans.x,&x);
    Fp4_mul(&tmp,&tmp,&lambda);
    Fp4_sub(&t_ans.y,&tmp,&P->y);
    
    EFp4_set(ANS,&t_ans);
    
    Fp4_clear(&x);
    Fp4_clear(&lambda);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&t_ans);
}
void EFp4_ECD_Sparse(struct EFp4 *ANS, struct EFp4 *P){
    if(P->infity==TRUE){
        EFp4_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp4_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp4_set_infity(ANS);
        return;
    }
    
    struct Fp4 x,y,lambda,tmp,a_betainv;
    struct EFp4 t_ans;
    Fp4_init(&x);
    Fp4_init(&lambda);
    Fp4_init(&tmp);
    Fp4_init(&y);
    Fp4_init(&a_betainv);
    
    EFp4_init(&t_ans);
    
    
    Fp4_mul(&x,&P->x,&P->x);
    Fp4_add(&tmp,&x,&x);
    Fp4_add(&x,&tmp,&x);
    
    Fp4_mul_betainv(&a_betainv);
    Fp4_add(&x, &x, &a_betainv);
    
    // Fp_add_mpz(&x.x0.x0,&x.x0.x0,a_x);
    Fp4_add(&y,&P->y,&P->y);
    Fp4_div(&lambda,&x,&y);
    Fp4_mul(&tmp,&lambda,&lambda);
    Fp4_add(&x,&P->x,&P->x);
    Fp4_sub(&x,&tmp,&x);
    Fp4_sub(&tmp,&P->x,&x);
    Fp4_set(&t_ans.x,&x);
    Fp4_mul(&tmp,&tmp,&lambda);
    Fp4_sub(&t_ans.y,&tmp,&P->y);
    
    EFp4_set(ANS,&t_ans);
    
    Fp4_clear(&x);
    Fp4_clear(&lambda);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&t_ans);
}
void EFp4_ECD_Pseudo_Sparse(struct EFp4 *ANS, struct EFp4 *P){
    if(P->infity==TRUE){
        EFp4_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp4_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp4_set_infity(ANS);
        return;
    }
    
    struct Fp4 x,y,lambda,tmp,a_betainv;
    struct EFp4 t_ans;
    Fp4_init(&x);
    Fp4_init(&lambda);
    Fp4_init(&tmp);
    Fp4_init(&y);
    Fp4_init(&a_betainv);
    
    EFp4_init(&t_ans);
    
    
    Fp4_mul(&x,&P->x,&P->x);
    Fp4_add(&tmp,&x,&x);
    Fp4_add(&x,&tmp,&x);
    
    Fp4_mul_betainv(&a_betainv);
    Fp4_mul(&a_betainv, &a_betainv, &z_inv2);
    Fp4_add(&x, &x, &a_betainv);
    
    Fp4_add(&y,&P->y,&P->y);
    Fp4_div(&lambda,&x,&y);
    Fp4_mul(&tmp,&lambda,&lambda);
    Fp4_add(&x,&P->x,&P->x);
    Fp4_sub(&x,&tmp,&x);
    Fp4_sub(&tmp,&P->x,&x);
    Fp4_set(&t_ans.x,&x);
    Fp4_mul(&tmp,&tmp,&lambda);
    Fp4_sub(&t_ans.y,&tmp,&P->y);
    
    EFp4_set(ANS,&t_ans);
    
    Fp4_clear(&x);
    Fp4_clear(&lambda);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&t_ans);
}
void EFp4_ECA(struct EFp4 *ANS, struct EFp4 *P1, struct EFp4 *P2){
    if(P2->infity==TRUE){//if P2==inf
        EFp4_set(ANS,P1);
        return;
    }
    else if(P1->infity==TRUE){//if P1==inf
        EFp4_set(ANS,P2);
        return;
    }
    else if(Fp4_cmp(&P1->x,&P2->x)==0&&Fp4_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp4_set_infity(ANS);
        return;
    }
    else if(EFp4_cmp(P1,P2)==0){ // P=Q
        EFp4_ECD(ANS,P1);
        return;
    }
    
    struct Fp4 x,y,lambda,tmp;
    struct EFp4 t_ans;
    
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&lambda);
    Fp4_init(&tmp);
    EFp4_init(&t_ans);
    
    Fp4_sub(&x,&P2->x,&P1->x);
    Fp4_sub(&y,&P2->y,&P1->y);
    Fp4_div(&lambda,&y,&x);
    Fp4_mul(&tmp,&lambda,&lambda);
    Fp4_add(&x,&P1->x,&P2->x);
    Fp4_sub(&x,&tmp,&x);
    Fp4_sub(&tmp,&P1->x,&x);
    Fp4_set(&t_ans.x,&x);
    Fp4_mul(&tmp,&tmp,&lambda);
    Fp4_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp4_set(ANS,&t_ans);
    
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&lambda);
    Fp4_clear(&tmp);
    EFp4_clear(&t_ans);
}
int EFp4_cmp(struct EFp4 *A,struct EFp4 *B){
    if(Fp4_cmp(&A->x,&B->x)==0 && Fp4_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}


//-----------------------------------------------------------------------------------------
// #pragma mark EFp8 methods
void EFp8_init(struct EFp8 *A){
    Fp8_init(&A->x);
    Fp8_init(&A->y);
    A->infity=FALSE;
}
void EFp8_set(struct EFp8 *A,struct EFp8 *B){
    Fp8_set(&A->x,&B->x);
    Fp8_set(&A->y,&B->y);
    A->infity=B->infity;
}
void EFp8_set_infity(struct EFp8 *A){
    Fp8_set_ui(&A->x,0);
    Fp8_set_ui(&A->y,0);
    A->infity=TRUE;
}
void EFp8_clear(struct EFp8 *A){
    Fp8_clear(&A->x);
    Fp8_clear(&A->y);
}
void EFp8_printf(struct EFp8 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd",A->x.x0.x0.x0.x0,A->x.x0.x0.x1.x0,A->x.x0.x1.x0.x0,A->x.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd)\n",A->x.x1.x0.x0.x0,A->x.x1.x0.x1.x0,A->x.x1.x1.x0.x0,A->x.x1.x1.x1.x0);
    
    gmp_printf("(%Zd,%Zd,%Zd,%Zd",A->y.x0.x0.x0.x0,A->y.x0.x0.x1.x0,A->y.x0.x1.x0.x0,A->y.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd)\n",A->y.x1.x0.x0.x0,A->y.x1.x0.x1.x0,A->y.x1.x1.x0.x0,A->y.x1.x1.x1.x0);
    
}
void EFp8_SCM_BIN(struct EFp8 *ANS,struct EFp8 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp8 Q,R;
    EFp8_init(&Q);
    EFp8_set(&Q,P);
    EFp8_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp8_ECD(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp8_ECA(&Q,&Q,P);
        }
    }
    EFp8_set(ANS,&Q);
    
    EFp8_clear(&Q);
    EFp8_clear(&R);
    return;
}
void EFp8_ECD(struct EFp8 *ANS, struct EFp8 *P){
    if(P->infity==TRUE){
        EFp8_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp8_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp8_set_infity(ANS);
        return;
    }
    
    struct Fp8 x,y,lambda,tmp;
    struct EFp8 t_ans;
    Fp8_init(&x);
    Fp8_init(&lambda);
    Fp8_init(&tmp);
    Fp8_init(&y);
    EFp8_init(&t_ans);
    
    Fp8_mul(&x,&P->x,&P->x);
    Fp8_add(&tmp,&x,&x);
    Fp8_add(&x,&tmp,&x);
    Fp8_add(&y,&P->y,&P->y);
    Fp8_div(&lambda,&x,&y);
    Fp8_mul(&tmp,&lambda,&lambda);
    Fp8_add(&x,&P->x,&P->x);
    Fp8_sub(&x,&tmp,&x);
    Fp8_sub(&tmp,&P->x,&x);
    Fp8_set(&t_ans.x,&x);
    Fp8_mul(&tmp,&tmp,&lambda);
    Fp8_sub(&t_ans.y,&tmp,&P->y);
    
    EFp8_set(ANS,&t_ans);
    
    Fp8_clear(&x);
    Fp8_clear(&lambda);
    Fp8_clear(&y);
    Fp8_clear(&tmp);
    EFp8_clear(&t_ans);
}
void EFp8_ECA(struct EFp8 *ANS, struct EFp8 *P1, struct EFp8 *P2){
    if(P2->infity==TRUE){//if P2==inf
        EFp8_set(ANS,P1);
        return;
    }
    else if(P1->infity==TRUE){//if P1==inf
        EFp8_set(ANS,P2);
        return;
    }
    else if(Fp8_cmp(&P1->x,&P2->x)==0&&Fp8_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp8_set_infity(ANS);
        return;
    }
    else if(EFp8_cmp(P1,P2)==0){ // P=Q
        EFp8_ECD(ANS,P1);
        return;
    }
    
    struct Fp8 x,y,lambda,tmp;
    struct EFp8 t_ans;
    
    Fp8_init(&x);
    Fp8_init(&y);
    Fp8_init(&lambda);
    Fp8_init(&tmp);
    EFp8_init(&t_ans);
    
    Fp8_sub(&x,&P2->x,&P1->x);
    Fp8_sub(&y,&P2->y,&P1->y);
    Fp8_div(&lambda,&y,&x);
    Fp8_mul(&tmp,&lambda,&lambda);
    Fp8_add(&x,&P1->x,&P2->x);
    Fp8_sub(&x,&tmp,&x);
    Fp8_sub(&tmp,&P1->x,&x);
    Fp8_set(&t_ans.x,&x);
    Fp8_mul(&tmp,&tmp,&lambda);
    Fp8_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp8_set(ANS,&t_ans);
    
    Fp8_clear(&x);
    Fp8_clear(&y);
    Fp8_clear(&lambda);
    Fp8_clear(&tmp);
    EFp8_clear(&t_ans);
}
int EFp8_cmp(struct EFp8 *A,struct EFp8 *B){
    if(Fp8_cmp(&A->x,&B->x)==0 && Fp8_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}

//-----------------------------------------------------------------------------------------
// #pragma mark EFp16 methods
void EFp16_init(struct EFp16 *A){
    Fp16_init(&A->x);
    Fp16_init(&A->y);
    A->infity=FALSE;
}
void EFp16_set(struct EFp16 *A,struct EFp16 *B){
    Fp16_set(&A->x,&B->x);
    Fp16_set(&A->y,&B->y);
    A->infity=B->infity;
}
void EFp16_set_infity(struct EFp16 *A){
    Fp16_set_ui(&A->x,0);
    Fp16_set_ui(&A->y,0);
    A->infity=TRUE;
}
void EFp16_set_EFp(struct EFp16 *A,struct EFp *B){
    Fp16_set_ui(&A->x,0);
    Fp16_set_ui(&A->y,0);
    
    Fp_set(&A->x.x0.x0.x0.x0,&B->x);
    Fp_set(&A->y.x0.x0.x0.x0,&B->y);
    A->infity=B->infity;
}
void EFp16_clear(struct EFp16 *A){
    Fp16_clear(&A->x);
    Fp16_clear(&A->y);
}
void EFp16_printf(struct EFp16 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd,\n",A->x.x0.x0.x0.x0.x0,A->x.x0.x0.x0.x1.x0,A->x.x0.x0.x1.x0.x0,A->x.x0.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x.x0.x1.x0.x0.x0,A->x.x0.x1.x0.x1.x0,A->x.x0.x1.x1.x0.x0,A->x.x0.x1.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x.x1.x0.x0.x0.x0,A->x.x1.x0.x0.x1.x0,A->x.x1.x0.x1.x0.x0,A->x.x1.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd)\n",A->x.x1.x1.x0.x0.x0,A->x.x1.x1.x0.x1.x0,A->x.x1.x1.x1.x0.x0,A->x.x1.x1.x1.x1.x0);
    
    gmp_printf("(%Zd,%Zd,%Zd,%Zd,\n",A->y.x0.x0.x0.x0.x0,A->y.x0.x0.x0.x1.x0,A->y.x0.x0.x1.x0.x0,A->y.x0.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->y.x0.x1.x0.x0.x0,A->y.x0.x1.x0.x1.x0,A->y.x0.x1.x1.x0.x0,A->y.x0.x1.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->y.x1.x0.x0.x0.x0,A->y.x1.x0.x0.x1.x0,A->y.x1.x0.x1.x0.x0,A->y.x1.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd)\n",A->y.x1.x1.x0.x0.x0,A->y.x1.x1.x0.x1.x0,A->y.x1.x1.x1.x0.x0,A->y.x1.x1.x1.x1.x0);
}
void EFp16_SCM_BIN(struct EFp16 *ANS,struct EFp16 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp16 Q,R;
    EFp16_init(&Q);
    EFp16_set(&Q,P);
    EFp16_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp16_ECD(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp16_ECA(&Q,&Q,P);
        }
    }
    EFp16_set(ANS,&Q);
    
    EFp16_clear(&Q);
    EFp16_clear(&R);
    return;
}
void EFp16_ECD(struct EFp16 *ANS, struct EFp16 *P){
    if(P->infity==TRUE){
        EFp16_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp16_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp16_set_infity(ANS);
        return;
    }
    
    struct Fp16 x,y,lambda,tmp;
    struct EFp16 t_ans;
    Fp16_init(&x);
    Fp16_init(&lambda);
    Fp16_init(&tmp);
    Fp16_init(&y);
    EFp16_init(&t_ans);
    
    Fp16_mul(&x,&P->x,&P->x);
    Fp16_add(&tmp,&x,&x);
    Fp16_add(&x,&tmp,&x);
    Fp_add_mpz(&x.x0.x0.x0.x0,&x.x0.x0.x0.x0,a_x);
    Fp16_add(&y,&P->y,&P->y);
    Fp16_div(&lambda,&x,&y);
    Fp16_mul(&tmp,&lambda,&lambda);
    Fp16_add(&x,&P->x,&P->x);
    Fp16_sub(&x,&tmp,&x);
    Fp16_sub(&tmp,&P->x,&x);
    Fp16_set(&t_ans.x,&x);
    Fp16_mul(&tmp,&tmp,&lambda);
    Fp16_sub(&t_ans.y,&tmp,&P->y);
    
    EFp16_set(ANS,&t_ans);
    
    Fp16_clear(&x);
    Fp16_clear(&lambda);
    Fp16_clear(&y);
    Fp16_clear(&tmp);
    EFp16_clear(&t_ans);
}
void EFp16_ECA(struct EFp16 *ANS, struct EFp16 *P1, struct EFp16 *P2){
    if(P2->infity==TRUE){//if P2==inf
        EFp16_set(ANS,P1);
        return;
    }
    else if(P1->infity==TRUE){//if P1==inf
        EFp16_set(ANS,P2);
        return;
    }
    else if(Fp16_cmp(&P1->x,&P2->x)==0&&Fp16_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp16_set_infity(ANS);
        return;
    }
    else if(EFp16_cmp(P1,P2)==0){ // P=Q
        EFp16_ECD(ANS,P1);
        return;
    }
    
    struct Fp16 x,y,lambda,tmp;
    struct EFp16 t_ans;
    
    Fp16_init(&x);
    Fp16_init(&y);
    Fp16_init(&lambda);
    Fp16_init(&tmp);
    EFp16_init(&t_ans);
    
    Fp16_sub(&x,&P2->x,&P1->x);
    Fp16_sub(&y,&P2->y,&P1->y);
    Fp16_div(&lambda,&y,&x);
    Fp16_mul(&tmp,&lambda,&lambda);
    Fp16_add(&x,&P1->x,&P2->x);
    Fp16_sub(&x,&tmp,&x);
    Fp16_sub(&tmp,&P1->x,&x);
    Fp16_set(&t_ans.x,&x);
    Fp16_mul(&tmp,&tmp,&lambda);
    Fp16_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp16_set(ANS,&t_ans);
    
    Fp16_clear(&x);
    Fp16_clear(&y);
    Fp16_clear(&lambda);
    Fp16_clear(&tmp);
    EFp16_clear(&t_ans);
}
int EFp16_cmp(struct EFp16 *A,struct EFp16 *B){
    if(Fp16_cmp(&A->x,&B->x)==0 && Fp16_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}
void EFp16_random_set(struct EFp16 *ANS){
    struct EFp16 P,ans_temp;
    EFp16_init(&P);
    EFp16_init(&ans_temp);
    EFp16_init(&P);
    
    struct Fp16 x,a,tmp16;
    Fp16_init(&a);
    Fp16_init(&x);
    Fp16_init(&tmp16);
    
    //t16=a^16+b^16=((t^2-2p)^2-2p^2)^2-2p^4)^2-2p^8
    mpz_t t2,p2,p22,p4,p8,tmp1,tmp2,t16;
    mpz_init(t2);
    mpz_init(p2);
    mpz_init(p22);
    mpz_init(p4);
    mpz_init(p8);
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(t16);
    
    mpz_pow_ui(tmp1,trace_t,2);//t^2
    mpz_mul_ui(p2,PRIME_P,2);//2p
    mpz_sub(t2,tmp1,p2); //t2=t^2-2p
    mpz_pow_ui(t2,t2,2);//t2=(t^2-2p)^2
    
    mpz_pow_ui(tmp1,PRIME_P,2); //p^2
    mpz_mul_ui(p22,tmp1,2);//2p^2
    mpz_sub(tmp1,t2,p22); // (t^2-2p)^2-2p^2
    mpz_pow_ui(tmp2,tmp1,2);//tmp2=((t^2-2p)^2-2p^2)^2
    
    mpz_pow_ui(tmp1,PRIME_P,4); //p^4
    mpz_mul_ui(p4,tmp1,2);//2p^4
    mpz_sub(tmp1,tmp2,p4); //(((t^2-2p)^2-2p^2)^2-2p^4)
    mpz_pow_ui(tmp2,tmp1,2);//tmp2=(((t^2-2p)^2-2p^2)^2-2p^4)^2
    
    mpz_pow_ui(tmp1,PRIME_P,8); //p^8
    mpz_mul_ui(p8,tmp1,2);//2p^8
    mpz_sub(t16,tmp2,p8);
    
    
    mpz_t r16_div_r2,sEFp_16;
    mpz_init(r16_div_r2);
    mpz_init(sEFp_16);
    mpz_pow_ui(tmp1,PRIME_P,16);
    mpz_add_ui(tmp1,tmp1,1);
    mpz_sub(sEFp_16,tmp1,t16);
    
    mpz_pow_ui(tmp1,order_r,2);
    
    //  printf("r^2 divisible %d\n",(int)mpz_divisible_p(sEFp_16,tmp1));
    mpz_tdiv_q(r16_div_r2,sEFp_16,order_r);
    mpz_tdiv_q(r16_div_r2,r16_div_r2,order_r);
    do{
        Fp16_random(&x);
        Fp16_mul(&a,&x,&x);
        Fp16_mul(&a,&a,&x);//x^3
        //        Fp16_mul_mpz(&tmp16,&x, a_x); //ax
        Fp16_add(&a, &a, &x);//x^3+ax
    }while(Fp16_legendre(&a)!=1);
    Fp16_sqrt(&P.y,&a);
    Fp16_set(&P.x,&x);
    EFp16_SCM_BIN(ANS,&P,r16_div_r2);//R
    
    EFp16_clear(&P);
    Fp16_clear(&a);
    Fp16_clear(&x);
    mpz_clear(t2);
    mpz_clear(p2);
    mpz_clear(p22);
    mpz_clear(p4);
    mpz_clear(p8);
    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(t16);
}
void EFp16_random_set_G1(struct EFp16 *ANS){
    struct EFp P;
    EFp_init(&P);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_tdiv_q(exp,EFp_total,order_r);
    EFp_random_set(&P);
    EFp_SCM_BIN(&P,&P,exp);
    Fp16_set_ui(&ANS->x,0);
    Fp_set(&ANS->x.x0.x0.x0.x0,&P.x);
    Fp16_set_ui(&ANS->y,0);
    Fp_set(&ANS->y.x0.x0.x0.x0,&P.y);
    ANS->infity=0;
    
    EFp_clear(&P);
    mpz_clear(exp);
}
void EFp16_random_set_G2(struct EFp16 *ANS){
    struct EFp16 P,P_frobenius,tmp_EFp16;
    EFp16_init(&P);
    EFp16_init(&P_frobenius);
    EFp16_init(&tmp_EFp16);
    
    EFp16_random_set(&P);
    
    EFp16_Frobenius_map(&P_frobenius,&P);
    Fp16_neg(&tmp_EFp16.y,&P.y);
    Fp16_set(&tmp_EFp16.x,&P.x);
    
    EFp16_ECA(&tmp_EFp16,&tmp_EFp16,&P_frobenius);
    //    EFp16_printf(&tmp_EFp16);
    EFp16_set(ANS,&tmp_EFp16);
    
    EFp16_clear(&P);
    EFp16_clear(&P_frobenius);
    EFp16_clear(&tmp_EFp16);
}
//-----------------------------------------------------------------------------------------
void EFp16_to_EFp4_map(struct EFp4 *ANS,struct EFp16 *A){
    Fp4_set_ui(&ANS->x,0);
    Fp4_set_ui(&ANS->y,0);
    Fp4_set(&ANS->x,&A->x.x0.x1);
    Fp4_set(&ANS->y,&A->y.x1.x1);
    ANS->infity=A->infity;
}

void EFp4_to_EFp16_map(struct EFp16 *ANS,struct EFp4 *A){
    Fp16_set_ui(&ANS->x,0);
    Fp16_set_ui(&ANS->y,0);
    Fp4_set(&ANS->x.x0.x1,&A->x);
    Fp4_set(&ANS->y.x1.x1,&A->y);
    ANS->infity=A->infity;
}
//-----------------------------------------------------------------------------------------
float timedifference_msec(struct timeval t0, struct timeval t1){
    return (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
}
void Fp16_Frobenius_map_old(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp16_pow(&tmp,A,PRIME_P);
    Fp16_set(ANS,&tmp);
    Fp16_clear(&tmp);
}
void EFp16_Frobenius_map(struct EFp16 *ANS,struct EFp16 *A){
    struct EFp16 tmp;
    EFp16_init(&tmp);
    
    Fp16_frobenius_map(&tmp.x,&A->x);
    Fp16_frobenius_map(&tmp.y,&A->y);
    
    EFp16_set(ANS,&tmp);
    
    EFp16_clear(&tmp);
}

void EFp4_Skew_Frobenius_p(struct EFp4 *ANS, struct EFp4 *Qt)
{
    struct EFp4 tmp_ans;
    EFp4_init(&tmp_ans);
    
    struct Fp4 Qt_x, Qt_y;
    Fp4_init(&Qt_x);
    Fp4_init(&Qt_y);
    Fp4_set(&tmp_ans.x, &Qt->x);
    Fp4_set(&tmp_ans.y, &Qt->y);
    
    
    Fp_mul(&Qt_x.x0.x0, &tmp_ans.x.x0.x1,&p_M_C_C8);
    Fp_mul(&Qt_x.x0.x1,&tmp_ans.x.x0.x0,&p_C8);
    Fp_mul(&Qt_x.x1.x0, &tmp_ans.x.x1.x1, &p_M_C_C4_C8);
    Fp_mul(&Qt_x.x1.x1, &tmp_ans.x.x1.x0, &p_C4_C16);
    
    Fp_mul(&Qt_y.x0.x0, &tmp_ans.y.x1.x1,&p_M_C_C_C4_C8_C16);
    Fp_mul(&Qt_y.x0.x1,&tmp_ans.y.x1.x0,&p_C4_C8_C16);
    Fp_mul(&Qt_y.x1.x0, &tmp_ans.y.x0.x0, &p_C_C8_C16);
    Fp_mul(&Qt_y.x1.x1, &tmp_ans.y.x0.x1, &p_M_C_C8_C16);
    
    Fp4_set(&tmp_ans.x, &Qt_x);
    Fp4_set(&tmp_ans.y, &Qt_y);
    
    EFp4_set(ANS,&tmp_ans);
    EFp4_clear(&tmp_ans);
    Fp4_clear(&Qt_x);
    Fp4_clear(&Qt_y);
}

void EFp4_Skew_Frobenius_p2(struct EFp4 *ANS, struct EFp4 *Qt){
    struct EFp4 ANS_tmp;
    EFp4_init(&ANS_tmp);
    EFp4_set(&ANS_tmp, Qt);
    
    Fp_mul(&ANS_tmp.x.x0.x0, &ANS_tmp.x.x0.x0,&p2_C8);
    Fp_mul(&ANS_tmp.x.x0.x1, &ANS_tmp.x.x0.x1,&p2_C8);
    Fp_mul(&ANS_tmp.x.x1.x0, &ANS_tmp.x.x1.x0,&p2_C4_C8);
    Fp_mul(&ANS_tmp.x.x1.x1, &ANS_tmp.x.x1.x1,&p2_C4_C8);
    
    Fp_mul(&ANS_tmp.y.x0.x0, &Qt->y.x0.x1,&p2_C_C8_C16);
    Fp_mul(&ANS_tmp.y.x0.x1, &Qt->y.x0.x0,&p2_C8_C16);
    Fp_mul(&ANS_tmp.y.x1.x0, &Qt->y.x1.x1,&p2_C_C4_C8_C16);
    Fp_mul(&ANS_tmp.y.x1.x1, &Qt->y.x1.x0,&p2_C4_C8_C16);
    
    EFp4_set(ANS, &ANS_tmp);
    
    EFp4_clear(&ANS_tmp);
}

void EFp4_Skew_Frobenius_p3(struct EFp4 *ANS, struct EFp4 *Qt){
    struct EFp4 tmp_ans;
    EFp4_init(&tmp_ans);
    
    Fp_mul(&tmp_ans.x.x0.x0, &Qt->x.x0.x1, &p3_M_C_C8);
    Fp_mul(&tmp_ans.x.x0.x1, &Qt->x.x0.x0, &p3_C8);
    Fp_mul(&tmp_ans.x.x1.x0, &Qt->x.x1.x1, &p3_M_C_C4_C8);
    Fp_mul(&tmp_ans.x.x1.x1, &Qt->x.x1.x0, &p3_C4_C8);
    
    Fp_mul(&tmp_ans.y.x0.x0, &Qt->y.x1.x0, &p3_C_C4_C8_C16);
    Fp_mul(&tmp_ans.y.x0.x1, &Qt->y.x1.x1, &p3_M_C_C4_C8_C16);
    Fp_mul(&tmp_ans.y.x1.x0, &Qt->y.x0.x1, &p3_M_C_C8_C16);
    Fp_mul(&tmp_ans.y.x1.x1, &Qt->y.x0.x0, &p3_C8_C16);
    
    EFp4_set(ANS,&tmp_ans);
    ANS->infity = Qt->infity;
    EFp4_clear(&tmp_ans);
}

void EFp4_Skew_Frobenius_p4(struct EFp4 *ANS,struct EFp4 *P){
    //x
    Fp_mul(&ANS->x.x0.x0, &P->x.x0.x0, &p4_C8);
    Fp_mul(&ANS->x.x0.x1, &P->x.x0.x1, &p4_C8);
    Fp_mul(&ANS->x.x1.x0, &P->x.x1.x0, &p4_C4_C8);
    Fp_mul(&ANS->x.x1.x1, &P->x.x1.x1, &p4_C4_C8);
    //y
    Fp_mul(&ANS->y.x0.x0, &P->y.x0.x0, &p4_C8_C16);
    Fp_mul(&ANS->y.x0.x1, &P->y.x0.x1, &p4_C8_C16);
    Fp_mul(&ANS->y.x1.x0, &P->y.x1.x0, &p4_C4_C8_C16);
    Fp_mul(&ANS->y.x1.x1, &P->y.x1.x1, &p4_C4_C8_C16);
    
    ANS->infity=P->infity;
}

void EFp4_Skew_Frobenius_p5(struct EFp4 *ANS,struct EFp4 *P){
    struct EFp4 tmp_ans;
    EFp4_init(&tmp_ans);
    
    //a4-a7
    Fp_mul(&tmp_ans.x.x0.x0, &P->x.x0.x1, &p5_M_C8);
    Fp_mul(&tmp_ans.x.x0.x1, &P->x.x0.x0, &p5_C8);
    Fp_mul(&tmp_ans.x.x1.x0, &P->x.x1.x1, &p5_M_C_C4_C8);
    Fp_mul(&tmp_ans.x.x1.x1, &P->x.x1.x0, &p5_C4_C16);
    
    //a12-a15
    Fp_mul(&tmp_ans.y.x0.x0,  &P->y.x1.x1, &p5_M_C_C_C4_C8_C16);
    Fp_mul(&tmp_ans.y.x0.x1,  &P->y.x1.x0, &p5_C4_C8_C16);
    Fp_mul(&tmp_ans.y.x1.x0,  &P->y.x0.x0, &p5_C_C8_C16);
    Fp_mul(&tmp_ans.y.x1.x1,  &P->y.x0.x1, &p5_M_C_C8_C16);
    
    EFp4_set(ANS,&tmp_ans);
    ANS->infity = P->infity;
    EFp4_clear(&tmp_ans);
}

void EFp4_Skew_Frobenius_p6(struct EFp4 *ANS,struct EFp4 *P){
    struct EFp4 tmp_ans;
    EFp4_init(&tmp_ans);
    //a4-a7
    Fp_mul(&tmp_ans.x.x0.x0, &P->x.x0.x0, &p6_C8);
    Fp_mul(&tmp_ans.x.x0.x1, &P->x.x0.x1, &p6_C8);
    Fp_mul(&tmp_ans.x.x1.x0, &P->x.x1.x0, &p6_M_C8);
    Fp_mul(&tmp_ans.x.x1.x1, &P->x.x1.x1, &p6_M_C8);
    
    //a12-a15
    Fp_mul(&tmp_ans.y.x0.x0, &P->y.x0.x1, &p6_C_C8_C16);
    Fp_mul(&tmp_ans.y.x0.x1, &P->y.x0.x0, &p6_C8_C16);
    Fp_mul(&tmp_ans.y.x1.x0, &P->y.x1.x1, &p6_M_C_C8_C16);
    Fp_mul(&tmp_ans.y.x1.x1, &P->y.x1.x0, &p6_M_C8_C16);
    
    EFp4_set(ANS,&tmp_ans);
    ANS->infity = P->infity;
    EFp4_clear(&tmp_ans);
}

void EFp4_Skew_Frobenius_p7(struct EFp4 *ANS,struct EFp4 *P){
    struct EFp4 tmp_ans;
    EFp4_init(&tmp_ans);
    //5-8
    Fp_mul(&tmp_ans.x.x0.x0, &P->x.x0.x1, &p7_M_C_C8); //-cc8 a5
    Fp_mul(&tmp_ans.x.x0.x1, &P->x.x0.x0, &p7_C8);//c8 a4
    Fp_mul(&tmp_ans.x.x1.x0, &P->x.x1.x1, &p7_M_C_C4_C8); // -c8c4c a7
    Fp_mul(&tmp_ans.x.x1.x1, &P->x.x1.x0, &p7_C8_C4); //c8c4
    //from 9-16
    Fp_mul(&tmp_ans.y.x0.x0, &P->y.x1.x0, &p7_C_C4_C8_C16);//cc4c8c16
    Fp_mul(&tmp_ans.y.x0.x1, &P->y.x1.x1, &p7_M_C_C4_C8_C16);//-cc4c8c16
    Fp_mul(&tmp_ans.y.x1.x0, &P->y.x0.x1, &p7_M_C_C8_C16);//-cc8c16
    Fp_mul(&tmp_ans.y.x1.x1, &P->y.x0.x0, &p7_C8_C16);//c8c16
    
    EFp4_set(ANS,&tmp_ans);
    ANS->infity = P->infity;
    EFp4_clear(&tmp_ans);
}

//-----------------------------------------------------------------------------------------
void KSS_16_parameters(void){
    
    mpz_t tmp1,tmp2,two;
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(two);
    
    
    //set p,r
    mpz_t p_tmp,r_tmp,t_tmp;
    mpz_t xpow2,xpow4,xpow5,xpow6,xpow8,xpow9,xpow10;
    //    mpz_t tmp1,tmp2;
    
    mpz_init(p_tmp);
    mpz_init(r_tmp);
    mpz_init(t_tmp);
    mpz_init(xpow2);
    mpz_init(xpow4);
    mpz_init(xpow5);
    mpz_init(xpow6);
    mpz_init(xpow8);
    mpz_init(xpow9);
    mpz_init(xpow10);
    mpz_init(tmp1);
    mpz_init(tmp2);
    
    mpz_mul(xpow2,X,X);
    mpz_mul(xpow4,xpow2,xpow2);
    mpz_mul(xpow5,xpow4,X);
    mpz_mul(xpow6,xpow5,X);
    mpz_mul(xpow8,xpow6,xpow2);
    mpz_mul(xpow9,xpow8,X);
    mpz_mul(xpow10,xpow9,X);
    
    //t=1/35(2x^5+41x+35)
    mpz_mul_ui(tmp1,X,41);
    mpz_add_ui(tmp1,tmp1,35);
    mpz_mul_ui(tmp2,xpow5,2);
    mpz_add(t_tmp,tmp1,tmp2);
    
    mpz_div_ui(trace_t,t_tmp,35);
    
    //    gmp_printf ("trace = %Zd\n",trace_t);
    //r=x^8+48x^4+625
    mpz_mul_ui(tmp1,xpow4,48);
    mpz_add_ui(r_tmp,xpow8,625);
    mpz_add(tmp2,tmp1,r_tmp);
    mpz_tdiv_q_ui(order_r,tmp2,61250);
    //     mpz_tdiv_q_ui(order_r,tmp2,49);
    //     mpz_set(order_r,tmp2);
    //    gmp_printf ("order = %Zd\n",order_r);
    // mpz_set(r,r_tmp);
    
    //p=1/980(x^10+2x^9+5x^8+48x^6+152x^5+240x^4+625x^2+2398x+3125)
    mpz_mul_ui(tmp1,xpow9,2);
    mpz_add(p_tmp,tmp1,xpow10);
    mpz_mul_ui(tmp1,xpow8,5);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,xpow6,48);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,xpow5,152);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,xpow4,240);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,xpow2,625);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,X,2398);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_add_ui(p_tmp,p_tmp,3125);
    
    mpz_div_ui(PRIME_P,p_tmp,980);
    // mpz_set(p,p_tmp);
    
    mpz_t mod, p1;
    mpz_init(mod);
    mpz_init(p1);
    
//    mpz_pow_ui(p1,PRIME_P,3);
//    mpz_mul_ui(p1,p1,3);
//    mpz_sub_ui(p1,p1,5);
    
//    int l = (int)mpz_divisible_p(p1,order_r);
//    int k = (int)mpz_divisible_ui_p(p1,8);
//    gmp_printf("\n\n p^4-1/8 = %d\n",k);
    //
    //        mpz_set_ui(mod,8);
    //
    //        mpz_mod(mod,p1,order_r);
    //        gmp_printf("\n == %Zd\n",mod);
    //
    //    gmp_printf("p:%Zd\n",prime);
    
    mpz_add_ui(order_EFp,PRIME_P,1);
    mpz_sub(order_EFp,order_EFp,trace_t);
    
    if(mpz_probab_prime_p(PRIME_P,25)==0){
        gmp_printf("p: %Zd\n",PRIME_P);
        printf("not  prime number!\n");
        exit(0);
    }
    
    //    if (mpz_divisible_p(order_EFp,order_r) == 0) {
    //        printf("Not Divisible \n");
    //    }
    //    else{
    //        printf("Divisible \n");
    //    }
    
    struct EFp P,ANS;
    int legendre;
    struct Fp rhs,tmp_ax,x;
    mpz_init(tmp_a);
    Fp_init(&rhs);
    Fp_init(&tmp_ax);
    EFp_init(&P);
    EFp_init(&ANS);
    Fp_init(&x);
    mpz_init(tmp_a);
    mpz_set_ui(tmp_a,0);
    
    for(;;){
        mpz_add_ui(tmp_a,tmp_a,1);
        Fp_set_ui(&x,1);
        legendre=0;
        while(legendre !=1){
            mpz_powm_ui(rhs.x0,x.x0,3,PRIME_P);
            //            gmp_printf("tmp %Zd\n",tmp_a);
            mpz_mul(tmp_ax.x0,x.x0,tmp_a);
            Fp_add(&rhs, &rhs, &tmp_ax);
            if((legendre = mpz_legendre(rhs.x0,PRIME_P))==1){
                //gmp_printf("a in while = %Zd\n",rhs.x0);
                Fp_printf(&rhs);
                Fp_sqrt(&P.y,&rhs);
                Fp_set(&P.x,&x);
                EFp_SCM_BIN(&ANS,&P,order_EFp);
                //                printf("SCM  ==\n");
                //                EFp_printf(&ANS);
                if(ANS.infity == TRUE){
                    mpz_set(a_x,tmp_a);
                    // mpz_clear(tmp_a);
                    Fp_clear(&rhs);
                    Fp_clear(&x);
                    EFp_clear(&P);
                    EFp_clear(&ANS);
                    return;
                }
            }
            Fp_add_ui(&x,&x,1);
        }
    }
    return;
}
//-----------------------------------------------------------------------------------------
void EFp4_SCM_ML(struct EFp4 *RES, struct EFp4 *P,mpz_t scalar){
    int i,length;
    length= (int)mpz_sizeinbase(scalar,2);
    char r_binary[length];
    mpz_get_str(r_binary,2,scalar);
    struct EFp4 T0,T1;
    EFp4_init(&T0);
    
    EFp4_set_infity(&T0);
    EFp4_init(&T1);
    EFp4_set(&T1, P);
    //    EFp_ECD(&T1,&T1);
    for(i=0;r_binary[i]!='\0';i++){
        if(r_binary[i]=='1'){
            EFp4_ECA(&T0,&T0,&T1);
            EFp4_ECD(&T1,&T1);
        }
        else if(r_binary[i]=='0'){
            EFp4_ECA(&T1,&T0,&T1);
            EFp4_ECD(&T0,&T0);
        }
    }
    EFp4_set(RES,&T0);
    EFp4_clear(&T0);
    EFp4_clear(&T1);
    return;
}
void EFp4_SCM_WIN(struct EFp4 *ANS, struct EFp4 *P, mpz_t scalar){
    int i,length;
    length= (int)mpz_sizeinbase(scalar,2);
    char r_binary[length];
    mpz_get_str(r_binary,2,scalar);
    
    int win_size = 4;
    int rem = length % win_size;
    int mul3 = length - rem;
    
    
    struct EFp4 R,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,ANS_temp;
    EFp4_init(&R);
    EFp4_init(&R2);
    EFp4_init(&R3);
    EFp4_init(&R4);
    EFp4_init(&R5);
    EFp4_init(&R6);
    EFp4_init(&R7);
    EFp4_init(&R8);
    EFp4_init(&R9);
    EFp4_init(&R10);
    EFp4_init(&R11);
    EFp4_init(&R12);
    EFp4_init(&R13);
    EFp4_init(&R14);
    EFp4_init(&R15);
    EFp4_init(&ANS_temp);
    
    EFp4_set_infity(&ANS_temp);
    EFp4_set(&R, P);
    EFp4_ECA(&R2, &R, &R);
    EFp4_ECA(&R3, &R2, &R);
    EFp4_ECA(&R4, &R3, &R);
    EFp4_ECA(&R5, &R4, &R);
    EFp4_ECA(&R6, &R5, &R);
    EFp4_ECA(&R7, &R6, &R);
    EFp4_ECA(&R8, &R7, &R);
    EFp4_ECA(&R9, &R8, &R);
    EFp4_ECA(&R10, &R9, &R);
    EFp4_ECA(&R11, &R10, &R);
    EFp4_ECA(&R12, &R11, &R);
    EFp4_ECA(&R13, &R12, &R);
    EFp4_ECA(&R14, &R13, &R);
    EFp4_ECA(&R15, &R14, &R);
    
    for(i=0; i< mul3;i=i+4){
        
        EFp4_ECD(&ANS_temp,&ANS_temp);
        EFp4_ECD(&ANS_temp,&ANS_temp);
        EFp4_ECD(&ANS_temp,&ANS_temp);
        EFp4_ECD(&ANS_temp,&ANS_temp);
        
        if (i+2 < length)
        {
            if(r_binary[i] =='1' && r_binary[i+1] =='1' && r_binary[i+2] =='1' && r_binary[i+3] =='1'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R15);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='1' && r_binary[i+2] =='1' && r_binary[i+3] =='0'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R14);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='1' && r_binary[i+2] =='0' && r_binary[i+3] =='1'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R13);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='1' && r_binary[i+2] =='0' && r_binary[i+3] =='0'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R12);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='0' && r_binary[i+2] =='1' && r_binary[i+3] =='1'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R11);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='0' && r_binary[i+2] =='1' && r_binary[i+3] =='0'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R10);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='0' && r_binary[i+2] =='0' && r_binary[i+3] =='1'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R9);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='0' && r_binary[i+2] =='0' && r_binary[i+3] =='0'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R8);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='1' && r_binary[i+2] =='1' && r_binary[i+3] =='1'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R7);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='1' && r_binary[i+2] =='1' && r_binary[i+3] =='0'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R6);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='1' && r_binary[i+2] =='0' && r_binary[i+3] =='1'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R5);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='1' && r_binary[i+2] =='0' && r_binary[i+3] =='0'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R4);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='0' && r_binary[i+2] =='1' && r_binary[i+3] =='1'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R3);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='0' && r_binary[i+2] =='1' && r_binary[i+3] =='0'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R2);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='0' && r_binary[i+2] =='0' && r_binary[i+3] =='1'){
                EFp4_ECA(&ANS_temp, &ANS_temp, &R);
            }
            else{
            }
        }
    }
    
    for(i = mul3;r_binary[i]!='\0';i++){
        EFp4_ECD(&ANS_temp,&ANS_temp);
        if(r_binary[i]=='1'){
            EFp4_ECA(&ANS_temp,&ANS_temp,P);
        }
    }
    EFp4_set(ANS,&ANS_temp);
    EFp4_clear(&ANS_temp);
    EFp4_clear(&R);
    EFp4_clear(&R2);
    EFp4_clear(&R3);
    EFp4_clear(&R4);
    EFp4_clear(&R5);
    EFp4_clear(&R6);
    EFp4_clear(&R7);
    return;
}

void EFp16_SCM_WIN(struct EFp16 *ANS, struct EFp16 *P, mpz_t scalar){
    int i,length;
    length= (int)mpz_sizeinbase(scalar,2);
    char r_binary[length];
    mpz_get_str(r_binary,2,scalar);
    
    int win_size = 4;
    int rem = length % win_size;
    int mul3 = length - rem;
    
    
    struct EFp16 R,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,ANS_temp;
    EFp16_init(&R);
    EFp16_init(&R2);
    EFp16_init(&R3);
    EFp16_init(&R4);
    EFp16_init(&R5);
    EFp16_init(&R6);
    EFp16_init(&R7);
    EFp16_init(&R8);
    EFp16_init(&R9);
    EFp16_init(&R10);
    EFp16_init(&R11);
    EFp16_init(&R12);
    EFp16_init(&R13);
    EFp16_init(&R14);
    EFp16_init(&R15);
    EFp16_init(&ANS_temp);
    
    EFp16_set_infity(&ANS_temp);
    EFp16_set(&R, P);
    EFp16_ECA(&R2, &R, &R);
    EFp16_ECA(&R3, &R2, &R);
    EFp16_ECA(&R4, &R3, &R);
    EFp16_ECA(&R5, &R4, &R);
    EFp16_ECA(&R6, &R5, &R);
    EFp16_ECA(&R7, &R6, &R);
    EFp16_ECA(&R8, &R7, &R);
    EFp16_ECA(&R9, &R8, &R);
    EFp16_ECA(&R10, &R9, &R);
    EFp16_ECA(&R11, &R10, &R);
    EFp16_ECA(&R12, &R11, &R);
    EFp16_ECA(&R13, &R12, &R);
    EFp16_ECA(&R14, &R13, &R);
    EFp16_ECA(&R15, &R14, &R);
    
    for(i=0; i< mul3;i=i+4){
        
        EFp16_ECD(&ANS_temp,&ANS_temp);
        EFp16_ECD(&ANS_temp,&ANS_temp);
        EFp16_ECD(&ANS_temp,&ANS_temp);
        EFp16_ECD(&ANS_temp,&ANS_temp);
        
        if (i+2 < length)
        {
            if(r_binary[i] =='1' && r_binary[i+1] =='1' && r_binary[i+2] =='1' && r_binary[i+3] =='1'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R15);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='1' && r_binary[i+2] =='1' && r_binary[i+3] =='0'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R14);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='1' && r_binary[i+2] =='0' && r_binary[i+3] =='1'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R13);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='1' && r_binary[i+2] =='0' && r_binary[i+3] =='0'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R12);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='0' && r_binary[i+2] =='1' && r_binary[i+3] =='1'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R11);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='0' && r_binary[i+2] =='1' && r_binary[i+3] =='0'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R10);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='0' && r_binary[i+2] =='0' && r_binary[i+3] =='1'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R9);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='0' && r_binary[i+2] =='0' && r_binary[i+3] =='0'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R8);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='1' && r_binary[i+2] =='1' && r_binary[i+3] =='1'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R7);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='1' && r_binary[i+2] =='1' && r_binary[i+3] =='0'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R6);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='1' && r_binary[i+2] =='0' && r_binary[i+3] =='1'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R5);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='1' && r_binary[i+2] =='0' && r_binary[i+3] =='0'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R4);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='0' && r_binary[i+2] =='1' && r_binary[i+3] =='1'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R3);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='0' && r_binary[i+2] =='1' && r_binary[i+3] =='0'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R2);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='0' && r_binary[i+2] =='0' && r_binary[i+3] =='1'){
                EFp16_ECA(&ANS_temp, &ANS_temp, &R);
            }
            else{
            }
        }
    }
    
    for(i = mul3;r_binary[i]!='\0';i++){
        EFp16_ECD(&ANS_temp,&ANS_temp);
        if(r_binary[i]=='1'){
            EFp16_ECA(&ANS_temp,&ANS_temp,P);
        }
    }
    EFp16_set(ANS,&ANS_temp);
    EFp16_clear(&ANS_temp);
    EFp16_clear(&R);
    EFp16_clear(&R2);
    EFp16_clear(&R3);
    EFp16_clear(&R4);
    EFp16_clear(&R5);
    EFp16_clear(&R6);
    EFp16_clear(&R7);
    return;
}

void EFp16_SCM_ML(struct EFp16 *ANS, struct EFp16 *P,mpz_t scalar){
    count_eca_ML_map = 0;
    count_ecd_ML_map =0;
    int i,length;
    length= (int)mpz_sizeinbase(scalar,2);
    char r_binary[length];
    mpz_get_str(r_binary,2,scalar);
    struct EFp16 T0,T1;
    EFp16_init(&T0);
    
    EFp16_set_infity(&T0);
    EFp16_init(&T1);
    EFp16_set(&T1, P);
    //    EFp_ECD(&T1,&T1);
    for(i=0;r_binary[i]!='\0';i++){
        if(r_binary[i]=='1'){
            EFp16_ECA(&T0,&T0,&T1);
            EFp16_ECD(&T1,&T1);
            count_eca_ML_map ++;
            count_ecd_ML_map ++;
        }
        else if(r_binary[i]=='0'){
            EFp16_ECA(&T1,&T0,&T1);
            EFp16_ECD(&T0,&T0);
            count_eca_ML_map ++;
            count_ecd_ML_map ++;
        }
    }
    EFp16_set(ANS,&T0);
    EFp16_clear(&T0);
    EFp16_clear(&T1);
    return;
}


void Check_SCM() {
    
    struct EFp4 ans_temp_EFp4,twisted_Q;
    EFp4_init(&ans_temp_EFp4);
    EFp4_init(&twisted_Q);
    
    struct EFp16 temp_EFp16, QinG2;
    EFp16_init(&temp_EFp16);
    EFp16_init(&QinG2);
    
    EFp16_random_set_G2(&QinG2);
    EFp16_printf(&QinG2);
    
    EFp16_to_EFp4_map(&twisted_Q, &QinG2);
    EFp4_printf(&twisted_Q);
    count_eca_new_pre=0; count_ecd_new_pre=0;
    
    //===================== SM test============
    mpz_t s;
    mpz_init(s);
    
    
    
    FILE *fp;
    fp=fopen("KSS16_128_JICCE.csv","a");
    
    if(fp == NULL){
        printf("Couldn't open file\n");
    }gmp_fprintf(fp,"S,HW,s_size,NEW_ECA,NEW_ECD,BIN_ECA,BIN_ECD,BIN_ECA_MAP,BIN_ECD_MAP,WIN_ECA,WIN_ECD,WIN_ECA_MAP,WIN_ECD_MAP,NAF_ECA,NAF_ECD,NAF_ECA_MAP,NAF_ECD_MAP,ML_ECA,ML_ECD,ML_ECA_MAP,ML_ECD_MAP,NEW_time,BIN_time,BIN_time_map,WIN_time,WIN_time_map,NAF_time,NAF_time_map,ML_time,ML_time_map\n");
    
    
    int k = 0;
    for (k = 0; k < 100; k++) {
        
        mpz_random(s,10);
        mpz_mod(s,s,order_r);
        int length = (int)mpz_sizeinbase(s,2);
        
        struct timeval t0;
        struct timeval t1;
        float elapsed_new = 0.0;
        float elapsed_BIN;
        float elapsed_BIN_map;
        float elapsed_WIN;
        float elapsed_WIN_map;
        float elapsed_NAF = 0.0;
        float elapsed_NAF_map = 0.0;
        float elapsed_ML;
        float elapsed_ML_map;
        
        
        printf("\n Binary method with mapping\n");
        gettimeofday(&t0, 0);
        EFp4_SCM_BIN_Sparse(&ans_temp_EFp4, &twisted_Q, s);
        gettimeofday(&t1, 0);
        EFp4_to_EFp16_map(&temp_EFp16, &ans_temp_EFp4);
        elapsed_BIN_map = timedifference_msec(t0, t1);
        printf("Old count ECA = %d, ECD = %d, Elasped old = %.2f\n",count_eca_BIN_map,count_ecd_BIN_map,elapsed_BIN_map);
        EFp16_printf(&temp_EFp16);
        
        printf("\n BIN method mapping\n");
        gettimeofday(&t0, 0);
        EFp16_SCM_BIN(&temp_EFp16, &QinG2, s);
        gettimeofday(&t1, 0);
        elapsed_BIN = timedifference_msec(t0, t1);
        printf("Count ECA = %d, ECD = %d, Elasped old = %.2f\n",count_eca_BIN,count_ecd_BIN,elapsed_BIN);
        EFp16_printf(&temp_EFp16);
        
        printf("\n  Window method map\n");
        gettimeofday(&t0, 0);
        EFp4_SCM_WIN(&ans_temp_EFp4, &twisted_Q, s);
        gettimeofday(&t1, 0);
        elapsed_WIN_map = timedifference_msec(t0, t1);
        printf("Old count ECA = %d, ECD = %d, Elasped old = %.2f\n",count_eca_WIN_map,count_ecd_WIN_map,elapsed_WIN_map);
        EFp16_printf(&temp_EFp16);
        
        printf("\n WIN method without mapping\n");
        gettimeofday(&t0, 0);
        EFp16_SCM_WIN(&temp_EFp16, &QinG2, s);
        gettimeofday(&t1, 0);
        elapsed_WIN = timedifference_msec(t0, t1);
        printf("Count ECA = %d, ECD = %d, Elasped old = %.2f\n",count_eca_WIN,count_ecd_WIN,elapsed_WIN);
        EFp16_printf(&temp_EFp16);
        
        
        printf("\n ML method wihtout mapping\n");
        gettimeofday(&t0, 0);
        EFp16_SCM_ML(&temp_EFp16, &QinG2, s);
        gettimeofday(&t1, 0);
        elapsed_ML = timedifference_msec(t0, t1);
        printf("Count ECA = %d, ECD = %d, Elasped old = %.2f\n",count_eca_ML,count_ecd_ML,elapsed_ML);
        EFp16_printf(&temp_EFp16);
        
        
        printf("\n ML method with mapping\n");
        gettimeofday(&t0, 0);
        EFp4_SCM_ML(&ans_temp_EFp4, &twisted_Q, s);
        gettimeofday(&t1, 0);
        EFp4_to_EFp16_map(&temp_EFp16, &ans_temp_EFp4);
        elapsed_ML_map = timedifference_msec(t0, t1);
        printf("Count ECA = %d, ECD = %d, Elasped old = %.2f\n",count_eca_ML_map,count_ecd_ML_map,elapsed_ML_map);
        EFp16_printf(&temp_EFp16);
        
        
        if(fp == NULL){
            printf("Couldn't open file\n");
        }
        int r_p = count_eca_new + count_eca_new_pre;
        //        gmp_fprintf(fp,"%Zd,%ld,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n",s, mpz_popcount(s),length,r_p, count_ecd_new, count_eca_BIN, count_ecd_BIN,count_eca_WIN,count_ecd_WIN,count_eca_NAF,count_ecd_NAF, elapsed_new,elapsed_BIN,elapsed_WIN,elapsed_NAF);
        
        gmp_fprintf(fp,"%Zd,%ld,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n\n\n",s, mpz_popcount(s),length,r_p, count_ecd_new, count_eca_BIN, count_ecd_BIN,count_eca_BIN_map,count_ecd_BIN_map,count_eca_WIN,count_ecd_WIN,count_eca_WIN_map,count_ecd_WIN_map,count_eca_NAF,count_ecd_NAF,count_eca_NAF_map,count_ecd_NAF_map,count_eca_ML,count_ecd_ML,count_eca_ML_map,count_ecd_ML_map,elapsed_new,elapsed_BIN,elapsed_BIN_map,elapsed_WIN,elapsed_WIN_map,elapsed_NAF,elapsed_NAF_map,elapsed_ML,elapsed_ML_map);
    }
    fclose(fp);
    
    EFp4_clear(&ans_temp_EFp4);
    EFp4_clear(&twisted_Q);
    EFp16_clear(&temp_EFp16);
    mpz_clear(s);
    
}

//-----------------------------------------------------------------------------------------
void Miller_algo(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q, mpz_t loop){
    struct Fp16 l_sum,v_sum;
    Fp16_init(&l_sum);
    Fp16_init(&v_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    Fp_set_ui(&v_sum.x0.x0.x0.x0,1);
    
    struct EFp16 T;
    EFp16_init(&T);
    EFp16_set(&T,P);
    
    
    struct Fp16 ltt,ltp,v2t,vtp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    Fp16_init(&v2t);
    Fp16_init(&vtp);
    
    int i;
    struct Fp16 tmp1;
    Fp16_init(&tmp1);
    //    Fp16_init(&lambda);
    int r_bit;//bit数
    
    r_bit= (int)mpz_sizeinbase(loop,2);
    
    for(i=r_bit-2;i>=0;i--){
        Fp16_mul(&l_sum,&l_sum,&l_sum);
        Fp16_mul(&v_sum,&v_sum,&v_sum);
        
        ltt_q(&ltt,&T,Q);
        Fp16_mul(&l_sum,&l_sum,&ltt);
        
        EFp16_ECD(&T,&T);
        v2t_q(&v2t,&T,Q);
        Fp16_mul(&v_sum,&v_sum,&v2t);
        
        if(mpz_tstbit(loop,i)==1){
            ltp_q(&ltp,&T,P,Q);
            Fp16_mul(&l_sum,&l_sum,&ltp);
            
            EFp16_ECA(&T,&T,P);
            vtp_q(&vtp,&T,Q);
            Fp16_mul(&v_sum,&v_sum,&vtp);
        }
    }
    
    
    // EFp16_printf(&T);
    Fp16_div(ANS,&l_sum,&v_sum);
    
    Fp16_clear(&l_sum);
    Fp16_clear(&v_sum);
    EFp16_clear(&T);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    Fp16_clear(&v2t);
    Fp16_clear(&vtp);
    Fp16_clear(&tmp1);
}
void Optimal_Miller(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q, mpz_t loop){
    
    struct EFp16 T,EFp16_tmp;
    EFp16_init(&T);
    EFp16_init(&EFp16_tmp);
    
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct Fp16 Px_neg;
    Fp16_init(&Px_neg);
    Fp16_neg(&Px_neg,&P->x);//TODO Why neg Px?
    
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    struct EFp16 Q_neg;
    EFp16_init(&Q_neg);
    Fp16_neg(&Q_neg.y,&Q->y);
    Fp16_set(&Q_neg.x,&Q->x);
    
    if(X_bit_binary[x_bit]==-1){
        printf("if \n %d\n",x_bit);
        EFp16_set(&T,&Q_neg);
    }else{
        printf("else \n %d\n",x_bit);
        EFp16_set(&T,Q);
    }
    int i;
    for(i=x_bit-1;i>=0;i--){
        switch (X_bit_binary[i]){
            case 0:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                Fp16_mul(&l_sum,&l_sum,&ltt);
                break;
                
            case 1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                
                DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
            case -1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                
                DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                ADD_LINE(&ltp,&T,&T,&Q_neg,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
        }
    }
    
    //  EFp16_SCM_BIN(&EFp_tmp,Q,prime);
    struct EFp4 Q_bar,EFp4_tmp;
    EFp4_init(&Q_bar);
    EFp4_init(&EFp4_tmp);
    EFp16_to_EFp4_map(&Q_bar, Q);
    EFp4_Skew_Frobenius_p(&EFp4_tmp,&Q_bar);
    EFp4_to_EFp16_map(&EFp16_tmp, &EFp4_tmp);
    //  EFp16_frobenius_map(&EFp_tmp, Q);
    
    
    ltp_q(&ltp,&T,&EFp16_tmp,P);
    Fp16_mul(&l_sum,&l_sum,&ltp);
    
    
    struct Fp16 tmp_f;
    Fp16_init(&tmp_f);
    
    
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    //    Fp16_pow(&l_sum, &tmp_f,prime);
    Fp16_frobenius_map(&l_sum, &tmp_f);
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    
    ltt_q(&ltt,Q,P);
    
    
    Fp16_mul(&l_sum,&tmp_f,&ltt);
    Fp16_set(ANS,&l_sum);
    
    
    Fp16_clear(&l_sum);
    EFp16_clear(&T);
    EFp16_clear(&EFp16_tmp);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    EFp4_clear(&Q_bar);
    EFp4_clear(&EFp4_tmp);
}

//void pre_calc_vector_final_exp(void)
//{
//    printf("\n\n VECTOR PRE_CALC \n");
//    struct timeval t0;
//    struct timeval t1;
//    gettimeofday(&t0, 0);
//
//
//
//    mpz_init(B);
//    mpz_init(A);
//
//    mpz_add_ui(B,X,1);
//    mpz_mul(B,B,B);
//    mpz_add_ui(B,B,4);
//
//
//    mpz_pow_ui(A,X,3);
//    mpz_mul(A,A,B);
//    mpz_add_ui(A,A,56);
//    gmp_printf("A = %Zd\n",A);
//    gmp_printf("B = %Zd\n",B);
//
//
//    mpz_t x2,x3,x4;
//    mpz_inits(x2,x3,x4,(mpz_ptr) NULL);
//    mpz_mul(x2,X,X);
//    mpz_mul(x3,x2,X);
//    mpz_mul(x4,x2,x2);
//
//    mpz_t m00, m11, m22, m33, m44, m55, m66, m77, tmp1;
//    mpz_inits(m00, m11, m22, m33, m44, m55, m66, m77,tmp1, (mpz_ptr) NULL);
//
//    //m00:=2*u^3*A+55*u^2*B;
//    mpz_mul(m00,x3,A);
//    mpz_mul_ui(m00,m00,2);
//    mpz_mul_ui(tmp1,x2,55);
//    mpz_mul(tmp1,tmp1,B);
//    mpz_add(m00,m00,tmp1);
//    //    gmp_printf("m00 = %Zd\n",m00);
//
//    //m11:=-4*u^2*A-75*u*B;
//    mpz_mul(m11,x2,A);
//    mpz_mul_si(m11,m11,-4);
//    mpz_mul(tmp1,X,B);
//    mpz_mul_ui(tmp1,tmp1,75);
//    mpz_sub(m11,m11,tmp1);
//    //    gmp_printf("m11 = %Zd\n",m11);
//    //m22:=-2*u*A-125*B;
//    mpz_mul(m22,X,A);
//    mpz_mul_si(m22,m22,-2);
//    mpz_mul_ui(tmp1,B,125);
//    mpz_sub(m22,m22,tmp1);
//    //    gmp_printf("m22 = %Zd\n",m22);
//    //m33:=-u^4*A-24*u^3*B+196;
//    mpz_mul(tmp1,x4,A);
//    mpz_neg(m33,tmp1);
//    mpz_mul(tmp1,x3,B);
//    mpz_mul_si(tmp1,tmp1,-24);
//    mpz_add(m33,m33,tmp1);
//    mpz_add_ui(m33,m33,196);
//    //    gmp_printf("m33 = %Zd\n",m33);
//
//    //m44:=u^3*A+10*u^2*B;
//    mpz_mul(m44,x3,A);
//    mpz_mul(tmp1,x2,B);
//    mpz_mul_ui(tmp1,tmp1,10);
//    mpz_add(m44,m44,tmp1);
//    //    gmp_printf("m44 = %Zd\n",m44);
//    //m55:=3*u^2*A+100*u*B;
//    mpz_mul(m55,x2,A);
//    mpz_mul_ui(m55,m55,3);
//    mpz_mul(tmp1,X,B);
//    mpz_mul_ui(tmp1,tmp1,100);
//    mpz_add(m55,m55,tmp1);
//    //    gmp_printf("m44 = %Zd\n",m55);
//
//    //m66:=-11*u*A-250*B;
//    mpz_mul(m66,X,A);
//    mpz_mul_si(m66,m66,-11);
//    mpz_mul_si(tmp1,B,-250);
//    mpz_add(m66,m66,tmp1);
//    //    gmp_printf("m66 = %Zd\n",m66);
//    //m77:=7*A;
//    mpz_mul_ui(m77,A,7);
//    //    gmp_printf("m77 = %Zd\n",m77);
//    gettimeofday(&t1, 0);
//    double elapsed = timedifference_msec(t0, t1);
//    printf("FINAL PREC CALC VECTOR ms: %f [ms]\n", elapsed);
//}

void final_exp_hard(struct Fp16 *ANS,struct Fp16 *fd)
{
    mpz_t temp_exp;
    mpz_init(temp_exp);
    
    struct Fp16 t,t0, t1, t2, t3, t4, t5, t6,t7, t8, t9, t10, t11, t12, tmp16, tmp16_1;
    Fp16_init(&t0);
    Fp16_init(&t1);
    Fp16_init(&t2);
    Fp16_init(&t3);
    Fp16_init(&t4);
    Fp16_init(&t5);
    Fp16_init(&t6);
    Fp16_init(&t7);
    Fp16_init(&t8);
    Fp16_init(&t9);
    Fp16_init(&t10);
    Fp16_init(&t11);
    Fp16_init(&t12);
    Fp16_init(&tmp16);
    Fp16_init(&tmp16_1);
    Fp16_init(&t);
    
    
    struct Fp16 t00, t01, t02, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27,t28, t29,t30, t31, t32,t33, t37, s0, s1, s2, s3;
    Fp16_init(&t00);
    Fp16_init(&t01);
    Fp16_init(&t02);
    Fp16_init(&t13);
    Fp16_init(&t14);
    Fp16_init(&t15);
    Fp16_init(&t16);
    Fp16_init(&t17);
    Fp16_init(&t18);
    Fp16_init(&t19);
    Fp16_init(&t20);
    Fp16_init(&t21);
    Fp16_init(&t22);
    Fp16_init(&t23);
    Fp16_init(&t24);
    Fp16_init(&t25);
    Fp16_init(&t26);
    Fp16_init(&t27);
    Fp16_init(&t28);
    Fp16_init(&t29);
    Fp16_init(&t30);
    Fp16_init(&t31);
    Fp16_init(&t32);
    Fp16_init(&t33);
    Fp16_init(&t37);
    Fp16_init(&s0);
    Fp16_init(&s1);
    Fp16_init(&s2);
    Fp16_init(&s3);
    
    Fp16_mul(&t0, fd, fd); //fd=M
    Fp16_mul(&t1, &t0, &t0);
    Fp16_pow_mother_parameter(&t2,fd);
    Fp16_mul(&t2, &t2, fd);
    Fp16_pow_mother_parameter(&t3,&t2);
    Fp16_mul(&t3, &t3, &t2);
    Fp16_mul(&t4, &t3, &t1);
    
    Fp16_pow_mother_parameter(&t5,&t4);
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&t6, &t4,temp_exp);
    mpz_set_ui(temp_exp,8);
    Fp16_pow(&t7, &t1,temp_exp);
    mpz_set_ui(temp_exp,2);
    Fp16_pow(&t8, &t7, temp_exp);
    Fp8_set(&tmp16_1.x0, &t1.x0);
    Fp8_neg(&tmp16_1.x1, &t1.x1);
    Fp16_mul(&t9,&t7,&tmp16_1);
    Fp16_mul(&t10, &t9, &t9);
    Fp16_pow_mother_parameter(&t11, &t5);
    Fp16_pow_mother_parameter(&t12, &t11);
    //t01:=t12*t10;
    Fp16_mul(&t01, &t12, &t10);
    
//    Fp16_printf(&t01);
//
//    mpz_set_str(temp_exp,"4567203436159856725189380193589211078618796570314944",10);
//    //    Fp16_pow(&tmp16, fd, A);
//    Fp16_pow(&tmp16, fd, temp_exp);
//    //    Fp16_invert(&tmp16, &tmp16);
//    Fp8_set(&tmp16_1.x0, &tmp16.x0);
//    Fp8_neg(&tmp16_1.x1, &tmp16.x1);
//    printf("IF t01 eq M^A %d\n",Fp16_cmp(&t01, &tmp16_1));
//    Fp16_printf(&tmp16_1);
    
    
    Fp16_pow_mother_parameter(&t14, &t01);
    Fp16_mul(&tmp16, &t14, &t14);
    
    Fp8_set(&t13.x0, &tmp16.x0);
    Fp8_neg(&t13.x1, &tmp16.x1);
    
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&t00, &t6, temp_exp);
    Fp16_pow(&t15, &t00, temp_exp);
    
    Fp8_set(&tmp16.x0, &t15.x0);
    Fp8_neg(&tmp16.x1, &t15.x1);
    Fp16_mul(&t0, &t13, &tmp16);
    
    Fp16_mul(&t16, &t0, &t0);
    mpz_set_ui(temp_exp,4);
    Fp16_pow(&t17, &t13, temp_exp);
    Fp16_mul(&t18, &t17, &t14);
    
    //t2:=t16*t18;
    Fp16_mul(&t2, &t16, &t18);
    
    //t19:=t14^(u);
    
    Fp16_pow_mother_parameter(&t19, &t14);
    //t20:=t19^(u);
    Fp16_pow_mother_parameter(&t20, &t19);
    //t21:=t20^(u);
    Fp16_pow_mother_parameter(&t21, &t20);
    //t22:=t19^2;
    Fp16_mul(&t22, &t19, &t19);
    
    //t23:=t5^5;
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&t23, &t5, temp_exp);
    //t24:=t23^5;
    Fp16_pow(&t24, &t23, temp_exp);
    //t25:=t24^3;
    mpz_set_ui(temp_exp,3);
    Fp16_pow(&t25, &t24, temp_exp);
    //t26:=t24*t25;
    Fp16_mul(&t26, &t24, &t25);
    //t27:=t22^2;
    Fp16_mul(&t27, &t22, &t22);
    //t37:=(t27*t25)^(-1);
    Fp16_mul(&tmp16, &t27, &t25);
    
    // Fp16_invert(&t37, &tmp16);
    Fp8_set(&t37.x0, &tmp16.x0);
    Fp8_neg(&t37.x1, &tmp16.x1);
    
    
    //t28:=t27*t19^(-1);
    Fp8_set(&tmp16.x0, &t19.x0);
    Fp8_neg(&tmp16.x1, &t19.x1);
    Fp16_mul(&t28, &t27, &tmp16);
    
    //t3:=t28*t26;
    Fp16_mul(&t3, &t28, &t26);
    
    //t29:=t11^5;
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&t29, &t11, temp_exp);
    //t30:=t29^2;
    Fp16_mul(&t30, &t29, &t29);
    //t4:=t20*t30;
    Fp16_mul(&t4, &t20, &t30);
    
    //s0:=t20^2;
    Fp16_mul(&s0, &t20, &t20);
    //s1:=t30^5;
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&s1, &t30, temp_exp);
    //s2:=s1*t29;
    Fp16_mul(&s2, &s1, &t29);
    //s3:=s0*s2;
    Fp16_mul(&s3, &s0, &s2);
    
    //t31:=t12^(24);
    mpz_set_ui(temp_exp,24);
    Fp16_pow(&t31, &t12, temp_exp);
    
    //t5:=t21^(-1)*t31^(-1);
    Fp8_set(&tmp16_1.x0, &t31.x0);
    Fp8_neg(&tmp16_1.x1, &t31.x1);
    
    Fp8_set(&tmp16.x0, &t21.x0);
    Fp8_neg(&tmp16.x1, &t21.x1);
    
    Fp16_mul(&t5, &tmp16, &tmp16_1);
    
    //t6:=t8^3*t1;
    mpz_set_ui(temp_exp,3);
    Fp16_pow(&tmp16, &t8, temp_exp);
    Fp16_mul(&t6, &tmp16,&t1);
    
    //t7:=t5*t6;
    Fp16_mul(&t7, &t5, &t6);
    
    //t8:=t01^7;
    mpz_set_ui(temp_exp,7);
    Fp16_pow(&t8, &t01, temp_exp);
    
    //t32:=t37^p*t7^(p^3);
    Fp16_frobenius_map(&tmp16, &t37);
    Fp16_frobenius_map_p3(&tmp16_1, &t7);
    Fp16_mul(&t32, &tmp16_1, &tmp16);
    
    //t32:=t32*t3^(p^5);
    Fp16_frobenius_map_p5(&tmp16, &t3);
    Fp16_mul(&t32, &t32, &tmp16);
    
    //t32:=t32*t8^(p^7);
    Fp16_frobenius_map_p7(&tmp16, &t8);
    Fp16_mul(&t32, &t32, &tmp16);
    
    //t33:=t0^(p^2)*t2^(p^6);
    Fp16_frobenius_map_p2(&tmp16, &t0);
    Fp16_frobenius_map_p6(&tmp16_1, &t2);
    Fp16_mul(&t33, &tmp16, &tmp16_1);
    // t7=M^m33 and t2=M^m66.
    
    //t:=t32*t33;
    Fp16_mul(&t, &t33, &t32);
    
    //t:=t*t4^(p^4);
    Fp16_frobenius_map_p4(&tmp16, &t4);
    Fp16_mul(&t, &t, &tmp16);
    
    //t:=t*s3;
    Fp16_mul(&t, &t, &s3);
    
    Fp16_set(ANS, &t);
    
    
//    mpz_clears(m00, m11, m22, m33, m44, m55, m66, m77,tmp1, (mpz_ptr) NULL);
//    mpz_clears(x2,x3,(mpz_ptr) NULL);
    
    Fp16_clear(&t0);
    Fp16_clear(&t1);
    Fp16_clear(&t2);
    Fp16_clear(&t3);
    Fp16_clear(&t4);
    Fp16_clear(&t5);
    Fp16_clear(&t6);
    Fp16_clear(&t7);
    Fp16_clear(&t8);
    Fp16_clear(&t9);
    Fp16_clear(&t10);
    Fp16_clear(&t11);
    Fp16_clear(&t12);
    Fp16_clear(&tmp16);
    Fp16_clear(&tmp16_1);
    Fp16_clear(&t);
    Fp16_clear(&t00);
    Fp16_clear(&t01);
    Fp16_clear(&t02);
    Fp16_clear(&t13);
    Fp16_clear(&t14);
    Fp16_clear(&t15);
    Fp16_clear(&t16);
    Fp16_clear(&t17);
    Fp16_clear(&t18);
    Fp16_clear(&t19);
    Fp16_clear(&t20);
    Fp16_clear(&t21);
    Fp16_clear(&t22);
    Fp16_clear(&t23);
    Fp16_clear(&t24);
    Fp16_clear(&t25);
    Fp16_clear(&t26);
    Fp16_clear(&t27);
    Fp16_clear(&t28);
    Fp16_clear(&t29);
    Fp16_clear(&t30);
    Fp16_clear(&t31);
    Fp16_clear(&t32);
    Fp16_clear(&t33);
    Fp16_clear(&t37);
    Fp16_clear(&s0);
    Fp16_clear(&s1);
    Fp16_clear(&s2);
    Fp16_clear(&s3);
}
void Final_Exp(struct Fp16 *ANS,struct Fp16 *A){
    
    struct Fp16 M,A_p8;
    Fp16_init(&M);
    Fp16_init(&A_p8);
    
    printf("\n\n f^p^8-1 \n");
    struct timeval t0;
    struct timeval t1;
    gettimeofday(&t0, 0);
     Init_mpz_Cost(&mpz_cost);
    Fp16_frobenius_map_p8(&A_p8,A);
    Fp16_div(&M, &A_p8, A);
    Print_mpz_Cost(&mpz_cost,"Soft Part:");
    gettimeofday(&t1, 0);
    double elapsed_0 = timedifference_msec(t0, t1);
    printf("FINAL EXP EASY PART: %f [ms]\n", elapsed_0);
    
    
    printf("\n\n M^P^8+1 \n");
    gettimeofday(&t0, 0);
     Init_mpz_Cost(&mpz_cost);
    final_exp_hard(ANS, &M);
    Print_mpz_Cost(&mpz_cost,"Hard Part:");
    //    mpz_t temp_exp;
    //    mpz_init(temp_exp);
    //    mpz_pow_ui(temp_exp,PRIME_P,8);
    //    mpz_add_ui(temp_exp,temp_exp,1);
    //    mpz_tdiv_q(temp_exp,temp_exp,order_r);
    //    Fp16_pow(ANS, &M, temp_exp);
    gettimeofday(&t1, 0);
    double elapsed = timedifference_msec(t0, t1);
    printf("FINAL EXP HARD PART ms: %f [ms]\n", elapsed);
    
    Fp16_clear(&M);
    Fp16_clear(&A_p8);
    
    ////expo:=(14*u^3*(p^8+1)) div (125*r);
    //    mpz_t expo, temp_exp;
    //    mpz_init(expo);
    //    mpz_init(temp_exp);
    //    mpz_pow_ui(expo,X,3);
    ////    mpz_mul_ui(expo,expo,857500);
    //     mpz_mul_ui(expo,expo,14);
    //    mpz_pow_ui(temp_exp,PRIME_P,8);
    //    mpz_add_ui(temp_exp,temp_exp,1);
    //    mpz_mul(expo,expo,temp_exp);
    //    mpz_mul_ui(temp_exp,order_r,125);
    //    mpz_tdiv_q(expo,expo,temp_exp);
    
    //    gmp_printf("exp = %Zd\n",expo);
    //    gmp_printf("p^8+1/r = %Zd\n\n",p8p1dr);
    //    Fp16_pow(ANS, &M, expo);
    
    
}
void Tate_Pairing(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2){
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    
    Miller_algo(&t_ans,G2,G1,order_r);
    Final_Exp(&t_ans,&t_ans);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&t_ans);
}
void Ate_Pairing(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2){
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    
    mpz_t tm1;
    mpz_init(tm1);
    mpz_sub_ui(tm1,trace_t,1);
    
    Miller_algo(&t_ans,G2,G1,tm1);
    Final_Exp(&t_ans,&t_ans);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&t_ans);
}
void Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2){
    struct Fp16 Miller_X, t_ans;
    Fp16_init(&Miller_X);
    Fp16_init(&t_ans);
    
    Optimal_Miller(&Miller_X,G1,G2,X);
    
    Final_Exp(&t_ans,&Miller_X);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&Miller_X);
    Fp16_clear(&t_ans);
}
void ADD_LINE(struct Fp16 *l_ANS,struct EFp16 *T_ANS,struct EFp16 *T,struct EFp16 *P,struct EFp16 *Q,struct Fp16 *Qx_neg){
    struct Fp16 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&tmp4);
    Fp16_init(&lambda);
    Fp16_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp16 x,y,tmp;
    Fp16_init(&x);
    Fp16_init(&y);
    Fp16_init(&tmp);
    
    struct EFp16 x3_tmp;
    EFp16_init(&x3_tmp);
    struct Fp16 A,B,C,D,E,F;
    Fp16_init(&A);
    Fp16_init(&B);
    Fp16_init(&C);
    Fp16_init(&D);
    Fp16_init(&E);
    Fp16_init(&F);
    
    
    Fp16_sub(&A,&P->x,&T->x);//xt-xp
    Fp16_sub(&B,&P->y,&T->y);//yt-yp
    Fp16_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp16_add(&D,&T->x,&P->x);
    Fp16_mul(&tmp1,&C,&C);
    Fp16_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp16_mul(&tmp2,&C,&T->x);
    Fp16_sub(&E,&tmp2,&T->y);
    
    Fp16_mul(&tmp3,&C,&x3_tmp.x);
    Fp16_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp16_set(&l_tmp,&Q->y);
    
    Fp16_add(&l_tmp,&l_tmp,&E);
    
    Fp16_mul(&F,&C,Qx_neg);
    Fp16_add(&l_tmp,&l_tmp,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp16_set(T_ANS,&x3_tmp);
    
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&tmp4);
    Fp16_clear(&lambda);
    Fp16_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp16_clear(&x);
    Fp16_clear(&y);
    Fp16_clear(&tmp);
    EFp16_clear(&x3_tmp);
    Fp16_clear(&A);
    Fp16_clear(&B);
    Fp16_clear(&C);
    Fp16_clear(&D);
    Fp16_clear(&E);
    Fp16_clear(&F);
}
//TODO *Q = *P, *Qx_neg = * Px_neg
void DBL_LINE(struct Fp16 *l_ANS,struct EFp16 *T_ANS,struct EFp16 *T,struct EFp16 *P,struct Fp16 *Px_neg){
    struct Fp16 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&tmp4);
    Fp16_init(&lambda);
    Fp16_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp16 x,y,tmp;
    Fp16_init(&x);
    Fp16_init(&y);
    Fp16_init(&tmp);
    
    struct EFp16 x3_tmp;
    EFp16_init(&x3_tmp);
    struct Fp16 A,B,C,D,E,F;
    Fp16_init(&A);
    Fp16_init(&B);
    Fp16_init(&C);
    Fp16_init(&D);
    Fp16_init(&E);
    Fp16_init(&F);
    
    
    Fp16_add(&A,&T->y,&T->y);//2y
    Fp16_mul(&B,&T->x,&T->x);//x^2
    Fp16_mul_ui(&B,&B,3);//3x^2
    Fp_add_mpz(&B.x0.x0.x0.x0,&B.x0.x0.x0.x0,a_x);//lambda=3x^2+a
    Fp16_div(&C,&B,&A);//lambda=3x^2+a/2y
    
    Fp16_add(&D,&T->x,&T->x); //D=2x
    Fp16_mul(&tmp1,&C,&C);// lamda^2
    Fp16_sub(&x3_tmp.x,&tmp1,&D); //x3.x=lamda^2-x2t
    
    Fp16_mul(&tmp2,&C,&T->x);//xt.lamda
    Fp16_sub(&E,&tmp2,&T->y);//xt.lamda-yt
    
    Fp16_mul(&tmp3,&C,&x3_tmp.x); //x3.lamda
    Fp16_sub(&x3_tmp.y,&E,&tmp3); // x3.y = xt.lamda-yt - x3.lamda
    
    Fp16_set(&l_tmp,&P->y);
    // Fp_set_ui(&l_tmp.x0.x0.x0,1);
    
    Fp16_add(&l_tmp,&l_tmp,&E);
    
    Fp16_mul(&F,&C,Px_neg);
    Fp16_add(&l_tmp,&l_tmp,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp16_set(T_ANS,&x3_tmp);
    
    if(T->infity==TRUE){
        EFp16_set(T_ANS,T);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp16_cmp_mpz(&T->y,cmp)==0){//P.y==0
        EFp16_set_infity(T_ANS);
        return;
    }
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&tmp4);
    Fp16_clear(&lambda);
    Fp16_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp16_clear(&x);
    Fp16_clear(&y);
    Fp16_clear(&tmp);
    EFp16_clear(&x3_tmp);
    Fp16_clear(&A);
    Fp16_clear(&B);
    Fp16_clear(&C);
    Fp16_clear(&D);
    Fp16_clear(&E);
    Fp16_clear(&F);
    mpz_clear(cmp);
}
//-------------------
void ltt_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q){
    struct Fp16 tmp1,tmp2,tmp3,lambda,ltt;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&lambda);
    Fp16_init(&ltt);
    
    Fp16_mul(&tmp1,&T->x,&T->x);//xt^2
    Fp16_add(&tmp2,&tmp1,&tmp1);
    Fp16_add(&tmp1,&tmp2,&tmp1);//3xt^2
    Fp_add_mpz(&tmp1.x0.x0.x0.x0,&tmp1.x0.x0.x0.x0,a_x);//TODO
    Fp16_add(&tmp2,&T->y,&T->y);//2yt
    
    Fp16_div(&lambda,&tmp1,&tmp2);//lambda=3xt^2+a/2yt
    Fp16_sub(&tmp3,&Q->x,&T->x);//tmp3=xq-xt
    Fp16_mul(&tmp3,&tmp3,&lambda);//tmp3=lambda(xq-xt)
    
    Fp16_sub(&ltt,&Q->y,&T->y);//yq-yt
    Fp16_sub(&ltt,&ltt,&tmp3);//ltt=yq-yt-lambda(xq-xt)
    
    Fp16_set(ANS,&ltt);
    
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&lambda);
    Fp16_clear(&ltt);
}
void v2t_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q){
    struct Fp16 v2t;
    Fp16_init(&v2t);
    
    Fp16_sub(&v2t,&Q->x,&T->x);//v2t=xq-xt
    Fp16_set(ANS,&v2t);
    
    Fp16_clear(&v2t);
}
void ltp_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *P,struct EFp16 *Q){
    struct Fp16 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&tmp4);
    Fp16_init(&lambda);
    Fp16_init(&ltp);
    
    if((Fp16_cmp(&T->x,&P->x))==0&&(Fp16_cmp(&T->y,&P->y))!=0){//xt==xp&&yt!=yp
        Fp16_sub(&ltp,&Q->x,&T->x);
        Fp16_set(ANS,&ltp);
        
        return;
    }
    
    Fp16_sub(&tmp1,&T->x,&P->x);//xt-xp
    Fp16_sub(&tmp2,&T->y,&P->y);//yt-yp
    Fp16_div(&lambda,&tmp2,&tmp1);//lambda=(yt-tp)/(xt-xp)
    
    Fp16_sub(&tmp3,&Q->x,&T->x);//tmp3=(xq-xt)
    Fp16_mul(&tmp4,&tmp3,&lambda);//tmp4=lambda(xq-xt)
    
    Fp16_sub(&ltp,&Q->y,&T->y);//ltp=yq-yt
    Fp16_sub(&ltp,&ltp,&tmp4);//ltp=yq-yt-lambda(xq-xt)
    
    Fp16_set(ANS,&ltp);
    
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&tmp4);
    Fp16_clear(&lambda);
    Fp16_clear(&ltp);
}
void vtp_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q){
    struct Fp16 vtp;
    Fp16_init(&vtp);
    if(T->infity==1){//if T is infity
        Fp16_set_ui(ANS,0);
        Fp_set_ui(&ANS->x0.x0.x0.x0,1);
        return;
    }
    
    Fp16_sub(&vtp,&Q->x,&T->x);
    Fp16_set(ANS,&vtp);
    
    Fp16_clear(&vtp);
}

#pragma mark Pseudo Sparse
//P=Q_map, Q =P_map
void Pseudo_type1_ADD_LINE(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *P,struct EFp4 *Q,struct Fp4 *L){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    
    
    Fp4_sub(&A,&P->x,&T->x);//xt-xp
    Fp4_sub(&B,&P->y,&T->y);//yt-yp
    Fp4_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp4_add(&D,&T->x,&P->x);
    Fp4_mul(&tmp1,&C,&C);
    Fp4_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp4_mul(&tmp2,&C,&T->x);
    Fp4_sub(&E,&tmp2,&T->y);
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);
    Fp4_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp_set_ui(&l_tmp.x0.x0.x0.x0,1);
    
    Fp4_mul(&l_tmp.x1.x1,&E,L);
    
    Fp4_neg(&l_tmp.x1.x0,&C);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp4_set(T_ANS,&x3_tmp);
    // if((Fp4_cmp(&T->x,&Q->x))==0&&(Fp4_cmp(&T->y,&Q->y))!=0){//xt==xp&&yt!=yp
    //     Fp4_sub(&ltp,&P->x,&T->x);
    //     Fp4_set(&l_ANS->x0.x0,&ltp);
    // }
    // if(T->infity==TRUE){//if P2==inf
    //     EFp4_set(T_ANS,P);
    //     return;
    // }
    // else if(P->infity==TRUE){//if P1==inf
    //     EFp4_set(T_ANS,T);
    //     return;
    // }
    // else if(Fp4_cmp(&T->x,&P->x)==0&&Fp4_cmp(&T->y,&P->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
    //     EFp4_set_infity(T_ANS);
    //     return;
    // }
    // else if(EFp4_cmp(T,P)==0){ // P=P
    //     EFp4_ECD(T_ANS,T);
    //     return;
    // }
    
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
}
void Pseudo_type1_DBL_LINE(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *Q,struct Fp4 *L){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltt;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltt);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    
//    mpz_t ac_inv;
//    mpz_init(ac_inv);
    
    Fp4_add(&A,&T->y,&T->y);//2yt
    Fp4_mul(&B,&T->x,&T->x);//xt^2
    Fp4_mul_ui(&B,&B,3);//3xt^2
    Fp4_mul_betainv(&tmp);
    Fp4_mul(&tmp, &tmp, &z_inv2);
    Fp4_add(&B, &B, &tmp);
    
    Fp4_div(&C,&B,&A);//
    
    Fp4_add(&D,&T->x,&T->x);//D=2xt
    Fp4_mul(&tmp1,&C,&C);//C^2
    Fp4_sub(&x3_tmp.x,&tmp1,&D);//x3.x = C^2-D
    
    Fp4_mul(&tmp2,&C,&T->x); // C*xt
    Fp4_sub(&E,&tmp2,&T->y); //E=C*xt-yt
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);//C*x3.x
    Fp4_sub(&x3_tmp.y,&E,&tmp3);//x3.y = E-C*x3.x
    
    Fp_set_ui(&l_tmp.x0.x0.x0.x0,1);
    
    Fp4_mul(&l_tmp.x1.x1,&E,L);//l_tmp=
    
    Fp4_neg(&l_tmp.x1.x0,&C);
    Fp16_set(l_ANS,&l_tmp);
    
    EFp4_set(T_ANS,&x3_tmp);
    
    // if(T->infity==TRUE){
    //     EFp4_set(T_ANS,T);
    //     return;
    // }
    // mpz_t cmp;
    // mpz_init(cmp);
    // mpz_set_ui(cmp,0);
    // if(Fp4_cmp_mpz(&T->y,cmp)==0){//P.y==0
    //     EFp4_set_infity(T_ANS);
    //     return;
    // }
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltt);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
}
void Pseudo_type1_mul(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    
    struct Fp16 a,b;
    Fp16_init(&a);
    Fp16_init(&b);
    
    struct Fp4 a0,a1,a2,a3,b2,b3,C_0,C_1,C_2,C_3,T0,T1,T2,T3,T4,a2a3,tmp2,tmp3;
    Fp4_init(&a0);
    Fp4_init(&a1);
    Fp4_init(&a2);
    Fp4_init(&a3);
    Fp4_init(&b2);
    Fp4_init(&b3);
    
    Fp4_init(&C_0);
    Fp4_init(&C_1);
    Fp4_init(&C_2);
    Fp4_init(&C_3);
    
    Fp4_init(&T0);
    Fp4_init(&T1);
    Fp4_init(&T2);
    Fp4_init(&T3);
    Fp4_init(&T4);
    
    Fp4_init(&a2a3);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    
    
    Fp16_set(&a, A);
    Fp16_set(&b, B);
    
    Fp4_set(&a0, &a.x0.x0);
    Fp4_set(&a1, &a.x0.x1);
    Fp4_set(&a2, &a.x1.x0);
    Fp4_set(&a3, &a.x1.x1);
    
    Fp4_set(&b2, &b.x1.x0);
    Fp4_set(&b3, &b.x1.x1);
    
    
    Fp4_add(&a2a3, &a2, &a3);   /**< (a2+a3) */
    Fp4_add(&T4, &b2, &b3);     /**< t4=(b2+b3) */
    
    Fp4_mul(&T1, &a2, &b2);     /**< t1=(a2*b2) */
    Fp4_mul(&T0, &a3, &b3);     /**< t0=(a3*b3) */
    //    Fp4_mul_v(&T0, &tmp2);/**< t0=(a3*b3)*beta */
    
    Fp4_mul(&tmp2, &a2a3, &T4);
    Fp4_sub(&tmp3, &tmp2, &T1);
    Fp4_sub(&C_0, &tmp3, &T0);
    Fp4_mul_v(&tmp3, &C_0);
    Fp4_add(&C_0, &tmp3, &a0);
    
    Fp4_mul_v(&tmp2, &T0);
    Fp4_add(&C_1, &T1, &tmp2);
    Fp4_add(&C_1, &C_1, &a1);
    
    Fp4_mul(&T3, &a0, &b2);
    Fp4_mul(&T2, &a1, &b3);
    Fp4_mul_v(&tmp2, &T2);
    Fp4_add(&C_2, &T3, &tmp2);
    Fp4_add(&C_2, &C_2, &a2);
    
    
    Fp4_add(&tmp2, &a0, &a1);
    Fp4_mul(&tmp3, &tmp2, &T4);
    Fp4_sub(&tmp2, &tmp3, &T3);
    Fp4_sub(&C_3, &tmp2, &T2);
    Fp4_add(&C_3, &C_3, &a3);
    
    Fp4_set(&a.x0.x0,&C_0);
    Fp4_set(&a.x0.x1,&C_1);
    Fp4_set(&a.x1.x0,&C_2);
    Fp4_set(&a.x1.x1,&C_3);
    
    Fp16_set(ANS, &a);
    
    
    Fp16_clear(&a);
    Fp16_clear(&b);
    Fp4_clear(&a0);
    Fp4_clear(&a1);
    Fp4_clear(&a2);
    Fp4_clear(&a3);
    Fp4_clear(&b2);
    Fp4_clear(&b3);
    Fp4_clear(&C_0);
    Fp4_clear(&C_1);
    Fp4_clear(&C_2);
    Fp4_clear(&C_3);
    Fp4_clear(&T0);
    Fp4_clear(&T1);
    Fp4_clear(&T2);
    Fp4_clear(&T3);
    Fp4_clear(&T4);
}

void Pseudo_type1_Optimal_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){//Q:G2,P:G1
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct EFp4 T,P_map,Q_map,EFp4_tmp;
    EFp4_init(&T);
    EFp4_init(&P_map);
    EFp4_init(&Q_map);
    EFp4_init(&EFp4_tmp);
    
    Init_mpz_Cost(&mpz_cost);
    struct Fp4 L,xy,xy_2,y_inv,tmp,y_tmp;
    Fp4_init(&L);
    Fp4_init(&xy);
    Fp4_init(&xy_2);
    Fp4_init(&y_inv);
    Fp4_init(&tmp);
    Fp4_init(&y_tmp);
    Fp4_init(&z_inv2);
    
    Fp4_invert(&y_inv,&P->y);//yp^-1
    Fp4_mul(&xy,&P->x,&y_inv);//xp.yp^-1
    
    Fp4_mul(&xy_2,&xy,&xy);//xy_2 = xp^2.yp^-2
    Fp4_mul(&P_map.x,&xy_2,&P->x);//P.x= xp^3.yp^-2
    Fp4_set(&P_map.y,&P_map.x);
    
    Fp4_mul(&y_tmp,&xy_2,&xy);// xp^2.yp^-2 * xp.yp^-1 = xp^3.yp^-3
    Fp4_mul(&Q_map.y,&y_tmp,&Q->y); //Q_map.y = yQ'.xp^3.yp^-3
    Fp4_mul(&Q_map.x,&xy_2,&Q->x); //Q_map.x = xQ'.xp^2.yp^-2
    
    Fp4_invert(&L,&P_map.y); // L =yp_bar^-1
    Fp4_set(&z_inv2,&xy_2);
    Fp4_mul(&z_inv2, &z_inv2, &z_inv2);
     Print_mpz_Cost(&mpz_cost,"Pre-calc Miller loop");
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    int i;
    
    struct EFp4 Q_neg;
    EFp4_init(&Q_neg);
    Fp4_neg(&Q_neg.y,&Q_map.y);
    Fp4_set(&Q_neg.x,&Q_map.x);
    
    if(X_bit_binary[x_bit]==-1){
        EFp4_set(&T,&Q_neg);
    }else{
        EFp4_set(&T,&Q_map);
    }
    
     Init_mpz_Cost(&mpz_cost);
    for(i=x_bit-1;i>=0;i--){
        switch (X_bit_binary[i]){
            case 0:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                Pseudo_type1_DBL_LINE(&ltt,&T,&T,&P_map,&L);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
                // Fp16_printf(&l_sum);
                break;
            case 1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                Pseudo_type1_DBL_LINE(&ltt,&T,&T,&P_map,&L);
                Pseudo_type1_ADD_LINE(&ltp,&T,&T,&Q_map,&P_map,&L);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
                break;
            case -1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                Pseudo_type1_DBL_LINE(&ltt,&T,&T,&P_map,&L);
                Pseudo_type1_ADD_LINE(&ltp,&T,&T,&Q_neg,&P_map,&L);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
                break;
        }
    }
    
    EFp4_Skew_Frobenius_p(&EFp4_tmp, &Q_map);
    Pseudo_type1_ADD_LINE(&ltp,&T,&T,&EFp4_tmp,&P_map,&L);
    Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
    
    struct Fp16 tmp_f;
    Fp16_init(&tmp_f);
    
    Fp16_frobenius_map_p3(&tmp_f, &l_sum);
    
    Pseudo_type1_DBL_LINE(&ltt,&T,&Q_map,&P_map,&L);
    Pseudo_type1_mul(&l_sum,&tmp_f,&ltt);
    Fp16_set(ANS,&l_sum);
    Print_mpz_Cost(&mpz_cost," Miller loop:");
    EFp4_clear(&Q_neg);
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    EFp4_clear(&P_map);
    EFp4_clear(&Q_map);
    EFp4_clear(&EFp4_tmp);
    Fp4_clear(&L);
    Fp4_clear(&xy);
    Fp4_clear(&xy_2);
    Fp4_clear(&y_inv);
    Fp4_clear(&tmp);
    Fp4_clear(&y_tmp);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    Fp16_clear(&tmp_f);
}

void Pseudo_Sparse_Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2){
    struct EFp4 G1_EFp4,G2_EFp4;
    EFp4_init(&G1_EFp4);
    EFp4_init(&G2_EFp4);
    
    struct Fp16 ltp,Miller_X,t_ans;
    Fp16_init(&ltp);
    Fp16_init(&Miller_X);
    Fp16_init(&t_ans);
    
    struct timeval t0;
    struct timeval t1;
    
    EFp4_set_EFp(&G1_EFp4,G1);
    EFp16_to_EFp4_map(&G2_EFp4,G2);
    printf("\n\nPseudo Sparse Miller algo\n");
    gettimeofday(&t0, 0);
    Big_M = 0; Small_m =0; Big_add=0; Sqr=0; fp16_pow =0;
    Pseudo_type1_Optimal_Miller(&Miller_X,&G1_EFp4,&G2_EFp4,X);
    gettimeofday(&t1, 0);
    double elapsed = timedifference_msec(t0,t1);
    printf("Miller Time : %f [ms]\n", elapsed);
    printf("Fp M =%d, Fp SQR = %d, Basis m = %d, ADD = %d, SUB= %d\n",Big_M, Sqr, Small_m,Big_add,fp16_pow);
    Big_M = 0; Small_m =0; Big_add=0; Sqr=0; fp16_pow =0;
    Final_Exp(&t_ans,&Miller_X);
    printf("FINAL EXP: Fp M =%d, Fp SQR = %d, Basis m = %d, ADD = %d, SUB= %d\n",Big_M, Sqr, Small_m,Big_add,fp16_pow);
    Fp16_set(ANS, &t_ans);
    
    Fp16_clear(&ltp);
    Fp16_clear(&Miller_X);
    Fp16_clear(&t_ans);
}

void Pseudo_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){//Q:G2,P:G1
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct EFp4 T,P_map,Q_map,EFp4_tmp;
    EFp4_init(&T);
    EFp4_init(&P_map);
    EFp4_init(&Q_map);
    EFp4_init(&EFp4_tmp);
    
    struct Fp4 L,xy,xy_2,y_inv,tmp,y_tmp;
    Fp4_init(&L);
    Fp4_init(&xy);
    Fp4_init(&xy_2);
    Fp4_init(&y_inv);
    Fp4_init(&tmp);
    Fp4_init(&y_tmp);
    Fp4_init(&z_inv2);
    
    Fp4_invert(&y_inv,&P->y);//yp^-1
    Fp4_mul(&xy,&P->x,&y_inv);//xp.yp^-1
    
    Fp4_mul(&xy_2,&xy,&xy);//xy_2 = xp^2.yp^-2
    Fp4_mul(&P_map.x,&xy_2,&P->x);//P.x= xp^3.yp^-2
    Fp4_set(&P_map.y,&P_map.x);
    
    Fp4_mul(&y_tmp,&xy_2,&xy);// xp^2.yp^-2 * xp.yp^-1 = xp^3.yp^-3
    Fp4_mul(&Q_map.y,&y_tmp,&Q->y); //Q_map.y = yQ'.xp^3.yp^-3
    Fp4_mul(&Q_map.x,&xy_2,&Q->x); //Q_map.x = xQ'.xp^2.yp^-2
    
    EFp4_set(&T,&Q_map);
    Fp4_invert(&L,&P_map.y); // L =yp_bar^-1
    Fp4_set(&z_inv2,&xy_2);
    Fp4_mul(&z_inv2, &z_inv2, &z_inv2);
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    
    //    rational_point_check(&Q_map);
    //    Pseudo_type1_DBL_LINE(&ltt,&T,&T,&P_map,&L);
    //    rational_point_check(&T);
    //    Pseudo_type1_ADD_LINE(&ltp,&T,&T,&Q_map,&P_map,&L);
    //    Fp16_printf(&ltp);
    //    rational_point_check(&T);
    
    int i;
    int r_bit;
    r_bit= (int)mpz_sizeinbase(loop,2);
    
    for(i=r_bit-2;i>=0;i--){
        if(mpz_tstbit(loop,i)==1){
            //            printf("\n%d",i);
            Fp16_mul(&l_sum,&l_sum,&l_sum);
            Pseudo_type1_DBL_LINE(&ltt,&T,&T,&P_map,&L);
            Pseudo_type1_ADD_LINE(&ltp,&T,&T,&Q_map,&P_map,&L);
            //            rational_point_check(&T);
            Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
            Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
            //            Fp16_mul(&l_sum,&l_sum,&ltt);
            //            Fp16_mul(&l_sum,&l_sum,&ltp);
        }else{
            Fp16_mul(&l_sum,&l_sum,&l_sum);
            Pseudo_type1_DBL_LINE(&ltt,&T,&T,&P_map,&L);
            //            rational_point_check(&T);
            Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
            //            Fp16_mul(&l_sum,&l_sum,&ltt);
        }
    }
    // EFp4_printf(&T);
    Fp16_set(ANS,&l_sum);
    
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    EFp4_clear(&P_map);
    EFp4_clear(&Q_map);
    EFp4_clear(&EFp4_tmp);
    Fp4_clear(&L);
    Fp4_clear(&xy);
    Fp4_clear(&xy_2);
    Fp4_clear(&y_inv);
    Fp4_clear(&tmp);
    Fp4_clear(&y_tmp);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
}

void Pseudo_Sparse_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2){
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    struct EFp4 EFp4_G1,EFp4_G2;
    EFp4_init(&EFp4_G1);
    EFp4_init(&EFp4_G2);
    mpz_t tm1;
    mpz_init(tm1);
    mpz_sub_ui(tm1,trace_t,1);
    
    EFp4_set_EFp(&EFp4_G1,G1);//get P of FP in Fp4
    EFp16_to_EFp4_map(&EFp4_G2,G2); //non isomorphic map of Q to Q'
    Pseudo_Miller(&t_ans,&EFp4_G1,&EFp4_G2,tm1);
    
    Final_Exp(&t_ans,&t_ans);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&t_ans);
    mpz_clear(tm1);
}

void Sparse_type1_ADD_LINE(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *P,struct EFp4 *Q,struct Fp4 *Qx_neg){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E,F;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    Fp4_init(&F);
    
    Fp4_sub(&A,&P->x,&T->x);//xt-xp
    Fp4_sub(&B,&P->y,&T->y);//yt-yp
    Fp4_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp4_add(&D,&T->x,&P->x);
    Fp4_mul(&tmp1,&C,&C);
    Fp4_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp4_mul(&tmp2,&C,&T->x);
    Fp4_sub(&E,&tmp2,&T->y);
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);
    Fp4_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp4_set(&l_tmp.x0.x0,&Q->y);
    
    Fp4_set(&l_tmp.x1.x1,&E);
    
    Fp4_mul(&F,&C,Qx_neg);
    Fp4_set(&l_tmp.x1.x0,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp4_set(T_ANS,&x3_tmp);
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
    Fp4_clear(&F);
}

void Sparse_type1_DBL_LINE(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *Q,struct Fp4 *Qx_neg){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E,F;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    Fp4_init(&F);
    
    
    Fp4_add(&A,&T->y,&T->y);//xt-xp
    Fp4_mul(&B,&T->x,&T->x);
    Fp4_mul_ui(&B,&B,3);
    struct Fp4 ac_inv;
    Fp4_init(&ac_inv);
    Fp4_mul_betainv(&ac_inv);
    Fp4_add(&B,&B,&ac_inv);
    Fp4_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp4_add(&D,&T->x,&T->x);
    Fp4_mul(&tmp1,&C,&C);
    Fp4_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp4_mul(&tmp2,&C,&T->x);
    Fp4_sub(&E,&tmp2,&T->y);
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);
    Fp4_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp4_set(&l_tmp.x0.x0,&Q->y);
    
    Fp4_set(&l_tmp.x1.x1,&E);
    
    Fp4_mul(&F,&C,Qx_neg);
    Fp4_set(&l_tmp.x1.x0,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp4_set(T_ANS,&x3_tmp);
    
    if(T->infity==TRUE){
        EFp4_set(T_ANS,T);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp4_cmp_mpz(&T->y,cmp)==0){//P.y==0
        EFp4_set_infity(T_ANS);
        return;
    }
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
    Fp4_clear(&F);
    mpz_clear(cmp);
}

void Sparse_type1_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct EFp4 T,EFp_tmp;
    EFp4_init(&T);
    EFp4_init(&EFp_tmp);
    
    struct Fp4 Px_neg;
    Fp4_init(&Px_neg);
    Fp4_neg(&Px_neg,&P->x);
    
    EFp4_set(&T,Q);
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    int i;
    int r_bit;
    r_bit= (int)mpz_sizeinbase(loop,2);
    
    for(i=r_bit-2;i>=0;i--){
        if(mpz_tstbit(loop,i)==1){
            Fp16_mul(&l_sum,&l_sum,&l_sum);
            Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
            Sparse_type1_ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
            Fp16_mul(&l_sum,&l_sum,&ltt);
            Fp16_mul(&l_sum,&l_sum,&ltp);
        }else{
            Fp16_mul(&l_sum,&l_sum,&l_sum);
            Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
            Fp16_mul(&l_sum,&l_sum,&ltt);
        }
    }
    Fp16_set(ANS,&l_sum);
    
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    EFp4_clear(&EFp_tmp);
    Fp4_clear(&Px_neg);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
}

void Sparse_type1_Optimal_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct EFp4 T,EFp4_tmp;
    EFp4_init(&T);
    EFp4_init(&EFp4_tmp);
    
    mpz_t p3;
    mpz_init(p3);
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    int i;
    
    struct Fp4 Px_neg;
    Fp4_init(&Px_neg);
    
    Fp4_neg(&Px_neg,&P->x);
    
    struct EFp4 Q_neg;
    EFp4_init(&Q_neg);
    Fp4_neg(&Q_neg.y,&Q->y);
    Fp4_set(&Q_neg.x,&Q->x);
    
    if(X_bit_binary[x_bit]==-1){
        EFp4_set(&T,&Q_neg);
    }else{
        EFp4_set(&T,Q);
    }
    for(i=x_bit-1;i>=0;i--){
        switch (X_bit_binary[i]){
            case 0:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                Fp16_mul(&l_sum,&l_sum,&ltt);
                break;
                
            case 1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                
                Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                Sparse_type1_ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
            case -1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                
                Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                Sparse_type1_ADD_LINE(&ltp,&T,&T,&Q_neg,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
        }
    }
    
    EFp4_Skew_Frobenius_p(&EFp4_tmp,Q);
    //    EFp4_SCM_BIN(&EFp4_tmp, &Q_map, prime);
    
    Sparse_type1_ADD_LINE(&ltp,&T,&T,&EFp4_tmp,P,&Px_neg);
    Fp16_mul(&l_sum,&l_sum,&ltp);
    
    
    struct Fp16 tmp_f;
    Fp16_init(&tmp_f);
    
    
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    //    Fp16_pow(&l_sum, &tmp_f,prime);
    Fp16_frobenius_map(&l_sum, &tmp_f);
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    
    
    //ltt_q(&ltt,Q,P);
    Sparse_type1_DBL_LINE(&ltt,&T,Q,P,&Px_neg);
    
    Fp16_mul(&l_sum,&tmp_f,&ltt);
    Fp16_set(ANS,&l_sum);
    
    
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    EFp4_clear(&EFp4_tmp);
    
}
void EFp4_SCM_BIN_Sparse(struct EFp4 *ANS,struct EFp4 *P,mpz_t j){
    int i,length;
    length = (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp4 Q,R;
    EFp4_init(&Q);
    EFp4_set(&Q,P);
    EFp4_init(&R);
    for(i=1;i < length; i++){
        EFp4_ECD_Sparse(&Q,&Q);
        scm_cost.loop_ecd++;
        if(j_binary[i]=='1'){
            EFp4_ECA(&Q,&Q,P);
            scm_cost.loop_eca++;
        }
    }
//    Print_scm_Cost(&scm_cost, "Normal SCM");
    EFp4_set(ANS,&Q);
    
    EFp4_clear(&Q);
    EFp4_clear(&R);
    return;
}
void EFp4_SCM_BIN_Pseudo_Sparse(struct EFp4 *ANS,struct EFp4 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp4 Q,R;
    EFp4_init(&Q);
    EFp4_set(&Q,P);
    EFp4_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp4_ECD_Pseudo_Sparse(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp4_ECA(&Q,&Q,P);
        }
    }
    EFp4_set(ANS,&Q);
    
    EFp4_clear(&Q);
    EFp4_clear(&R);
    return;
}



void Sparse_Ate_Pairing(struct Fp16 *ANS,struct EFp4 *G1,struct EFp4 *G2){
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    
    mpz_t tm1;
    mpz_init(tm1);
    mpz_sub_ui(tm1,trace_t,1);
    
    Sparse_type1_Miller(&t_ans,G1,G2,tm1);
    Final_Exp(&t_ans,&t_ans);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&t_ans);
    mpz_clear(tm1);
}


void Sparse_Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2){
    struct EFp4 G1_EFp4,G2_EFp4;
    EFp4_init(&G1_EFp4);
    EFp4_init(&G2_EFp4);
    
    struct Fp16 ltp,Miller_X,t_ans;
    Fp16_init(&ltp);
    Fp16_init(&Miller_X);
    Fp16_init(&t_ans);
    
    EFp4_set_EFp(&G1_EFp4,G1);
    EFp16_to_EFp4_map(&G2_EFp4,G2);
    
    Sparse_type1_Optimal_Miller(&Miller_X,&G1_EFp4,&G2_EFp4,X);
    Final_Exp(&t_ans,&Miller_X);
    Fp16_set(ANS, &t_ans);
    
    Fp16_clear(&ltp);
    Fp16_clear(&Miller_X);
    Fp16_clear(&t_ans);
}

void check_Pairing(void){
    struct EFp P_EFp, R_EFp;
    EFp_init(&P_EFp);
    EFp_init(&R_EFp);
    
    
    struct EFp16 P_Fp16,Q_Fp16,Q_G2,R_Fp16,S_Fp16;
    EFp16_init(&P_Fp16);
    EFp16_init(&Q_Fp16);
    EFp16_init(&R_Fp16);
    EFp16_init(&S_Fp16);
    EFp16_init(&Q_G2);
    
    struct Fp16 ans_Fp16,tmp1_Fp16,tmp2_Fp16,tmp3_Fp16;
    Fp16_init(&ans_Fp16);
    Fp16_init(&tmp1_Fp16);
    Fp16_init(&tmp2_Fp16);
    Fp16_init(&tmp3_Fp16);
    
    struct EFp4 P_EFp4, Q_EFp4, R_EFp4, S_EFp4;
    EFp4_init(&P_EFp4);
    EFp4_init(&Q_EFp4);
    EFp4_init(&R_EFp4);
    EFp4_init(&S_EFp4);
    
    mpz_t a,b,ab;
    mpz_init(a);
    mpz_init(b);
    mpz_init(ab);
    
    mpz_set_ui(a,31);
    mpz_set_ui(b,11);
    mpz_mul(ab,a,b);
    
    //    mpz_t x,y;
    //    mpz_inits(x,y,(mpz_ptr) NULL);
    //    mpz_set_str(x,"585634432000126707887057629201458798521445673169852167601690963969803976546284829414253996667593644070", 10);
    //    mpz_set_str(y,"354142610898165644039571248952651466071208533406755516950532183093299777682028181848353083070437124285", 10);
    //    Fp_set_mpz(&P_EFp.x,x);
    //    Fp_set_mpz(&P_EFp.y,y);
    //    EFp16_set_EFp(&P_Fp16,&P_EFp);
    //    //    EFp16_printf(&P_Fp16);
    //    mpz_clears(x,y,(mpz_ptr) NULL);
    //
    //    mpz_t x1, x2, x3, x4, y1, y2, y3, y4;
    //    mpz_inits(x1,x2,x3,x4,y1,y2,y3,y4,(mpz_ptr) NULL);
    //
    //    mpz_set_str(x1,"244491741495785299367612042330502507448971522000836606635913541960426621971155507113603357863094662406", 10);
    //    mpz_set_str(x2,"218139782565182912336439776658495825407576298090161716969444130503315739115994708740875588793005720529", 10);
    //    mpz_set_str(x3,"489284606975189850057489668560782153703917092806842888347849447445917065112045634272024913067442383327", 10);
    //    mpz_set_str(x4,"246161136867400494929650830056723287620786052741538082500453356439354042813046547782023695596667138627", 10);
    //    mpz_set_str(y1,"48629828925217256502074899893821067838880538260652902679139733711894991914937403050096788961116072841", 10);
    //    mpz_set_str(y2,"338343177397929794197126001730091007944155731601968256694856657191699538746843500939799293603025576582", 10);
    //    mpz_set_str(y3,"38952231663050097715674495759233974069804273632029475897530255418880509349197207272848037344062136227", 10);
    //    mpz_set_str(y4,"250604271284793810746824678011335939079727656816914850518650242122895946469531880085329357599729742073", 10);
    //    Fp_set_mpz(&Q_Fp16.x.x0.x1.x0.x0,x1);
    //    Fp_set_mpz(&Q_Fp16.x.x0.x1.x0.x1,x2);
    //    Fp_set_mpz(&Q_Fp16.x.x0.x1.x1.x0,x3);
    //    Fp_set_mpz(&Q_Fp16.x.x0.x1.x1.x1,x4);
    //
    //    Fp_set_mpz(&Q_Fp16.y.x1.x1.x0.x0,y1);
    //    Fp_set_mpz(&Q_Fp16.y.x1.x1.x0.x1,y2);
    //    Fp_set_mpz(&Q_Fp16.y.x1.x1.x1.x0,y3);
    //    Fp_set_mpz(&Q_Fp16.y.x1.x1.x1.x1,y4);
    
    
    //    printf("G1=");
    //    EFp16_printf(&P_Fp16);
    //    printf("G2=");s
    //    EFp16_printf(&Q_Fp16);
    //---------------------------------------------------
    //    EFp_random_set(&P_EFp);
    //    EFp16_set_EFp(&P_Fp16,&P_EFp);
    //    EFp16_random_set(&Q_Fp16);
    //    printf("\nTate Pairing\n");
    //
    //    printf("G1=");
    //    EFp16_printf(&P_Fp16);
    //    printf("G2=");
    //    EFp16_printf(&Q_Fp16);
    //
    //    Tate_Pairing(&tmp1_Fp16,&P_Fp16,&Q_Fp16);
    //
    //    Fp16_pow(&tmp1_Fp16,&tmp1_Fp16,ab);
    //    printf("\nf^ab=");
    //    Fp16_printf(&tmp1_Fp16);
    //
    //    EFp16_SCM_BIN(&R_Fp16,&P_Fp16,a);
    //    EFp16_SCM_BIN(&S_Fp16,&Q_Fp16,b);
    //
    //    Tate_Pairing(&tmp2_Fp16,&R_Fp16,&S_Fp16);
    //
    //    printf("\nf'  =");
    //    Fp16_printf(&tmp2_Fp16);
    //    //----------------------------------------------------
    //    printf("\nAte Pairing\n");
    //
    //    printf("\n\n Ate Pairing\n");
    //    EFp_random_set(&P_EFp);
    //    printf("G1=");
    //    EFp_printf(&P_EFp);
    //    EFp16_set_EFp(&P_Fp16, &P_EFp);
    //    printf("\nG2=");
    //    EFp16_random_set_G2(&Q_Fp16);
    //    EFp16_printf(&Q_Fp16);
    //
    //    Ate_Pairing(&tmp1_Fp16,&P_Fp16,&Q_Fp16);
    //
    //    Fp16_pow(&tmp1_Fp16,&tmp1_Fp16,ab);
    //    printf("\nf^ab=");
    //    Fp16_printf(&tmp1_Fp16);
    //
    //    EFp16_SCM_BIN(&R_Fp16,&P_Fp16,a);
    //    EFp16_SCM_BIN(&S_Fp16,&Q_Fp16,b);
    //
    //    Ate_Pairing(&tmp2_Fp16,&R_Fp16,&S_Fp16);
    //
    //    printf("\nf'  =");
    //    Fp16_printf(&tmp2_Fp16);
    //----------------------------------------------------
    //    printf("\nAte Pairing Sparse\n");
    //
    //    EFp4_set_EFp(&P_EFp4,&P_EFp);
    //    EFp16_to_EFp4_map(&Q_EFp4,&Q_Fp16);
    //
    //     printf("G1=");
    //     EFp4_printf(&P_EFp4);
    //     printf("G2=");
    //     EFp4_printf(&Q_EFp4);
    //
    //    Sparse_Ate_Pairing(&tmp1_Fp16,&P_EFp4,&Q_EFp4);
    //    //
    //    Fp16_pow(&tmp1_Fp16,&tmp1_Fp16,ab);
    //    printf("\nf^ab=");
    //    Fp16_printf(&tmp1_Fp16);
    //
    //    EFp4_SCM_BIN(&R_EFp4,&P_EFp4,a);
    //    EFp4_SCM_BIN_Sparse(&S_EFp4,&Q_EFp4,b);
    //
    //    Sparse_Ate_Pairing(&tmp2_Fp16,&R_EFp4,&S_EFp4);
    //
    //    printf("f'  =");
    //    Fp16_printf(&tmp2_Fp16);
    
    // ----------------------------------------------------
    //    printf("\nPseudo Sparse Ate Pairing isomorphic Twist\n");
    //    printf("G1=");
    //    EFp_printf(&P_EFp);
    //    printf("G2=");
    //    EFp16_printf(&Q_Fp16);
    //
    //    Pseudo_Sparse_Ate_Pairing(&tmp1_Fp16,&P_EFp,&Q_Fp16);
    //    Fp16_pow(&tmp1_Fp16,&tmp1_Fp16,ab);
    //    printf("\nf^ab=");
    //    Fp16_printf(&tmp1_Fp16);
    //
    //    EFp_SCM_BIN(&R_EFp,&P_EFp,a);
    //    EFp16_SCM_BIN(&S_Fp16,&Q_Fp16,b);
    //
    //    Pseudo_Sparse_Ate_Pairing(&tmp2_Fp16,&R_EFp,&S_Fp16);
    //
    //    printf("\nf'  =");
    //    Fp16_printf(&tmp2_Fp16);
    // ----------------------------------------------------
    
    //    printf("\n\nOptimal Ate Pairing\n");
    //    EFp_random_set(&P_EFp);
    //    printf("G1=");
    //    EFp_printf(&P_EFp);
    //    EFp16_set_EFp(&P_Fp16, &P_EFp);
    //    printf("\nG2=");
    //    EFp16_random_set_G2(&Q_Fp16);
    //    EFp16_printf(&Q_Fp16);
    //
    //
    //    Optimal_Ate_Pairing(&tmp1_Fp16,&P_Fp16,&Q_Fp16);
    //    Fp16_pow(&tmp1_Fp16,&tmp1_Fp16,ab);
    //    printf("\nf^ab=");
    //    Fp16_printf(&tmp1_Fp16);
    //
    //    EFp16_SCM_BIN(&R_Fp16,&P_Fp16,a);
    //    EFp16_SCM_BIN(&S_Fp16,&Q_Fp16,b);
    //
    //    Optimal_Ate_Pairing(&tmp2_Fp16,&R_Fp16,&S_Fp16);
    //
    //    printf("\nf'  =");
    //    Fp16_printf(&tmp2_Fp16);
    //----------------------------------------------------
//    printf("\n\n Optimal Ate Pairing Sparse\n");
//    EFp_random_set(&P_EFp);
//    printf("G1=");
//    EFp_printf(&P_EFp);
//    printf("\nG2=");
//    EFp16_random_set_G2(&Q_Fp16);
//    EFp16_printf(&Q_Fp16);
//
//    Sparse_Optimal_Ate_Pairing(&tmp1_Fp16,&P_EFp,&Q_Fp16);
//
//    Fp16_pow(&tmp1_Fp16,&tmp1_Fp16,ab);
//    printf("\nf^ab=");
//    Fp16_printf(&tmp1_Fp16);
//
//    EFp_SCM_BIN(&R_EFp,&P_EFp,a);
//    EFp16_SCM_BIN(&S_Fp16,&Q_Fp16,b);
//
//    Sparse_Optimal_Ate_Pairing(&tmp2_Fp16,&R_EFp,&S_Fp16);
//
//    printf("\nf'  =");
//    Fp16_printf(&tmp2_Fp16);
    //----------------------------------------------------
    EFp_random_set(&P_EFp);
    printf("G1=");
    EFp_printf(&P_EFp);
    printf("\nG2=");
    EFp16_random_set_G2(&Q_Fp16);
    EFp16_printf(&Q_Fp16);


    printf("\n\nPseudo Sparse Optimal Ate Pairing \n");
    Big_M = 0; Small_m =0; Big_add=0; Sqr=0; fp16_pow =0;
    Pseudo_Sparse_Optimal_Ate_Pairing(&tmp1_Fp16,&P_EFp,&Q_Fp16);
    printf("Fp SQR in Fp16 mul = Fp16 pow %d, %d, %d, %d, %d\n",Big_M, Sqr, Small_m,Big_add,fp16_pow);
    Fp16_pow(&tmp1_Fp16,&tmp1_Fp16,ab);
    printf("\nf^ab=");
    Fp16_printf(&tmp1_Fp16);

    EFp_SCM_BIN(&R_EFp,&P_EFp,a);
    EFp16_SCM_BIN(&S_Fp16,&Q_Fp16,b);
    Pseudo_Sparse_Optimal_Ate_Pairing(&tmp2_Fp16,&R_EFp,&S_Fp16);
    printf("\nf'  =");
    Fp16_printf(&tmp2_Fp16);
    
    
    printf("Bilinearity Check: \n");
    if (Fp16_cmp(&tmp2_Fp16, &tmp1_Fp16) == 0) {
        printf("Success \n");
    }
    else{
        printf("Failure\n");
    }
    
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(ab);
    
    EFp_clear(&P_EFp);
    EFp_clear(&R_EFp);
    
    EFp16_clear(&P_Fp16);
    EFp16_clear(&Q_Fp16);
    EFp16_clear(&R_Fp16);
    EFp16_clear(&S_Fp16);
    
    Fp16_clear(&ans_Fp16);
    Fp16_clear(&tmp1_Fp16);
    Fp16_clear(&tmp2_Fp16);
    Fp16_clear(&tmp3_Fp16);
    
    EFp4_clear(&P_EFp4);
    EFp4_clear(&Q_EFp4);
    EFp4_clear(&R_EFp4);
    EFp4_clear(&S_EFp4);
}
void rational_point_check(struct EFp4 *A){
    struct EFp4 SS,QQ,TT;
    EFp4_init(&SS);
    EFp4_init(&QQ);
    EFp4_init(&TT);
    EFp4_set(&QQ,A);
    Fp4_mul(&SS.y,&QQ.y,&QQ.y);
    Fp4_mul(&SS.x,&QQ.x,&QQ.x);
    Fp4_mul(&SS.x,&SS.x,&QQ.x);
    
    Fp4_mul_betainv(&TT.y);
    Fp4_mul(&TT.x,&z_inv2,&QQ.x);
    // Fp4_mul(&TT.x,&TT.y,&TT.x);
    Fp4_mul(&TT.x,&TT.x,&TT.y);
    Fp4_add(&SS.x,&SS.x,&TT.x);
    if(Fp4_cmp(&SS.x,&SS.y)!=0){
        printf("\nNot Rational point\n");
        EFp4_printf(&SS);
        EFp4_printf(A);
    }
}

//-----------------------------------------------------------------------------------------
#pragma mark G1 SCM methods
void JSF(int **binary,mpz_t S[2],int *loop_length){
    int i,j;
    unsigned long int u;
    mpz_t mod_2,mod_4,mod_8;
    mpz_init(mod_2);
    mpz_init(mod_4);
    mpz_init(mod_8);
    
    mpz_t k[2];
    mpz_init(k[0]);
    mpz_init(k[1]);
    //set
    j=0;
    mpz_set(k[0],S[0]);
    mpz_set(k[1],S[1]);
    
    while(mpz_cmp_ui(k[0],0)>0 || mpz_cmp_ui(k[1],0)>0){
        for(i=0; i<2; i++){
            mpz_mod_ui(mod_2,k[i],2);
            if(mpz_cmp_ui(mod_2,0)==0){
                u=0;
            }else{
                mpz_mod_ui(mod_4,k[i],4);
                u=mpz_get_ui(mod_4);
                if(u==3){
                    u=-1;
                }
                mpz_mod_ui(mod_8,k[i],8);
                mpz_mod_ui(mod_4,k[1-i],4);
                if((mpz_cmp_ui(mod_8,3)==0 || mpz_cmp_ui(mod_8,5)==0) && mpz_cmp_ui(mod_4,2)==0){
                    u=-u;
                }
            }
            binary[i][j]=(int)u;
        }
        for(i=0; i<2; i++){
            u=binary[i][j];
            switch (u){
                case 1:
                    mpz_sub_ui(k[i],k[i],1);
                    break;
                case -1:
                    mpz_add_ui(k[i],k[i],1);
                    break;
                default:
                    break;
            }
            mpz_tdiv_q_ui(k[i],k[i],2);
        }
        j=j+1;
    }
    *loop_length=j-1;
    
    mpz_clear(mod_2);
    mpz_clear(mod_4);
    mpz_clear(mod_8);
    mpz_clear(k[0]);
    mpz_clear(k[1]);
}
void G1_SCM_normal(struct EFp16 *ANS,struct EFp16 *P, mpz_t scalar){
    struct EFp tmp_P;
    EFp_init(&tmp_P);
    
    //set tmp_P
    Fp_set(&tmp_P.x,&P->x.x0.x0.x0.x0);
    Fp_set(&tmp_P.y,&P->y.x0.x0.x0.x0);
    tmp_P.infity=P->infity;
    
    //SCM
    EFp_SCM_BIN(&tmp_P,&tmp_P,scalar);
    
    //set ANS
    Fp16_set_ui(&ANS->x,0);
    Fp_set(&ANS->x.x0.x0.x0.x0,&tmp_P.x);
    Fp16_set_ui(&ANS->y,0);
    Fp_set(&ANS->y.x0.x0.x0.x0,&tmp_P.y);
    ANS->infity=tmp_P.infity;
    
    EFp_clear(&tmp_P);
}
void G1_SCM_2split(struct EFp16 *ANS,struct EFp16 *P,mpz_t scalar){
    //s=s0+s1[x^4
    int i,length_s[2],loop_length;struct EFp next_tmp_P,tmp_P,tmp_P_x4,skew_P,skew_7P;
    EFp_init(&next_tmp_P);
    EFp_init(&tmp_P);
    EFp_init(&skew_P);
    EFp_init(&skew_7P);
    EFp_init(&tmp_P_x4);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    struct EFp table[4];
    for(i=0; i<4; i++){
        EFp_init(&table[i]);
    }
    
    //set tmp_P
    EFp_set_EFp16(&tmp_P,P);
    //set tmp_P_x4
    EFp_Skew_Frobenius_p4(&skew_P,&tmp_P);
    EFp_ECD(&skew_7P,&skew_P);              //[2]skew_P
    EFp_ECA(&skew_7P,&skew_7P,&skew_P);     //[3]skew_P
    EFp_ECD(&skew_7P,&skew_7P);             //[6]skew_P
    EFp_ECA(&skew_7P,&skew_7P,&skew_P);     //[7]skew_P
    EFp_init(&skew_P);
    EFp_ECD(&skew_P,&tmp_P);                //[2]tmp_P
    EFp_ECA(&skew_P,&skew_P,&tmp_P);        //[3]tmp_P
    EFp_ECD(&skew_P,&skew_P);               //[6]tmp_P
    EFp_ECD(&skew_P,&skew_P);               //[12]tmp_P
    EFp_ECD(&skew_P,&skew_P);               //[24]tmp_P
    EFp_ECA(&tmp_P_x4,&skew_P,&skew_7P);
    EFp_neg(&tmp_P_x4,&tmp_P_x4);
    
    //set table
    table[0].infity=1;                        //00
    EFp_set(&table[1],&tmp_P);                //01
    EFp_set(&table[2],&tmp_P_x4);            //10
    EFp_ECA(&table[3],&tmp_P,&tmp_P_x4);    //11
    
    //s0,s1
    mpz_pow_ui(buf,X,4);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);
    
    //binary
    loop_length=0;
    for(i=0; i<2; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        printf("G2 length_s%d:%d\n",i,length_s[i]);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    printf("\n");
    //set binary
    char binary_s[2][loop_length+1];
    char str[5],*e;
    int binary[loop_length+1];
    for(i=0; i<2; i++){
        if(length_s[i]==loop_length){
            mpz_get_str(binary_s[i],2,s[i]);
        }else{
            char binary_buf[loop_length+1];
            mpz_get_str(binary_buf,2,s[i]);
            memset(binary_s[i],'0',sizeof(binary_s[i]));
            memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
        }
    }
    for(i=0; i<loop_length; i++){
        sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
        binary[i]=(int)strtol(str,&e,2);
    }
    EFp_set(&next_tmp_P,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        EFp_ECD(&next_tmp_P,&next_tmp_P);
        EFp_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
    }
    
    EFp16_set_EFp(ANS,&next_tmp_P);
    
    mpz_clear(buf);
    EFp_clear(&next_tmp_P);
    EFp_clear(&tmp_P);
    EFp_clear(&skew_P);
    EFp_clear(&skew_7P);
    EFp_clear(&tmp_P_x4);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    for(i=0; i<4; i++){
        EFp_clear(&table[i]);
    }
}
void G1_SCM_2split_JSF(struct EFp16 *ANS,struct EFp16 *P, mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;struct EFp next_tmp_P,tmp_P,tmp_P_neg,tmp_P_x4,tmp_P_x4_neg,skew_P,skew_7P;
    EFp_init(&next_tmp_P);
    EFp_init(&tmp_P);
    EFp_init(&tmp_P_neg);
    EFp_init(&skew_P);
    EFp_init(&skew_7P);
    EFp_init(&tmp_P_x4);
    EFp_init(&tmp_P_x4_neg);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    struct EFp table[9];
    for(i=0; i<9; i++){
        EFp_init(&table[i]);
    }
    
    //set tmp_P
    EFp_set_EFp16(&tmp_P,P);
    EFp_neg(&tmp_P_neg,&tmp_P);
    //set tmp_P_x4
    EFp_Skew_Frobenius_p4(&skew_P,&tmp_P);
    EFp_ECD(&skew_7P,&skew_P);              //[2]skew_P
    EFp_ECA(&skew_7P,&skew_7P,&skew_P);     //[3]skew_P
    EFp_ECD(&skew_7P,&skew_7P);             //[6]skew_P
    EFp_ECA(&skew_7P,&skew_7P,&skew_P);     //[7]skew_P
    EFp_init(&skew_P);
    EFp_ECD(&skew_P,&tmp_P);                //[2]tmp_P
    EFp_ECA(&skew_P,&skew_P,&tmp_P);        //[3]tmp_P
    EFp_ECD(&skew_P,&skew_P);               //[6]tmp_P
    EFp_ECD(&skew_P,&skew_P);               //[12]tmp_P
    EFp_ECD(&skew_P,&skew_P);               //[24]tmp_P
    EFp_ECA(&tmp_P_x4,&skew_P,&skew_7P);
    EFp_neg(&tmp_P_x4,&tmp_P_x4);
    EFp_neg(&tmp_P_x4_neg,&tmp_P_x4);
    
    //set table
    table[0].infity=1;                                 //00
    EFp_set(&table[1],&tmp_P);                        //01
    EFp_set(&table[2],&tmp_P_x4);                    //10
    EFp_ECA(&table[3],&tmp_P_x4,&tmp_P);            //11
    EFp_set(&table[4],&tmp_P_neg);                    //0-1
    EFp_set(&table[5],&tmp_P_x4_neg);                //-10
    EFp_ECA(&table[6],&tmp_P_x4_neg,&tmp_P_neg);    //-1-1
    EFp_ECA(&table[7],&tmp_P_x4,&tmp_P_neg);        //1-1
    EFp_ECA(&table[8],&tmp_P_x4_neg,&tmp_P);        //-11
    
    //s0,s1
    mpz_pow_ui(buf,X,4);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);
    
    //get loop_length
    loop_length=0;
    for(i=0; i<2; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    //JSF
    int JSF_length;
    int JSF_binary[2][loop_length+1];
    //    char check[5];
    for(i=0; i<loop_length; i++){
        JSF_binary[0][i]=0;
        JSF_binary[1][i]=0;
    }
    int *JSF_pointer[2];
    JSF_pointer[0]=JSF_binary[0];
    JSF_pointer[1]=JSF_binary[1];
    JSF(JSF_pointer,s,&JSF_length);
    int binary[JSF_length+1];
    
    for(i=JSF_length; i>=0; i--){
        if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0)         binary[i]=0;
        else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1)     binary[i]=1;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0)     binary[i]=2;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)    binary[i]=3;
        else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)    binary[i]=4;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)    binary[i]=5;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)    binary[i]=6;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)    binary[i]=7;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)    binary[i]=8;
    }
    EFp_set(&next_tmp_P,&table[binary[JSF_length]]);
    //SCM
    for(i=JSF_length-1; i>=0; i--){
        EFp_ECD(&next_tmp_P,&next_tmp_P);
        EFp_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
    }
    EFp16_set_EFp(ANS,&next_tmp_P);
    
    mpz_clear(buf);
    EFp_clear(&next_tmp_P);
    EFp_clear(&tmp_P);
    EFp_clear(&tmp_P_neg);
    EFp_clear(&skew_P);
    EFp_clear(&skew_7P);
    EFp_clear(&tmp_P_x4);
    EFp_clear(&tmp_P_x4_neg);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    for(i=0; i<9; i++){
        EFp_clear(&table[i]);
    }
}
void check_G1_SCM(){
    struct timeval t0,t1;
    struct EFp P_EFp;
    EFp_init(&P_EFp);
    struct EFp16 P_EFp16,result,test1,test2,test3;
    EFp16_init(&P_EFp16);
    EFp16_init(&result);
    EFp16_init(&test1);
    EFp16_init(&test2);
    EFp16_init(&test3);
    mpz_t scalar;
    mpz_init(scalar);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    printf("\n-----------------------------------------------------------------------------------------------------------------\ntest\n\n");
    //select scalar
    mpz_urandomm(scalar,state,order_r);    //S1
    printf("scalar\n");
    gmp_printf("%Zd",scalar);
    printf("\n\n");
    
    //set P_EFp & EFP16
    EFp16_random_set_G1(&P_EFp16);
    EFp_set_EFp16(&P_EFp,&P_EFp16);
    
    printf("P\n");
    EFp16_printf(&P_EFp16);
    printf("\n");
    
    //currect result
    EFp16_SCM_BIN(&result,&P_EFp16,scalar);
    printf("[scalar]P\n");
    EFp16_printf(&result);
    printf("\n");
    
    //test1
    printf("--------------NORMAL G1 SCM---------------\n");
    gettimeofday(&t0,NULL);
    G1_SCM_normal(&test1,&P_EFp16,scalar);
    gettimeofday(&t1,NULL);
    printf("[scalar]P\n");
    EFp16_printf(&test1);
    printf("NORMAL G1 SCM : %.2f[ms]\n\n",timedifference_msec(t0,t1));
    
    //test2
    printf("--------------2SPLIT G1 SCM---------------\n");
    gettimeofday(&t0,NULL);
    G1_SCM_2split(&test2,&P_EFp16,scalar);
    gettimeofday(&t1,NULL);
    printf("[scalar]P\n");
    EFp16_printf(&test2);
    printf("2SPLIT G1 SCM : %.2f[ms]\n\n",timedifference_msec(t0,t1));
    
    //test3
    printf("-----------2SPLIT G1 SCM (JSF) ------------\n");
    gettimeofday(&t0,NULL);
    G1_SCM_2split_JSF(&test3,&P_EFp16,scalar);
    gettimeofday(&t1,NULL);
    printf("[scalar]P\n");
    EFp16_printf(&test3);
    printf("2SPLIT G1 SCM (JSF) : %.2f[ms]\n\n",timedifference_msec(t0,t1));
    
    if(Fp16_cmp(&test1.x,&test2.x)==0 && Fp16_cmp(&test1.y,&test2.y)==0
       && Fp16_cmp(&test1.x,&test3.x)==0 && Fp16_cmp(&test1.y,&test3.y)==0){
        printf("test success\n");
    }
    
    mpz_clear(scalar);
    EFp_clear(&P_EFp);
    EFp16_clear(&P_EFp16);
    EFp16_clear(&result);
    EFp16_clear(&test1);
    EFp16_clear(&test2);
    EFp16_clear(&test3);
}
//-----------------------------------------------------------------------------------------
#pragma mark G2 SCM methods
void G2_SCM_normal(struct EFp16 *ANS,struct EFp16 *P, mpz_t scalar){
    struct EFp4 twisted_P;
    EFp4_init(&twisted_P);
    
    EFp16_to_EFp4_map(&twisted_P,P);
    Init_scm_cost(&scm_cost);
    EFp4_SCM_BIN_Sparse(&twisted_P,&twisted_P,scalar);
    Print_scm_Cost(&scm_cost, "Normal SCM");
    EFp4_to_EFp16_map(ANS,&twisted_P);
    
    EFp4_clear(&twisted_P);
}
void G2_SCM_2split(struct EFp16 *ANS,struct EFp16 *P, mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    struct EFp4 next_tmp_P,tmp_P,tmp_P_x4,skew_P,skew_7P;
    EFp4_init(&next_tmp_P);
    EFp4_init(&tmp_P);
    EFp4_init(&skew_P);
    EFp4_init(&skew_7P);
    EFp4_init(&tmp_P_x4);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    struct EFp4 table[4];
    for(i=0; i<4; i++){
        EFp4_init(&table[i]);
    }
    
    //set tmp_P
    EFp16_to_EFp4_map(&tmp_P,P);
    //set tmp_P_x4
    EFp4_Skew_Frobenius_p4(&skew_P,&tmp_P);
//    Skew_Frobenius_map(&next_tmp_P, &tmp_P);
//    Skew_Frobenius_map(&skew_P, &next_tmp_P);
//    Skew_Frobenius_map(&next_tmp_P, &skew_P);
//    Skew_Frobenius_map(&skew_P, &next_tmp_P);
    
//    EFp4_ECD_Sparse(&skew_7P,&skew_P);              //[2]skew_P
//    EFp4_ECA(&skew_7P,&skew_7P,&skew_P);     //[3]skew_P
//    EFp4_ECD_Sparse(&skew_7P,&skew_7P);             //[6]skew_P
//    EFp4_ECA(&skew_7P,&skew_7P,&skew_P);     //[7]skew_P
//    EFp4_init(&skew_P);
//    EFp4_ECD_Sparse(&skew_P,&tmp_P);                //[2]tmp_P
//    EFp4_ECA(&skew_P,&skew_P,&tmp_P);        //[3]tmp_P
//    EFp4_ECD_Sparse(&skew_P,&skew_P);               //[6]tmp_P
//    EFp4_ECD_Sparse(&skew_P,&skew_P);               //[12]tmp_P
//    EFp4_ECD_Sparse(&skew_P,&skew_P);               //[24]tmp_P
//    EFp4_ECA(&tmp_P_x4,&skew_P,&skew_7P);
//    EFp4_neg(&tmp_P_x4,&tmp_P_x4);
    
    Init_scm_cost(&scm_cost);
     mpz_t tmp_scalar; mpz_init(tmp_scalar);
    mpz_set_ui(tmp_scalar,7);
    EFp4_SCM_BIN_Sparse(&skew_7P, &skew_P, tmp_scalar);         //[7]p^4
    mpz_set_ui(tmp_scalar,24);
    EFp4_SCM_BIN_Sparse(&skew_P, &tmp_P, tmp_scalar);         //[24]Q'
    scm_cost.pre_eca = scm_cost.loop_eca;
    scm_cost.pre_ecd = scm_cost.loop_ecd;
    
    EFp4_ECA(&tmp_P_x4, &skew_7P, &skew_P);
    scm_cost.pre_eca++;
    EFp4_neg(&tmp_P_x4, &tmp_P_x4);                         //-([7]p^4+[24]Q);
   
    
    
    //set table
    table[0].infity=1;                        //00
    EFp4_set(&table[1],&tmp_P);                //01
    EFp4_set(&table[2],&tmp_P_x4);            //10
    EFp4_ECA(&table[3],&tmp_P,&tmp_P_x4);    //11
    scm_cost.pre_eca++;
    
    //s0,s1
    mpz_pow_ui(buf,X,4);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);
    
    //binary
    loop_length=0;
    for(i=0; i<2; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        printf("G2 length_s%d:%d\n",i,length_s[i]);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    printf("\n");
    //set binary
    char binary_s[2][loop_length+1];
    char str[5],*e;
    int binary[loop_length+1];
    for(i=0; i<2; i++){
        if(length_s[i]==loop_length){
            mpz_get_str(binary_s[i],2,s[i]);
        }else{
            char binary_buf[loop_length+1];
            mpz_get_str(binary_buf,2,s[i]);
            memset(binary_s[i],'0',sizeof(binary_s[i]));
            memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
        }
    }
    for(i=0; i<loop_length; i++){
        sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
        binary[i]=(int)strtol(str,&e,2);
    }
    EFp4_set(&next_tmp_P,&table[binary[0]]);
    
    //SCM
    scm_cost.loop_eca=0;
    scm_cost.loop_ecd=0;
    for(i=1; i<loop_length; i++){
        EFp4_ECD_Sparse(&next_tmp_P,&next_tmp_P);
        EFp4_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
        scm_cost.loop_eca++;
        scm_cost.loop_ecd++;
    }
    
    EFp4_to_EFp16_map(ANS,&next_tmp_P);
    Print_scm_Cost(&scm_cost, "2 split");
    
    mpz_clear(buf);
    EFp4_clear(&next_tmp_P);
    EFp4_clear(&tmp_P);
    EFp4_clear(&skew_P);
    EFp4_clear(&skew_7P);
    EFp4_clear(&tmp_P_x4);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    for(i=0; i<4; i++){
        EFp4_clear(&table[i]);
    }
}
void G2_SCM_2split_JSF(struct EFp16 *ANS,struct EFp16 *P, mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;struct EFp4 next_tmp_P,tmp_P,tmp_P_neg,tmp_P_x4,tmp_P_x4_neg,skew_P,skew_7P;
    EFp4_init(&next_tmp_P);
    EFp4_init(&tmp_P);
    EFp4_init(&tmp_P_neg);
    EFp4_init(&skew_P);
    EFp4_init(&skew_7P);
    EFp4_init(&tmp_P_x4);
    EFp4_init(&tmp_P_x4_neg);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    struct EFp4 table[9];
    for(i=0; i<9; i++){
        EFp4_init(&table[i]);
    }
    
    Init_scm_cost(&scm_cost);
    //set tmp_P
    EFp16_to_EFp4_map(&tmp_P,P);
    EFp4_neg(&tmp_P_neg,&tmp_P);
    //set tmp_P_x4
    EFp4_Skew_Frobenius_p4(&skew_P,&tmp_P);
    EFp4_ECD_Sparse(&skew_7P,&skew_P);              //[2]skew_P
    EFp4_ECA(&skew_7P,&skew_7P,&skew_P);     //[3]skew_P
    EFp4_ECD_Sparse(&skew_7P,&skew_7P);             //[6]skew_P
    EFp4_ECA(&skew_7P,&skew_7P,&skew_P);     //[7]skew_P
    EFp4_init(&skew_P);
    EFp4_ECD_Sparse(&skew_P,&tmp_P);                //[2]tmp_P
    EFp4_ECA(&skew_P,&skew_P,&tmp_P);        //[3]tmp_P
    EFp4_ECD_Sparse(&skew_P,&skew_P);               //[6]tmp_P
    EFp4_ECD_Sparse(&skew_P,&skew_P);               //[12]tmp_P
    EFp4_ECD_Sparse(&skew_P,&skew_P);               //[24]tmp_P
    EFp4_ECA(&tmp_P_x4,&skew_P,&skew_7P);
    EFp4_neg(&tmp_P_x4,&tmp_P_x4);
    EFp4_neg(&tmp_P_x4_neg,&tmp_P_x4);
    scm_cost.pre_ecd = scm_cost.pre_ecd + 6;
    scm_cost.pre_eca = scm_cost.pre_eca + 4;
    
    //set table
    table[0].infity=1;                                 //00
    EFp4_set(&table[1],&tmp_P);                        //01
    EFp4_set(&table[2],&tmp_P_x4);                    //10
    EFp4_ECA(&table[3],&tmp_P_x4,&tmp_P);            //11
    EFp4_set(&table[4],&tmp_P_neg);                    //0-1
    EFp4_set(&table[5],&tmp_P_x4_neg);                //-10
    EFp4_ECA(&table[6],&tmp_P_x4_neg,&tmp_P_neg);    //-1-1
    EFp4_ECA(&table[7],&tmp_P_x4,&tmp_P_neg);        //1-1
    EFp4_ECA(&table[8],&tmp_P_x4_neg,&tmp_P);        //-11
     scm_cost.pre_eca = scm_cost.pre_eca + 4;
    //s0,s1
    mpz_pow_ui(buf,X,4);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);
    
    //get loop_length
    loop_length=0;
    for(i=0; i<2; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    //JSF
    int JSF_length;
    int JSF_binary[2][loop_length+1];
    //    char check[5];
    for(i=0; i<loop_length; i++){
        JSF_binary[0][i]=0;
        JSF_binary[1][i]=0;
    }
    int *JSF_pointer[2];
    JSF_pointer[0]=JSF_binary[0];
    JSF_pointer[1]=JSF_binary[1];
    JSF(JSF_pointer,s,&JSF_length);
    int binary[JSF_length+1];
    
    for(i=JSF_length; i>=0; i--){
        if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0)         binary[i]=0;
        else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1)     binary[i]=1;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0)     binary[i]=2;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)    binary[i]=3;
        else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)    binary[i]=4;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)    binary[i]=5;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)    binary[i]=6;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)    binary[i]=7;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)    binary[i]=8;
    }
    EFp4_set(&next_tmp_P,&table[binary[JSF_length]]);
    //SCM
    scm_cost.loop_eca = 0;
    scm_cost.loop_ecd = 0;
    for(i=JSF_length-1; i>=0; i--){
        EFp4_ECD_Sparse(&next_tmp_P,&next_tmp_P);
        EFp4_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
        scm_cost.loop_ecd++;
        scm_cost.loop_eca++;
    }
    EFp4_to_EFp16_map(ANS,&next_tmp_P);
    Print_scm_Cost(&scm_cost, "2 split JFS");
    
    mpz_clear(buf);
    EFp4_clear(&next_tmp_P);
    EFp4_clear(&tmp_P);
    EFp4_clear(&tmp_P_neg);
    EFp4_clear(&skew_P);
    EFp4_clear(&skew_7P);
    EFp4_clear(&tmp_P_x4);
    EFp4_clear(&tmp_P_x4_neg);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    for(i=0; i<9; i++){
        EFp4_clear(&table[i]);
    }
}

void G2_SCM_8split(struct EFp16 *ANS,struct EFp16 *Q,mpz_t S){
    //s=s0+s1[x]+s2[x^2]+s3[x^3]+s4[x^4]+s5[x^5]+s6[x^6]+s7[x^7]
    int i,length_s[8],loop_length;
    struct EFp4 next_twisted_Q,twisted_Q_x[8];
    EFp4_init(&next_twisted_Q);
    for(i=0; i<8; i++){
        EFp4_init(&twisted_Q_x[i]);
    }
    struct EFp4 table0to3[16],table4to7[16];
    for(i=0; i<16; i++){
        EFp4_init(&table0to3[i]);
        EFp4_init(&table4to7[i]);
    }
    
    struct EFp4 tmp_EFp4, skew_p, skew_p2, skew_p3, skew_p4, skew_p5, skew_p6, skew_p7;
    EFp4_init(&tmp_EFp4);
    EFp4_init(&skew_p);
    EFp4_init(&skew_p2);
    EFp4_init(&skew_p3);
    EFp4_init(&skew_p4);
    EFp4_init(&skew_p5);
    EFp4_init(&skew_p6);
    EFp4_init(&skew_p7);
    
    //set twisted_Q
    EFp16_to_EFp4_map(&twisted_Q_x[0],Q);
    EFp4_Skew_Frobenius_p (&skew_p,  &twisted_Q_x[0]);
    EFp4_Skew_Frobenius_p2(&skew_p2, &twisted_Q_x[0]);
    EFp4_Skew_Frobenius_p3(&skew_p3, &twisted_Q_x[0]);
    EFp4_Skew_Frobenius_p4(&skew_p4, &twisted_Q_x[0]);
    EFp4_Skew_Frobenius_p5(&skew_p5, &twisted_Q_x[0]);
    EFp4_Skew_Frobenius_p6(&skew_p6, &twisted_Q_x[0]);
    EFp4_Skew_Frobenius_p7(&skew_p7, &twisted_Q_x[0]);
    
    Init_scm_cost(&scm_cost);
    
    mpz_t tmp_scalar; mpz_init(tmp_scalar);
    //u \equiv 2p^5-p
    EFp4_ECD_Sparse(&twisted_Q_x[1], &skew_p5);            scm_cost.pre_ecd++;           //[2]p^5
    EFp4_neg(&tmp_EFp4, &skew_p);                           //-p
    EFp4_ECA(&twisted_Q_x[1], &twisted_Q_x[1], &tmp_EFp4); scm_cost.pre_eca++; //2p^5-p
    
   
    
    
    //u^2 \equiv -(4p^6+3p^2);
     mpz_set_ui(tmp_scalar,4);
//    EFp4_ECD_Sparse(&twisted_Q_x[2], &skew_p6);             //[2]p^6
//    EFp4_ECD_Sparse(&twisted_Q_x[2], &twisted_Q_x[2]);      //[4]p^6
    scm_cost.loop_eca=0;
    scm_cost.loop_ecd = 0;
    EFp4_SCM_BIN_Sparse(&twisted_Q_x[2], &skew_p6, tmp_scalar);
    scm_cost.pre_eca =  scm_cost.pre_eca + scm_cost.loop_eca;
    scm_cost.pre_ecd =  scm_cost.pre_ecd + scm_cost.loop_ecd;
    
    EFp4_ECD_Sparse(&tmp_EFp4, &skew_p2);                    scm_cost.pre_ecd++;              //[2]p^2
    EFp4_ECA(&tmp_EFp4, &tmp_EFp4, &skew_p2);                scm_cost.pre_eca++;              //[3]p^2
    EFp4_ECA(&twisted_Q_x[2], &twisted_Q_x[2], &tmp_EFp4);   scm_cost.pre_eca++;           //[4]p^6+[3]p^2;
    EFp4_neg(&twisted_Q_x[2], &twisted_Q_x[2]);              //-([4]p^6+[3]p^2);
    scm_cost.pre_eca++;
  
    
    //u^3 \equiv 11p^3-2p^7;
    mpz_set_ui(tmp_scalar,11);
    EFp4_SCM_BIN_Sparse(&twisted_Q_x[3], &skew_p3, tmp_scalar);         //[11]p^3
    EFp4_ECD_Sparse(&tmp_EFp4, &skew_p7);                               //[2]p^7
    EFp4_neg(&tmp_EFp4, &tmp_EFp4);
    EFp4_ECA(&twisted_Q_x[3], &twisted_Q_x[3], &tmp_EFp4);
    
    //u^4 \equiv -(7p^4+24);
    mpz_set_ui(tmp_scalar,7);
    EFp4_SCM_BIN_Sparse(&twisted_Q_x[4], &skew_p4, tmp_scalar);         //[7]p^4
    mpz_set_ui(tmp_scalar,24);
    EFp4_SCM_BIN_Sparse(&tmp_EFp4, &twisted_Q_x[0], tmp_scalar);        //[24]Q'
    EFp4_ECA(&twisted_Q_x[4], &twisted_Q_x[4], &tmp_EFp4);
    EFp4_neg(&twisted_Q_x[4], &twisted_Q_x[4]);                         //-([7]p^4+[24]Q);
    
     //u^5 \equiv 38p-41p^5);
    mpz_set_ui(tmp_scalar,38);
    EFp4_SCM_BIN_Sparse(&twisted_Q_x[5], &skew_p, tmp_scalar);         //[38]p
    mpz_set_ui(tmp_scalar,41);
    EFp4_SCM_BIN_Sparse(&tmp_EFp4, &skew_p5, tmp_scalar);                //[41]p^5
    EFp4_neg(&tmp_EFp4, &tmp_EFp4);
    EFp4_ECA(&twisted_Q_x[5], &twisted_Q_x [5], &tmp_EFp4);
    
    //u^6 \equiv 117p^6+44p^2);
    mpz_set_ui(tmp_scalar,117);
    EFp4_SCM_BIN_Sparse(&twisted_Q_x[6], &skew_p6, tmp_scalar);         //[117]p^6
    mpz_set_ui(tmp_scalar,44);
    EFp4_SCM_BIN_Sparse(&tmp_EFp4, &skew_p2, tmp_scalar);                //[44]p^2
    EFp4_ECA(&twisted_Q_x[6], &twisted_Q_x[6], &tmp_EFp4);
    
    //u^7 \equiv -278p^3-29p^7;
    mpz_set_ui(tmp_scalar,278);
    EFp4_SCM_BIN_Sparse(&twisted_Q_x[7], &skew_p3, tmp_scalar);         //[278]p^3
    mpz_set_ui(tmp_scalar,29);
    EFp4_SCM_BIN_Sparse(&tmp_EFp4, &skew_p7, tmp_scalar);                //[29]p^7
    EFp4_ECA(&twisted_Q_x[7], &twisted_Q_x[7], &tmp_EFp4);
    EFp4_neg(&twisted_Q_x[7], &twisted_Q_x[7]);
    
    
    //set table
    table0to3[0].infity=1;                                      //0000
    EFp4_set(&table0to3[1],&twisted_Q_x[0]);                    //0001
    EFp4_set(&table0to3[2],&twisted_Q_x[1]);                    //0010
    EFp4_ECA(&table0to3[3],&twisted_Q_x[0],&twisted_Q_x[1]);    //0011
    EFp4_set(&table0to3[4],&twisted_Q_x[2]);                    //0100
    EFp4_ECA(&table0to3[5],&twisted_Q_x[2],&twisted_Q_x[0]);    //0101
    EFp4_ECA(&table0to3[6],&twisted_Q_x[2],&twisted_Q_x[1]);    //0110
    EFp4_ECA(&table0to3[7],&table0to3[6],&twisted_Q_x[0]);      //0111
    EFp4_set(&table0to3[8],&twisted_Q_x[3]);                    //1000
    EFp4_ECA(&table0to3[9],&twisted_Q_x[3],&twisted_Q_x[0]);    //1001
    EFp4_ECA(&table0to3[10],&twisted_Q_x[3],&twisted_Q_x[1]);   //1010
    EFp4_ECA(&table0to3[11],&table0to3[10],&twisted_Q_x[0]);    //1011
    EFp4_ECA(&table0to3[12],&twisted_Q_x[3],&twisted_Q_x[2]);   //1100
    EFp4_ECA(&table0to3[13],&table0to3[12],&twisted_Q_x[0]);    //1101
    EFp4_ECA(&table0to3[14],&table0to3[12],&twisted_Q_x[1]);    //1110
    EFp4_ECA(&table0to3[15],&table0to3[14],&twisted_Q_x[0]);    //1111
    
    table4to7[0].infity=1;                                      //0000
    EFp4_set(&table4to7[1],&twisted_Q_x[4]);                    //0001
    EFp4_set(&table4to7[2],&twisted_Q_x[5]);                    //0010
    EFp4_ECA(&table4to7[3],&twisted_Q_x[4],&twisted_Q_x[5]);    //0011
    EFp4_set(&table4to7[4],&twisted_Q_x[6]);                    //0100
    EFp4_ECA(&table4to7[5],&twisted_Q_x[6],&twisted_Q_x[4]);    //0101
    EFp4_ECA(&table4to7[6],&twisted_Q_x[6],&twisted_Q_x[5]);    //0110
    EFp4_ECA(&table4to7[7],&table4to7[6],&twisted_Q_x[4]);      //0111
    EFp4_set(&table4to7[8],&twisted_Q_x[7]);                    //1000
    EFp4_ECA(&table4to7[9],&twisted_Q_x[7],&twisted_Q_x[4]);    //1001
    EFp4_ECA(&table4to7[10],&twisted_Q_x[7],&twisted_Q_x[5]);   //1010
    EFp4_ECA(&table4to7[11],&table4to7[10],&twisted_Q_x[4]);    //1011
    EFp4_ECA(&table4to7[12],&twisted_Q_x[7],&twisted_Q_x[6]);   //1100
    EFp4_ECA(&table4to7[13],&table4to7[12],&twisted_Q_x[4]);    //1101
    EFp4_ECA(&table4to7[14],&table4to7[12],&twisted_Q_x[5]);    //1110
    EFp4_ECA(&table4to7[15],&table4to7[14],&twisted_Q_x[4]);    //1111
    
    mpz_t A,B,C,D,s[8],x_1,x_2,x_4;
    mpz_init(A);
    mpz_init(B);
    mpz_init(C);
    mpz_init(D);
    mpz_init(x_1);
    mpz_init(x_2);
    mpz_init(x_4);
    for(i=0; i<8; i++){
        mpz_init(s[i]);
    }
    
    mpz_set(x_1,X);
    mpz_mul(x_2,x_1,x_1);
    mpz_mul(x_4,x_2,x_2);
    
    //s0,s1,s2,s3,s4,s5,s6,s7
    mpz_tdiv_qr(B,A,S,x_4);
    mpz_tdiv_qr(D,C,A,x_2);
    mpz_tdiv_qr(s[1],s[0],C,x_1);
    mpz_tdiv_qr(s[3],s[2],D,x_1);
    mpz_tdiv_qr(D,C,B,x_2);
    mpz_tdiv_qr(s[5],s[4],C,x_1);
    mpz_tdiv_qr(s[7],s[6],D,x_1);
    
    //binary
    loop_length=0;
    for(i=0; i<8; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        printf("G2 length_s%d:%d\n",i,length_s[i]);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    //set binary
    char binary_s[8][loop_length+1];
    char str[5],*e;
    long binary0to3[loop_length],binary4to7[loop_length];
    for(i=0; i<8; i++){
        if(length_s[i]==loop_length){
            mpz_get_str(binary_s[i],2,s[i]);
        }
        else{
            char binary_buf[loop_length+1];
            mpz_get_str(binary_buf,2,s[i]);
            memset(binary_s[i],'0',sizeof(binary_s[i]));
            memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
        }
    }
    
    for(i=0; i<loop_length; i++){
        sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
        binary0to3[i]=strtol(str,&e,2);
        sprintf(str,"%c%c%c%c",binary_s[7][i],binary_s[6][i],binary_s[5][i],binary_s[4][i]);
        binary4to7[i]=strtol(str,&e,2);
    }
    EFp4_ECA(&next_twisted_Q,&table0to3[binary0to3[0]],&table4to7[binary4to7[0]]);
    
    
    for(i=1; i<loop_length; i++){
        EFp4_ECD_Sparse(&next_twisted_Q,&next_twisted_Q);
        EFp4_ECA(&next_twisted_Q,&next_twisted_Q,&table0to3[binary0to3[i]]);
        EFp4_ECA(&next_twisted_Q,&next_twisted_Q,&table4to7[binary4to7[i]]);
    }
    EFp4_to_EFp16_map(ANS,&next_twisted_Q);
    ANS->infity=next_twisted_Q.infity;
    
    mpz_clear(A);
    mpz_clear(B);
    mpz_clear(C);
    mpz_clear(D);
    mpz_clear(x_1);
    mpz_clear(x_2);
    mpz_clear(x_4);
    EFp4_clear(&next_twisted_Q);
    for(i=0; i<8; i++){
        EFp4_clear(&twisted_Q_x[i]);
    }
    for(i=0; i<8; i++){
        mpz_clear(s[i]);
    }
    for(i=0; i<16; i++){
        EFp4_clear(&table0to3[i]);
        EFp4_clear(&table4to7[i]);
    }
    
    mpz_clear(tmp_scalar);
    EFp4_clear(&skew_p);
    EFp4_clear(&skew_p2);
    EFp4_clear(&skew_p3);
    EFp4_clear(&skew_p4);
    EFp4_clear(&skew_p5);
    EFp4_clear(&skew_p6);
    EFp4_clear(&skew_p7);
    EFp4_clear(&tmp_EFp4);
}

void G2_SCM_4split(struct EFp16 *ANS,struct EFp16 *Q,mpz_t S){
    //s=s0+s1[x^2]+s2[x^4]+s3[x^6]
    long i,length_s[4],loop_length;
    struct EFp4 next_twisted_Q,twisted_Q,twisted_Q_2x,twisted_Q_4x,twisted_Q_6x;
    EFp4_init(&next_twisted_Q);
    EFp4_init(&twisted_Q);
    EFp4_init(&twisted_Q_2x);
    EFp4_init(&twisted_Q_4x);
    EFp4_init(&twisted_Q_6x);
    
    struct EFp4 tmp_EFp4, skew_p2, skew_p4, skew_p6;
    EFp4_init(&tmp_EFp4);
    EFp4_init(&skew_p2);
    EFp4_init(&skew_p4);
    EFp4_init(&skew_p6);
    
    //set twisted_Q
    EFp16_to_EFp4_map(&twisted_Q,Q);                        //twisted_Q
    EFp4_Skew_Frobenius_p2(&skew_p2,&twisted_Q);            //twisted_Q_2x
    EFp4_Skew_Frobenius_p4(&skew_p4,&twisted_Q);            //twisted_Q_4x
    EFp4_Skew_Frobenius_p6(&skew_p6,&twisted_Q);            //twisted_Q_6x
    
    mpz_t tmp_scalar; mpz_init(tmp_scalar);
    
    //u^2 \equiv -(4p^6+3p^2);
    mpz_set_ui(tmp_scalar,4);
    EFp4_SCM_BIN_Sparse(&twisted_Q_2x, &skew_p6, tmp_scalar);
    EFp4_ECD_Sparse(&tmp_EFp4, &skew_p2);                    //[2]p^2
    EFp4_ECA(&tmp_EFp4, &tmp_EFp4, &skew_p2);                //[3]p^2
    EFp4_ECA(&twisted_Q_2x, &twisted_Q_2x, &tmp_EFp4);       //[4]p^6+[3]p^2;
    EFp4_neg(&twisted_Q_2x, &twisted_Q_2x);                  //-([4]p^6+[3]p^2);

    //u^4 \equiv -(7p^4+24);
    mpz_set_ui(tmp_scalar,7);
    EFp4_SCM_BIN_Sparse(&twisted_Q_4x, &skew_p4, tmp_scalar);         //[7]p^4
    mpz_set_ui(tmp_scalar,24);
    EFp4_SCM_BIN_Sparse(&tmp_EFp4, &twisted_Q, tmp_scalar);        //[24]Q'
    EFp4_ECA(&twisted_Q_4x, &twisted_Q_4x, &tmp_EFp4);
    EFp4_neg(&twisted_Q_4x, &twisted_Q_4x);                         //-([7]p^4+[24]Q);
    
    //u^6 \equiv 117p^6+44p^2);
    mpz_set_ui(tmp_scalar,117);
    EFp4_SCM_BIN_Sparse(&twisted_Q_6x, &skew_p6, tmp_scalar);         //[117]p^6
    mpz_set_ui(tmp_scalar,44);
    EFp4_SCM_BIN_Sparse(&tmp_EFp4, &skew_p2, tmp_scalar);             //[44]p^2
    EFp4_ECA(&twisted_Q_6x, &twisted_Q_6x, &tmp_EFp4);
    
    //table
    struct EFp4 table[16];
    for(i=0; i<16; i++){
        EFp4_init(&table[i]);
    }
    //set table
    table[0].infity=1;                                   //0000
    EFp4_set(&table[1],&twisted_Q);                      //0001
    EFp4_set(&table[2],&twisted_Q_2x);                   //0010
    EFp4_ECA(&table[3],&twisted_Q_2x,&twisted_Q);        //0011
    EFp4_set(&table[4],&twisted_Q_4x);                   //0100
    EFp4_ECA(&table[5],&twisted_Q_4x,&twisted_Q);        //0101
    EFp4_ECA(&table[6],&twisted_Q_4x,&twisted_Q_2x);     //0110
    EFp4_ECA(&table[7],&table[6],&twisted_Q);            //0111
    EFp4_set(&table[8],&twisted_Q_6x);                   //1000
    EFp4_ECA(&table[9],&twisted_Q_6x,&twisted_Q);        //1001
    EFp4_ECA(&table[10],&twisted_Q_6x,&twisted_Q_2x);    //1010
    EFp4_ECA(&table[11],&twisted_Q_6x,&table[3]);        //1011
    EFp4_ECA(&table[12],&twisted_Q_6x,&twisted_Q_4x);    //1100
    EFp4_ECA(&table[13],&table[12],&twisted_Q);          //1101
    EFp4_ECA(&table[14],&table[12],&twisted_Q_2x);       //1110
    EFp4_ECA(&table[15],&table[14],&twisted_Q);          //1111
    
    
    mpz_t A,B,s[4],x_4,x_2;
    mpz_init(A);
    mpz_init(B);
    mpz_init(x_2);
    mpz_init(x_4);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //s0,s1,s2,s3
    mpz_mul(x_2,X,X);
    mpz_mul(x_4,x_2,x_2);
    mpz_tdiv_qr(B,A,S,x_4);
    mpz_tdiv_qr(s[1],s[0],A,x_2);
    mpz_tdiv_qr(s[3],s[2],B,x_2);
    
    //binary
    loop_length=0;
    for(i=0; i<4; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        printf("G2 length_s%ld:%ld\n",i,length_s[i]);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    printf("\n");
    //set binary
    char binary_s[4][loop_length+1];
    char str[5],*e;
    long binary[loop_length+1];
    for(i=0; i<4; i++){
        if(length_s[i]==loop_length){
            mpz_get_str(binary_s[i],2,s[i]);
        }
        else{
            char binary_buf[loop_length+1];
            mpz_get_str(binary_buf,2,s[i]);
            memset(binary_s[i],'0',sizeof(binary_s[i]));
            memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
        }
    }
    for(i=0; i<loop_length; i++){
        sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
        binary[i]=strtol(str,&e,2);
    }
    
    EFp4_set(&next_twisted_Q,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        EFp4_ECD_Sparse(&next_twisted_Q,&next_twisted_Q);
        EFp4_ECA(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
    }
    
    EFp4_to_EFp16_map(ANS,&next_twisted_Q);
    ANS->infity=next_twisted_Q.infity;
    
    EFp4_clear(&next_twisted_Q);
    EFp4_clear(&twisted_Q);
    EFp4_clear(&twisted_Q_2x);
    EFp4_clear(&twisted_Q_4x);
    EFp4_init(&twisted_Q_6x);
    mpz_clear(x_2);
    mpz_clear(x_4);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
    for(i=0; i<16; i++){
        EFp4_clear(&table[i]);
    }
    EFp4_clear(&tmp_EFp4);
    EFp4_clear(&skew_p2);
    EFp4_clear(&skew_p4);
    EFp4_clear(&skew_p6);
    mpz_clear(tmp_scalar);
}


void check_G2_SCM(){
    struct timeval t0,t1;
    struct EFp4 P_EFp4;
    EFp4_init(&P_EFp4);
    struct EFp16 P_EFp16,result,test1,test2,test3;
    EFp16_init(&P_EFp16);
    EFp16_init(&result);
    EFp16_init(&test1);
    EFp16_init(&test2);
    EFp16_init(&test3);
    mpz_t scalar;
    mpz_init(scalar);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    printf("\n-----------------------------------------------------------------------------------------------------------------\ntest\n\n");
    //select scalar
    mpz_urandomm(scalar,state,order_r);    //S1
    printf("scalar\n");
    gmp_printf("%Zd",scalar);
    printf("\n\n");
    
    //set P_EFp & EFP16
    EFp16_random_set_G2(&P_EFp16);
    printf("Q\n");
    EFp16_printf(&P_EFp16);
    printf("\n\n");
    
    //currect result
    EFp16_SCM_BIN(&result,&P_EFp16,scalar);
    printf("[scalar]P\n");
    EFp16_printf(&result);
    printf("\n");
    
    //test1
    printf("--------------NORMAL G2 SCM---------------\n");
    gettimeofday(&t0,NULL);
    G2_SCM_normal(&test1,&P_EFp16,scalar);
    gettimeofday(&t1,NULL);
    printf("[scalar]P\n");
    EFp16_printf(&test1);
    printf("NORMAL G2 SCM : %.2f[ms]\n\n",timedifference_msec(t0,t1));
    
    //test2
    printf("--------------2SPLIT G2 SCM---------------\n");
    gettimeofday(&t0,NULL);
    G2_SCM_2split(&test2,&P_EFp16,scalar);
    gettimeofday(&t1,NULL);
    printf("[scalar]P\n");
    EFp16_printf(&test2);
    printf("2SPLIT G2 SCM : %.2f[ms]\n\n",timedifference_msec(t0,t1));
    
    printf("--------------4 SPLIT G2 SCM---------------\n");
    gettimeofday(&t0,NULL);
    G2_SCM_4split(&test2,&P_EFp16,scalar);
    gettimeofday(&t1,NULL);
    printf("[scalar]P\n");
    EFp16_printf(&test2);
    printf("4 SPLIT G2 SCM : %.2f[ms]\n\n",timedifference_msec(t0,t1));
    
    
    printf("--------------8 SPLIT G2 SCM---------------\n");
    gettimeofday(&t0,NULL);
    G2_SCM_8split(&test2,&P_EFp16,scalar);
    gettimeofday(&t1,NULL);
    printf("[scalar]P\n");
    EFp16_printf(&test2);
    printf("8 SPLIT G2 SCM : %.2f[ms]\n\n",timedifference_msec(t0,t1));
    
    
    //test3
    printf("-----------2SPLIT G2 SCM (JSF) ------------\n");
    gettimeofday(&t0,NULL);
    G2_SCM_2split_JSF(&test3,&P_EFp16,scalar);
    gettimeofday(&t1,NULL);
    printf("[scalar]P\n");
    EFp16_printf(&test3);
    printf("2SPLIT G1 SCM (JSF): %.2f[ms]\n\n",timedifference_msec(t0,t1));
    
    if(Fp16_cmp(&test1.x,&test2.x)==0 && Fp16_cmp(&test1.y,&test2.y)==0
       && Fp16_cmp(&test1.x,&test3.x)==0 && Fp16_cmp(&test1.y,&test3.y)==0){
        printf("test success\n");
    }
    
    mpz_clear(scalar);
    EFp4_clear(&P_EFp4);
    EFp16_clear(&P_EFp16);
    EFp16_clear(&result);
    EFp16_clear(&test1);
    EFp16_clear(&test2);
    EFp16_clear(&test3);
}
void check_G2_order(){
    printf("---------------------------------------------------------------\n");
    struct EFp4 P_EFp4,test;
    EFp4_init(&P_EFp4);
    EFp4_init(&test);
    struct Fp4 beta_inv,x_tmp1,x_tmp2,beta1;
    Fp4_init(&beta_inv);
    Fp4_init(&x_tmp1);
    Fp4_init(&x_tmp2);
    Fp4_init(&beta1);
    mpz_t f4,buf,order1,order2,order3,order4,mod_r;
    mpz_init(f4);
    mpz_init(buf);
    mpz_init(order1);
    mpz_init(order2);
    mpz_init(order3);
    mpz_init(order4);
    mpz_init(mod_r);
    
    
    mpz_mul(buf,trace4,trace4);     //t^2
    mpz_mul_ui(f4,prime4,4);        //4p
    mpz_sub(f4,f4,buf);             //4p-t^2
    if(mpz_perfect_square_p(f4)!=0){
        mpz_sqrt(f4,f4);
        printf("f4 : ");
        mpz_out_str(stdout,10,f4);
        printf("\n");
    }else{
        printf("failed\n");
    }
    mpz_add_ui(order1,prime4,1);
    mpz_sub(order1,order1,f4);          //order1 : p^4+1-f4
    mpz_add_ui(order2,prime4,1);
    mpz_add(order2,order2,f4);          //order2 : p^4+1+f4
    
    mpz_add_ui(order3,prime4,1);
    mpz_sub(order3,order3,trace4);      //order3 : p^4+1-t4
    mpz_add_ui(order4,prime4,1);
    mpz_add(order4,order4,trace4);      //order4 : p^4+1+t4
    
    printf("order1 (p^4+1-f4) : ");
    mpz_out_str(stdout,10,order1);
    printf("\n");
    printf("order2 (p^4+1+f4) : ");
    mpz_out_str(stdout,10,order2);
    printf("\n");
    printf("order3 (p^4+1-t4) : ");
    mpz_out_str(stdout,10,order3);
    printf("\n");
    printf("order4 (p^4+1+t4) : ");
    mpz_out_str(stdout,10,order4);
    printf("\n");
    
    mpz_mod(mod_r,order1,order_r);
    printf("order1 (p^4+1-f4) mod r result : ");
    mpz_out_str(stdout,10,mod_r);
    printf("\n");
    mpz_mod(mod_r,order2,order_r);
    printf("order2 (p^4+1+f4) mod r result : ");
    mpz_out_str(stdout,10,mod_r);
    printf("\n");
    mpz_mod(mod_r,order3,order_r);
    printf("order3 (p^4+1-t4) mod r result : ");
    mpz_out_str(stdout,10,mod_r);
    printf("\n");
    mpz_mod(mod_r,order4,order_r);
    printf("order4 (p^4+1+t4) mod r result : ");
    mpz_out_str(stdout,10,mod_r);
    printf("\n\n");
    
    
    //d=4 twist order = order2
    
    Fp_set_ui(&beta1.x1.x0,1);
    Fp4_invert(&beta_inv,&beta1);
    while(1){
        Fp4_random(&P_EFp4.x);
        Fp4_mul(&x_tmp1,&P_EFp4.x,&P_EFp4.x);
        Fp4_mul(&x_tmp1,&x_tmp1,&P_EFp4.x);
        Fp4_mul_mpz(&x_tmp2,&P_EFp4.x,a_x);
        Fp4_mul(&x_tmp2,&x_tmp2,&beta_inv);
        Fp4_add(&x_tmp1,&x_tmp1,&x_tmp2);
        if(Fp4_legendre(&x_tmp1)==1){
            Fp4_sqrt(&P_EFp4.y,&x_tmp1);
            break;
        }
    }
    printf("P in E':y^2=x^3+a_x*x*β\n");
    EFp4_printf(&P_EFp4);
    printf("\n");
    
    printf("result\n");
    EFp4_init(&test);
    EFp4_SCM_BIN_Sparse(&test,&P_EFp4,order1);
    if(test.infity==1){
        printf("order1 (p^4+1-f4)\n");
    }
    EFp4_init(&test);
    EFp4_SCM_BIN_Sparse(&test,&P_EFp4,order2);
    if(test.infity==1){
        printf("order2 (p^4+1+f4)\n");
    }
    EFp4_init(&test);
    EFp4_SCM_BIN_Sparse(&test,&P_EFp4,order3);
    if(test.infity==1){
        printf("order3 (p^4+1-t4)\n");
    }
    EFp4_init(&test);
    EFp4_SCM_BIN_Sparse(&test,&P_EFp4,order4);
    if(test.infity==1){
        printf("order4 (p^4+1+t4)\n");
    }
    
    
    mpz_clear(f4);
    mpz_clear(buf);
    mpz_clear(order1);
    mpz_clear(order2);
    mpz_clear(order3);
    mpz_clear(order4);
    mpz_clear(mod_r);
    Fp4_clear(&beta_inv);
    Fp4_clear(&beta1);
    Fp4_clear(&x_tmp1);
    Fp4_clear(&x_tmp2);
    EFp4_clear(&P_EFp4);
    EFp4_clear(&test);
}
//-----------------------------------------------------------------------------------------
#pragma mark G3 EXP methods
void G3_EXP_2split(struct Fp16 *ANS,struct Fp16 *A, mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    struct Fp16 next_tmp_A,tmp_A,tmp_A_x4,skew_A,skew_7P;
    Fp16_init(&next_tmp_A);
    Fp16_init(&tmp_A);
    Fp16_init(&skew_A);
    Fp16_init(&skew_7P);
    Fp16_init(&tmp_A_x4);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    struct Fp16 table[4];
    for(i=0; i<4; i++){
        Fp16_init(&table[i]);
    }
    
    //set tmp_A
    Fp16_set(&tmp_A,A);
    //set tmp_A_x4
    Fp16_frobenius_map_p4(&skew_A,&tmp_A);
    Fp16_sqr(&skew_7P,&skew_A);                     //[2]skew_A
    Fp16_mul(&skew_7P,&skew_7P,&skew_A);            //[3]skew_A
    Fp16_sqr(&skew_7P,&skew_7P);                    //[6]skew_A
    Fp16_mul(&skew_7P,&skew_7P,&skew_A);            //[7]skew_A
    Fp16_init(&skew_A);
    Fp16_sqr(&skew_A,&tmp_A);                       //[2]tmp_A
    Fp16_mul(&skew_A,&skew_A,&tmp_A);               //[3]tmp_A
    Fp16_sqr(&skew_A,&skew_A);                      //[6]tmp_A
    Fp16_sqr(&skew_A,&skew_A);                      //[12]tmp_A
    Fp16_sqr(&skew_A,&skew_A);                      //[24]tmp_A
    Fp16_mul(&tmp_A_x4,&skew_A,&skew_7P);
    Fp16_frobenius_map_p8(&tmp_A_x4,&tmp_A_x4);
    
    //set table
    Fp_set_ui(&table[0].x0.x0.x0.x0,1);        //00
    Fp16_set(&table[1],&tmp_A);                //01
    Fp16_set(&table[2],&tmp_A_x4);            //10
    Fp16_mul(&table[3],&tmp_A,&tmp_A_x4);    //11
    
    //s0,s1
    mpz_pow_ui(buf,X,4);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);
    
    //binary
    loop_length=0;
    for(i=0; i<2; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        printf("G2 length_s%d:%d\n",i,length_s[i]);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    printf("\n");
    //set binary
    char binary_s[2][loop_length+1];
    char str[5],*e;
    int binary[loop_length+1];
    for(i=0; i<2; i++){
        if(length_s[i]==loop_length){
            mpz_get_str(binary_s[i],2,s[i]);
        }else{
            char binary_buf[loop_length+1];
            mpz_get_str(binary_buf,2,s[i]);
            memset(binary_s[i],'0',sizeof(binary_s[i]));
            memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
        }
    }
    for(i=0; i<loop_length; i++){
        sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
        binary[i] = (int)strtol(str,&e,2);
    }
    Fp16_set(&next_tmp_A,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        Fp16_sqr(&next_tmp_A,&next_tmp_A);
        Fp16_mul(&next_tmp_A,&next_tmp_A,&table[binary[i]]);
    }
    
    Fp16_set(ANS,&next_tmp_A);
    
    mpz_clear(buf);
    Fp16_clear(&next_tmp_A);
    Fp16_clear(&tmp_A);
    Fp16_clear(&skew_A);
    Fp16_clear(&skew_7P);
    Fp16_clear(&tmp_A_x4);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    for(i=0; i<4; i++){
        Fp16_clear(&table[i]);
    }
}
void G3_EXP_2split_JSF(struct Fp16 *ANS,struct Fp16 *A, mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    struct Fp16 next_tmp_A,tmp_A,tmp_A_neg,tmp_A_x4,tmp_A_x4_neg,skew_A,skew_7P;
    Fp16_init(&next_tmp_A);
    Fp16_init(&tmp_A);
    Fp16_init(&tmp_A_neg);
    Fp16_init(&skew_A);
    Fp16_init(&skew_7P);
    Fp16_init(&tmp_A_x4);
    Fp16_init(&tmp_A_x4_neg);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    struct Fp16 table[9];
    for(i=0; i<9; i++){
        Fp16_init(&table[i]);
    }
    
    //set tmp_A
    Fp16_set(&tmp_A,A);
    Fp16_frobenius_map_p8(&tmp_A_neg,&tmp_A);
    //set tmp_A_x4
    Fp16_frobenius_map_p4(&skew_A,&tmp_A);
    Fp16_sqr(&skew_7P,&skew_A);                     //[2]skew_A
    Fp16_mul(&skew_7P,&skew_7P,&skew_A);            //[3]skew_A
    Fp16_sqr(&skew_7P,&skew_7P);                    //[6]skew_A
    Fp16_mul(&skew_7P,&skew_7P,&skew_A);            //[7]skew_A
    Fp16_init(&skew_A);
    Fp16_sqr(&skew_A,&tmp_A);                       //[2]tmp_A
    Fp16_mul(&skew_A,&skew_A,&tmp_A);               //[3]tmp_A
    Fp16_sqr(&skew_A,&skew_A);                      //[6]tmp_A
    Fp16_sqr(&skew_A,&skew_A);                      //[12]tmp_A
    Fp16_sqr(&skew_A,&skew_A);                      //[24]tmp_A
    Fp16_mul(&tmp_A_x4,&skew_A,&skew_7P);
    Fp16_frobenius_map_p8(&tmp_A_x4,&tmp_A_x4);
    Fp16_frobenius_map_p8(&tmp_A_x4_neg,&tmp_A_x4);
    
    //set table
    Fp_set_ui(&table[0].x0.x0.x0.x0,1);                //00
    Fp16_set(&table[1],&tmp_A);                        //01
    Fp16_set(&table[2],&tmp_A_x4);                    //10
    Fp16_mul(&table[3],&tmp_A_x4,&tmp_A);            //11
    Fp16_set(&table[4],&tmp_A_neg);                    //0-1
    Fp16_set(&table[5],&tmp_A_x4_neg);                //-10
    Fp16_mul(&table[6],&tmp_A_x4_neg,&tmp_A_neg);    //-1-1
    Fp16_mul(&table[7],&tmp_A_x4,&tmp_A_neg);        //1-1
    Fp16_mul(&table[8],&tmp_A_x4_neg,&tmp_A);        //-11
    
    //s0,s1
    mpz_pow_ui(buf,X,4);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);
    
    //get loop_length
    loop_length=0;
    for(i=0; i<2; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    //JSF
    int JSF_length;
    int JSF_binary[2][loop_length+1];
    //    char check[5];
    for(i=0; i<loop_length; i++){
        JSF_binary[0][i]=0;
        JSF_binary[1][i]=0;
    }
    int *JSF_pointer[2];
    JSF_pointer[0]=JSF_binary[0];
    JSF_pointer[1]=JSF_binary[1];
    JSF(JSF_pointer,s,&JSF_length);
    int binary[JSF_length+1];
    
    for(i=JSF_length; i>=0; i--){
        if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0)         binary[i]=0;
        else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1)     binary[i]=1;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0)     binary[i]=2;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)    binary[i]=3;
        else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)    binary[i]=4;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)    binary[i]=5;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)    binary[i]=6;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)    binary[i]=7;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)    binary[i]=8;
    }
    Fp16_set(&next_tmp_A,&table[binary[JSF_length]]);
    //SCM
    for(i=JSF_length-1; i>=0; i--){
        Fp16_sqr(&next_tmp_A,&next_tmp_A);
        Fp16_mul(&next_tmp_A,&next_tmp_A,&table[binary[i]]);
    }
    Fp16_set(ANS,&next_tmp_A);
    
    mpz_clear(buf);
    Fp16_clear(&next_tmp_A);
    Fp16_clear(&tmp_A);
    Fp16_clear(&tmp_A_neg);
    Fp16_clear(&skew_A);
    Fp16_clear(&skew_7P);
    Fp16_clear(&tmp_A_x4);
    Fp16_clear(&tmp_A_x4_neg);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    for(i=0; i<9; i++){
        Fp16_clear(&table[i]);
    }
}
void check_G3_EXP(){
    printf("\n-----------------------------------------------------------------------------------------------------------------\ntest\n\n");
    struct timeval t0,t1;
    struct EFp16 P,Q,s1_P,s2_Q;
    EFp16_init(&P);
    EFp16_init(&Q);
    EFp16_init(&s1_P);
    EFp16_init(&s2_Q);
    struct EFp tmp_P;
    EFp_init(&tmp_P);
    struct Fp16 Z1,Z2,test1,test2,test3;
    Fp16_init(&Z1);
    Fp16_init(&Z2);
    Fp16_init(&test1);
    Fp16_init(&test2);
    Fp16_init(&test3);
    mpz_t s1,s2,s12;
    mpz_init(s1);
    mpz_init(s2);
    mpz_init(s12);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,order_r);     //s1
    mpz_urandomm(s2,state,order_r);     //s2
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order_r);           //s12
    
    EFp16_random_set_G1(&P);
    G1_SCM_2split(&s1_P,&P,s1);
    EFp16_random_set_G2(&Q);
    G2_SCM_2split(&s2_Q,&Q,s2);
    
    EFp_set_EFp16(&tmp_P,&s1_P);
    Pseudo_Sparse_Optimal_Ate_Pairing(&Z1,&tmp_P,&s2_Q);
    EFp_set_EFp16(&tmp_P,&P);
    Pseudo_Sparse_Optimal_Ate_Pairing(&Z2,&tmp_P,&Q);
    
    //test1
    printf("--------------NORMAL G3 SCM---------------\n");
    gettimeofday(&t0,NULL);
    Fp16_pow(&test1,&Z2,s12);
    gettimeofday(&t1,NULL);
    printf("Z^s12\n");
    Fp16_printf(&test1);
    printf("NORMAL G3 SCM : %.2f[ms]\n\n",timedifference_msec(t0,t1));
    
    //test2
    printf("--------------2SPLIT G3 SCM---------------\n");
    gettimeofday(&t0,NULL);
    G3_EXP_2split(&test2,&Z2,s12);
    gettimeofday(&t1,NULL);
    printf("Z^s12\n");
    Fp16_printf(&test2);
    printf("2SPLIT G3 SCM : %.2f[ms]\n\n",timedifference_msec(t0,t1));
    
    //test3
    printf("-----------2SPLIT G3 SCM (JSF)-------------\n");
    gettimeofday(&t0,NULL);
    G3_EXP_2split_JSF(&test3,&Z2,s12);
    gettimeofday(&t1,NULL);
    printf("Z^s12\n");
    Fp16_printf(&test3);
    printf("2SPLIT G3 SCM (JSF) : %.2f[ms]\n\n",timedifference_msec(t0,t1));
    
    
    printf("Bilinearity Check : ");
    if(Fp16_cmp(&Z1,&test1)==0 && Fp16_cmp(&Z1,&test2)==0 && Fp16_cmp(&Z1,&test3)==0){
        printf("success\n");
    }else{
        printf("failed\n");
    }
    
    
    EFp16_clear(&P);
    EFp16_clear(&Q);
    EFp16_clear(&s1_P);
    EFp16_clear(&s2_Q);
    EFp_clear(&tmp_P);
    Fp16_clear(&Z1);
    Fp16_clear(&Z2);
    Fp16_clear(&test1);
    Fp16_clear(&test2);
    Fp16_clear(&test3);
    mpz_clear(s1);
    mpz_clear(s2);
    mpz_clear(s12);
}

