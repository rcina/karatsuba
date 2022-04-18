/*
 * author: Robert Cina
 * created: April 18, 2022
 * description: This program implements the Karatsuba multiplication algorithm
 * using The GNU Multiple Precision Arithmetic Library (gmp.h)
 * 
 */
#include <stdio.h>
#include <gmp.h>
#include <math.h>

/*
 * getSize calculates the number of digits in num
 * @param {num} The number to calculate the number of digits of
 * @returns The count of the number of digits in num
 */
int getSize(mpz_t num){
    int count = 0;
    long int base = 10;
    mpz_t zero;
    mpz_t quotient, divisor;
    mpz_init(zero);
    mpz_set_si(zero,0);
    mpz_init_set_si(divisor, base);
    mpz_init_set(quotient,num);

    while(mpz_cmp(quotient, zero) > 0) {
        count++;
        mpz_tdiv_q(quotient,quotient,divisor);
    }
    return count;
}

/*
 * karatsuba performs the multiplication of two n-digit numbers
 * using the Karatsuba multiplication algorithm
 * @param {result} holds the product or result of the multiplication
 * @param {x} the first number or operand to multiply
 * @param {y} the second number or operanad to multiply
 */

void karatsuba(mpz_t result, mpz_t x, mpz_t y) {
    unsigned long int size; /* holds the max number of digits of x or y */
    unsigned long int n;  /* divides number of digits in half */
    unsigned long base = 10; /* represents the base of the decimal number system */

    mpz_t a; /* reprsents the upper half of the integer x */
    mpz_t b; /* represents the lower half of the integer x */
    mpz_t c; /* represents the upper half of the integer y */
    mpz_t d; /* represents the lower half of the integer y */
    mpz_t quotient; /* holds the result of division */
    mpz_t remainder;
    mpz_t power; /* holds the value of base raised to the power n */
    mpz_t ac, bd,abcd,mul_result; /* variables hold intermediate results of karatsuba algorithm */
    mpz_t sum_ab; /* holds sum of a + b */
    mpz_t sum_cd; /* holds sum of c + d */
    mpz_t sum; /* holds sum of intermediate calculation */

    /* term1 - term 4 hold results of intermediate calculations */
    mpz_t term1;
    mpz_t term2;
    mpz_t term3;
    mpz_t term4;

    mpz_init(a);
    mpz_init(b);
    mpz_init(c);
    mpz_init(d);
    mpz_init(quotient);
    mpz_init(remainder);
    mpz_init(power);
    mpz_init(ac);
    mpz_init(bd);
    mpz_init(abcd);
    mpz_init(mul_result);
    mpz_init(sum_ab);
    mpz_init(sum_cd);
    mpz_init(sum);
    mpz_init(term1);
    mpz_init(term2);
    mpz_init(term3);
    mpz_init(term4);
    
    
    size = (getSize(x) >= getSize(y)) ? getSize(x): getSize(y);
    n = (unsigned long int) ceil(size/2.0);

    mpz_ui_pow_ui(power, base, n);

    mpz_fdiv_q(quotient,x, power);
    mpz_set(a, quotient);
    mpz_fdiv_r(remainder,x, power);
    mpz_set(b, remainder);

    mpz_fdiv_q(quotient, y, power);
    mpz_set(c, quotient);
    mpz_fdiv_r(remainder, y, power);
    mpz_set(d, remainder);
    
    /* if size of problem is equal to one return the result of the multiplication */
    if(n < 2){
        mpz_mul(mul_result,x,y);
        mpz_add(result, result,mul_result);
        return;
    }
    /* otherwise recursively call karatsuba on a*c */
    karatsuba(ac,a,c);
    /* recursively call karatsuba on b*d */
    karatsuba(bd, b,d);
    /* recursively call karatsuba on (a+b) * (c+d) */
    mpz_add(sum_ab,a,b);
    mpz_add(sum_cd, c,d);
    karatsuba(abcd, sum_ab, sum_cd);
    
    /* calculate abcd term after recursive calls*/
    mpz_sub(abcd, abcd, ac);
    mpz_sub(abcd, abcd, bd);

    /* combine terms and sum karatsuba result of ac*10^(2*n) + abcd*10^n + bd */
    mpz_ui_pow_ui(term1, base,2*n);
    mpz_mul(term3, term1, ac);

    mpz_ui_pow_ui(term2, base,n);
    mpz_mul(term4, term2, abcd);

    mpz_add(sum, term3, term4);
    mpz_add(result, sum, bd);
}

int main(){
    mpz_t x,y, result;
    mpz_init(x);
    mpz_init(y);
    mpz_init(result);
    mpz_set_ui(result, 0);
    
    /* 64 digit numbers for test input */
    /* x=3141592653589793238462643383279502884197169399375105820974944592 */
    /* y=2718281828459045235360287471352662497757247093699959574966967627 */
    printf("Enter first number: ");
    gmp_scanf("%Zd",x);

    printf("Enter second number: ");
    gmp_scanf("%Zd",y);
    karatsuba(result, x,y);
    gmp_printf("%Zd x %Zd =\n %Zd\n", x,y,result);
    
    return 0;
}
