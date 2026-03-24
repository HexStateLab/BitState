/*
 * ═══════════════════════════════════════════════════════════════════════════════
 * BIGINT LIBRARY — 4096-bit Arbitrary Precision Integer Arithmetic
 * ═══════════════════════════════════════════════════════════════════════════════
 * C port of bigint.asm. Faithful translation of the NASM logic.
 */

#include "bigint.h"
#include <string.h>

/* ─── Core Operations ───────────────────────────────────────────────────────── */

void bigint_clear(BigInt *a)
{
    memset(a->limbs, 0, BIGINT_BYTES);
}

void bigint_copy(BigInt *dst, const BigInt *src)
{
    memcpy(dst->limbs, src->limbs, BIGINT_BYTES);
}

int bigint_is_zero(const BigInt *a)
{
    for (int i = 0; i < BIGINT_LIMBS; i++) {
        if (a->limbs[i] != 0) return 0;
    }
    return 1;
}

int bigint_cmp(const BigInt *a, const BigInt *b)
{
    for (int i = BIGINT_LIMBS - 1; i >= 0; i--) {
        if (a->limbs[i] > b->limbs[i]) return 1;
        if (a->limbs[i] < b->limbs[i]) return -1;
    }
    return 0;
}

/* ─── Arithmetic ────────────────────────────────────────────────────────────── */

void bigint_add(BigInt *result, const BigInt *a, const BigInt *b)
{
    __uint128_t carry = 0;
    for (int i = 0; i < BIGINT_LIMBS; i++) {
        __uint128_t sum = (__uint128_t)a->limbs[i] + b->limbs[i] + carry;
        result->limbs[i] = (uint64_t)sum;
        carry = sum >> 64;
    }
}

void bigint_sub(BigInt *result, const BigInt *a, const BigInt *b)
{
    __uint128_t borrow = 0;
    for (int i = 0; i < BIGINT_LIMBS; i++) {
        __uint128_t diff = (__uint128_t)a->limbs[i] - b->limbs[i] - borrow;
        result->limbs[i] = (uint64_t)diff;
        borrow = (diff >> 64) & 1;  /* Borrow if underflow */
    }
}

void bigint_mul(BigInt *result, const BigInt *a, const BigInt *b)
{
    uint64_t temp[BIGINT_LIMBS * 2];
    memset(temp, 0, sizeof(temp));

    for (int i = 0; i < BIGINT_LIMBS; i++) {
        if (a->limbs[i] == 0) continue;
        uint64_t carry = 0;
        for (int j = 0; j < BIGINT_LIMBS; j++) {
            __uint128_t prod = (__uint128_t)a->limbs[i] * b->limbs[j]
                             + temp[i + j] + carry;
            temp[i + j] = (uint64_t)prod;
            carry = (uint64_t)(prod >> 64);
        }
        if (i + BIGINT_LIMBS < BIGINT_LIMBS * 2) {
            temp[i + BIGINT_LIMBS] += carry;
        }
    }

    /* Truncate to BIGINT_LIMBS (lower half) */
    memcpy(result->limbs, temp, BIGINT_BYTES);
}

/* ─── Bit Operations ────────────────────────────────────────────────────────── */

void bigint_shl1(BigInt *a)
{
    uint64_t carry = 0;
    for (int i = 0; i < BIGINT_LIMBS; i++) {
        uint64_t new_carry = a->limbs[i] >> 63;
        a->limbs[i] = (a->limbs[i] << 1) | carry;
        carry = new_carry;
    }
}

void bigint_shr1(BigInt *a)
{
    uint64_t carry = 0;
    for (int i = BIGINT_LIMBS - 1; i >= 0; i--) {
        uint64_t new_carry = a->limbs[i] & 1;
        a->limbs[i] = (a->limbs[i] >> 1) | (carry << 63);
        carry = new_carry;
    }
}

int bigint_get_bit(const BigInt *a, uint32_t bit_index)
{
    uint32_t limb = bit_index / 64;
    uint32_t bit  = bit_index % 64;
    if (limb >= BIGINT_LIMBS) return 0;
    return (a->limbs[limb] >> bit) & 1;
}

void bigint_set_bit(BigInt *a, uint32_t bit_index)
{
    uint32_t limb = bit_index / 64;
    uint32_t bit  = bit_index % 64;
    if (limb >= BIGINT_LIMBS) return;
    a->limbs[limb] |= (1ULL << bit);
}

void bigint_clr_bit(BigInt *a, uint32_t bit_index)
{
    uint32_t limb = bit_index / 64;
    uint32_t bit  = bit_index % 64;
    if (limb >= BIGINT_LIMBS) return;
    a->limbs[limb] &= ~(1ULL << bit);
}

uint32_t bigint_bitlen(const BigInt *a)
{
    for (int i = BIGINT_LIMBS - 1; i >= 0; i--) {
        if (a->limbs[i] != 0) {
            /* Find highest set bit in this limb */
            uint32_t bit = 63;
            uint64_t v = a->limbs[i];
            while (!(v & (1ULL << bit))) bit--;
            return (uint32_t)i * 64 + bit + 1;
        }
    }
    return 0;
}

/* ─── Conversion ────────────────────────────────────────────────────────────── */

void bigint_set_u64(BigInt *a, uint64_t val)
{
    bigint_clear(a);
    a->limbs[0] = val;
}

uint64_t bigint_to_u64(const BigInt *a)
{
    return a->limbs[0];
}

/* ─── Division / Modulo (Long Division, Bit-by-Bit) ─────────────────────────── */

void bigint_div_mod(const BigInt *dividend, const BigInt *divisor,
                    BigInt *quotient, BigInt *remainder)
{
    bigint_clear(quotient);
    bigint_clear(remainder);

    uint32_t bits = bigint_bitlen(dividend);
    if (bits == 0) return;  /* dividend is zero */

    for (int i = (int)bits - 1; i >= 0; i--) {
        /* remainder <<= 1 */
        bigint_shl1(remainder);

        /* remainder.bit[0] = dividend.bit[i] */
        if (bigint_get_bit(dividend, (uint32_t)i)) {
            bigint_set_bit(remainder, 0);
        }

        /* if remainder >= divisor */
        if (bigint_cmp(remainder, divisor) >= 0) {
            BigInt temp;
            bigint_sub(&temp, remainder, divisor);
            bigint_copy(remainder, &temp);
            bigint_set_bit(quotient, (uint32_t)i);
        }
    }
}

/* ─── GCD (Euclidean) ───────────────────────────────────────────────────────── */

void bigint_gcd(BigInt *result, const BigInt *a, const BigInt *b)
{
    BigInt x, y, q, r;
    bigint_copy(&x, a);
    bigint_copy(&y, b);

    while (!bigint_is_zero(&y)) {
        bigint_div_mod(&x, &y, &q, &r);
        bigint_copy(&x, &y);
        bigint_copy(&y, &r);
    }

    bigint_copy(result, &x);
}

/* ─── Modular Exponentiation (Left-to-Right Binary) ─────────────────────────── */

void bigint_pow_mod(BigInt *result, const BigInt *base,
                    const BigInt *exp, const BigInt *mod)
{
    BigInt b, temp, q;
    bigint_copy(&b, base);
    bigint_set_u64(result, 1);

    uint32_t bits = bigint_bitlen(exp);

    for (uint32_t i = 0; i < bits; i++) {
        if (bigint_get_bit(exp, i)) {
            /* result = (result * base) % mod */
            bigint_mul(&temp, result, &b);
            bigint_div_mod(&temp, mod, &q, result);
        }
        /* base = (base * base) % mod */
        bigint_mul(&temp, &b, &b);
        bigint_div_mod(&temp, mod, &q, &b);
    }
}

/* ─── String Conversion ─────────────────────────────────────────────────────── */

int bigint_from_decimal(BigInt *a, const char *str)
{
    bigint_clear(a);
    if (!str || !*str) return -1;

    /* Skip leading whitespace */
    while (*str == ' ' || *str == '\t') str++;
    if (!*str) return -1;

    BigInt ten, digit_bi, temp;
    bigint_set_u64(&ten, 10);

    for (const char *p = str; *p; p++) {
        if (*p < '0' || *p > '9') {
            if (*p == '_' || *p == ' ') continue;  /* Allow digit separators */
            return -1;
        }
        uint64_t d = (uint64_t)(*p - '0');

        /* a = a * 10 + d */
        bigint_mul(&temp, a, &ten);
        bigint_copy(a, &temp);

        bigint_set_u64(&digit_bi, d);
        bigint_add(&temp, a, &digit_bi);
        bigint_copy(a, &temp);
    }

    return 0;
}

void bigint_to_decimal(char *buf, size_t bufsize, const BigInt *a)
{
    if (bufsize == 0) return;

    if (bigint_is_zero(a)) {
        buf[0] = '0';
        buf[1] = '\0';
        return;
    }

    /* Extract digits by repeated division by 10 */
    char tmp[1240];  /* 4096 bits ≈ 1233 decimal digits max */
    int pos = 0;

    BigInt val, quot, rem, ten;
    bigint_copy(&val, a);
    bigint_set_u64(&ten, 10);

    while (!bigint_is_zero(&val) && pos < 1239) {
        bigint_div_mod(&val, &ten, &quot, &rem);
        tmp[pos++] = '0' + (char)bigint_to_u64(&rem);
        bigint_copy(&val, &quot);
    }

    /* Reverse into output buffer */
    size_t len = (size_t)pos;
    if (len >= bufsize) len = bufsize - 1;
    for (size_t i = 0; i < len; i++) {
        buf[i] = tmp[pos - 1 - (int)i];
    }
    buf[len] = '\0';
}
