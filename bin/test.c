#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "ldpc_codes.h"
#include "ldpc_encoder.h"
#include "ldpc_decoder.h"

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

/* Test vectors for generator matrix.
 * This is the CRC32 of every uint32_t in the generator matrix,
 * in the byte-swapped endianness that init_generator gives you.
 */
const uint32_t test_vector_g[] = {
    [LDPC_CODE_N128_K64]    = 0xDC64D486,
    [LDPC_CODE_N256_K128]   = 0xD78B5564,
    [LDPC_CODE_N512_K256]   = 0x6AF9EC6A,
    [LDPC_CODE_N1280_K1024] = 0x452FE118,
    [LDPC_CODE_N1536_K1024] = 0xBCCBA8D0,
    [LDPC_CODE_N2048_K1024] = 0x1597B6F6,
};

/* Test vectors for parity check matrix.
 * This is the CRC32 of every uint32_t in the parity check matrix.
 */
const uint32_t test_vector_h[] = {
    [LDPC_CODE_N128_K64]    = 0x4FDF9E5A,
    [LDPC_CODE_N256_K128]   = 0x588971F8,
    [LDPC_CODE_N512_K256]   = 0x33BDB5C2,
    [LDPC_CODE_N1280_K1024] = 0x90224F9A,
    [LDPC_CODE_N1536_K1024] = 0x0A8EFA1C,
    [LDPC_CODE_N2048_K1024] = 0x2CD11363,
};

/* Test vectors for encoded message, where message is (0, 1, 2, 3, ..., k/8].
 * This is the CRC32 for the codeword.
 */
const uint32_t test_vector_txcode[] = {
    [LDPC_CODE_N128_K64]    = 0x07279866,
    [LDPC_CODE_N256_K128]   = 0x964F9176,
    [LDPC_CODE_N512_K256]   = 0x441CE45D,
    [LDPC_CODE_N1280_K1024] = 0x99AE48D8,
    [LDPC_CODE_N1536_K1024] = 0x3BA467B3,
    [LDPC_CODE_N2048_K1024] = 0xC7253610,
};

/* Test vectors for sparse parity check matrix.
 * Each is the CRC32 for {ci, cs, vi, vs}.
 */
const uint32_t test_vector_sparse_h[][4] = {
    [LDPC_CODE_N128_K64]    = {0xEFA94DF3, 0xAACD61F7, 0x2A434B32, 0xC2C357FD},
    [LDPC_CODE_N256_K128]   = {0x16D7B220, 0x45C89470, 0x8C0B1BB5, 0xB52CD823},
    [LDPC_CODE_N512_K256]   = {0xBF7FE8DE, 0x91ECF7A7, 0x9F3B08E5, 0xADF27DAC},
    [LDPC_CODE_N1280_K1024] = {0x0E90A077, 0x513AA9FC, 0xC04FC69E, 0x3015B990},
    [LDPC_CODE_N1536_K1024] = {0x79393640, 0x88868300, 0x5B47BD80, 0xC1270F59},
    [LDPC_CODE_N2048_K1024] = {0xA9E48635, 0x07AB0B57, 0x551B70B1, 0x09153CF1},
};

/* Simple CRC-32 implementation for comparing to test vectors. */
uint32_t crc32(uint8_t* data, size_t length)
{
    size_t i, j;
    uint32_t crc=0xFFFFFFFF, mask;
    for(i=0; i<length; i++) {
        crc ^= data[i];
        for(j=0; j<8; j++) {
            mask = -(crc & 1);
            crc = (crc >> 1) ^ (0xEDB88320 & mask);
        }
    }
    return ~crc;
}

bool test_code(enum ldpc_code code)
{
    /* Code parameters. */
    int n=0, k=0, p=0, m=0, b=0, s=0;

    /* Counters and timers and so forth for the test function. */
    int i;
    bool ok, overall_ok = true;
    struct timespec t_start, t_end;
    double time_taken;

    /* Pointers for memory we'll allocate to use the codes. */
    uint32_t *g, *h;
    size_t size_g, size_h;
    uint16_t *ci, *cs, *vi, *vs;
    size_t size_ci, size_cs, size_vi, size_vs;
    uint8_t *txdata, *txcode_small, *txcode_fast, *rxcode, *rxcode_llr;
    uint8_t *rxdata, *bf_wa;
    float *rxllrs, *mp_wa;
    size_t size_bf_wa, size_mp_wa;

    (void)rxdata; (void)bf_wa; (void)rxllrs; (void)mp_wa;
    (void)size_bf_wa; (void)size_mp_wa;

    /* Time this test. */
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t_start);

    ldpc_codes_get_params(code, &n, &k, &p, &m, &b, &s);

    printf("************************************************************\n");
    printf("* (%d, %d) code\n", n, k);
    printf("* Parameters: " KYEL "n=%d k=%d p=%d m=%d b=%d s=%d" KNRM "\n",
           n, k, p, m, b, s);
    printf("* ----------------------------------------------------------\n");

    /* Initialise the generator matrix for fast encoding. */
    printf("* Check generator matrix:                      ");
    size_g = ldpc_codes_size_generator(code);
    g = malloc(size_g);
    ldpc_codes_init_generator(code, g);
    if(crc32((uint8_t*)g, size_g) == test_vector_g[code]) {
        printf(KGRN "OK" KNRM "\n");
    } else {
        printf(KRED "FAIL" KNRM "\n");
        overall_ok = false;
    }

    /* Make up some data to transmit. */
    txdata = malloc(k/8);
    for(i=0; i<k/8; i++) {
        txdata[i] = ~i;
    }

    /* Perform a slow encode. */
    printf("* Check small encoder:                         ");
    txcode_small = malloc(n/8);
    ldpc_encode_small(code, txdata, txcode_small);
    if(crc32(txcode_small, n/8) == test_vector_txcode[code]) {
        printf(KGRN "OK" KNRM "\n");
    } else {
        printf(KRED "FAIL" KNRM "\n");
        overall_ok = false;
    }

    /* Perform a fast encode. */
    printf("* Check fast encoder:                          ");
    txcode_fast  = malloc(n/8);
    ldpc_encode_fast(code, g, txdata, txcode_fast);
    if(crc32(txcode_fast, n/8) == test_vector_txcode[code]) {
        printf(KGRN "OK" KNRM "\n");
    } else {
        printf(KRED "FAIL" KNRM "\n");
        overall_ok = false;
    }

    /* Check the encodings are equal. */
    printf("* Check fast vs small encoder:                 ");
    ok = true;
    for(i=0; i<n/8; i++) {
        if(txcode_small[i] != txcode_fast[i]) {
            ok = false;
            break;
        }
    }
    if(ok) {
        printf(KGRN "OK" KNRM "\n");
    } else {
        printf(KRED "FAIL" KNRM "\n");
        overall_ok = false;
    }

    /* Initialise the parity check matrix. */
    printf("* Check parity matrix:                         ");
    size_h = ldpc_codes_size_paritycheck(code);
    h = malloc(size_h);
    ldpc_codes_init_paritycheck(code, h);
    if(crc32((uint8_t*)h, size_h) == test_vector_h[code]) {
        printf(KGRN "OK" KNRM "\n");
    } else {
        printf(KRED "FAIL" KNRM "\n");
        overall_ok = false;
    }

    /* Initialise the sparse parity check matrix. */
    printf("* Check sparse parity matrix:                  ");
    ldpc_codes_size_sparse_paritycheck(code, &size_ci, &size_cs,
                                             &size_vi, &size_vs);
    ci = malloc(size_ci);
    cs = malloc(size_cs);
    vi = malloc(size_vi);
    vs = malloc(size_vs);
    ldpc_codes_init_sparse_paritycheck(code, ci, cs, vi, vs);
    if(   crc32((uint8_t*)ci, size_ci/2) == test_vector_sparse_h[code][0]
       && crc32((uint8_t*)cs, size_cs/2) == test_vector_sparse_h[code][1]
       && crc32((uint8_t*)vi, size_vi/2) == test_vector_sparse_h[code][2]
       && crc32((uint8_t*)vs, size_vs/2) == test_vector_sparse_h[code][3])
    {
        printf(KGRN "OK" KNRM "\n");
    } else {
        printf(KRED "FAIL" KNRM "\n");
        overall_ok = false;
    }

    /* Copy txcode into rxcode, flip the first transmitted bit. */
    rxcode = malloc(n/8);
    for(i=0; i<n/8; i++) {
        rxcode[i] = txcode_fast[i];
    }
    rxcode[0] ^= (1<<7);

    /* Generate RX LLRs, check we can round-trip. */
    printf("* Check round-tripping hard info to LLRs:      ");
    rxcode_llr = malloc(n/8);
    rxllrs = malloc(ldpc_decode_size_llrs(code));
    ldpc_decode_hard_to_llrs(code, rxcode, rxllrs);
    ldpc_decode_llrs_to_hard(code, rxllrs, rxcode_llr);
    ok = true;
    for(i=0; i<n/8; i++) {
        if(rxcode[i] != rxcode_llr[i]) {
            ok = false;
            break;
        }
    }
    if(ok) {
        printf(KGRN "OK" KNRM "\n");
    } else {
        printf(KRED "FAIL" KNRM "\n");
        overall_ok = false;
    }

    /* Bit flip decoder. */
    printf("* Check bit-flipping decoder:                  ");
    rxdata = malloc(ldpc_decode_size_out(code));
    bf_wa = malloc(ldpc_decode_size_bf_wa(code));
    ldpc_decode_bf(code, ci, cs, vi, vs, rxcode, rxdata, bf_wa);
    ok = true;
    for(i=0; i<k/8; i++) {
        if(rxdata[i] != txdata[i]) {
            ok = false;
            break;
        }
    }
    if(ok) {
        printf(KGRN "OK" KNRM "\n");
    } else {
        printf(KRED "FAIL" KNRM "\n");
        overall_ok = false;
    }

    /* Message passing decoder. */
    printf("* Check message-passing decoder:               ");
    mp_wa = malloc(ldpc_decode_size_mp_wa(code));
    ldpc_decode_mp(code, ci, cs, vi, vs, rxllrs, rxdata, mp_wa);
    ok = true;
    for(i=0; i<k/8; i++) {
        if(rxdata[i] != txdata[i]) {
            ok = false;
            break;
        }
    }
    if(ok) {
        printf(KGRN "OK" KNRM "\n");
    } else {
        printf(KRED "FAIL" KNRM "\n");
        overall_ok = false;
    }


    /* Free all that memory we allocated. */
    free(g);
    free(txdata);
    free(txcode_small);
    free(txcode_fast);
    free(h);
    free(ci);
    free(cs);
    free(vi);
    free(vs);
    free(rxcode);
    free(rxcode_llr);
    free(rxllrs);
    free(bf_wa);
    free(rxdata);
    free(mp_wa);

    /* Report time taken to run this test. */
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t_end);
    time_taken = (double)t_end.tv_sec + (double)t_end.tv_nsec * 1e-9 -
                 (double)t_start.tv_sec - (double)t_start.tv_nsec * 1e-9;
    printf("* ----------------------------------------------------------\n");
    printf("* Time taken: " KYEL " %.1fms" KNRM "\n", time_taken * 1000.0f);
    printf("************************************************************\n\n");

    return overall_ok;
}

int main()
{
    enum ldpc_code codes[] = {
        LDPC_CODE_N128_K64,
        LDPC_CODE_N256_K128,
        LDPC_CODE_N512_K256,
        LDPC_CODE_N1280_K1024,
        LDPC_CODE_N1536_K1024,
        LDPC_CODE_N2048_K1024
    };

    size_t i;

    bool ok = true;

    printf("\n");

    for(i=0; i<sizeof(codes)/sizeof(codes[0]); i++) {
        ok &= test_code(codes[i]);
    }

    if(ok) {
        printf(KGRN "All tests passed.\n" KNRM);
        return 0;
    } else {
        printf(KRED "Test failure.\n" KNRM);
        return 1;
    }
}
