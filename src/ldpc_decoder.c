/*
 * LDPC decoding functions
 * Copyright 2016 Adam Greig
 */

#include <math.h>
#include <float.h>
#include <string.h>
#include "ldpc_decoder.h"

static inline float sign(float x) {
    return (x>0.0f) - (x<0.0f);
}

/*
 * Erasure decoding algorithm to pre-process punctured codes before using the
 * bit flipping algorithm.
 * Since bit flipping can't handle erasures, we instead try and decode them.
 * The basic idea is:
 *  * For each erased bit `a`:
 *      * For each check `i` that `a` is associated with:
 *          * If `a` is the only erasure that `i` is associated with,
 *            then compute the parity of `i`, and cast a vote for the
 *            value of `a` that would give even parity
 *          * Otherwise ignore `i`
 *      * If there is a majority vote, set `a` to the winning value,
 *        and mark it no longer erased. Otherwise, leave it erased.
 * This is based on the paper
 * Novel multi-Gbps bit-flipping decoders for punctured LDPC codes,
 * by Archonta, Kanistras and Paliouras, MOCAST 2016
 *
 * ci, cs, vi, and vs must all have been initialised appropriately.
 * output must be (n+p)/8 long, with the first n/8 bytes already set to the
 * received hard information, and the punctured bits in it will be updated.
 * working must be (n+p) long. (we could use p/8 but since the bit flipping
 *                              decoder needs n+p working area, we lose
 *                              nothing and gain speed by using one byte
 *                              per bit here).
 */
void ldpc_decode_erasures(enum ldpc_code code,
                          uint16_t* ci, uint16_t* cs,
                          uint16_t* vi, uint16_t* vs,
                          uint8_t* output,
                          uint8_t* working,
                          uint16_t* iters_run)
{
    int n, k, p, i, a, b, i_b, a_i, iters, bits_fixed=0;
    uint8_t *erasures = working;
    uint8_t parity;
    int8_t votes;
    bool only_one_erasure;
    const int max_iters = 16;

    /* Check ci, cs, vi, and vs have actually been initialised. */
    if(ci == NULL || cs == NULL || vi == NULL || vs == NULL) {
        return;
    }

    ldpc_codes_get_params(code, &n, &k, &p, NULL, NULL, NULL);

    /* Set all the punctured bits to be erased and have value 0 (arbitrary) */
    memset(erasures, 0, n);
    memset(erasures+n, 1, p);
    memset(output + n/8, 0x00, p/8);

    /* Run until either we run out of iterations or all bits get fixed. */
    for(iters=0; iters<max_iters && bits_fixed<p; iters++) {
        /* For each punctured bit */
        for(a=n; a<n+p; a++) {
            /* Skip bits that are no longer marked as erased */
            if(!erasures[a]) {
                continue;
            }

            /* Track votes for 0 (negative) or 1 (positive). */
            votes = 0;

            /* For each check this bit is associated with */
            for(a_i=vs[a]; a_i<vs[a+1]; a_i++) {
                i = vi[a_i];
                parity = 0;

                /* See what the check parity is, and quit without voting if
                 * the parity check has another erasure already.
                 */
                only_one_erasure = true;
                for(i_b=cs[i]; i_b<cs[i+1]; i_b++) {
                    b = ci[i_b];

                    /* Skip the punctured bit we're looking at */
                    if(a == b) {
                        continue;
                    }

                    /* If we see another erasure, stop */
                    if(erasures[b]) {
                        only_one_erasure = false;
                        break;
                    }

                    /* Add up the parity for this check */
                    parity += (output[b/8] >> (7-(b%8))) & 1;

                }

                /* Cast a vote if we only have one erasure.
                 * If all the bits except this one add to odd parity,
                 * we vote for this one to be 1 (to get even parity). And vv.
                 */
                if(only_one_erasure) {
                    if((parity & 1) == 1) {
                        votes++;
                    } else {
                        votes--;
                    }
                }
            }

            /* If we had a majority vote one way or the other, great!
             * Set ourselves to the majority vote value and clear our erasure.
             */
            if(votes != 0) {
                erasures[a] = 0;
                bits_fixed++;

                if(votes > 0) {
                    output[a/8] |=  (1<<(7-(a%8)));
                } else {
                    output[a/8] &= ~(1<<(7-(a%8)));
                }
            }
        }
    }

    *iters_run = iters;
}

bool ldpc_decode_bf(enum ldpc_code code,
                    uint16_t* ci, uint16_t* cs, uint16_t* vi, uint16_t* vs,
                    const uint8_t* input, uint8_t* output, uint8_t* working,
                    uint16_t* iters_run)
{
    int n, k, p, i, a, i_a, max_violations, iters;
    const int max_iters = 20;

    *iters_run = 0;

    /* Use working area to store parity-violation counts per bit */
    uint8_t* violations = working;

    if(code == LDPC_CODE_NONE) {
        return false;
    }

    ldpc_codes_get_params(code, &n, &k, &p, NULL, NULL, NULL);

    /* Copy input to codeword space */
    memcpy(output, input, n/8);

    /* If the code is punctured, first try and fix erasures */
    if(p > 0) {
        ldpc_decode_erasures(code, ci, cs, vi, vs, output, working, iters_run);
    }

    /* Run the bit flipping algorithm */
    for(iters=0; iters<max_iters; iters++) {
        /* Clear the violations counters */
        memset(violations, 0, n+p);
        max_violations = 0;

        /* For each parity check, work out the parity sum */
        for(i=0; i<(n-k+p); i++) {
            int parity = 0;
            /* For each bit, update the parity sum */
            for(i_a=cs[i]; i_a<cs[i+1]; i_a++) {
                a = ci[i_a];
                parity += (output[a/8] >> (7-(a%8))) & 1;
            }

            /* If the check has odd parity, add one violation to each
             * variable node involved in the check */
            if((parity & 1) != 0) {
                for(i_a=cs[i]; i_a<cs[i+1]; i_a++) {
                    a = ci[i_a];
                    violations[a]++;
                }
            }
        }

        /* Find the maximum number of violations */
        for(a=0; a<(n+p); a++) {
            if(violations[a] > max_violations) {
                max_violations = violations[a];
            }
        }

        if(max_violations == 0) {
            /* No violations means valid codeword */
            *iters_run += iters;
            return true;
        } else {
            /* Otherwise flip all bits that had the maximum violations */
            for(a=0; a<(n+p); a++) {
                if(violations[a] == max_violations) {
                    output[a/8] ^= 1<<(7-(a%8));
                }
            }
        }
    }

    /* If we didn't successfully decode after max_iters, fail here. */
    *iters_run += iters;
    return false;
}

size_t ldpc_decode_size_bf_wa(enum ldpc_code code)
{
    int n, p;
    ldpc_codes_get_params(code, &n, NULL, &p, NULL, NULL, NULL);
    return n+p;
}

bool ldpc_decode_mp(enum ldpc_code code,
                    uint16_t* ci, uint16_t* cs,
                    uint16_t* vi, uint16_t* vs,
                    const float* llrs, uint8_t* output, float* working,
                    uint16_t* iters_run)
{
    /* Code parameters */
    int n, k, p, s;
    /* Variable and check node numbers and lookup indices */
    int a, b, i_a, i_b, j_b, i, j, a_i, a_j, b_j;
    /* Iteration counting */
    int iters;
    const int max_iters = 20;
    /* Message storage */
    float * u, * v;
    /* Output computation/checking */
    float llr_a;
    int parity;
    bool parity_ok;
    float prev_v_ai;
    float fabs_v_bj;
    float sgnprod;
    float minacc;

    *iters_run = 0;

    if(code == LDPC_CODE_NONE) {
        return false;
    }

    ldpc_codes_get_params(code, &n, &k, &p, NULL, NULL, &s);

    /* Split up working area.
     * u(i->a) holds messages from checks to variables,
     * v(a->i) holds messages from variables to checks.
     */
    u = working;
    v = working + s;

    /* Initialise u(i->a) to 0 */
    memset(u, 0, sizeof(float) * s);

    /* Initialise v(a->i) to 0 */
    memset(v, 0, sizeof(float) * s);

    /* Run the message passing algorithm */
    for(iters=0; iters<max_iters; iters++) {
        /* Keep track of whether the overall parity is met.
         * Will be set to false as soon as one invalid parity equation is
         * encountered.
         */
        parity_ok = true;
        memset(output, 0, (n+p)/8);

        /* Update variable nodes' messages to check nodes.
         * For each variable node, for each check node connected to it,
         * initialise this message v(a->i) to the LLR (or 0 for punctured bits)
         * and then add on all of the incoming messages not from the current
         * check node.
         * Additionally we sum all u(i->a) for this a to marginalise this
         * variable node and see if the hard decoding gives a valid codeword,
         * which is our signal to stop iterating.
         */
        for(a=0; a<n+p; a++) {
            llr_a = a < n ? llrs[a] : 0.0f;

            /* For each check node i connected to variable node a */
            for(a_i=vs[a]; a_i<vs[a+1]; a_i++) {
                i = vi[a_i];
                prev_v_ai = v[a_i];
                v[a_i] = a < n ? llrs[a] : 0.0f;

                /* For each check node j connected to variable node a */
                for(a_j=vs[a]; a_j<vs[a+1]; a_j++) {
                    j = vi[a_j];

                    /* We need to find where the incoming messages u(j->a)
                     * are stored in u. This means going through every variable
                     * node connected to check node j, and seeing if it's equal
                     * to a, and if so, using that message.
                     * This loop could be replaced by another index table the
                     * same size as ci, which might save time if this section
                     * proves to be slow.
                     */
                    for(j_b=cs[j]; j_b<cs[j+1]; j_b++) {
                        b = ci[j_b];
                        if(a == b) {
                            /* Sum up just the incoming messages not from i
                             * for v(a->i)
                             */
                            if(j != i) {
                                v[a_i] += u[j_b];
                            }

                            /* We sum up all the incoming messages for llr_a */
                            llr_a += u[j_b];

                            /* As soon as we've found our a, we can stop
                             * looking.
                             */
                            break;
                        }
                    }
                }

                /* Our min sum correction trick will be to zero any messages
                 * that have changed sign, as per Savin 2009:
                 * http://arxiv.org/abs/0803.1090v2
                 */
                if(prev_v_ai != 0.0f && sign(v[a_i]) != sign(prev_v_ai)) {
                    v[a_i] = 0.0f;
                }
            }

            /* Hard decode llr_a to determine this output bit */
            output[a/8] |= (llr_a <= 0.0f) << (7 - (a%8));

        }

        /* Update check nodes' messages to variable nodes.
         * For each check node, for each variable node connected to it,
         * initialise the message u(i->a) to FLT_MAX and then find the
         * minimum of all incoming messages as well as the product of all their
         * signs. Additionally we use this loop to keep track of the parity
         * sum for this check node under hard decoding, and use that to see if
         * the overall message has been decoded OK.
         */
        for(i=0; i<n-k+p; i++) {
            parity = 0;

            /* For each variable node a connected to check node i */
            for(i_a=cs[i]; i_a<cs[i+1]; i_a++) {
                a = ci[i_a];
                sgnprod = 1.0f;
                minacc = FLT_MAX;

                /* For each variable node b connected to check node i */
                for(i_b=cs[i]; i_b<cs[i+1]; i_b++) {
                    b = ci[i_b];

                    /* Don't process the message from the variable node we're
                     * currently updating.
                     */
                    if(b == a) {
                        continue;
                    }

                    /* We need to find where the incoming messages v(b->i)
                     * are stored in v. As with the u(j->a) messages, we need
                     * to go through each check node j associated with variable
                     * node b, and if j==i we can use the message.
                     * This loop could also be replaced by another index table
                     * the same size as vi, which might save time.
                     */
                    for(b_j=vs[b]; b_j<vs[b+1]; b_j++) {
                        j = vi[b_j];
                        if(i == j) {
                            sgnprod *= sign(v[b_j]);
                            fabs_v_bj = fabs(v[b_j]);
                            if(fabs_v_bj < minacc) {
                                minacc = fabs_v_bj;
                            }

                            /* As soon as we find the node we're looking for,
                             * we can stop looking.
                             */
                            break;
                        }
                    }
                }
                u[i_a] = sgnprod * minacc;

                /* Work out this node's parity */
                parity += (output[a/8] >> (7 - (a%8))) & 1;
            }

            /* Odd parity is bad parity */
            if((parity & 1) == 1) {
                parity_ok = false;
            }
        }

        /* If every parity check was satisfied, we're done. */
        if(parity_ok) {
            *iters_run = iters;
            return true;
        }
    }

    /* If we got here, we ran out of iterations :( */
    *iters_run = iters;
    return false;
}

size_t ldpc_decode_size_mp_wa(enum ldpc_code code)
{
    int s;
    ldpc_codes_get_params(code, NULL, NULL, NULL, NULL, NULL, &s);
    return 2*s*sizeof(float);
}

size_t ldpc_decode_size_out(enum ldpc_code code)
{
    int n, p;
    ldpc_codes_get_params(code, &n, NULL, &p, NULL, NULL, NULL);
    return (n+p)/8;
}

void ldpc_decode_hard_to_llrs_ber(enum ldpc_code code, const uint8_t* input,
                                  float* llrs, float ber)
{
    int i, n=0;
    float logber;

    if(code == LDPC_CODE_NONE) {
        return;
    }

    ldpc_codes_get_params(code, &n, NULL, NULL, NULL, NULL, NULL);
    logber = logf(ber);
    for(i=0; i<n; i++) {
        if((input[i/8] >> (7-(i%8))) & 1) {
            llrs[i] = logber;
        } else {
            llrs[i] = -logber;
        }
    }
}

void ldpc_decode_hard_to_llrs(enum ldpc_code code, const uint8_t* input,
                           float* llrs)
{
    ldpc_decode_hard_to_llrs_ber(code, input, llrs, 0.05f);
}

void ldpc_decode_llrs_to_hard(enum ldpc_code code, const float* llrs,
                              uint8_t* output)
{
    int i, n=0;

    if(code == LDPC_CODE_NONE) {
        return;
    }

    ldpc_codes_get_params(code, &n, NULL, NULL, NULL, NULL, NULL);

    memset(output, 0, n/8);

    for(i=0; i<n; i++) {
        output[i/8] |= (llrs[i] <= 0.0f) << (7 - (i%8));
    }
}

size_t ldpc_decode_size_llrs(enum ldpc_code code)
{
    int n;
    ldpc_codes_get_params(code, &n, NULL, NULL, NULL, NULL, NULL);
    return sizeof(float) * n;
}
