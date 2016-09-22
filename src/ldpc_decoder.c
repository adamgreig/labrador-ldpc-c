/*
 * LDPC decoding functions
 * Copyright 2016 Adam Greig
 */

#include <math.h>
#include <float.h>
#include <string.h>
#include "ldpc_decoder.h"

#include <stdio.h>

static inline float sign(float x) {
    return (x>0.0f) - (x<0.0f);
}

bool ldpc_decode_bf(enum ldpc_code code,
                    uint16_t* ci, uint16_t* cs,
                    const uint8_t* input, uint8_t* output, uint8_t* working)
{
    int n, k, p, i, a, i_a, max_violations, iters;
    const int max_iters = 20;

    uint8_t * codeword, * violations;

    if(code == LDPC_CODE_NONE) {
        return false;
    }

    ldpc_codes_get_params(code, &n, &k, &p, NULL, NULL, NULL);

    /* Split up the working area */
    codeword = working;
    violations = working + (n + p)/8;

    /* Copy input to codeword space */
    memcpy(codeword, input, n/8);
    /* Set all the punctured bits to 0 */
    memset(codeword + n/8, 0, p/8);

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
                parity += (codeword[a/8] >> (7-(a%8))) & 1;
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
            memcpy(output, codeword, k/8);
            return true;
        } else {
            /* Otherwise flip all bits that had the maximum violations */
            for(a=0; a<(n+p); a++) {
                if(violations[a] == max_violations) {
                    codeword[a/8] ^= 1<<(7-(a%8));
                }
            }
        }
    }

    /* If we didn't successfully decode in max_iters, fail here. */
    return false;
}

size_t ldpc_decode_size_bf_wa(enum ldpc_code code)
{
    int n;
    ldpc_codes_get_params(code, &n, NULL, NULL, NULL, NULL, NULL);
    return (9*n)/8;
}

bool ldpc_decode_mp(enum ldpc_code code,
                    uint16_t* ci, uint16_t* cs,
                    uint16_t* vi, uint16_t* vs,
                    const float* llrs, uint8_t* output, float* working)
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
            return true;
        }
    }

    /* If we got here, we ran out of iterations :( */
    return false;
}

size_t ldpc_decode_size_mp_wa(enum ldpc_code code)
{
    int s;
    ldpc_codes_get_params(code, NULL, NULL, NULL, NULL, NULL, &s);
    return 2*s*sizeof(float);
}

size_t ldpc_decode_size_mp_out(enum ldpc_code code)
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
