#ifndef LDPC_DECODER_H
#define LDPC_DECODER_H
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include "ldpc_codes.h"

/* Decode received input into output using bit flipping algorithm.
 * This algorithm is very quick, uses little memory, and only requires
 * hard information, but is around 1dB less capable at decoding than
 * the message passing algorithm with hard information.
 *
 * (ci, cs, vi, vs) must all have been initialised by
 * ldpc_codes_init_sparse_paritycheck for the appropriate code, except for the
 * non-punctured codes (n=128, n=256, n=512), where vi and vs are unused and
 * may be NULL.
 *
 * input must be n/8 bytes long and each bit is a hard decision.
 * output must be (n+p)/8 bytes long and is written with the decoded codeword,
 * so the user data is present in the first k/8 bytes.
 * working must be n + p bytes long and is used as a scratch working area:
 *
 * Code             Length of working area
 * (128, 64)        128
 * (256, 128)       256
 * (512, 256)       512
 * (1280, 1024)     1408
 * (1536, 1024)     1792
 * (2048, 1024)     2560
 *
 * The required size of the output and the working area are available as
 * LDPC_SIZE_BF_WA(CODE) and LDPC_SIZE_OUT(CODE), or from the
 * ldpc_decode_size_bf_wa(code) and ldpc_decode_size_out(code) functions.
 *
 * Calling this function requires at least 100 bytes of stack.
 *
 * Returns true on decoding success, false otherwise. Even a failed decode
 * may have corrected some of the bit errors, but the result is not a valid
 * codeword.
 */
bool ldpc_decode_bf(enum ldpc_code code,
                    uint16_t* ci, uint16_t* cs, uint16_t* vi, uint16_t* vs,
                    const uint8_t* input, uint8_t* output, uint8_t* working);

/* Find the size (in BYTES) required for the working area of the BF algorithm.
 * This is the same as described in the associated comment,
 * and is n + p.
 * The same information is available statically from the LDPC_SIZE_BF_WA macro
 * in ldpc_sizes.h.
 */
size_t ldpc_decode_size_bf_wa(enum ldpc_code code);

/* Decode LLRs into data using message passing algorithm.
 * This algorithm is slower and ideally requires soft information,
 * but decodes very close to optimal. If you don't have soft information
 * but do have the channel BER, you can use ldpc_decode_hard_to_llrs_ber to
 * go from hard information (bytes from a receiver) to soft information.
 * If you don't have that, you can use ldpc_decode_hard_to_llrs to generate
 * arbitrary LLRs from the hard information.
 *
 * (ci, cs, vi, vs) must all have been initialised by
 * ldpc_codes_init_sparse_paritycheck for the appropriate code.
 * llrs must be n floats long, where positive numbers are more likely to be 0.
 * output must be (n+p)/8 bytes long, of which the first k/8 bytes will be set
 *     to the original transmitted message (then followed by parity bits).
 * working must be 2*s long:
 *
 * Code             Length of output        Length of working area
 * (128, 64)        16                      1024
 * (256, 128)       32                      2048
 * (512, 256)       64                      4096
 * (1280, 1024)     176                     9984
 * (1536, 1024)     224                     11776
 * (2048, 1024)     320                     15360
 *
 * The byte sizes are statically available from LDPC_SIZE_MP_WA(CODE) and
 * LDPC_SIZE_OUT(CODE) macros in ldpc_sizes.h, or from
 * ldpc_decode_size_mp_wa(code) and ldpc_decode_size_out(code).
 *
 * Calling this function uses at least 100 bytes of stack.
 *
 * Returns true on decoding success, false otherwise.
 */
bool ldpc_decode_mp(enum ldpc_code code,
                    uint16_t* ci, uint16_t* cs,
                    uint16_t* vi, uint16_t* vs,
                    const float* llrs, uint8_t* output, float* working);

/* Find the size (in BYTES) required for the working area of the MP algorithm.
 * This is the same as described in the associated comment, and is
 * 2*s*sizeof(float). The same information is available statically from the
 * LDPC_SIZE_MP_WA macro in ldpc_sizes.h.
 */
size_t ldpc_decode_size_mp_wa(enum ldpc_code code);

/* Find the size (in BYTES) required for the output of the decoders.
 * This is the same as described in the associated comments, and is
 * (n+p)/8. The same information is available statically from the
 * LDPC_SIZE_OUT macro in ldpc_sizes.h.
 */
size_t ldpc_decode_size_out(enum ldpc_code code);

/* Create approximate LLRs using just the channel BER and the received data.
 * Can be used to feed the message passing algorithm soft-ish information.
 * input must be n/8 bytes long, llrs must be n floats long.
 */
void ldpc_decode_hard_to_llrs_ber(enum ldpc_code code, const uint8_t* input,
                                  float* llrs, float ber);

/* Create hard LLRs from hard received data.
 * input must be n/8 bytes long, llrs must be n floats long.
 */
void ldpc_decode_hard_to_llrs(enum ldpc_code code, const uint8_t* input,
                              float* llrs);

/* Create hard information from received LLRs.
 * llrs must be n floats long, output must be n/8 bytes long.
 */
void ldpc_decode_llrs_to_hard(enum ldpc_code code, const float* llrs,
                              uint8_t* output);

/* Find the size (in BYTES) required to store the LLRs for the given code.
 * This is sizeof(float)*n.
 * The same information is available statically from the LDPC_SIZE_LLRS macro
 * in ldpc_sizes.h.
 */
size_t ldpc_decode_size_llrs(enum ldpc_code code);

#endif
