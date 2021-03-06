#ifndef LDPC_CODES_H
#define LDPC_CODES_H
#include <stddef.h>
#include <stdint.h>

/* Available LDPC codes.
 * n is the block length (bits transmitted over the air),
 * k is the data length (number of user bits encoded).
 * LDPC_CODE_NONE causes functions to all no-op.
 */
enum ldpc_code {
    LDPC_CODE_NONE,
    LDPC_CODE_N128_K64,
    LDPC_CODE_N256_K128,
    LDPC_CODE_N512_K256,
    LDPC_CODE_N1280_K1024,
    LDPC_CODE_N1536_K1024,
    LDPC_CODE_N2048_K1024,
};

/* Get parameters corresponding to a given code.
 * n: code length (block length; number of transmitted bits)
 * k: code dimension (information length; number of user bits)
 * p: punctured checks (number of parity bits not transmitted)
 * m: sub-matrix size (used in code definition)
 * b: circulant block size (used in generator construction)
 * s: sum of the H matrix, i.e. total number of parity check edges
 * Specify NULL if you don't care about a certain parameter.
 */
void ldpc_codes_get_params(enum ldpc_code code,
                           int* n, int* k, int* p, int* m, int* b, int* s);

/* Fill h with the appropriate parity check matrix,
 * densely packed (32 columns per word).
 * Required size of h is (n+p)*(n-k+p)/8 bytes:
 *  (128, 64)       1024 bytes      256 words       uint32_t[64][4]
 *  (256, 128)      4096 bytes      1024 words      uint32_t[128][8]
 *  (512, 256)      16384 bytes     4096 words      uint32_t[256][16]
 *  (1280, 1024)    67584 bytes     16896 words     uint32_t[384][44]
 *  (1536, 1024)    172032 bytes    43008 words     uint32_t[768][56]
 *  (2048, 1024)    491520 bytes    122880 words    uint32_t[1536][80]
 * Note the larger codes are punctured so the parity check matrix may be larger
 * than the usual (n-k, n) size.
 */
void ldpc_codes_init_paritycheck(enum ldpc_code code, uint32_t* h);

/* Find the required size (in BYTES) for a given parity check matrix.
 * This reflects the information in the comment for init_paritycheck.
 * The same information is available statically via the LDPC_SIZE_H
 * macro in ldpc_sizes.h.
 * Returns (n+p)*(n-k+p)/8 for the given code.
 */
size_t ldpc_codes_size_paritycheck(enum ldpc_code code);

/* Fill sparse representations of parity check matrix.
 * This representation has two 1d-lists, ci and vi, one representing the
 * non-zero indices along each row (check nodes), and the other along each
 * column (variable nodes). They allow iterating through the parity matrix
 * connections either from row to column (check to variable node) or vice
 * versa very efficiently.
 *
 * To index into the lists, each has a list of starting points cs and vs.
 * The length of cs is equal to the number of parity check equations (before
 * puncturing) plus 1, while the length of vs is equal to the number of
 * variable nodes plus 1 (or the codeword length, n, + 1).
 * The lengths of each sub-list (the degree of that node) is implicit from the
 * starting point of the next sub-list. The final entry in cs and vs is set to
 * one after the end of ci/vi.
 *
 * All of these must have been allocated before calling this function,
 * with the following lengths (multiply by 2 for size in bytes).
 * Code             ci, vi        cs          vs
 * (n, k, p, s)        s     n-k+p+1       n+p+1
 * (128, 64)         512          65         129
 * (256, 128)       1024         129         257
 * (512, 256)       2048         257         513
 * (1280, 1024)     4992         385        1409
 * (1536, 1024)     5888         769        1793
 * (2048, 1024)     7680        1537        2561
 *
 * If you only want to use the bitflipping decoder, you can initialise only
 * the row results, using ldpc_codes_init_sparse_paritycheck_rows.
 *
 */
void ldpc_codes_init_sparse_paritycheck(enum ldpc_code code,
                                        uint16_t* ci, uint16_t* cs,
                                        uint16_t* vi, uint16_t* vs);
void ldpc_codes_init_sparse_paritycheck_rows(enum ldpc_code code,
                                             uint16_t* ci, uint16_t* cs);

/* Find the sizes (in BYTES) for ci, cs, vi, and vs, as used in
 * init_sparse_paritycheck. These sizes reflect the information in the comment.
 * The same information is available statically via the LDPC_SIZE_CI,
 * LDPC_SIZE_CS, LDPC_SIZE_VI, and LDPC_SIZE_VS macros in ldpc_sizes.h.
 * For their sum, you can also use LDPC_SIZE_SPARSE_H.
 * Returns sizeof(uint16_t) * (s, n-k+p+1, n+p+1) for the given code.
 */
void ldpc_codes_size_sparse_paritycheck(enum ldpc_code code,
                                        size_t* ci, size_t* cs,
                                        size_t* vi, size_t* vs);

/* Gets a pointer to the relevant constants for compact generator matrices.
 * Also sets n and k (code size) and b (circulant block size).
 * The returned pointer is to a uint32_t[k/b][(n-k)/32],
 * each row has (n-k) bits ((n-k)/32 words),
 * there are k/b rows, which each represent b rows of the actual generator
 * matrix, which has k rows.
 */
const uint32_t * ldpc_codes_get_compact_generator(enum ldpc_code code,
                                                  int* n, int* k, int* b);

/* Initialise a generator matrix expanded from the compact circulant form.
 * This allows quicker encoding at the cost of more memory usage.
 * Note this will only initialise the parity part of G, and not the identity
 * matrix, since all supported codes are systematic. This matches what's
 * expected by the non-compact encoder.
 * Required size of g is k*(n-k)/8 bytes:
 *  (128, 64)       512 bytes       128 words       uint32_t[64][2]
 *  (256, 128)      2048 bytes      512 words       uint32_t[128][4]
 *  (512, 256)      8192 bytes      2048 words      uint32_t[256][8]
 *  (1280, 1024)    32768 bytes     8192 words      uint32_t[1024][8]
 *  (1536, 1024)    65536 bytes     16384 words     uint32_t[1024][16]
 *  (2048, 1024)    131072 bytes    32768 words     uint32_t[1024][32]
 */
void ldpc_codes_init_generator(enum ldpc_code code, uint32_t* g);

/* Find the size (in BYTES) required for the generator matrix g.
 * This is the same as the information in the comment for init_generator,
 * and is equal to k*(n-k)/8.
 * The same information is available statically from the LDPC_SIZE_G macro
 * in ldpc_sizes.h.
 */
size_t ldpc_codes_size_generator(enum ldpc_code code);

#endif
