#ifndef LDPC_CODE_SIZES_H
#define LDPC_CODE_SIZES_H

/* Look up N from code name */
#define LDPC_SIZE_N_LDPC_CODE_NONE          (0)
#define LDPC_SIZE_N_LDPC_CODE_N128_K64      (128)
#define LDPC_SIZE_N_LDPC_CODE_N256_K128     (256)
#define LDPC_SIZE_N_LDPC_CODE_N512_K256     (512)
#define LDPC_SIZE_N_LDPC_CODE_N1280_K1024   (1280)
#define LDPC_SIZE_N_LDPC_CODE_N1536_K1024   (1536)
#define LDPC_SIZE_N_LDPC_CODE_N2048_K1024   (2048)
#define LDPC_SIZE_N(CODE) LDPC_SIZE_N_(CODE)
#define LDPC_SIZE_N_(CODE) LDPC_SIZE_N_##CODE

/* Look up K from code name */
#define LDPC_SIZE_K_LDPC_CODE_NONE          (0)
#define LDPC_SIZE_K_LDPC_CODE_N128_K64      (64)
#define LDPC_SIZE_K_LDPC_CODE_N256_K128     (128)
#define LDPC_SIZE_K_LDPC_CODE_N512_K256     (256)
#define LDPC_SIZE_K_LDPC_CODE_N1280_K1024   (1024)
#define LDPC_SIZE_K_LDPC_CODE_N1536_K1024   (1024)
#define LDPC_SIZE_K_LDPC_CODE_N2048_K1024   (1024)
#define LDPC_SIZE_K(CODE) LDPC_SIZE_K_(CODE)
#define LDPC_SIZE_K_(CODE) LDPC_SIZE_K_##CODE

/* Look up P (number of punctured checks) from code names */
#define LDPC_SIZE_P_LDPC_CODE_NONE          (0)
#define LDPC_SIZE_P_LDPC_CODE_N128_K64      (0)
#define LDPC_SIZE_P_LDPC_CODE_N256_K128     (0)
#define LDPC_SIZE_P_LDPC_CODE_N512_K256     (0)
#define LDPC_SIZE_P_LDPC_CODE_N1280_K1024   (128)
#define LDPC_SIZE_P_LDPC_CODE_N1536_K1024   (256)
#define LDPC_SIZE_P_LDPC_CODE_N2048_K1024   (512)
#define LDPC_SIZE_P(CODE) LDPC_SIZE_P_(CODE)
#define LDPC_SIZE_P_(CODE) LDPC_SIZE_P_##CODE

/* Look up S (sum of H) from code names */
#define LDPC_SIZE_S_LDPC_CODE_NONE          (0)
#define LDPC_SIZE_S_LDPC_CODE_N128_K64      (512)
#define LDPC_SIZE_S_LDPC_CODE_N256_K128     (1024)
#define LDPC_SIZE_S_LDPC_CODE_N512_K256     (2048)
#define LDPC_SIZE_S_LDPC_CODE_N1280_K1024   (4992)
#define LDPC_SIZE_S_LDPC_CODE_N1536_K1024   (5888)
#define LDPC_SIZE_S_LDPC_CODE_N2048_K1024   (7680)
#define LDPC_SIZE_S(CODE) LDPC_SIZE_S_(CODE)
#define LDPC_SIZE_S_(CODE) LDPC_SIZE_S_##CODE

/* Parity check matrix size in bytes (= ((n+p) * (n-k+p)) / 8) */
#define LDPC_SIZE_H(CODE)                   \
    ((LDPC_SIZE_N(CODE) + LDPC_SIZE_P(CODE)) * \
     (LDPC_SIZE_N(CODE) + LDPC_SIZE_P(CODE) - LDPC_SIZE_K(CODE)) / 8)

/* Sparse parity CI and VI size in bytes (= s uint16_t) */
#define LDPC_SIZE_CI_VI(CODE)               \
    (sizeof(uint16_t) * LDPC_SIZE_S(CODE))
#define LDPC_SIZE_CI(CODE) LDPC_SIZE_CI_VI(CODE)
#define LDPC_SIZE_VI(CODE) LDPC_SIZE_CI_VI(CODE)

/* Sparse parity CS size in bytes (= n - k + p + 1 uint16_t) */
#define LDPC_SIZE_CS(CODE)                  \
    (sizeof(uint16_t) *                     \
     (LDPC_SIZE_N(CODE) - LDPC_SIZE_K(CODE) + LDPC_SIZE_P(CODE) + 1))

/* Sparse parity VS size in bytes (= n + p + 1 uint16_t)*/
#define LDPC_SIZE_VS(CODE)                  \
    (sizeof(uint16_t) *                     \
     (LDPC_SIZE_N(CODE) + LDPC_SIZE_P(CODE) + 1))

/* Sparse parity overall size in bytes */
#define LDPC_SIZE_SPARSE_H(CODE)            \
    (LDPC_SIZE_CI(CODE) + LDPC_SIZE_CS(CODE)\
     + LDPC_SIZE_VI(CODE) + LDPC_SIZE_VS(CODE))

/* Expanded generator matrix sizes in bytes, for fast encoder.
 * There are k*(n-k)/8 parity bits (excluding systematic identity matrix).
 */
#define LDPC_SIZE_G(CODE)                   \
    (LDPC_SIZE_K(CODE) * (LDPC_SIZE_N(CODE) - LDPC_SIZE_K(CODE)) / 8)

/* Required working area for the bitflipping algorithm, in bytes. (=n/8 + n) */
#define LDPC_SIZE_BF_WA(CODE)               \
    ((LDPC_SIZE_N(CODE) * 9) / 8)

/* Size of LLRs for message passing decoder in bytes (=n floats) */
#define LDPC_SIZE_MP_LLRS(CODE)             \
    (LDPC_SIZE_N(CODE) * sizeof(float))

/* Size of working area message passing decoder in bytes (=2s floats) */
#define LDPC_SIZE_MP_WA(CODE)               \
    (2 * sizeof(float) * LDPC_SIZE_S(CODE))

/* Size of output for message passing decoder in bytes (=(n+p)/8) */
#define LDPC_SIZE_MP_OUT(CODE)              \
    ((LDPC_SIZE_N(CODE) + LDPC_SIZE_P(CODE)) / 8)

/* Add up the required size for each type of encoder, including the size
 * needed to store the resulting output.
 */
#define LDPC_SIZE_TX_SMALL(CODE)    (LDPC_SIZE_N(CODE) / 8)
#define LDPC_SIZE_TX_FAST(CODE)     (LDPC_SIZE_G(CODE) + LDPC_SIZE_N(CODE) / 8)

/* Add up the required size for each type of decoder, including the size
 * needed to store the input and resulting output.
 */
#define LDPC_SIZE_RX_BF(CODE)       \
    (LDPC_SIZE_SPARSE_H(CODE) + LDPC_SIZE_BF_WA(CODE) + LDPC_SIZE_N(CODE)/8 + \
     LDPC_SIZE_K(CODE)/8)
#define LDPC_SIZE_RX_MP(CODE)       \
    (LDPC_SIZE_SPARSE_H(CODE) + LDPC_SIZE_N(CODE)/8 + LDPC_SIZE_MP_LLRS(CODE) \
     + LDPC_SIZE_MP_WA(CODE) + LDPC_SIZE_MP_OUT(CODE))

/* Find the size needed, parameterised over the encoder, decoder,
 * TX code, and RX code, to store all the relevant expanded codes, working
 * areas, and the outputs of the encoders and decoders.
 */
#define LDPC_SIZE(ENCODER, TXCODE, DECODER, RXCODE)     \
    (LDPC_SIZE_TX_##ENCODER(TXCODE) + LDPC_SIZE_RX_##DECODER(RXCODE))

#endif
