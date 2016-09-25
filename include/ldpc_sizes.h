#ifndef LDPC_CODE_SIZES_H
#define LDPC_CODE_SIZES_H

/* Look up N from code name */
#define LDPC_PARAM_N_LDPC_CODE_NONE          (0)
#define LDPC_PARAM_N_LDPC_CODE_N128_K64      (128)
#define LDPC_PARAM_N_LDPC_CODE_N256_K128     (256)
#define LDPC_PARAM_N_LDPC_CODE_N512_K256     (512)
#define LDPC_PARAM_N_LDPC_CODE_N1280_K1024   (1280)
#define LDPC_PARAM_N_LDPC_CODE_N1536_K1024   (1536)
#define LDPC_PARAM_N_LDPC_CODE_N2048_K1024   (2048)
#define LDPC_PARAM_N(CODE) LDPC_PARAM_N_(CODE)
#define LDPC_PARAM_N_(CODE) LDPC_PARAM_N_##CODE

/* Look up K from code name */
#define LDPC_PARAM_K_LDPC_CODE_NONE          (0)
#define LDPC_PARAM_K_LDPC_CODE_N128_K64      (64)
#define LDPC_PARAM_K_LDPC_CODE_N256_K128     (128)
#define LDPC_PARAM_K_LDPC_CODE_N512_K256     (256)
#define LDPC_PARAM_K_LDPC_CODE_N1280_K1024   (1024)
#define LDPC_PARAM_K_LDPC_CODE_N1536_K1024   (1024)
#define LDPC_PARAM_K_LDPC_CODE_N2048_K1024   (1024)
#define LDPC_PARAM_K(CODE) LDPC_PARAM_K_(CODE)
#define LDPC_PARAM_K_(CODE) LDPC_PARAM_K_##CODE

/* Look up P (number of punctured checks) from code names */
#define LDPC_PARAM_P_LDPC_CODE_NONE          (0)
#define LDPC_PARAM_P_LDPC_CODE_N128_K64      (0)
#define LDPC_PARAM_P_LDPC_CODE_N256_K128     (0)
#define LDPC_PARAM_P_LDPC_CODE_N512_K256     (0)
#define LDPC_PARAM_P_LDPC_CODE_N1280_K1024   (128)
#define LDPC_PARAM_P_LDPC_CODE_N1536_K1024   (256)
#define LDPC_PARAM_P_LDPC_CODE_N2048_K1024   (512)
#define LDPC_PARAM_P(CODE) LDPC_PARAM_P_(CODE)
#define LDPC_PARAM_P_(CODE) LDPC_PARAM_P_##CODE

/* Look up S (sum of H) from code names */
#define LDPC_PARAM_S_LDPC_CODE_NONE          (0)
#define LDPC_PARAM_S_LDPC_CODE_N128_K64      (512)
#define LDPC_PARAM_S_LDPC_CODE_N256_K128     (1024)
#define LDPC_PARAM_S_LDPC_CODE_N512_K256     (2048)
#define LDPC_PARAM_S_LDPC_CODE_N1280_K1024   (4992)
#define LDPC_PARAM_S_LDPC_CODE_N1536_K1024   (5888)
#define LDPC_PARAM_S_LDPC_CODE_N2048_K1024   (7680)
#define LDPC_PARAM_S(CODE) LDPC_PARAM_S_(CODE)
#define LDPC_PARAM_S_(CODE) LDPC_PARAM_S_##CODE

/* Parity check matrix size in bytes (= ((n+p) * (n-k+p)) / 8) */
#define LDPC_SIZE_H(CODE)                   \
    ((LDPC_PARAM_N(CODE) + LDPC_PARAM_P(CODE)) * \
     (LDPC_PARAM_N(CODE) + LDPC_PARAM_P(CODE) - LDPC_PARAM_K(CODE)) / 8)

/* Parity check matrix length in uint32_t */
#define LDPC_LENGTH_H(CODE) (LDPC_SIZE_H(CODE) / 4)

/* Sparse parity CI and VI length (in uint16_t) and size (in bytes) */
#define LDPC_LENGTH_CI(CODE) (LDPC_PARAM_S(CODE))
#define LDPC_LENGTH_VI(CODE) (LDPC_PARAM_S(CODE))
#define LDPC_SIZE_CI(CODE)   (LDPC_LENGTH_CI(CODE) * sizeof(uint16_t))
#define LDPC_SIZE_VI(CODE)   (LDPC_LENGTH_VI(CODE) * sizeof(uint16_t))

/* Sparse parity CS length (in uint16_t) and size (in bytes) */
#define LDPC_LENGTH_CS(CODE) (  LDPC_PARAM_N(CODE) \
                              - LDPC_PARAM_K(CODE) \
                              + LDPC_PARAM_P(CODE) \
                              + 1)
#define LDPC_SIZE_CS(CODE)   (LDPC_LENGTH_CS(CODE) * sizeof(uint16_t))

/* Sparse parity VS length (in uint16_t) and size (in bytes) */
#define LDPC_LENGTH_VS(CODE) (LDPC_PARAM_N(CODE) + LDPC_PARAM_P(CODE) + 1)
#define LDPC_SIZE_VS(CODE)   (LDPC_LENGTH_VS(CODE) * sizeof(uint16_t))

/* Sparse parity overall size in bytes */
#define LDPC_SIZE_SPARSE_H(CODE) (  LDPC_SIZE_CI(CODE) \
                                  + LDPC_SIZE_CS(CODE) \
                                  + LDPC_SIZE_VI(CODE) \
                                  + LDPC_SIZE_VS(CODE))

/* Expanded generator matrix sizes in bytes, for fast encoder.
 * There are k*(n-k)/8 parity bits (excluding systematic identity matrix).
 */
#define LDPC_SIZE_G(CODE) (LDPC_PARAM_K(CODE) \
                           * (LDPC_PARAM_N(CODE) - LDPC_PARAM_K(CODE)) / 8)

/* Expanded generator matrix length in uint32_t */
#define LDPC_LENGTH_G(CODE) (LDPC_SIZE_G(CODE) / 4)

/* Required working area for the bitflipping algorithm, in bytes.
 * (=(n+p)/8 + n + p)
 */
#define LDPC_SIZE_BF_WA(CODE) (((LDPC_PARAM_N(CODE) + LDPC_PARAM_P(CODE)) \
                               * 9) / 8)
#define LDPC_LENGTH_BF_WA(CODE) LDPC_SIZE_BF_WA(CODE)

/* MP decoder LLR length (in floats) and size (in bytes) */
#define LDPC_LENGTH_MP_LLRS(CODE) (LDPC_PARAM_N(CODE))
#define LDPC_SIZE_MP_LLRS(CODE)   (LDPC_LENGTH_MP_LLRS(CODE) * sizeof(float))

/* MP decoder working area length (in floats) and size (in bytes) */
#define LDPC_LENGTH_MP_WA(CODE) (2 * LDPC_PARAM_S(CODE))
#define LDPC_SIZE_MP_WA(CODE)   (LDPC_LENGTH_MP_WA(CODE) * sizeof(float))

/* Size of output for message passing decoder in bytes (=(n+p)/8) */
#define LDPC_SIZE_MP_OUT(CODE) ((LDPC_PARAM_N(CODE) + LDPC_PARAM_P(CODE)) / 8)
#define LDPC_LENGTH_MP_OUT(CODE) LDPC_SIZE_MP_OUT(CODE)

/* Add up the required size for each type of encoder, including the size
 * needed to store the resulting output.
 */
#define LDPC_SIZE_TX_SMALL(CODE)  (LDPC_PARAM_N(CODE) / 8)
#define LDPC_SIZE_TX_FAST(CODE)   (LDPC_SIZE_G(CODE) + LDPC_PARAM_N(CODE) / 8)

/* Add up the required size for each type of decoder, including the size
 * needed to store the input and resulting output.
 */
#define LDPC_SIZE_RX_BF(CODE)    (  LDPC_SIZE_SPARSE_CI(CODE) \
                                  + LDPC_SIZE_SPARSE_CS(CODE) \
                                  + LDPC_SIZE_BF_WA(CODE)     \
                                  + LDPC_PARAM_N(CODE)/8      \
                                  + LDPC_PARAM_K(CODE)/8)

#define LDPC_SIZE_RX_MP(CODE)    (  LDPC_SIZE_SPARSE_H(CODE) \
                                  + LDPC_SIZE_MP_LLRS(CODE)  \
                                  + LDPC_SIZE_MP_WA(CODE)    \
                                  + LDPC_PARAM_N(CODE)/8     \
                                  + LDPC_SIZE_MP_OUT(CODE))

/* Find the size needed, parameterised over the encoder, decoder,
 * TX code, and RX code, to store all the relevant expanded codes, working
 * areas, and the outputs of the encoders and decoders.
 */
#define LDPC_SIZE(ENCODER, TXCODE, DECODER, RXCODE)     \
    (LDPC_SIZE_TX_##ENCODER(TXCODE) + LDPC_SIZE_RX_##DECODER(RXCODE))

#endif
