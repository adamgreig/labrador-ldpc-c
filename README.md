# Labrador LDPC Library

This library handles [LDPC][1] encoding and decoding. It is designed to be 
usable on embedded platforms as well as in phone, desktop, or server 
applications.

Various LDPC codes of different lengths and rates are available. Functions are 
available to encode and decode, both when RAM and CPU is abundant and (to some 
extent) when it is not. No dynamic allocations are performed so the user is 
entirely in control of the memory profile.

The codes used are specified by the [CCSDS][2] and are recommended for space 
telecommand and telemetry applications.

[1]: https://en.wikipedia.org/wiki/Low-density_parity-check_code
[2]: https://public.ccsds.org/

## Quick Summary

This summary contains the basic information on getting started. The subsequent 
sections have more technical detail.

#### Pick a Code

You can choose from different code sizes: 8, 16, 32, or 128 bytes of user 
information. For the 128 byte packets, you can choose code rate 1/2, 2/3, or 
4/5; for the other sizes only rate 1/2 is available. The rate means how much 
coding overhead/redundancy is added when transmitting. Longer codes have higher 
performance but take longer to transmit so increase latency and might be more 
susceptible to burst errors. The codes are chosen by the `ldpc_code` enum, for 
example `LDPC_CODE_N1280_K1024` for a 128-byte (1024-bit) rate 4/5 code.

For this example we'll use the `LDPC_CODE_N1280_K1024`, but you'll need to 
change the memory sizes for whichever code you're using (all detailed in 
`ldpc_codes.h`).

#### Initialisation

To store the data-to-send and the codeword-to-transmit:
```c
uint8_t txdata[128];
uint8_t txcode[160];
```

No further initialisation needed for the slower encoder. For quicker but more 
RAM-using encoding:
```c
uint32_t g[1024*8];
ldpc_codes_init_generator(LDPC_CODE_N1280_K1024, g);
```

To store hard information for a decoder (only if you don't have soft 
information):
```c
uint8_t rxcode[160];
```

To store soft information for a decoder (needed even if you only have hard 
information):
```c
float rxllrs[1280];
```

Initialising the decoder:
```c
uint32_t h[384*44];
uint16_t ci[4992], vi[4992], cs[385], vs[1409];
float workingarea[9984];

ldpc_codes_init_paritycheck(LDPC_CODE_N1280_K1024, h);
ldpc_codes_init_sparse_paritycheck(LDPC_CODE_N1280_K1024, h, ci, cs, vi, vs);
```

To store the decoded output:
```c
uint8_t rxdata[176];
```

#### Encoding

Write into `txdata`.

Slow but low-memory encoding:
```c
ldpc_encode_small(LDPC_CODE_N1280_K1024, txdata, txcode);
```

Quick encoding:
```c
ldpc_encode_fast(LDPC_CODE_N1280_K1024, g, txdata, txcode);
```

Now transmit `txcode` out the radio.

#### Decoding

If you have hard information in `rxcode`:
```c
ldpc_decode_hard_llrs(LDPC_CODE_K1280_N1024, rxcode, rxllrs);
```

Then, or if you already have soft information (more-positive means 
more-likely-to-be-zero):
```c
ldpc_decode_mp(LDPC_CODE_K1280_N1024, ci, cs, vi, vs,
               rxllrs, rxdata, workingarea);
```

If it returns `true`, a codeword was found, which almost certainly means valid 
data was recovered into the first *n* bytes of `rxdata`. If it returns `false`, 
it tried its best, and you might find the output useful, but it's not a valid 
codeword.


## The Codes

Various LDPC codes of different lengths and rates are available. Currently we 
support six codes. They're specified as (n, k), where *n* is the number of bits 
transmitted (the codeword length), and *k* is the number of information bits 
encoded. Divide *k* by 8 to get the packet size in bytes. Divide *n* by *k* to 
get the code rate.

All the codes use a protograph-plus-circulant construction which has various 
nice features; primarily that it is possible to represent both the generator 
matrix and the parity check matrix in a very compact form which is good for 
embedded platforms.

The short rate 1/2 codes come from the telecommand channel coding experimental 
LDPC standard, 231.1-O-1, currently in review to be included in the main 
standard. The longer codes are available in rates 1/2, 2/3, and 4/5, and come 
from the telemetry standard 131.0-B-2. All the standards are freely available 
from the CCSDS website.

* `LDPC_CODE_N128_K64`, a short code useful for very small packets or very 
  limited receive hardware, with the poorest performance of the lot;
* `LDPC_CODE_N256_K128`, a somewhat longer code perhaps useful for short 
  packets;
* `LDPC_CODE_N512_K256`, the longest code from the TC spec;
* `LDPC_CODE_N1280_K1024`, a high-rate code for situations where less error 
  coding is required;
* `LDPC_CODE_N1536_K1024`, a medium-rate code;
* `LDPC_CODE_N2048_K1024`, a low-rate code for longer packets.

It is straightforward to add support for the 4096-bit codes in the standard; 
they're only omitted to keep the code neater as no current use is foreseen. 
Adding the 16384-bit codes is also possible but will take a little extra work 
in improving the generator-matrix script.

## Encoding

In a mathematical sense, encoding takes some user data bits *x* and multiplies 
it with the code's *generator matrix*, *G*, to get a codeword *c=xG*.

To ease code space requirements, all the generator matrices are stored in a 
compact form afforded us by the nature of the code. These constants are either 
specified in the standard (for the TC codes) or computed via the 
`scripts/ccsds_131.py` script.

It's possible to encode a codeword directly from this form, using the 
`ldpc_encode_small` function. The RAM footprint is small, but it takes around 
two or three orders of magnitude more time to encode. Might not be an issue for 
many use cases with low-speed radios.

For faster encoding, the compact generator matrix can be expanded in RAM via 
the `ldpc_codes_init_generator` and then the encoder function 
`ldpc_encode_fast` can be used. The amount of RAM needed depends on the code, 
see `ldpc_codes.h` for details, but it's between 512 bytes and 131kB.

Both encoders take a `uint8_t[k/8]` and write to a `uint8_t[n/8]`.

## Decoding

Decoding takes information from the radio and attempts to find the codeword 
that fits the received data best. Sometimes this is not possible and decoding 
fails, but often we can find the original codeword and thus the original 
transmitted information. In the best case, the radio tells us how close to 
being a 0 or 1 each bit was (called soft information), but often we are only 
given its best guess (hard information).

Decoding relies on the *parity check matrix*, *H*. This is specified in the 
standards and is represented in broadly the same form in the source file, as it 
is a very compact encoding. To turn this compact encoding into an actual parity 
check matrix, use `ldpc_codes_init_paritycheck`. This requires a reasonable 
amount of RAM, see `ldpc_codes.h` for details. Both decoders require further 
processing of this matrix into a sparse representation, provided by 
`ldpc_codes_init_sparse_paritycheck`. The sparse representation lists the 
indices that are non-zero along each row and along each column, and separately 
has a list of the starting position in this list for each row and column. In 
the future it may be possible to omit the full generation step, which would 
save a lot of RAM on embedded platforms, but it is not currently implemented.

There are two available decoders, the bit-flipping decoder `ldpc_decode_bf` 
(fast but bad) and the message-passing decoder `ldpc_decode_mp` (slow but 
good). The difference between them amounts to 2dB of SNR (in other words, you'd 
need 1.6 times as much transmit power to make up the difference in decoders), 
but around a factor of 50 time difference.

The bit flipping decoder is fast and only uses hard information (in other 
words, you know only whether each bit is more like a 0 or a 1). The basic idea 
is to find all the bits that are involved in the most parity checks which are 
not satisfied, and try flipping them until it fixes things. It cannot operate 
on punctured codes, such as the 1024-bit codes, as they need to represent an 
"unknown" bit. A pre-processing erasure-decoding step could alleviate this 
restriction but has not been written.

The message passing decoder involves substantially more computation for each 
iteration and works by computing the probability of each bit given the other 
bits connected to it via the parity check matrix. It can handle soft 
information and as such is also suited for the punctured codes. The 
implementation here uses min-sum with self-correction [paper][3] which is 
reasonably efficient and performs almost as well as an "optimal" sum-product 
decoding.

[3]: https://arxiv.org/pdf/0803.1090.pdf

Possible improvements here include a better min-sum correction factor (for 
instance using normalised or offset min-sum may perform better if the correct 
normalisation factor for our various codes is computed), and extending the 
sparse representation to include an inverse lookup which would shave 20% off 
the runtime at the cost of twice the RAM. More drastic options include writing 
a fixed point decoder which could likely be very fast and take advantage of the 
STM32F4's DSP instruction set (since 8 bits per symbol appears to be plenty for 
this algorithm).


## Executables

Right now two test executables are included, `ber_trials` and 
`throughput_trials`, which run a number of encoder tests to determine BER 
against Eb/N0 and throughput in kbps against Eb/N0. To configure the code in 
use and the Eb/N0 points, edit the relevant `.c` and recompile.

The rest can then be used to plot scripts via `plot_ber.py` and 
`plot_throughput.py` in `scripts/` though this is still a semi-manual process.
