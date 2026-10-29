# mlkem-1k: size-optimized ML-KEM for x86-64

2026-10-28  Markku-Juhani O. Saarinen < mjos@iki.fi >

The goal is the smallest possible x86-64 machine code for a complete,
standards-compliant **ML-KEM-512** ([FIPS 203](https://doi.org/10.6028/NIST.FIPS.203)): key
generation, encapsulation and decapsulation, including all of SHA-3/SHAKE.

The open question was: Can ML-KEM-512 fit in **under 1000 bytes**? Other attempts
were at roughly 2000 bytes, including the initial effort with GPT 6.1 Sol.

None of this is meant for production. Speed, constant-time execution and
side-channel resistance are explicitly out of scope.

## Status

[`mlkem.S`](mlkem.S) is currently **986 bytes**, under 1000.
It passes all tests for ML-KEM-512, ML-KEM-768 and ML-KEM-1024. ML-KEM-768
and ML-KEM-1024 build to 981 bytes. ML-KEM-512 has two special cases under `#if`:
J and H share a length (`|z || c| = |ek| = 800`), and the CBD sampler picks
eta = 3 or 2 from the nonce. The other sets use eta = 2 throughout.

| Component | Functions | Bytes |
|-----------|-----------|------:|
| Keccak    | `keccak_f1600`, `keccak_absorb`, `xhash_out`, `hash_hdr`, `hash_G`, `hash_H` | 245 |
| Arithmetic mod q | `ntt_common` / `polymul_acc`, `ntt`, `invntt`, `dot` | 190 |
| Encoding  | `poly_encode` (= `poly_compress_encode`), `poly_decode`, `get_bits` | 100 |
| Sampling  | `sample_cbd` (also SampleNTT and ByteDecode_12), `cbd_vec` | 175 |
| ML-KEM (with K-PKE inlined) | `frame_enter`, `mlkem_keygen`, `mlkem_encaps`, `mlkem_decaps` | 276 |
| **Total** | | **986** |

Run `make size` for the current per-function numbers. History: 2284 bytes
for the first version, then 2272, 1454, 1217, 1138, 1102, 1077, 1070, 1067, 1061, 1041, 1039, 1029, 1022 and now 986.

## Rules

The score is the size of the `.text` section of `mlkem.S` built for
ML-KEM-512, as reported by `make size`. Object file headers don't count. A
submission must also:

* pass `make test`, or at least `./mlkem512_test`;
* implement all three functions with the C interface below;
* carry no data tables, no libc calls and no system calls (constants and
  tables inside `.text` count toward the size);
* be position independent and keep all state on the stack;
* follow the System V ABI at the three exported entry points (internal
  routines may use any convention).

Variable time, branches on secrets and slow code are all allowed. The current
code spends milliseconds on one key exchange, which is fine.

## Interface

```c
void mlkem_keygen(uint8_t ek[800], uint8_t dk[1632],
                  const uint8_t d[32], const uint8_t z[32]);   // Alg. 16
void mlkem_encaps(uint8_t ct[768], uint8_t k[32],
                  const uint8_t ek[800], const uint8_t m[32]);  // Alg. 17
void mlkem_decaps(uint8_t k[32], const uint8_t ct[768],
                  const uint8_t dk[1632]);                      // Alg. 18
```

These are the deterministic "internal" functions of FIPS 203. The caller
supplies the randomness `d`, `z` and `m`, which is also how the ACVP tests
work. Outputs must not overlap inputs. Sizes are shown for ML-KEM-512.
`-DMLKEM_K=3` or `-DMLKEM_K=4` selects ML-KEM-768 or ML-KEM-1024.

## Build and test

```
make test     # build and run the ML-KEM-512, -768 and -1024 harnesses
make size     # code size per parameter set, per-function sizes
make kat      # regenerate the test-vector block for mlkem_test.c
```

**NOTE:** `make kat` has dependencies that are not included in this repo.

It needs GCC (or Clang) and GNU binutils on x86-64 Linux. `make kat` also
needs Python 3 with `pycryptodome` and the sibling [py-acvp-pqc](../py-acvp-pqc)
directory. It takes about half a minute.

[`mlkem_test.c`](mlkem_test.c) is self-contained. It tests:

1. **ACVP KeyGen.** Two vectors per parameter set, checked against SHA3-256
   hashes of `ek` and `dk`.
2. **ACVP Encaps.** One vector, checking the shared key in full and a hash of
   the ciphertext.
3. **ACVP Decaps.** Two valid and two modified ciphertexts. The modified ones
   exercise implicit rejection.
4. **Round trip and rejection.** Single-bit flips at the boundaries of the
   ciphertext's two parts. The expected rejection key `J(z || c')` comes from
   an independent C Keccak inside the harness.
5. **Bulk.** 256 complete keygen, encaps, decaps and rejection cycles, folded
   into one digest that [`gen_kat.py`](gen_kat.py) computes with the Python
   reference model `fips203.py`. Few vectors miss rare-case bugs such as an
   off-by-one rounding constant in compression. This test catches them.

Every call is also checked for writes past the end of its outputs and for
modified inputs. Large outputs are stored as hashes, so the harness stays
small enough to read.

When changing the code, test with deliberately broken variants too. A size
trick that passes the ACVP vectors can still be wrong in a rare case, and the
bulk test exists for exactly that.

One class of bug no black-box test can catch: off-by-one errors inside
K-PKE.Decrypt, such as a wrong rounding constant in `Decompress_d`. For an
invalid ciphertext, decapsulation returns `J(z || c)` whatever message it
decrypts. For a valid one, the decryption noise margin absorbs a ±1 error.
Changes to the decryption arithmetic must therefore be shown exact, for
example by exhaustive checks over all inputs, not just by passing the
tests. The current `poly_decode`, `poly_compress_encode` and the 12-bit decode
were checked this way for every input they can receive.

The ACVP vectors also never contain non-canonical keys, meaning ek or dk
coefficients in 3329..4095. Changes to key decoding should be compared
against a previous version on random key bytes.

## Design

* **One prologue.** The three entry points start with `call frame_enter`.
  It saves the registers, maps the arguments, allocates one stack frame and
  calls back into the body. K-PKE is not a separate layer: encryption is
  inlined into `mlkem_encaps` and decryption into `mlkem_decaps`.
* **Decaps re-encrypts by calling encaps.** This recomputes `h = H(ek)`
  instead of reading `h` from `dk`. The two agree for every `dk` that
  passes the FIPS 203 input check, which is left to the caller. The
  ciphertext comparison (`repe cmpsq`) and the key selection branch, so
  decapsulation is not constant-time. Decryption needs `v - sᵀu`; instead
  of negating anything, decaps decodes `s` through `NTT⁻¹` with the scale
  factor `-128⁻¹` in place of `128⁻¹`, so the shared row product computes
  `v + (-s)ᵀu` directly.
* **Products in the normal domain.** All polynomial products use one
  negacyclic schoolbook multiply-accumulate (`polymul_acc`). `dot`
  computes `acc += Σ M[j]·v[j]` for all three row products: `A·s` in
  keygen, `Aᵀ·y` in encryption and `sᵀ·u` in decryption. The matrix A is
  never stored; each entry is sampled with SHAKE128 when needed.
* **The NTT only converts.** FIPS 203 defines `A^`, `t^` and `s^` in the
  NTT domain, so it is still needed: `NTT⁻¹` after `SampleNTT` and after
  decoding `t^` or `s^`, and `NTT` before encoding them in keygen. It is
  evaluated directly from its definition in O(n²), with a table of powers of
  17 pushed onto the stack. It shares its reduce-and-store tail with
  `polymul_acc`. There are no butterflies and no zeta tables.
* **Keccak.** One sponge serves H, G, J, PRF and XOF through a "mode" word
  holding the rate and padding byte. Round constants come from the LFSR.
  Rotation amounts and lane order come from the (x, y) walk. The hash
  front end has no frame of its own: the state sits right after the input
  (G and J hash in place), and the routine copies a fixed 32 bytes to an
  output pointer that the caller pushed. `hash_hdr` hashes from the
  header to RDI, which serves both G and J. CBD reads its PRF bits
  straight from the Keccak state, after the first block has been copied
  just below it.
* **Encoding and sampling.** One bit-serial reader (`get_bits`) handles
  every `ByteDecode_d`. One bit-serial writer (`rcr byte ptr [rdi], 1`)
  handles every `ByteEncode_d` and compresses each coefficient on the way:
  `floor((x·2^d + 1664) / q)`, one division, exact because q is odd; decompression multiplies by q and rounds. `SampleNTT`, `SamplePolyCBD` and the
  `ByteDecode_12` of `t^` and `s^` share one sampler, selected by SF and
  ZF. SampleNTT reads candidates straight from the Keccak state and
  rejects those `>= q`; the decode mode reads from the key and accepts
  every 12-bit value, so out-of-range values are reduced mod q, as FIPS 203
  requires when the caller skips the input checks.
* **Register contracts.** Internal routines return pointers where the next
  routine wants them (for example, start and end in RCX and RDI). Many
  routines fall through into the next one or share its `ret`. Each routine
  documents its contract in `mlkem.S`.
* **Stack.** Everything lives on the stack. Peak use is 7.2 kB
  (ML-KEM-512 encaps) to 14.4 kB (ML-KEM-1024 decaps). Each
  operation takes several milliseconds.

## Ideas for getting smaller

Successive passes saved 237, 79, 36, 25, 7, 3, 6, 20, 2, 10, 0, 7 and 36
bytes. The last pass came from register and frame plumbing: the hash front
end, R9 as `dot`'s scratch pointer, and R8 = q set once in the sampler.
Leads that were measured:

* **Mirrored multiplier.** `polymul_acc`'s 7-byte negacyclic sign fix
  becomes a 2-byte `adc` if the scratch polynomial holds (b, q - b). Both
  ways of producing that work, but cost more than they save: evaluating
  the inverse transform at 512 points (+6 net) or mirroring in `dot`
  (+17).
* **Keccak.** At 150 bytes, `keccak_f1600` is the largest routine, and its
  round loop closes with a 6-byte `jne rel32`. Precomputing the 24 round
  constants onto the stack shrinks the round body enough for a 2-byte
  jump, but the generator makes the routine 158 bytes. Zero slots that
  remove theta's top-row special case give 156. SSE2 (4-5 byte
  instructions, unaligned rows, no 64-bit rotate) and other theta, chi
  and iota forms were all counted larger.
* **Tried and measured, no gain:** one Keccak state shared across a
  vector's CBD polynomials (+18), a unified sampler and hash front end,
  merging keygen's flow with encaps, NTT-domain basemul, sharing a loop
  between SampleNTT and CBD (6-bit CBD samples straddle SHAKE256 blocks),
  and dozens of register and layout variants.
* **Built at 986, no gain:** pushing the first PRF block onto the stack
  instead of copying it, and counting `dot`'s rows down.

## Files

| File | Purpose |
|------|---------|
| [`mlkem.S`](mlkem.S) | the implementation (GAS, Intel syntax, run through the C preprocessor) |
| [`mlkem_test.c`](mlkem_test.c) | self-contained test harness with embedded test vectors |
| [`gen_kat.py`](gen_kat.py) | generates the test-vector block from the ACVP JSON files and `fips203.py` |
| [`Makefile`](Makefile) | `test`, `size`, `kat` and `clean` targets |

Test vectors come from NIST's
[ACVP-Server](https://github.com/usnistgov/ACVP-Server) JSON files, through
the local copy in `py-acvp-pqc/json-copy`.
