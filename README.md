# minipqc / mlkem: size-optimized ML-KEM for x86-64

The goal is the smallest possible x86-64 machine code for a complete,
standards-compliant **ML-KEM-512** ([FIPS 203](https://doi.org/10.6028/NIST.FIPS.203)): key
generation, encapsulation and decapsulation, including all of SHA-3/SHAKE.
It implements the deterministic internal algorithms (FIPS 203 Algs. 16–18);
random generation and input checking are left to the caller (see
[Interface](#interface)).

The open question was: Can ML-KEM-512
fit in **under 1000 bytes**? Other attempts are at roughly 2000 bytes.

None of this is meant for production. Speed, constant-time execution and
side-channel resistance are explicitly out of scope.

## Status

[`mlkem.S`](mlkem.S) is currently **897 bytes**, under 1000.
It passes all tests for ML-KEM-512, ML-KEM-768 and ML-KEM-1024. ML-KEM-768
and ML-KEM-1024 build to 899 and 900 bytes. ML-KEM-512 has four special cases
under `#if`: J and H share a length (`|z || c| = |ek| = 800`), the CBD
sampler picks eta = 3 or 2 from the nonce (the other sets use eta = 2
throughout), and encaps' two-row c1 loop flips BL with `neg bl`, and `invntt` runs in the
sampler's frame (about 340 bytes more stack, still below ML-KEM-1024's).

| Component | Functions | Bytes |
|-----------|-----------|------:|
| Keccak    | `keccak_f1600`, `keccak_absorb`, `xhash_out`, `hash_hdr`, `hash_G`, `hash_G_cbd`, `hash_H` | 234 |
| Arithmetic mod q | `ntt_common` / `polymul_acc`, `ntt`, `invntt`, `dot` / `kpke_rows` | 185 |
| Encoding and sampling | `poly_encode` (= `poly_compress_encode`), `poly_decode` / `get_bits`, `sample_cbd` (also SampleNTT and ByteDecode_12, including `.Lsn_accept`, placed after `dot`), `cbd_vec` | 233 |
| ML-KEM (with K-PKE inlined) | `frame_enter`, `mlkem_keygen`, `mlkem_encaps`, `mlkem_decaps` | 245 |
| **Total** | | **897** |

Run `make size` for the current per-function numbers. History: 2284 bytes
for the first version, then 2272, 1454, 1217, 1138, 1102, 1077, 1070, 1067, 1061, 1041, 1039, 1029, 1022, 986, 978, 953, 945, 943, 936, 930, 924, 922, 921, 920, 914, 907, 903, 899 and now 897.

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

**NOTE:** Some of the dependencies for `make kat` are missing from this repo.

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

The ACVP tests also check every call for writes past the end of its outputs
and for modified inputs. Large outputs are stored as hashes, so the harness stays
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
* **Decaps re-encrypts inside its own frame.** It calls into the
  encaps body (`.Lenc_body`) rather than `mlkem_encaps`, writing m' into
  the frame header and c' where its own `s` was. This recomputes `h = H(ek)`
  instead of reading `h` from `dk`. The two agree for every `dk` that
  passes the FIPS 203 input check, which is left to the caller. The
  ciphertext comparison (a fused `lodsb; xchg al, [rdi]; scasb; loope`
  loop that also copies c over c' up to the first difference) and the key selection branch, so
  decapsulation is not constant-time. Decryption needs `v - sᵀu`; instead
  of negating anything, decaps decodes `s` through `NTT⁻¹` with the scale
  factor `-128⁻¹` in place of `128⁻¹`, so the shared row product computes
  `v + (-s)ᵀu` directly.
* **One row loop.** keygen, encaps and decaps share `kpke_rows`, folded
  into the end of `dot`. keygen and encaps reach it through
  `hash_G_cbd`: a `call hash_G` placed just before `cbd_vec`, so G
  returns straight into the noise sampler, which ends in the row loop.
  keygen marks itself with EBX bit 31 (transposed A); decaps enters
  from `poly_decode` with one row.
* **32-bit coefficients.** Polynomials are 256 × uint32 (1 kB). The
  multiplier reads its operand straight from memory and several loads and
  stores get shorter (−10 bytes). Alone it would cost about 10 kB of
  stack; with decaps re-encrypting in its own frame the peak goes down.
* **Polynomial-end list.** `cbd_vec` pushes the vector's start, then it
  and `poly_decode` push the end address of every polynomial they write;
  the sampler writes each noise polynomial at the current top of the
  list, so nothing else tracks the output pointer. The row loop takes each row's
  accumulator from the top of that list instead of computing it, and pops
  the K − 1 unused entries when it is done (−6).
* **Products in the normal domain.** All polynomial products use one
  negacyclic schoolbook multiply-accumulate (`polymul_acc`). `dot`
  computes `acc += Σ M[j]·v[j]` for all three row products: `A·s` in
  keygen, `Aᵀ·y` in encryption and `sᵀ·u` in decryption. The matrix A is
  never stored; each entry is sampled with SHAKE128 when needed.
* **The NTT only converts.** FIPS 203 defines `A^`, `t^` and `s^` in the
  NTT domain, so it is still needed: `NTT⁻¹` after `SampleNTT` and after
  decoding `t^` or `s^`, and `NTT` before encoding them in keygen. It is
  evaluated directly from its definition in O(n²), with a table of powers of
  17 pushed onto the stack. Both directions use the same table of 17^(−e):
  the forward transform flips bits 1..7 of its output index, because
  o(127 − i) = −o(i) mod 256. The direction comes from `cmp rsi, r9`, since
  the inverse always runs on `dot`'s scratch polynomial and keygen's s and t
  lie below it. It shares its reduce-and-store tail with `polymul_acc`.
  There are no butterflies and no zeta tables.
* **Keccak.** One sponge serves H, G, J, PRF and XOF through a "mode" word in R11W
  holding the rate and padding byte. Round constants come from the LFSR.
  Rotation amounts and lane order come from the (x, y) walk.
  `keccak_absorb` takes no length: every input is placed so that it ends
  exactly where the Keccak state begins, and the loop runs until RSI
  reaches RDI. The round
  loop closes with a 2-byte `jne` onto a `test al, 0xEB` whose immediate
  plus the next opcode byte form a hidden `jmp` back to the round start; a
  build-time check (in a non-loaded section) fails if that distance moves. The hash
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
  `ByteDecode_12` of `t^` and `s^` share one sampler, selected by SF, OF
  and CF. SampleNTT reads candidates straight from the Keccak state and
  rejects those `>= q`; the decode mode reads from the key and accepts
  every 12-bit value, so out-of-range values are reduced mod q, as FIPS 203
  requires when the caller skips the input checks.
* **Register contracts.** Internal routines return pointers where the next
  routine wants them (for example, start and end in RCX and RDI). Many
  routines fall through into the next one or share its `ret`. Each routine
  documents its contract in `mlkem.S`.
* **Stack.** Everything lives on the stack. Peak use is about
  10.9 kB (ML-KEM-512) to 14.6 kB (ML-KEM-1024), nearly the same for all
  three operations. Each
  operation takes several milliseconds.

## Ideas for getting smaller

**Stack for size.** A pass with stack use allowed up to 32 kB found that
almost nothing gets smaller with more stack: decaps calling
`mlkem_encaps` in its own frame again is 0 (stack 24.9 kB), an NTT that
reads a pushed copy of its input is +6, an aligned frame is +3, and one
Keccak frame per call is +2. Its one win, the polynomial-end list above,
costs only 16 to 56 bytes of stack. Counted: the sampler's work area in
the body frame (+6), never reusing buffers (+2 to +20), one worst-case
frame (disp32 offsets), a mirrored multiplier operand (9 to 16 to save
5).

**Latest passes.** The polynomial-end list (−6) and then its extension
to the vector start, a fused compare-and-copy in decaps, and decaps
entering `poly_decode` mid-loop (−7). Rejected there: a list of starts
instead of ends (+2), qword compares (+4), `hash_H` reading ek from R10
(+1), `dot` leaving RSI = RBP (+2).

**Pass at 907 (−4).** Decaps' decode row uses BL = 255 − K, with a
brute-forced mode test `add al, 128 − K; add al, 128 + K` (SF ≠ OF for
decode, CF for SampleNTT); `dot` sets CL = 48 so the 12-bit decode
advances with `lea r10, [rsi+rcx*8]` instead of a 32-bit displacement;
and both producers hand the row loop BH = 0, so `kpke_rows` drops `mov
bh, 0`. Brute force found no one-instruction decode test and no shorter
K = 2 width choice. Rejected: two direct encoder calls for c1 (+1), a
32-bit accept test (0), `shld`/`rol` seed swaps (+2).

**Where it stands.** Single instructions, adjacent pairs and every
window of up to six instructions have been searched mechanically, and
every routine order measured. The one stack-for-size trade found
(`invntt` reusing `sample_cbd`'s frame, ML-KEM-512 only) is now adopted.
A pass at 903 also searched non-local pairs, about 1.6 million in all:
any two non-adjacent deletions; a 2- or 3-byte setter replaced by a 1-
or 2-byte instruction elsewhere, up to 40 instructions back or at a
routine boundary; re-registered push, pop and `xchg`; changed `reg,reg`
operands; and `mov r8, imm8` constants. The only screens that passed were
known-invalid (the rounding `adc`, `cdq`, keygen leaving RBX to the
caller). Hashing ek inside dk (+1) and a `rol`-based encoder phase (0)
were built and rejected.

**Pass at 899 (−4).** Whole routines were rewritten from the
specification instead of edited locally. The win came from the NTT's
contract rather than its algorithm: one power table built with
`imul eax, eax, 17` (exact for every constant, including −128⁻¹), the
forward direction reading it through o(127 − i) = −o(i), and the direction
taken from RSI vs R9. The transforms no longer touch RDI, so keygen drops
its `mov edi, imm32` and the push and pop around `ntt`. Built and
rejected: the unified table with EDI still as the direction flag (+2).
Counted and not worth building: polymul through a mirrored table of b
(+18), decaps computing J(z‖c) first (+11), keygen building ek in the
frame and hashing it in place (+17), mu folded into the v compression
(≥ 0), iota counting with ECX as the bit position (+1), theta fused into
the pi cycle (≥ 0), column parity accumulated during chi, out-of-place
pi, the pi cycle walked backwards (0 to +2), and Horner or butterfly
NTTs. Searched with no solution: ending the rho loop on flags alone with
CH and CL scaled by 64. Bluestein's NTT-as-polymul needs a 512th root of
unity mod q, which does not exist.

**Pass at 897 (−2).** Re-running the mechanical sweeps on the new NTT
contract (about 940,000 jobs) found two places where `xor edx, edx`
becomes a `cdq` one instruction earlier, while EAX is provably below
2^31: before the table loop's `imul` (EAX is c, or a remainder < q) and
before the output loop's `pop` (EAX is a quotient < 2^32 / q, an index
< 512 or a product < 2^24). After the `imul` or the `pop` it would be
wrong: decaps' 17c and a polymul slot can have bit 31 set. Ruled out: a
2-byte direction test. EAX, ECX, EDX and EBX hold loop indices, RSP and
RBP lie on the same side of both kinds of polynomial, and the flags are
reset every term. The one candidate the sweep found, `cmp esi, edi`
(896), changes keygen's dk depending on whether the caller's buffer is
static or on the stack. Also ruled out: a 2-byte setter for keygen's
c = 1. On the first call RAX is the vector start left by `dot`, and the
only zero register, RCX, must stay zero, so `inc eax` or `mov al, 1`
needs a fix that costs the byte it saves.
Counted: parity in RCX (+1), b in RDI (+1), `push 48; pop rcx` in dot
(+1), `not al` for the index flip (0).

Successive passes saved 237, 79, 36, 25, 7, 3, 6, 20, 2, 10, 0, 7, 36, 8,
25, 8, 2, 7, 6, 6, 2 and 1 bytes, and later passes 1, 6, 7, 4, 0, 4 and 2.
One of them was a superoptimizer-style search. A unicorn trace of full KEM runs gave per-instruction liveness
and about 120 real machine states per instruction; every window of 1–6
instructions was then tried against every shorter instruction sequence
(all sequences of up to 3 bytes, all valid 4-byte instructions, and a
restricted pool of 4–6 byte pairs and triples), about 5·10^10 runs in
all. It found one byte: `get_bits`' `cmp cl, 12; je` became `cmp edx,
ebx; ja`, valid in all three modes because decaps' EBX carries DEC_NEG in
its upper bits. A layout search showed every branch is already rel8 and
no routine order removes a jump. The pass before let `keccak_absorb`
consume the mode with `xchg eax, r11d` (−1), and moved `cbd_vec`'s
counter so its flags select CBD mode without `sahf` (−1). It also ran
four automated sweeps, which found nothing: deleting any one
instruction, any adjacent pair, each flag-only `cmp`/`test` with a
neighbour moved into its slot, and adjacent swaps plus a nearby
deletion. The pass before that moved the Keccak mode word
from R8W to R11W (−4). Keccak happens to hand R11 back unchanged, since
rho's junk-lane trick restores it, so R8 = q can be set once in
`frame_enter` and never changes, and the sampler and decoder use it
directly. Two missed encodings also went: the body pointer in RDX
instead of R11, and `rep movsb` for z (−2). The pass before moved SampleNTT's and the
12-bit decode's candidate loop inside `poly_decode`'s bit-reader loop
(−3; encoding and sampling code now interleave, so the table merges
them), hid all of `hash_G` inside the immediate of decaps' 10-byte
`mov r8, imm64` (−1), flipped encaps' K = 2 row counter with `neg bl`
(−1), and compared absorb pointers in 32 bits (−1). The pass before removed the absorb length (−4,
inputs end where the state begins), let `polymul_acc` push its starting
value straight from memory with the reduction ignoring the junk upper
half (−1), toggled BH for K = 2 (−1), and copied ek with `rep movsq`
(−1). The pass before switched to 32-bit coefficients
(−10) and paid for the stack by re-encrypting inside decaps' frame (+8).
The pass before changed the Keccak calling convention
(mode byte holds rate − 1, `keccak_absorb` takes length − 1 and loops
do-while), let two instructions double as another's immediate
(`and ecx, imm8` swallowing a `push rbx`, `mov r8d, imm32` swallowing
`hash_G`'s push), and passed G's output in RSI. An automated pass that
deleted each of the 421 instructions in turn found nothing removable
except one rounding `adc` that only the exhaustive check catches. The
pass before merged the three row loops into one
(`kpke_rows`, −10), made G return into the noise sampler (−6), and
shared tails between routines (keygen's z copy uses `xhash_out`'s copy,
`frame_enter` uses `ntt_common`'s epilogue). It also found an
algebraic one: CBD's `sbb eax, -1` became `sbb eax, ecx`, because the
extra sum of loop counters cancels through the negation between the two
halves. The pass before used three reusable tricks: `pushf`/`popf` so
CBD and SampleNTT share one `call keccak_absorb` and one `call
keccak_f1600`; a path that harmlessly runs through another path's code
(a 1-byte `cmp eax, imm32` swallowing a 4-byte instruction); and
`get_bits` placed inside `poly_decode`'s loop with a conditional return.
The pass before that came from register and frame plumbing: the hash front
end, R9 as `dot`'s scratch pointer, and R8 = q set once in the sampler.
Leads that were measured:

* **Mirrored multiplier.** `polymul_acc`'s 7-byte negacyclic sign fix
  becomes a 2-byte `adc` if the scratch polynomial holds (b, q - b). Both
  ways of producing that work, but cost more than they save: evaluating
  the inverse transform at 512 points (+6 net) or mirroring in `dot`
  (+17).
* **Keccak.** At 148 bytes, `keccak_f1600` is still the largest routine.
  The hidden-`jmp` trick already gives the round loop a 2-byte branch.
  Precomputed round constants (158), zero slots that remove theta's
  top-row special case (156), SSE2 (4-5 byte instructions, unaligned
  rows, no 64-bit rotate) and other theta, chi and iota forms were all
  larger.
* **Tried and measured, no gain:** one Keccak state shared across a
  vector's CBD polynomials (+18), a unified sampler and hash front end,
  merging keygen's flow with encaps, NTT-domain basemul, sharing a loop
  between SampleNTT and CBD (6-bit CBD samples straddle SHAKE256 blocks),
  and dozens of register and layout variants.
* **Built at 986, no gain:** pushing the first PRF block onto the stack
  instead of copying it, and counting `dot`'s rows down.
* **Stack for size.** `invntt` reusing `sample_cbd`'s frame (an
  `and ecx, imm32` swallowing `ntt_common`'s prologue) gives 944 bytes,
  but every operation needs 336 more bytes of stack (ML-KEM-1024 decaps
  15056). Not adopted. At 943 it fits only for ML-KEM-512 (942 bytes,
  K=2 stack +336). 32-bit coefficients with the old decaps call give 935
  bytes but about 24.4 kB of stack.
* **Found by the search, but wrong:** `cdq` for `xor edx, edx` before
  the reduction passes every test, but an ek whose t is all q − 1 makes a
  polymul sum reach 2.84·10^9 > 2^31 and the division faults. Others
  relied on stack addresses or fixed frame offsets.
* **Built at 922, no gain:** encoder rounding through R8 (+2 at best),
  CBD counting down (+1). Brute-forced with no solution: ending rho on
  flags alone, an LFSR round-end mask, a linear recurrence for the pi
  walk, and a K = 2 eta selection of 6 bytes or fewer. Counted: R8 as
  the SampleNTT threshold (+2), more global-constant registers (0 or
  worse), `frame_enter` rewrites (0 or worse). The K = 2-only stack trade
  now gives 921.
* **Built at 924, no gain:** R11 as a permanent header pointer (0).
  Brute-force searches found no 2-instruction mode test in `dot` that
  accepts decaps' BL and no 2-byte K = 2 eta selection. Counted: decaps
  running all K rows (+5), dividing by 4096 so ByteEncode_12 becomes
  Compress_12 (+2), hiding the encoder's raw path in keygen's immediate
  (no byte layout fits). The K = 2-only stack trade now gives 923.
* **Built at 930, no gain:** CBD nonces starting at −(K+1) (0), padding
  via `movzx` (failed), rho seeded from iota's EAX (breaks the hidden
  jump), `test ebx, ebx` instead of `pushf`/`popf` (0), and about ten
  sampler and keygen variants (0 to +4). A fresh delete-one-instruction
  sweep again found only the rounding `adc`. The K = 2-only stack trade
  now gives 929.
* **Built at 936, no gain:** RCX = 1 selecting polymul (+1, and RBP is
  taken). The `invntt`/`sample_cbd` frame trick now gives 935 for
  ML-KEM-512 only (K=2 stack 10792, overall peak unchanged) or 935 for
  all sets with peak 14888; not adopted. Counted: decaps flows with J
  first, J through `hash_H`, c' at header+32 (+3..+6); keygen reorders
  (+2..+7); no 2^k or b·2^8k with b < 256 is ≡ 1665 mod q, so mu can't
  be one bit or byte add.
* **Built at 943, no gain:** keeping z by shifting RBP (0, +32 bytes of
  stack). Not worth building: NTT-domain data in natural order (+10 net),
  a merged NTT/polymul inner loop (register conflicts), clean 12-bit
  decode values (+2 net), negating v via 2^d − y (inexact at the tie),
  12-bit decode moved into `poly_decode` (+10), several Keccak theta
  forms (+7, +9).
* **Built at 945, no gain:** the `dot` accumulator in RAX (+4), `hash_H`
  taking ek in RDX (0), raw ByteEncode_12 selected by `cmp cl, 12` (+2).
  Every jump is already rel8 and no routine order gives a new
  fall-through.
* **Built at 953:** Rejected: `rol`-based
  encoder bit phase (0), CBD loop before `cbd_vec` (+0..+3),
  `mov ebx, imm32` in encaps (+1). Worked out and not worth building:
  BMI1 `andn` in chi (−1, not baseline x86-64); dropping DEC_NEG,
  since Compress_1(−w) = Compress_1(w), saves 4 but needs −v (+7); a
  32-bit SampleNTT compare (+2 net); `get_bits` inlined in SampleNTT
  (rel8 range); `sub di, 512` (wrong when DI wraps).
* **Built at 978, no gain:** R9 set once in `frame_enter` (+4), z copied
  late (+5), a `lodsq`/`scasq` compare-and-copy in decaps (+2), only encaps
  swapping in `hash_G` (+1), SampleNTT falling into CBD's block copy (+1),
  Keccak preserving RSI (failed, fix ≥ 3 bytes), a forward NTT from the
  inverse table (no gain). Counted exactly and not worth building: a
  per-bit-refill CBD reader, the sampler's work area in the body frame
  (disp32), a rate-dependent hash copy length, and a `[b | q-b]`
  polynomial for the multiplier (≥ 8 bytes of loop control to save 5).

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
