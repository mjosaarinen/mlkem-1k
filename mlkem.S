//  mlkem.S
//  === A size-oriented ML-KEM (FIPS 203) for x86-64 (SysV ABI).
//
//  No precomputed tables or matrix storage; encoding and Keccak are shared
//  across all three FIPS 203 parameter sets.
//
//  NOT for production use: variable time (DIV, branches on secret data),
//  no wiping of secrets, no input validation (FIPS 203 Sect. 7.2/7.3
//  modulus / hash checks are left to the caller).
//
//  Build (GCC preprocesses .S files):
//      gcc -c -DMLKEM_K=2 mlkem.S          // ML-KEM-512  (default)
//      gcc -c -DMLKEM_K=3 mlkem.S          // ML-KEM-768
//      gcc -c -DMLKEM_K=4 mlkem.S          // ML-KEM-1024
//
//  C interface (the deterministic "internal" functions of FIPS 203):
//
//  void mlkem_keygen(uint8_t ek[EK], uint8_t dk[DK],
//                    const uint8_t d[32], const uint8_t z[32]);   // Alg. 16
//  void mlkem_encaps(uint8_t ct[CT], uint8_t k[32],
//                    const uint8_t ek[EK], const uint8_t m[32]);  // Alg. 17
//  void mlkem_decaps(uint8_t k[32], const uint8_t ct[CT],
//                    const uint8_t dk[DK]);                       // Alg. 18
//
//  Outputs must not overlap inputs. No static data is used: everything
//  lives on the stack (peak use is about 14.6 kB, ML-KEM-1024).
//
//  Structure: the exported functions share one prologue (frame_enter).
//  K-PKE.Encrypt is inlined into encaps and K-PKE.Decrypt into decaps;
//  decaps re-encrypts by calling the encaps body inside its own frame (no
//  second frame). All three share the row loop
//  "kpke_rows" / "dot" (keygen and encaps reach it through G and the noise
//  sampler, hash_G_cbd). All polynomial products are computed in the
//  normal domain by a negacyclic schoolbook multiplication; the O(n^2) NTT
//  evaluator only converts the NTT-domain values that FIPS 203 defines
//  (A^, t^, s^).
//  SampleNTT, SamplePolyCBD and ByteDecode_12 of t^ / s^ share a sampler.
//
//  Internal routines use custom register contracts, documented at each
//  routine; callers rely on them exactly (e.g. returned start / end
//  pointers in RCX / RDI). Polynomials are 256 x uint32 (1024 bytes): the
//  products read a, b and the accumulator with 32-bit operands (IMUL
//  takes a straight from memory); values are below 2^16.
//  Size tricks used throughout: push imm8 / pop reg for small constants
//  and push reg / pop reg for register copies, string instructions,
//  fall-through between routines, and values that routines are known to
//  leave behind (e.g. RCX = 0 or DL = 0 after Keccak, the flags
//  selecting the sampler mode), each noted where it is relied on. Two
//  registers are global: R8 = q for a whole operation (frame_enter sets
//  it; no routine changes it) and R11 = the Keccak mode (set by each hash
//  or sampler call; keccak_f1600 preserves it, keccak_absorb's padding
//  consumes it). Two
//  paths sharing a call keep their selecting flags with PUSHF / POPF; a
//  one-byte CMP EAX, imm32 skips a 4-byte instruction, one jump is hidden
//  inside another instruction (keccak_f1600, checked at build time), an
//  instruction ending a routine can take the next routine's first bytes as
//  its immediate (invntt's AND ECX; decaps' MOV R11, imm64 swallows all of
//  hash_G), a CALL followed by the next routine makes the callee return
//  into it (hash_G_cbd), SampleNTT's candidate loop runs inside
//  poly_decode's loop, and NEG BL; JNS runs a loop exactly twice from
//  BL = 0xFF (encaps' c1 loop for ML-KEM-512). cbd_vec and poly_decode
//  push the start of the vector and the end of every polynomial (the end
//  list): CBD writes at its top and dot takes each row's accumulator from
//  it. decaps compares c with c' while copying c over it (LODSB / XCHG /
//  SCASB / LOOPE). Loop counters are chosen for the values they leave
//  behind: cbd_vec counts BH down and poly_decode shifts a single BH bit
//  out, so both hand the row loop BH = 0; decaps' row number is the
//  BL = 0 its decoder wraps to. A per-term CL = 48 lets the 12-bit
//  decode advance R10 by 8 * RCX = 384 without a 32-bit displacement.

    .intel_syntax noprefix

//  --- parameter sets (FIPS 203, Table 2)

#ifndef MLKEM_K
#define MLKEM_K     2
#endif

#if     MLKEM_K == 2
#define ETA1        3
#define ETA2        2
#define DU          10
#define DV          4
#elif   MLKEM_K == 3
#define ETA1        2
#define ETA2        2
#define DU          10
#define DV          4
#elif   MLKEM_K == 4
#define ETA1        2
#define ETA2        2
#define DU          11
#define DV          5
#else
#error  "MLKEM_K must be 2, 3 or 4."
#endif

#define Q           3329
#define POLY        1024                         // bytes per polynomial
#define EK_BYTES    (384*MLKEM_K+32)
#define DK_BYTES    (768*MLKEM_K+96)
#define CT_BYTES    (32*(DU*MLKEM_K+DV))

//  Keccak sponge "modes": (rate in bytes - 1) | (domain padding byte << 8).
#define SHA3_256    0x0687                      // H: rate 136, pad 0x06
#define SHA3_512    0x0647                      // G: rate  72, pad 0x06
#define SHAKE128    0x1FA7                      // XOF: rate 168, pad 0x1F
#define SHAKE256    0x1F87                      // PRF, J: rate 136, pad 0x1F

//  Frame layout shared by the three bodies and dot. The frame starts with a
//  HDR-byte header (m || H(ek) || K || r in encaps, d || k || rho || sigma
//  in keygen; in decaps m' || H(ek) || K' || r for the encaps body, then
//  z || c for J) at RBP-HDR; RBP points just after it, at
//  the vector (s, y or u), then come the accumulators (e1[K] -> u[K] and
//  e2 -> v in encaps, v -> w in decaps, e[K] -> t[K] in keygen), then
//  dot's scratch (R9). HDR = 128 keeps all polynomials at 0 mod 8
//  (ntt_common needs bit 2 clear), the header within disp8 reach of RBP
//  and fits m || H(ek) || K || r, so G writes K || r right after its
//  input.
#define HDR         128
#define FUNC(f)     .type f, @function; f:
#define ENDF(f)     .size f, . - f

    .text

//  ========================================================================
//  Arithmetic mod q (the evaluator below the sampler, see ntt_common)
//  ========================================================================

//  The NTT is computed directly from its definition as an O(n^2)
//  evaluation (FIPS 203 eq. 4.8 / Alg. 9-10 give the same values):
//      fwd:  NTT[2i+b]  = sum_j f[2j+b] * z^(o(i)*j)
//      inv:  f[2j+b]    = 128^-1 * sum_i NTT[2i+b] * z^(-o(i)*j)
//  with z = 17 and o(i) = 2*BitRev7(i)+1. Exponents are taken mod 256 and
//  looked up in one table T[e] = c*z^-e for both directions (fwd: c = 1;
//  inv: c = 128^-1 = 3303, or -128^-1 for decaps' -s) that is pushed
//  onto the stack as P_k = c*17^(k+1) (IMUL by the immediate 17), so
//  T[e] = P_255-e sits at [rsp+8e]. The inverse uses e = o(i)*j as it
//  stands; the forward needs z^(+o(i)*j) and gets it from
//  o(127-i) = -o(i) (mod 256), i.e. by flipping bits 1..7 of its output
//  index before the bit reversal.
//  polymul_acc shares the output path: every output is pushed as an
//  unreduced qword slot whose low dword holds the sum (polymul pushes
//  acc[x] || acc[x+1] with one PUSH and ignores the high dword), and the
//  final loop pops, reduces the low dword (DIV R8D) and stores them to
//  RSI[255..0]. All entry points need R8 = q (frame_enter sets it for the
//  whole operation), preserve RSI, RBX, RBP, R8-R15 and return fully reduced
//  coefficients in [0, q). They clobber RAX, RCX and RDX (RDI is kept)
//  and use about 4 kB of stack.

//  poly_decode: Decompress_d(ByteDecode_d) of the ciphertext (decaps):
//  K polynomials with d = ECX (exactly; DU), then one with d = DV,
//  stored consecutively from RDI. BH = 64 >> K at entry: one bit that
//  SHL BH, 1 moves up after every polynomial end (DU below bit 7, DV at
//  bit 7, then BH = 0, which kpke_rows needs). The input is read from
//  RSI (kept) at the running bit position EDX, which must be 0 at entry;
//  values are below q (IMUL
//  by R8D = q). EBX bits 16+ must be nonzero (decaps' DEC_NEG: get_bits
//  takes EDX > EBX as the sampler). The entry is the per-coefficient
//  XOR EAX, EAX / INC BL / JNZ that ends every coefficient: with BL = 0xFF
//  it clears EAX, wraps BL to 0 (BL counts the 256 coefficients of each
//  polynomial) and falls through to the PUSH RDI that pushes every
//  polynomial end, so the first push is the start of the first
//  polynomial (the end list, see dot). Then falls through into
//  kpke_rows, where .Lpd_end pops the last end into R9 (dot's scratch),
//  with BL = 0 (wrapped by the last coefficient) and BH = 0: kpke_rows
//  makes BL = 255-K, one decode row. poly_decode sits after sample_cbd: its
//  loop is also SampleNTT's / decode's candidate loop (see sample_cbd).
//  Each pass first runs the sampler's block-end test CMP DH, 0xFF, which
//  never fires here (DH < 64: at most 12544 bits for K = 4).
//
//  get_bits (inside that loop): replaces the low CL (1..31) bits of EAX by
//  the CL bits at bit position DX of [RSI], little endian, and keeps EAX's
//  other bits (the callers clear them or rely on them); EDX advances by
//  CL. The 16-bit BT takes DX as a signed offset (SampleNTT reads the
//  state at negative offsets; EDX bits 16+ are ignored). Clobbers nothing
//  else. EDX > EBX selects the sampler, which continues at .Lsn_accept:
//  SampleNTT has DX >= -1600 + 12 (EDX >= 0xF9CC) and the 12-bit decode
//  EDX = 12(k+1), both above EBX = k < 256; poly_decode has EDX <= 12544
//  below EBX >= 2^16 (DEC_NEG) and continues with Decompress_d.

//  sample_cbd is the sampler body shared by three modes; the flags at
//  entry select one (no instruction before the JB changes them; PUSHF
//  keeps CF over the shared keccak_absorb call): SF != OF ByteDecode_12 of
//  a key polynomial, otherwise CF = 1 SampleNTT, SF = OF = CF = 0 (and
//  ZF = 0) SamplePolyCBD. All
//  modes need RCX < 256 (decode: RCX = 48, see below) and start with
//  EBX = 0.
//
//  The output polynomial is at R9 (kept) for SampleNTT and decode; CBD
//  writes at the top of the end list (see cbd_vec).
//
//  SampleNTT(rho || j || i) (Algorithm 7): RSI = rho (32 bytes), AL = j,
//  AH = i.
//  SHAKE128 is squeezed one block at a time; 12-bit candidates are read
//  from the state (1344 = 112 * 12 bits per block, so they never straddle
//  blocks) and accepted (CF = 1) below the threshold R11W = q.
//
//  Decode: R10 = RSI = a 384-byte ByteEncode_12 block (t^[j] in ek or
//  s^[j] in dk; decode branches off before the seed copy).
//  The same candidate loop reads 256 12-bit values from R10 with
//  R11W = SHAKE128 = 0x1FA7 > 4095 as the threshold, so every value is
//  accepted: 3329..4095 stay and invntt reduces them mod q (FIPS 203 uses
//  ByteDecode_12's output mod q, unlike SampleNTT's candidates). R10
//  advances by 8 * RCX = 384 (dot passes RCX = 48).
//
//  The candidate loop is poly_decode's loop: .Lsn_block sets DX and
//  CL = 12 and falls into it (poly_decode clears EAX after each
//  coefficient, not before);
//  get_bits (EDX > EBX here) jumps to .Lsn_accept (after dot), which stores
//  the candidate, loops back to .Lsn_candidate and finally leaves through
//  invntt. Both modes so return NTT^-1 of the sampled / decoded
//  polynomial in the normal domain at R9 (and RSI). R11W is the Keccak
//  mode and the threshold: the entry sets SHAKE128 (also decode's
//  threshold), SampleNTT then copies q from R8 (all modes need R8 = q).
//  get_bits keeps EAX bits 12+, which give invntt's c (see there): in
//  SampleNTT they are 0 (keccak_f1600 returns EAX = 127, later EAX is the
//  previous candidate); in decode they are dot's EBX bits 12+ (EAX = EBX
//  at entry; BH = j < 16, so bits 12..15 are 0).
//
//  cbd_vec: SamplePolyCBD_eta(PRF_eta(s, N)) (Algorithm 8 and Sect. 4.1)
//  for N = 0..n-1, n = BH >= 1 at entry (BH counts the polynomials left
//  down to 0, BL = N), stored consecutively
//  from RBP (the frame vector), whose start it pushes first (the bottom
//  of the end list, see dot). RSI = s (kept); BL starts at 0 (both
//  callers set EBX), RCX < 256. eta = eta1, except eta2 for N >= K when
//  BH starts at 2K+1 (encaps' e1 and e2; keygen has 2K and EBX bit 31):
//  only ML-KEM-512 has eta1 != eta2, and sample_cbd picks eta from
//  EAX = EBX (BH = polynomials left, BL = N). Each pass
//  copies N to AL, then INC BL (BL <= 2K) both makes sample_cbd save the
//  next N and leaves the CBD flags SF = ZF = OF = 0 (CF = 0: G's Keccak
//  ends with an equal CMP; at the loop's DEC BH the CF of CBD's last SBB
//  is 0); the CBD exit decrements BH, loops while BH > 0 and otherwise
//  jumps to .Lpd_end (R9 = the last end, then kpke_rows) with BL = 2K or
//  2K+1, BH = 0 and RCX = 0. Each
//  pass writes its polynomial at the top of the end list (the previous
//  end, or the vector start; read from the sampler frame at RBP+24) and
//  pushes the polynomial's end. Each coefficient runs two
//  halves of "negate, then subtract k + b_k for the next eta bits b_k,
//  k = eta..1" from s = R8D = q; with T = eta(eta+1)/2 that is
//  -(-s - T - x) - T - y = s + x - y. The value is
//  negative after the first half only (0 < s < 2^31), so its sign ends
//  the coefficient; it is stored as a dword: q + x - y, in
//  [q - eta, q + eta] (not reduced; polymul_acc accepts that, see there).
//
//  PRF_eta needs 64 * eta <= 192 bytes, more than the SHAKE256 rate (136)
//  for eta = 3: the first block is copied to RSP, just below the state,
//  and permuting the state puts the second block right after it.
//  SampleNTT and CBD share the keccak_absorb call and that permutation
//  call (SampleNTT's next block; it also runs the harmless block copy),
//  with PUSHF / POPF keeping the mode flags over them.
//
//  Frame: 136 bytes, then the Keccak state at RSP+136. keccak_absorb's
//  input ends where its state starts: the seed (32 bytes || 2-byte suffix)
//  is written at RSP+103; SampleNTT absorbs all 34 bytes with its state at
//  RSP+137 (hence 344 bytes), CBD's DEC RDI drops the suffix's second
//  byte (PRF input s || N) and puts its state at RSP+136, right after the
//  first-block copy.

//  hash_G_cbd: G (hash_G, output RSI), then cbd_vec (hash_G returns to
//  it) with RSI = the second half of G (sigma or r) as s; BH (the
//  number of polynomials) survives G. Shared by keygen and encaps, which continue with
//  kpke_rows (R10 = rho or ek).

FUNC(hash_G_cbd)
    call    hash_G                      // returns to cbd_vec
ENDF(hash_G_cbd)

FUNC(cbd_vec)
    push    rbp                         // start of the vector
.Lcv_poly:
    push    rsi                         // s
    mov     eax, ebx                    // N = BL, AH = BH
    inc     bl                          // the next N (sample_cbd saves it);
                                        // SF = ZF = OF = 0 (BL <= 2K), CF = 0
                                        // select CBD
ENDF(cbd_vec)

FUNC(sample_cbd)
    push    rbx
    push    0
    pop     rbx                         // EBX = 0 (flags kept)
    mov     r11w, SHAKE128              // decode threshold, SampleNTT mode
    enter   344, 0                      // 136 + 200 bytes (+1, SampleNTT)
    jl      .Lsn_decode                 // no flags changed since entry
    lea     rdi, [rsp+136-33]           // the seed ends at the state
    push    rdi
    mov     cl, 32                      // RCX < 256 at entry
    rep     movsb
    stosw                               // seed = s || N, or rho || j || i
    pop     rsi
    pushf                               // CF = 1: SampleNTT, kept over Keccak
    jb      .Lsn_absorb

    dec     rdi                         // s || N: the state starts at N's
                                        // successor (RSP+136)
#if ETA1 != ETA2
    //  eta = eta2 for encaps' e1 / e2 (N >= K, so BH = 2K+1-N <= K+1),
    //  else eta1 (keygen's EAX has bit 31 set): CF = 1 iff EAX < (K+2):0.
    cmp     eax, (MLKEM_K+2) << 8       // CF iff encaps' N >= K: EAX = (2K+1-N)
                                        // << 8 | N (keygen: bit 31 set)
    sbb     ebx, -ETA1                  // EBX = eta1 - CF (eta2 = eta1 - 1)
#else
    mov     bl, ETA1                    // eta1 = eta2
#endif
    mov     r11b, SHAKE256 & 0xFF       // same padding byte as SHAKE128
.Lsn_absorb:
    call    keccak_absorb               // CBD: SHAKE256 (rate 136), SampleNTT:
                                        // SHAKE128 (rate 168, EBX = 0 counts
                                        // accepted coefficients)
    mov     r11d, r8d                   // SampleNTT's threshold R11W = q
                                        // (unused by CBD)
    popf
    jb      .Lsn_block
    push    rdi
    pop     rsi                         // state
    push    rsp
    pop     rdi
    mov     cl, 136                     // RCX = 0 after Keccak
    rep     movsb                       // first block to RSP, RDI = state
    //  Shared with SampleNTT's next block, which arrives with ZF = 1 (CBD
    //  has ZF = 0 from entry: MOVS / PUSH / POP keep the flags).
.Lsn_permute:
    pushf
    call    keccak_f1600                // CBD: second block right after the
    popf                                // first
    je      .Lsn_block
    mov     rdi, [rbp+24]               // polynomial: top of the end list
    xor     esi, esi                    // bit position
    //  DL = 0 (Keccak's chi leaves it) counts the coefficients.
.Lcbd_coef:
    mov     eax, r8d                    // s = q
.Lcbd_half:
    neg     eax                         // -s, then s + T + sum(x)
    mov     ecx, ebx
.Lcbd_bit:
    bt      dword ptr [rsp], esi
    inc     esi                         // keeps CF
    sbb     eax, ecx                    // subtracts ECX + bit, SF = sign
    loop    .Lcbd_bit
    js      .Lcbd_half                  // < 0 after the x bits only
    stosd                               // q + sum(x) - sum(y)
    inc     dl                          // wrap after coefficient 255
    jnz     .Lcbd_coef
    leave                               // RCX = 0 (LOOP), RDI = end
    pop     rbx
    pop     rsi                         // s
    push    rdi                         // end of this polynomial
    dec     bh                          // BH = polynomials left
    jnz     .Lcv_poly
    jmp     .Lpd_end                    // R9 = last end, then kpke_rows

.Lsn_decode:
    lea     r10, [rsi+rcx*8]            // RSI = R10 = block, R10 = next block
                                        // (RCX = 48: dot)
    cdq                                 // EDX = 0: EAX = dot's EBX >= 0
    .byte   0x3d                        // cmp eax, imm32: skips the MOV DX;
                                        // R11W = 0x1FA7: accept all
    //  SampleNTT: Keccak returns RSI = state + 200; get_bits reads at the
    //  signed 16-bit offset DX, so EDX bits 16+ do not matter.
.Lsn_block:
    mov     dx, -1600                   // DX = -1600: bit 0 of the state
    mov     cl, 12                      // RCX = 0 after Keccak; get_bits keeps it
.Ldec_coef:
.Lsn_candidate:
    cmp     dh, 0xFF                    // SampleNTT: DX = -256 after 1344 bits
                                        // (decode, poly_decode: DH < 64)
    je      .Lsn_permute
get_bits:
    push    rcx
.Lgb_bit:
    bt      word ptr [rsi], dx
    rcr     eax, 1
    inc     edx
    loop    .Lgb_bit
    pop     rcx
    rol     eax, cl
    cmp     edx, ebx                    // EDX > EBX: sampler (see above)
    ja      .Lsn_accept                 // (SampleNTT / decode)
    imul    eax, r8d
    shr     eax, cl                     // carry = original product bit d-1
    adc     al, 0                       // Decompress_d rounding: no carry
                                        // out of AL for d = 4, 5, 10, 11
    stosd
ENDF(sample_cbd)
FUNC(poly_decode)
    xor     eax, eax                    // get_bits keeps the bits above d
    inc     bl                          // 256 coefficients from BL = 0
    jnz     .Ldec_coef
    push    rdi                         // end of this polynomial
    shl     bh, 1                       // a single bit from 1 << (6-K)
    jno     .Ldec_coef                  // below bit 7: same d
    mov     cl, DV
    jae     .Ldec_coef                  // bit 7: v; then BH = 0
.Lpd_end:
    pop     r9                          // = RDI: dot's scratch
    //  Fall through into kpke_rows.
ENDF(poly_decode)

//  dot: one row of the products of K-PKE.KeyGen, K-PKE.Encrypt and
//  K-PKE.Decrypt, in the normal domain, looping over the rows (see
//  kpke_rows, its only entry). [RSP] = start of the accumulator: the
//  producers (cbd_vec, poly_decode) push the start of the vector and then
//  the end of every polynomial they write, .Lpd_end pops the last one
//  into R9 (scratch, kept), so the stack then lists the starts of all
//  polynomials, the accumulators from the last one down on top; each row
//  pops its own. In all three operations exactly K entries remain, the
//  starts of vec[K-1] .. vec[0] (keygen 2K+1 pushed and K rows, encaps
//  2K+2 and K+1, decaps K+2 and 1), popped before the RET. vec[j] at
//  RBP+POLY*j, R10 = rho, ek or dk.
//  acc += sum_j M[j] * vec[j], j = 0..K-1 (BH), with
//      BL < K:  M[j] = NTT^-1(SampleNTT(rho || BL || j)) = A[j][BL],
//               rho = R10 (encaps rows)
//      EBX bit 31 set, BL = i < K: M[j] = NTT^-1(SampleNTT(rho || j || i))
//               = A[i][j] (keygen; the seed bytes are swapped)
//      BL = K or 255-K:
//               M[j] = +-NTT^-1(ByteDecode_12(R10 + 384j)), t (EBX bits
//               16+ = 0, encaps) or -s (EBX bits 16+ = DEC_NEG, decaps;
//               see sample_cbd and invntt); R10 then ends 384K bytes
//               further (rho, or ek inside dk).
//  All rows run through sample_cbd with the flags of ADD AL, 128-K / ADD
//  AL, 128+K, which keep AL (AL = BL, or AL = BH = j < K after the keygen
//  swap): SF != OF iff AL = K or 255-K (decode), CF = 1 for AL < K
//  (SampleNTT; the sum wraps).
//  dot keeps EBX bits 16+ (encaps keeps them 0, as frame_enter sets
//  them; keygen's bit 31 does not reach invntt, see sample_cbd). Needs
//  RCX < 256 (callers leave 0 or DV); each term sets CL = 48 for
//  sample_cbd's decode.
//  Returns BH = K (ML-KEM-512: BH = 0, see below), RDI = RSI = the last
//  accumulator and RCX = 0 (its RET is the routine's only exit). For K = 2 the term loop toggles BH with
//  XOR BH, 1: every producer leaves BH = 0 (cbd_vec counts it down,
//  poly_decode shifts its bit out) and every row returns it to 0.

//  kpke_rows: the rows of K-PKE.KeyGen / K-PKE.Encrypt / K-PKE.Decrypt,
//  one dot per row from the last accumulator down (each row pops its start
//  from the end list, see dot), while BL >= 0
//  after DEC BL. Entered from cbd_vec (via .Lpd_end) when it is done (BL = 2K or
//  2K+1 and BH = 0, so BL = K-1 .. 0 for keygen: t[i] = e[i] + sum_j A[i][j] s[j],
//  EBX bit 31 selecting the transposed seed; BL = K .. 0 for encaps:
//  v = e2 + t^T y first, dot decodes t from ek and leaves R10 at rho,
//  then u[i] = e1[i] + sum_j A[j][i] y[j]) or by poly_decode (BL = 0, so 255-K:
//  decaps' single decode row, w = v - s^T u). Returns as dot does, with
//  BL = 0xFF (keygen, encaps) and OF = 0 (DEC BL from 0).

FUNC(kpke_rows)
    sub     bl, MLKEM_K+1               // BH = 0 (XOR BH, 1 below keeps it)
ENDF(kpke_rows)

FUNC(dot)
.Ldot_row:                              // [RSP] = accumulator
#if MLKEM_K != 2
    mov     bh, 0
#endif
.Ldot_term:
    mov     cl, 48                      // RCX < 256: sample_cbd's decode
                                        // advances R10 by 8 * 48 = 384
    push    r10
    pop     rsi                         // rho, or the next 384-byte block
    mov     eax, ebx
    test    eax, eax                    // keygen: EBX bit 31 set
    jns     .Ldot_seed                  // encaps: seed i || j
    xchg    al, ah                      // keygen: seed j || i, AL = j < K
.Ldot_seed:
    add     al, 128-MLKEM_K             // decaps: BL = 255-K
    add     al, 128+MLKEM_K             // AL kept; decode: SF != OF, S: CF
    call    sample_cbd                  // R9 = A[j][i], t[j] or -s[j]
    pop     rsi
    push    rsi                         // accumulator
    push    rbp
    pop     rcx                         // a = vec[BH] (polymul_acc reads BH)
    call    polymul_acc
    xor     ecx, ecx
#if MLKEM_K == 2
    xor     bh, 1                       // j = 0, 1, then BH = 0 again
    jnz     .Ldot_term
#else
    inc     bh
    cmp     bh, MLKEM_K
    jb      .Ldot_term
#endif
    pop     rdi                         // RDI = RSI = accumulator
    dec     bl
    jns     .Ldot_row                   // the next row (kpke_rows)
    .rept   MLKEM_K
    pop     rax                         // the vector's K starts
    .endr
    ret
ENDF(dot)

//  The sampler's accept step (see sample_cbd), reached from get_bits with
//  EDX > EBX; it ends by falling through into invntt.

.Lsn_accept:
    cmp     ax, r11w                    // CF = 1: accept
    mov     [r9+rbx*4], eax             // a rejected value is overwritten later;
                                        // bits 16+ (DEC_NEG) are ignored: the
                                        // transform reads 16-bit inputs
    adc     bl, bh                      // BH = 0; CF after the 256th accept
    jnc     .Lsn_candidate
#if MLKEM_K != 2
    leave                               // RBX stays saved: ntt_common pops it
#endif
    //  Fall through to invntt (R9 = polynomial): SampleNTT returns
    //  A = NTT^-1(A^), decode t = NTT^-1(t^) or -s = -NTT^-1(s^).

//  ========================================================================
//  NTT and polynomial multiplication (ntt_common)
//  ========================================================================

//  invntt: R9 = polynomial, ntt: RSI = polynomial, transformed in place
//  (both return it in RSI, keep RDI and return RCX = 0, which keygen and
//  dot rely on; inside, RCX = 0 selects the transform and RSI < R9 the
//  forward direction: invntt has RSI = R9 (or R9 + 4, see .Lnt_body),
//  keygen's s and t lie below dot's scratch R9; ntt needs RCX = 0 on
//  entry, which keygen has). invntt's only entry is sample_cbd's fall-through,
//  which arrives with RCX = 12 and with sample_cbd's RBX still pushed: its
//  AND ECX, 0x53 clears RCX and takes ntt_common's PUSH RBX as the imm8,
//  so ntt_common's final POP RBX restores the caller's RBX for
//  sample_cbd (no extra stack). It sets AX = 3303, so
//  c = 3303 + 65536 * (EAX >> 16): sample_cbd leaves EAX >> 16 = 0
//  (c = 128^-1) except in decaps' decode rows, where it is DEC_NEG
//  (c = -128^-1; 17c = 0x8D72DB57 < 2^32, so the 32-bit IMUL is exact).
//  ntt is the forward entry; its only caller (keygen) passes EAX = c = 1.
//  Inputs may be up to 4095 (unreduced ByteDecode_12).
//  polymul_acc: RSI += a * R9 in Z_q[X]/(X^256+1) with a = RCX + POLY*BH
//  (dot passes RCX = RBP and BH = j: a = vec[j]); RSI must not alias a
//  or R9; RCX nonzero selects this mode. b = R9 <= q (b is
//  always a transform output); a and the accumulator may be up to q + 3
//  (CBD samples are left as q + v). The (unused) table is
//  still built from EAX, which needs bit 31 clear: the table loop's CDQ
//  (before the IMUL, which truncates and keeps EDX) must give EDX = 0 for
//  DIV, and later passes see remainders < q. Every transform c qualifies
//  (1, 3303, DEC_NEG << 16 | 3303 = 0x08520CE7), and dot calls polymul_acc
//  right after invntt, whose last DIV leaves EAX = quotient < 2^32 / q.
//  a, b and the accumulator are read as dwords
//  (clean values below 2^16); transform inputs as words.

FUNC(invntt)
    mov     ax, 3303                    // c = 128^-1
    push    r9
    pop     rsi                         // polynomial
#if MLKEM_K == 2
    //  Stack trade (ML-KEM-512 only, about +340 bytes of stack, peak still
    //  set by ML-KEM-1024): AND ECX, imm32 swallows all of ntt_common's
    //  prologue (PUSH RBX / PUSH RBP / PUSH RSP / POP RBP), so the transform
    //  runs in sample_cbd's still-open frame (no LEAVE above): RBP+8 holds
    //  the caller's RBX there too, and ntt_common's LEAVE / POP RBX unwind
    //  both.
    .byte   0x81, 0xe1                  // RCX = 12 & 0x5d545553 = 0
#else
    .byte   0x83, 0xe1                  // and ecx, imm8 0x53 = PUSH RBX below:
                                        // RCX = 12 & 0x53 = 0 (transform)
#endif
ENDF(invntt)

FUNC(ntt)
FUNC(ntt_common)
FUNC(polymul_acc)
    push    rbx                         // 0x53: invntt's imm8, skipped there
    push    rbp                         // RBP restores RSP after the table
    push    rsp
    pop     rbp
    xor     ebx, ebx                    // BL = output index x (and counter)
.Lnt_tab:
    cdq                                 // EDX = 0: EAX < 2^31 (below)
    imul    eax, eax, 17                // exact for every transform c
    div     r8d
    xchg    eax, edx
    push    rax
    inc     bl
    jnz     .Lnt_tab
    cdq                                 // EAX < q: EDX = 0 (polymul's i)
    mov     dh, [rbp+9]                 // caller's BH (polymul's j)
    //  x outputs are already pushed, so T[e] is at [rsp+8(e+x)].
.Lnt_out:
    jrcxz   .Lnt_body
    //  polymul: slot = acc[x] + sum_i a[i] * b[x-i] in Z_q[X]/(X^256+1);
    //  terms that wrap around use q - b[x-i]. With acc, a <= q+3 and
    //  b <= q the sum is below q+4 + 256*(q+3)*q < 2^32, so a 32-bit ADD
    //  suffices.
    //  EDX = 256j + i with i = 0 here: the inner loop wraps DL back to 0,
    //  so [RCX+RDX*2] = a[i].
    push    qword ptr [rsi+rbx*4]   // acc[x] (the high dword is junk)
.Lsb_term:
    mov     eax, ebx
    sub     al, dl                      // x - i mod 256, CF = wrapped
    mov     eax, [r9+rax*4]
    jnc     1f
    neg     eax
    add     eax, r8d                    // q - b
1:  imul    eax, [rcx+rdx*4]
    add     [rsp], eax
    inc     dl
    jnz     .Lsb_term
.Lnt_next:
    inc     bl
    jnz     .Lnt_out
.Lnt_pop:
    //  CDQ before the POP: EAX is the previous quotient (< 2^32 / q), or on
    //  entry the transform's table index (< 512) or polymul's last product
    //  (< 2^24), so EDX = 0; the slot itself may reach 2^32 (polymul).
    cdq
    pop     rax
    div     r8d                         // the slot's high dword is ignored
    dec     bl
    mov     [rsi+rbx*4], edx
    jnz     .Lnt_pop
.Lnt_leave:                             // also frame_enter's epilogue
    leave
    pop     rbx
    ret

    //  Transform output x = 2i+b: slot = sum over n = 2y of in[n+b] * T[e].
    //  RSI is offset by 4b: it toggles 256 times and ends unchanged. This
    //  needs bit 2 of RSI clear; all polynomials are at 0 mod 8 in
    //  frames whose base is 8-aligned (RSP is always a multiple of 8).
    //  inv: e = o(y)*j (A = input index), fwd: e = o(127-i)*y = -o(i)*y
    //  (A = output index), so the roles swap.
.Lnt_body:
    push    rcx                         // zero accumulator slot
.Lnt_in:
    mov     eax, ecx                    // inv: A = n, B = x
    mov     edx, ebx
    cmp     rsi, r9                     // inv: RSI = R9 or R9+4
    jae     1f
    xchg    eax, edx                    // fwd: A = x, B = n
    xor     al, 0xfe                    // A = 127 - i: o(127-i) = -o(i)
1:  mov     ah, 1                       // sentinel: EDX = B'<<9 | ? | o(A)
2:  shr     eax, 1
    rcl     edx, 1
    jnz     2b
    shld    eax, edx, 22                // B = B'/2
    mul     dl
    mov     ah, 0                       // e = o(A)*B mod 256
    add     eax, ebx
    movzx   edx, word ptr [rsi+rcx*4]
    imul    edx, [rsp+rax*8+8]
    add     [rsp], edx                  // 128 terms < 2^31
    add     cl, 2
    jnc     .Lnt_in                     // ends with RCX = 0 again
    xor     rsi, 4
    jmp     .Lnt_next
ENDF(polymul_acc)
ENDF(ntt_common)
ENDF(ntt)


//  ========================================================================
//  ML-KEM (Section 6), exported
//  ========================================================================

//  frame_enter: shared prologue/epilogue. An exported function starts with
//  "call frame_enter"; the rest of that function runs as a subroutine with
//  RBX and RBP saved (bodies set EBX themselves), RAX = 1st, RSI = 2nd and
//  R10 = 3rd argument (RDX holds the body's address: POP RDX / CALL RDX
//  are a byte shorter each than with R11), R8 = q for the whole operation
//  (ntt_common, the encoders and the samplers rely on it; nothing changes
//  R8), RBP = header + HDR (see the frame layout) for a FRAME-byte
//  scratch area that starts with the header. The 32 bytes at RCX (4th
//  argument: m in encaps, d in keygen, a dummy in decaps) are
//  copied to the header first, so bodies start with RDI = header + 32 and
//  RCX = 0. Bodies keep R12-R15.
//  The body's RET returns here, which restores everything (RBP = top from
//  the slot below the frame, then ntt_common's LEAVE / POP RBX / RET) and
//  returns to the caller. The frame base is 8-aligned (see ntt_common).

//  FRAME >= HDR + (2K+2)*POLY (vector, accumulators, dot's scratch) and
//  >= 32 + CT_BYTES + 200 (J's input z || c and Keccak state), 8-aligned.
//  decaps runs the encaps body in its own frame, so one frame is live.
#define FRAME       (HDR+(2*MLKEM_K+2)*POLY)

FUNC(frame_enter)
    mov     r10, rdx                    // caller-saved: R10 is never restored
    mov     r8d, Q                      // R8 = q for the whole operation
    pop     rdx                         // body address
    push    rbx
    xchg    rax, rdi                    // RAX = 1st argument
    enter   FRAME, 0                    // push rbp, RBP = top, allocate
    push    rsp
    pop     rdi                         // RDI = frame base (header)
    push    rbp                         // [base-8] = top
    push    rsi
    push    rcx
    pop     rsi
    push    32
    pop     rcx
    rep     movsb                       // 4th argument into the header
    pop     rsi
    lea     rbp, [rdi+HDR-32]           // RBP = vector, after the header
    call    rdx
    pop     rbp                         // RBP = top, then LEAVE restores
    jmp     .Lnt_leave                  // RSP, RBP and RBX (ntt_common's tail)
ENDF(frame_enter)


//  mlkem_encaps: ML-KEM.Encaps_internal(ek, m) (Algorithm 17) with
//  K-PKE.Encrypt (Algorithm 14) inlined, in the normal domain:
//  u[i] = e1[i] + sum_j A[j][i] * y[j] and v = e2 + sum_j t[j] * y[j] + mu,
//  so each accumulator starts as a sampled noise polynomial.
//  RDI = ciphertext, RSI = shared key K, RDX = ek, RCX = m.
//  .Lenc_body (also decaps' re-encryption, in decaps' frame): RAX = c, RSI =
//  K, R10 = ek, m in the header, RDI = header + 32, RCX < 256, any EBX.

//  Frame: header m (0), H(ek) (32), K || r (64, 96) = G(m || H(ek));
//  y[K], e1[K] -> u[K] and e2 -> v from RBP (the order of c = c1 || c2),
//  then the scratch polynomial.

    .globl  mlkem_encaps
FUNC(mlkem_encaps)
    call    frame_enter                 // copies m
.Lenc_body:
    push    rax                         // c, popped for c1
    push    rsi                         // K: G's output
    push    r10
    pop     rsi                         // ek
    call    hash_H                      // m || H(ek), RDI = header+64
    pop     rsi                         // K

    //  (K, r) = G(m || H(ek)): K goes to the output RSI; then
    //  y[j] = CBD_eta1(PRF(r, j)), e1[i] = CBD_eta2(PRF(r, K+i)) and
    //  e2 = CBD_eta2(PRF(r, 2K))
    //  then the rows v = e2 + t^T y and u[i] = e1[i] + sum_j A[j][i] y[j]
    xor     ebx, ebx
    mov     bh, 2*MLKEM_K+1             // y[K], e1[K], e2 from RBP
    call    hash_G_cbd                  // returns BH = K (0 if K = 2), RSI = u[0]

    //  c1[i] = ByteEncode_du(Compress_du(u[i])), i = 0..K-1 (BH)
    pop     rdi                         // c
.Lenc_c1:
    mov     cl, DU                      // RCX = 0 after dot, then DU
    call    poly_compress_encode        // RSI ends at the next polynomial
#if MLKEM_K == 2
    neg     bl                          // dot left BL = 0xFF: 0x01 (SF = 0),
    jns     .Lenc_c1                    // then 0xFF again: two passes
#else
    dec     bh
    jnz     .Lenc_c1
#endif
    //  v += Decompress_1(m), m at RBP-HDR. The row loop left BL = 0xFF
    //  and the c1 loop BH = 0 (poly_compress_encode keeps RBX), so EBX =
    //  255..0 indexes the coefficients.
.Lenc_mu:
    bt      dword ptr [rbp-HDR], ebx
    jnc     .Lenc_mu_next
    add     word ptr [rsi+rbx*4], (Q+1)/2
.Lenc_mu_next:
    dec     ebx
    jns     .Lenc_mu
    mov     cl, DV                      // ECX = DV (it was DU)
    //  Fall through to poly_compress_encode, which returns to frame_enter.
ENDF(mlkem_encaps)

//  poly_compress_encode: ByteEncode_d(Compress_d) of the polynomial at RSI
//  into RDI (32*d bytes), RCX = d (all 64 bits: LOOP counts the bits of
//  each coefficient from it). poly_encode is the same entry with RCX = 0
//  (keygen, after ntt): ByteEncode_12 without compression (d = 12 is
//  loaded per coefficient; the pushed RCX = 0 is restored). Compression
//  needs R8D = q (set by frame_enter). Inputs may exceed q (up to q + 1664 for
//  v + mu; adding q adds 2^d, which the packer drops). Compress_d(x) =
//  round(x * 2^d / q) is computed exactly as floor((x * 2^d + (q-1)/2) / q):
//  q is odd, so there are no ties. Only the low d bits of each value are
//  packed.
//  Each bit is rotated into the top of the output byte at RDI (RCR through
//  the carry: after 8 bits all of its old contents are gone); BL += 32
//  carries every 8 bits, then RDI moves to the next byte (SCASB). BH counts
//  256 coefficients. The output byte is read and rotated before it is
//  complete, so the output must not overlap the input (no caller does).
//  Returns RDI = end of output, RSI = end of the polynomial and RCX = d
//  (0 for poly_encode); keeps RBX, clobbers RAX and RDX.

FUNC(poly_compress_encode)
FUNC(poly_encode)
    push    rbx
    xor     ebx, ebx
.Lenc_coef:
    lodsd
    push    rcx
    jrcxz   .Lenc_raw                   // ByteEncode_12: bits 16+ unused
    shl     eax, cl                     // < 2^24
    add     eax, (Q-1)/2
    cdq
    div     r8d
.Lenc_bit:
    shr     eax, 1
    rcr     byte ptr [rdi], 1
    add     bl, 32
    jnc     .Lenc_next
    scasb                               // RDI += 1 (flags unused)
.Lenc_next:
    loop    .Lenc_bit
    pop     rcx
    dec     bh
    jnz     .Lenc_coef
.Lpe_ret:
    pop     rbx
    ret
.Lenc_raw:
    mov     cl, 12                      // RCX = 0 is pushed and restored
    jmp     .Lenc_bit
ENDF(poly_encode)
ENDF(poly_compress_encode)

//  ========================================================================
//  Hashing (G, H, J); mlkem_decaps and mlkem_keygen sit among these
//  routines for short jumps (decaps falls into hash_hdr, keygen ends with
//  a jump into xhash_out)
//  ========================================================================

//  keccak_absorb's input always ends where its state starts (RDI).
//  hash_hdr: hash the bytes from RBP-HDR up to RDI with mode R11W into the
//  output at [RSP], via xhash_out with the Keccak state at RDI, right
//  after the input (no return address of its own: mlkem_decaps falls in
//  for J, the state then lies in decaps' dead frame, with all of hash_G
//  as its MOV R11's immediate; hash_G falls in).
//  hash_G: G = SHA3-512 of the bytes from RBP-HDR up to RDI, with the
//  state at RDI (its 200 bytes reach into the vector, which is sampled
//  later: s or y), so G's 64 bytes are in place there; the first 32 are
//  also copied to the output pointer RSI (encaps passes its K output,
//  keygen its dk, whose first 32 bytes are written again later). It sets
//  R11W = SHA3_512 itself; its only caller is hash_G_cbd. Returns
//  RSI = RDI + 32 (the second half: sigma or r).
//  hash_H: H = SHA3-256 (32 bytes) of the EK_BYTES-byte ek at RSI into
//  RDI. ek is copied to RBP and the state follows the copy (s, dead in
//  keygen's tail, or y, sampled later in encaps); RCX < 256.
//  Both return RDI = output + 32, like xhash_out.

FUNC(hash_H)
    mov     r11w, SHA3_256
    push    rdi                         // output
    push    rbp
    pop     rdi
    mov     cl, EK_BYTES/8              // RCX < 256 in all callers
    rep     movsq                       // ek to RBP, the state right after it
    push    rbp
    pop     rsi
    jmp     xhash_out
ENDF(hash_H)


//  mlkem_decaps: ML-KEM.Decaps_internal(dk, c) (Algorithm 18) with
//  K-PKE.Decrypt (Algorithm 15) inlined, in the normal domain. Re-encryption
//  calls the encaps body (.Lenc_body) in this frame, which recomputes
//  h = H(ek) instead of reading it from dk; the two are equal for every dk that passes the FIPS 203
//  decapsulation key check. Variable time.
//  RDI = shared key K, RSI = ciphertext, RDX = dk.

//  DEC_NEG * 65536 + 3303 = -3303 = -128^-1 (mod q): see invntt and
//  sample_cbd.
#define DEC_NEG     0x852

//  Frame: u[K] from RBP and the accumulator w, as in encaps; then m' in
//  the header and the encaps body's frame, with c' at header+32 (its
//  H(ek) || K' || r slots and the dead y, all free when it encodes; m'
//  before it is still read by the mu step); then z (copied again from dk)
//  in the header, so z || c' is in place, and c is compared with c' and
//  copied over it byte by byte up to the first difference (LODSB / XCHG /
//  SCASB / LOOPE), then the rest with REP MOVSB: z || c for J.

    .globl  mlkem_decaps
FUNC(mlkem_decaps)
    push    rsi
    pop     rcx                         // c: any readable 32 bytes
    call    frame_enter
    push    rax                         // K: J's output, or popped at the end
    push    rsi                         // c, for the compare and the copy
    push    rax                         // K, for the encaps body

    //  u[i] = Decompress_du(ByteDecode_du(c1[i])) at RBP, as y in encaps
    push    rbp
    pop     rdi                         // RSI = c (frame_enter keeps it)
    //  then v = Decompress_dv(ByteDecode_dv(c2)) into the accumulator, and
    //  poly_decode falls into kpke_rows. BH = 64 >> K: K polynomials with
    //  DU before v (see poly_decode); BL = 0xFF and EDX = 0 for poly_decode's entry, then
    //  BL = 0 (wrapped), so BL = 255-K in kpke_rows, makes one row that
    //  decodes s from dk in dot.
    //  EBX bits 16+ = DEC_NEG make the invntt after sample_cbd's decode
    //  use c = -128^-1: dot adds (-s)^T u.
    mov     ebx, DEC_NEG << 16 | (64 >> MLKEM_K) << 8 | 0xFF
    xor     edx, edx
    mov     cl, DU                      // RCX = 0 after frame_enter
    call    poly_decode                 // u[K], v, then the row: RSI = acc

    //  acc = v - s^T u = w; m' = ByteEncode_1(Compress_1(w)) into the
    //  header, where the encaps body expects m.
    lea     rdi, [rbp-HDR]
    mov     cl, 1                       // RCX = 0 after dot
    call    poly_compress_encode        // RDI = header+32

    //  (K', c') = Encaps_internal(ek, m') in this frame: K' straight to the
    //  output, c' at header+32; on a mismatch K = J(z || c).
    pop     rsi                         // K
    push    rdi
    pop     rax                         // c' at header+32
    call    .Lenc_body                  // R10 = rho inside dk on return
    lea     rsi, [r10+64]               // z, after rho || H(ek)
    lea     rdi, [rbp-HDR]
    mov     cl, 32                      // RCX < 256 after the encaps body
    rep     movsb                       // z || c', RDI = c'
    //  (CT_BYTES = 768 for ML-KEM-512: MOV CH, 3)
    pop     rsi                         // c
#if CT_BYTES % 256
    mov     cx, CT_BYTES
#else
    mov     ch, CT_BYTES/256            // RCX = 0 after REP MOVSB
#endif
.Ldec_cmp:
    lodsb                               // compare c with c' and copy c over
    xchg    al, [rdi]                   // it up to the first difference
    scasb
    loope   .Ldec_cmp                   // ZF = 1: equal
    je      .Lpe_ret                    // equal: K = K' (encaps wrote it);
                                        // pops the K pointer into RBX
    rep     movsb                       // the rest: z || c
    //  K = J(z || c): the copy ended at RDI, the output K is at [RSP].
    //  MOV R11, imm64 sets R11W = SHAKE256 (keccak_absorb ignores R11's
    //  upper bits); the imm64's top six bytes are all of hash_G (its MOV
    //  R11W, SHA3_512 and PUSH RSI), skipped.
    .byte   0x49, 0xbb, SHAKE256 & 0xff, SHAKE256 >> 8
ENDF(mlkem_decaps)

FUNC(hash_G)
    mov     r11w, SHA3_512              // (with the PUSH: decaps' imm64)
    push    rsi                         // output, above the return address
ENDF(hash_G)

FUNC(hash_hdr)
    lea     rsi, [rbp-HDR]              // the input ends at the state (RDI)
    //  Fall through to xhash_out.
ENDF(hash_hdr)

//  xhash_out: one-shot SHA3 / SHAKE; the first 32 bytes of the first
//  block are copied to the output pointer at [RSP] (popped on return).
//  RDI = Keccak state (200 bytes of scratch; the output may be the state
//  itself), RSI = input, which ends at RDI (at least 1 byte), R11W = mode,
//  RCX < 256. Returns RSI = state + 32,
//  RDI = output + 32
//  and RCX = 0; preserves R8-R10 (not R11). mlkem_keygen jumps to .Lxo_copy to copy
//  z (RSI) to RDI and return.

FUNC(xhash_out)
    call    keccak_absorb               // RCX = 0
    push    rdi
    pop     rsi                         // RSI = state
    pop     rdi                         // RDI = output
.Lxo_copy:
    mov     cl, 32                      // RCX = 0 after Keccak
    rep     movsb
    ret
ENDF(xhash_out)

//  mlkem_keygen: ML-KEM.KeyGen_internal(d, z) (Algorithms 16 and 13).
//  RDI = ek, RSI = dk, RDX = d, RCX = z (swapped to RCX = d, RDX = z for
//  frame_enter, which leaves z in R10). Computes t = A*s + e in the normal domain
//  (A = NTT^-1(A^)); s^ = NTT(s) and t^ = NTT(t) are encoded. Placed
//  after xhash_out, whose final 32-byte copy it shares (z).

//  Header (at RBP-HDR): d || k (33 bytes), then rho || sigma = G(d || k);
//  s[K], then e[K] -> t[K] from RBP.
//  BH = B + n*KG_STEP after n encode steps (B = K left by the rows, 0 for
//  K = 2): below 0x80 for n < K, in [0x80, 0x100) for K <= n < 2K, carry
//  at n = 2K (K = 2, 3, 4).
#define KG_STEP     ((128+MLKEM_K-1)/MLKEM_K)

    .globl  mlkem_keygen
FUNC(mlkem_keygen)
    xchg    rcx, rdx                    // frame_enter copies d
    call    frame_enter
    push    r10                         // z, popped last
    push    rax                         // ek
    push    rsi                         // dk

    //  (rho, sigma) = G(d || k)
    mov     al, MLKEM_K
    stosb                               // d || k ends at rho
    push    rdi
    pop     r10                         // rho, where dot reads it
    //  G's output pointer RSI is dk: rho there is overwritten by s^ later

    //  then s[i] = CBD(PRF(sigma, i)), e[i] = CBD(PRF(sigma, K+i)) and the
    //  rows t[i] = e[i] + sum_j A[i][j] * s[j] (EBX bit 31: seed j || i)
    mov     ebx, 1 << 31 | 2*MLKEM_K << 8   // s[K], e[K] from RBP
    call    hash_G_cbd                  // input d || k: 33 bytes

    //  dk = ByteEncode_12(s^) || ByteEncode_12(t^) || rho || H(ek) || z,
    //  where the middle part ByteEncode_12(t^) || rho is ek.
    pop     rdi                         // dk
    push    rbp
    pop     rsi                         // s, then t
    //  The rows left BH = K (0 if K = 2) and OF = 0 (DEC BL from 0). BH steps by
    //  KG_STEP: the K-th step is the only one that crosses 0x7F (OF = 1),
    //  the 2K-th the only one that carries out (see KG_STEP).
.Lkg_enc:
    jno     .Lkg_cursor
    push    rdi                         // ek inside dk, after s^
.Lkg_cursor:
    push    1
    pop     rax                         // c = 1
    call    ntt                         // s^, t^ (RSI = polynomial < R9)
    call    poly_encode                 // RCX = 0 after ntt: d = 12
    add     bh, KG_STEP
    jnc     .Lkg_enc
    push    r10
    pop     rsi                         // rho
    mov     cl, 32                      // poly_encode left RCX = 0
    rep     movsb                       // ek inside dk is complete
    pop     rsi                         // ek inside dk
    pop     rdi                         // ek
    push    rdi
    mov     cl, EK_BYTES/8              // REP left RCX = 0
    rep     movsq                       // ek = copy from dk; RSI ends at
    push    rsi                         // the H(ek) field
    pop     rdi
    pop     rsi                         // ek
    call    hash_H                      // RDI ends at the z field
    pop     rsi                         // z
    jmp     .Lxo_copy                   // copy 32 bytes and return
ENDF(mlkem_keygen)

//  ========================================================================
//  Keccak
//  ========================================================================

//  keccak_absorb: RDI = state (200 bytes, gets zeroed), RSI = input, which
//  ends exactly at the state (no empty inputs; the length is RDI - RSI),
//  R11W = mode (upper bits ignored; keccak_f1600 keeps it over the
//  mid-input permutations), RCX < 256. EDX is not used (Keccak clobbers
//  DL and DH).
//  Absorbs, pads and permutes, so the first output block is ready in the
//  state. Preserves RDI and R8-R10; R11 is destroyed by the padding (XCHG:
//  no caller needs the mode afterwards).

FUNC(keccak_absorb)
    push    rdi
    mov     cl, 200                     // RCX < 256 in all callers
    xor     eax, eax
    rep     stosb                       // 200 bytes
    pop     rdi                         // RCX = 0: position in block
.Lab_byte:
    lodsb
    xor     [rdi+rcx], al
    inc     ecx
    cmp     cl, r11b                    // block full (CL > rate - 1)?
    jbe     .Lab_next
    push    rsi
    call    keccak_f1600
    pop     rsi                         // RCX = 0 after Keccak
.Lab_next:
    cmp     esi, edi                    // the input ends at the state: the low
    js      .Lab_byte                   // dwords' difference is negative while
                                        // RSI < RDI (RDI - RSI < 2^31, also
                                        // across a 2^32 boundary)
    xchg    eax, r11d                   // R11 is not needed after absorb
    xor     [rdi+rcx], ah               // domain + first padding bit
    mov     cl, al                      // RCX = rate - 1 (RCX < 256 already)
    xor     byte ptr [rdi+rcx], 0x80    // last padding bit
    //  Fall through to the permutation.
ENDF(keccak_absorb)

//  keccak_f1600: permute the 25-lane state at RDI (lane x + 5*y).
//  Round constants come from the LFSR x^8+x^6+x^5+x^4+1 in DH. Its
//  reciprocal right-shift form starts at 0x1d and uses 0x8e feedback: the
//  carry from SHR is the next iota bit. Each round builds its scratch
//  values with pushes; PUSH RBP / LEAVE frame them (about 500 bytes of
//  stack). Needs RCX < 256 (all callers have a small count there).
//  Clobbers RAX, DL, DH, RSI and flags (EDX bits 16+ are kept);
//  returns RCX = 0 (xhash_out, keccak_absorb and sample_cbd rely on it),
//  DL = 0 (the CBD coefficient count), RSI = RDI + 200 and EAX = 127
//  (SampleNTT relies on both: get_bits keeps EAX bits 12+, which must be
//  0). Preserves RDI, R8-R11 and all
//  callee-saved registers. R11 is rho's rotating register: the first
//  step parks it in lane 1 and the last step takes it back unchanged, so
//  it is preserved (keccak_absorb keeps its mode there, and SampleNTT its
//  threshold).

FUNC(keccak_f1600)
    mov     dh, 0x1d
.Lkf_round:
    push    rbp
    push    rsp
    pop     rbp
    //  theta: C[x] by a 5-deep running xor over lanes 24..0 (lanes
    //  24..20 start the columns, so they skip the xor)
    mov     cl, 25                      // RCX < 256 on entry, 0 later
.Lkf_col:
    mov     rax, [rdi+rcx*8-8]
    cmp     cl, 21
    jae     .Lkf_col_top
    xor     rax, [rsp+32]
.Lkf_col_top:
    push    rax
    loop    .Lkf_col
    //  [rsp + 8x] = C[x]. Lanes 24..0: D = C[x-1] ^ rol(C[x+1], 1),
    //  re-pushing C[x] keeps the window C[x+1] at [rsp], C[x-1] at +24.
    //  RAX = [rsp] on every entry (the last push, C[0] at first).
    mov     cl, 25
.Lkf_theta:
    rol     rax, 1
    xor     rax, [rsp+24]
    xor     [rdi+rcx*8-8], rax
    mov     rax, [rsp+32]
    push    rax
    loop    .Lkf_theta
    //  rho and pi: lane i = x + 5y moves to y + 5((2x + 3y) mod 5). EAX =
    //  2i and ESI = 2y, so EAX - ESI = 2x + 8y = 2x + 3y mod 5. Start one
    //  step early at (x, y) = (1, 1), whose successor is lane 1: that step
    //  stores junk there with rotation 0 and picks up A[1]; the last step
    //  overwrites the junk. CL = (t+1)t/2 rotation, CH = step count.
    mov     dl, 5                       // divisor here, chi's row count later
    push    2
    pop     rax                         // 2i = 2 for lane 6 (= 2x mod 5)
    push    rax
    pop     rsi                         // 2y = 2
.Lkf_rho:
    rol     r11, cl
    sub     eax, esi                    // 2x + 8y < 256, = 2x + 3y mod 5
    div     dl                          // AH = y' = (2x + 3y) mod 5, AL < 128
    //  TEST AL, 0xEB only sets flags (none are live here). The round loop's
    //  JNE enters at its imm8: EB C1 = "jmp .Lkf_round", whose rel8 is the
    //  first byte of the SHR (C1 E8 07), so the backward JNE fits in rel8.
    .byte   0xa8                        // test al, imm8
.Lkf_loop:
    .byte   0xeb                        // the imm8; "jmp .Lkf_round" from .Lkf_loop
    shr     eax, 7                      // 2y'
    //  Build-time check (a non-allocated section, no .text bytes): fails
    //  with "division by zero" unless EB C1 lands exactly on .Lkf_round.
    .pushsection .discard.mlkem_check, "", @progbits
    .byte   1 / (.Lkf_loop + 2 - .Lkf_round == 0x3f)
    .popsection
    xchg    eax, esi                    // ESI = 2y', EAX = 2y = 2x'
    lea     eax, [rax+rsi*4]
    add     eax, esi                    // 2i' = 2x' + 10y'
    xchg    [rdi+rax*4], r11
    inc     ch
    add     cl, ch
    cmp     ch, 25
    jne     .Lkf_rho
    //  chi: per row push A1 A0 A4 A3 A2 A1, then pop A[x+1] and read A[x+2]
    push    rdi
    pop     rsi                         // DL = 5 rows (set before rho)
.Lkf_chi_row:
    push    qword ptr [rsi+8]
    push    qword ptr [rsi]
    push    4
    pop     rcx
.Lkf_chi_push:
    push    qword ptr [rsi+rcx*8]
    loop    .Lkf_chi_push
    mov     cl, 5
.Lkf_chi:
    pop     rax
    not     rax
    and     rax, [rsp]
    xor     [rsi], rax
    lodsq
    loop    .Lkf_chi
    dec     dl
    jnz     .Lkf_chi_row
    //  iota: LFSR output bit j toggles lane bit 2^j - 1, j = 0..6
    xor     eax, eax                    // bit position 0, 1, 3, .., 63
    mov     cl, 7
.Lkf_iota:
    shr     dh, 1                       // emit old bit 0 into carry
    jnc     .Lkf_iota_skip
    xor     dh, 0x8e                    // reciprocal feedback if bit was 1
    btc     dword ptr [rdi], eax        // dword form reaches bit 63 too
.Lkf_iota_skip:
    stc
    adc     eax, eax                    // 2 * EAX + 1
    loop    .Lkf_iota
    leave
    cmp     dh, 0xe2                    // LFSR state after 24 rounds
    jne     .Lkf_loop                   // via the hidden jmp in rho
    ret
ENDF(keccak_f1600)

    .section .note.GNU-stack, "", @progbits
